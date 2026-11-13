load_pkg <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Required R package '%s' is missing; recreate workflow/envs/r.yaml", pkg), call. = FALSE)
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

load_pkg(c("vroom", "tibble", "glue", "dplyr", "tidyr", "purrr", "GenomicRanges", "stringr", "BiocParallel", "parallel"))

my_vroom <- function(file, na_append = NULL, delim = "\t", comment = "#", col_names = TRUE, ...) {
  if (!file.exists(file)) stop(sprintf("Annotation file does not exist: %s", file), call. = FALSE)
  if (!file.info(file)$size) {
    if (is.character(col_names)) return(tibble::as_tibble(empty_columns(col_names)))
    return(tibble::tibble())
  }
  # Read the header rather than assuming an AnnotSV release-specific width.
  dat <- vroom::vroom(
    file, delim = delim, col_types = vroom::cols(.default = vroom::col_character()),
    na = c("", "NA", "NaN", na_append), comment = comment, col_names = col_names,
    .name_repair = "check_unique", progress = FALSE, ...
  )
  problems <- vroom::problems(dat)
  if (nrow(problems)) {
    sample <- head(problems, 5L)
    details <- sprintf("row %s, column %s: expected %s; actual %s",
                       sample$row, sample$col, sample$expected, sample$actual)
    stop(sprintf("Malformed annotation table %s: %d parsing problems\n%s",
                 file, nrow(problems), paste(details, collapse = "\n")), call. = FALSE)
  }
  dat
}

read_annotsv_annotations <- function(file, sv_type) {
  # AnnotSV fields contain literal hashes (e.g. phenotype references) and quotes.
  # Treating # as a comment truncates valid rows and loses downstream columns.
  my_vroom(file, quote = "", comment = "") |>
    prepare_annotsv_annotations(sv_type = sv_type)
}

nothing <- function(x) x

create_grs <- function(string) as(string, "GRanges")

calculate_overlap_flag <- function(coord, sv) {
  if (!length(coord) || is.na(coord) || coord %in% c("", ".", "NA")) return(FALSE)
  gr <- coord |>
    stringr::str_replace_all("chr", "") |>
    stringr::str_split(";") |>
    unlist() |>
    create_grs()
  pcts_gr <- width(pintersect(gr, sv)) / width(gr)
  pcts_sv <- width(pintersect(gr, sv)) / width(sv)
  any(pcts_gr >= 0.9 & pcts_sv >= 0.8, na.rm = TRUE)
}

add_overlap_flag_annotsv <- function(tsv, sv_type, workers = detectCores()) {
  columns <- c("B_gain_coord", "B_loss_coord", "po_B_loss_allG_coord", "B_ins_coord", "po_B_gain_allG_coord", "B_inv_coord")
  flags <- c("flag_gain", "flag_loss", "flag_po_loss", "flag_ins", "flag_po_gain", "flag_inv")
  for (flag in flags) tsv[[flag]] <- rep(FALSE, nrow(tsv))
  if (!nrow(tsv)) return(tsv)
  starts <- suppressWarnings(as.numeric(tsv$SV_start))
  if (sv_type == "INS") {
    lengths <- suppressWarnings(as.numeric(tsv$SV_length))
    ends <- starts + lengths - 1
    valid <- is.finite(starts) & is.finite(lengths) & lengths > 0
    if (any(!valid)) warning("Insertion length/position unavailable: benign-overlap checks skipped for those rows", call. = FALSE)
  } else {
    ends <- suppressWarnings(as.numeric(tsv$SV_end))
    valid <- is.finite(starts) & is.finite(ends) & ends >= starts
    if (any(!valid)) stop("AnnotSV contains invalid SV_start/SV_end coordinates", call. = FALSE)
  }
  rows <- which(valid)
  if (!length(rows)) return(tsv)
  gr_sv <- GRanges(seqnames = sub("^chr", "", as.character(tsv$SV_chrom[rows])), ranges = IRanges(starts[rows], ends[rows]))
  # Explicit parameter avoids changing the process-wide BiocParallel registry.
  parallel <- if (workers > 1) MulticoreParam(workers = workers) else SerialParam()
  for (j in seq_along(columns)) {
    tsv[[flags[j]]][rows] <- unlist(bplapply(seq_along(rows), function(i) {
      calculate_overlap_flag(tsv[[columns[j]]][rows[i]], gr_sv[i])
    }, BPPARAM = parallel))
  }
  tsv
}

filter_vep <- function(v, maf_vep, t, missing_population_af = "keep") {
  missing_population_af <- match.arg(missing_population_af, c("keep", "exclude", "error"))
  if (!all(c("population_AF_max", "population_AF_status") %in% names(maf_vep))) {
    maf_vep <- prepare_vep_annotations(maf_vep, missing_population_af)
  }
  require_columns(maf_vep, c("vcf_id", "vcf_qual"), "VEP MAF")
  if (missing_population_af == "error" && any(is.na(maf_vep$population_AF_max))) {
    stop("VEP population frequency unknown; missing_population_af=error", call. = FALSE)
  }
  quality <- suppressWarnings(as.numeric(maf_vep$vcf_qual))
  # Preserve the existing thresholds: AF >=1%, quality cutoff currently disabled.
  excluded <- (!is.na(maf_vep$population_AF_max) & maf_vep$population_AF_max >= 0.01) |
    (!is.na(quality) & quality < 0)
  if (missing_population_af == "exclude") excluded <- excluded | is.na(maf_vep$population_AF_max)
  list(id_include = unique(na.omit(maf_vep$vcf_id)),
       id_exclude = unique(na.omit(maf_vep$vcf_id[excluded])))
}

filter_snpeff <- function(v, tsv_snpeff, t) {
  require_columns(tsv_snpeff, c("ID", "QUAL"), "SnpEff table")
  quality <- suppressWarnings(as.numeric(tsv_snpeff$QUAL))
  list(id_include = unique(na.omit(tsv_snpeff$ID)),
       id_exclude = unique(na.omit(tsv_snpeff$ID[!is.na(quality) & quality < 0])))
}

filter_annotsv <- function(v, tsv_annotsv, terms_relative) {
  include <- with(tsv_annotsv,
    grepl(terms_relative, GenCC_disease, ignore.case = TRUE) |
      grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) |
      grepl("4|5", ACMG_class) | suppressWarnings(as.numeric(AnnotSV_ranking_score)) >= 0.9)
  exclude <- with(tsv_annotsv, switch(v,
    BND = as.numeric(B_gain_AFmax) >= 0.05 & as.numeric(B_loss_AFmax) >= 0.05 & flag_gain & flag_loss,
    DEL = (as.numeric(B_loss_AFmax) >= 0.05 & flag_loss) | (!is.na(po_B_loss_allG_coord) & flag_po_loss),
    INS = (as.numeric(B_ins_AFmax) >= 0.05 & flag_ins) | (!is.na(po_B_gain_allG_coord) & flag_po_gain),
    INV = as.numeric(B_inv_AFmax) >= 0.05 & flag_inv,
    DUP = (as.numeric(B_gain_AFmax) >= 0.05 & flag_gain) | (!is.na(po_B_gain_allG_coord) & flag_po_gain)
  ))
  list(id_include = unique(na.omit(tsv_annotsv$ID[which(include)])),
       id_exclude = unique(na.omit(tsv_annotsv$ID[which(exclude)])))
}

merge_svid <- function(list_vep, list_snpeff, list_annotsv, stringent = TRUE) {
  if (stringent) {
    ids_include <- unique(c(list_vep$id_include, list_snpeff$id_include, list_annotsv$id_include))
    ids_exclude <- unique(c(list_vep$id_exclude, list_snpeff$id_exclude, list_annotsv$id_exclude))
    ids <- setdiff(ids_include, ids_exclude) |> union(list_annotsv$id_include) |> unique()
  } else {
    ids_vep <- setdiff(list_vep$id_include, list_vep$id_exclude)
    ids_snpeff <- setdiff(list_snpeff$id_include, list_snpeff$id_exclude)
    ids <- setdiff(c(ids_vep, ids_snpeff), list_annotsv$id_exclude) |> union(list_annotsv$id_include) |> unique()
  }
  ids
}

filter_svid <- function(caller, sv_type, id_tab = NULL, annotsv_result, vep_result,
                        snpeff_result, terms_relative, missing_population_af = "keep", ...) {
  missing_population_af <- match.arg(missing_population_af, c("keep", "exclude", "error"))
  empty <- list(id_include = NULL, id_exclude = NULL)
  list_vep <- if (!nrow(vep_result)) empty else filter_vep(sv_type, vep_result, caller, missing_population_af)
  list_snpeff <- if (!nrow(snpeff_result)) empty else filter_snpeff(sv_type, snpeff_result, caller)
  list_annotsv <- if (!nrow(annotsv_result)) empty else filter_annotsv(sv_type, annotsv_result, terms_relative)
  ids <- merge_svid(list_vep, list_snpeff, list_annotsv, ...)
  if (missing_population_af %in% c("exclude", "error")) {
    # Explicit strict policy applies to the final selection too, including IDs
    # missing entirely from the MAF. Any measured transcript AF counts as known.
    known <- if (nrow(vep_result)) vep_result$vcf_id[vep_result$population_AF_status == "observed"] else character()
    if (missing_population_af == "error" && length(setdiff(ids, known))) {
      stop("Population frequency unknown for selected IDs absent from the VEP MAF", call. = FALSE)
    }
    if (missing_population_af == "exclude") ids <- intersect(ids, known)
  }
  ids
}
