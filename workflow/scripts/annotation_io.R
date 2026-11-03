# Schema normalization and checked process execution for annotation filtering.
# These helpers use base R so error paths can be tested without annotation tools.

require_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(sprintf("%s is missing required columns: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(x)
}

empty_columns <- function(columns) {
  as.data.frame(setNames(rep(list(character()), length(columns)), columns), stringsAsFactors = FALSE)
}

add_optional_columns <- function(x, columns, label) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    if (nrow(x)) {
      warning(sprintf("%s: unavailable optional columns filled with NA: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
    }
    for (column in missing) x[[column]] <- rep(NA_character_, nrow(x))
  }
  x
}

prepare_vep_annotations <- function(x, missing_population_af = "keep") {
  policy <- match.arg(missing_population_af, c("keep", "exclude", "error"))
  required <- c("vcf_id", "Start_Position")
  if (!nrow(x) && !ncol(x)) x <- empty_columns(required)
  require_columns(x, required, "VEP MAF")
  frequency_columns <- c("AF", "EAS_AF", "ASN_AF", "gnomADe_AF", "gnomADe_EAS_AF",
                         "gnomADg_AF", "gnomADg_EAS_AF", "gnomAD_AF", "gnomAD_EAS_AF")
  present <- intersect(frequency_columns, names(x))
  maxima <- rep(NA_real_, nrow(x))
  if (length(present) && nrow(x)) {
    maxima <- vapply(seq_len(nrow(x)), function(i) {
      tokens <- unlist(strsplit(as.character(unlist(x[i, present, drop = FALSE])), "[,;|&]"))
      values <- suppressWarnings(as.numeric(tokens))
      if (any(is.finite(values) & (values < 0 | values > 1))) {
        stop("VEP population frequency outside [0,1]", call. = FALSE)
      }
      valid <- is.finite(values) & values >= 0 & values <= 1
      if (any(valid)) max(values[valid]) else NA_real_
    }, numeric(1))
  }
  x <- add_optional_columns(x, "vcf_qual", "VEP MAF")
  x$population_AF_max <- maxima
  x$population_AF_status <- ifelse(is.na(maxima), "unknown", "observed")
  unknown <- sum(is.na(maxima))
  if (unknown) {
    text <- sprintf("VEP population frequency unknown for %d/%d rows; missing_population_af=%s. Unknown is not AF=0.", unknown, nrow(x), policy)
    if (policy == "error") stop(text, call. = FALSE)
    warning(text, call. = FALSE)
  }
  x
}

prepare_snpeff_annotations <- function(x) {
  if (!nrow(x) && !ncol(x)) x <- empty_columns("ID")
  require_columns(x, "ID", "SnpEff table")
  add_optional_columns(x, "QUAL", "SnpEff table")
}

prepare_annotsv_annotations <- function(x, sv_type) {
  required <- c("ID", "SV_chrom", "SV_start")
  if (sv_type != "INS") required <- c(required, "SV_end")
  if (!nrow(x) && !ncol(x)) x <- empty_columns(required)
  require_columns(x, required, "AnnotSV table")
  optional <- c("SV_length", "GenCC_disease", "OMIM_phenotype", "ACMG_class", "AnnotSV_ranking_score",
                "B_gain_coord", "B_loss_coord", "po_B_loss_allG_coord", "B_ins_coord",
                "po_B_gain_allG_coord", "B_inv_coord", "B_gain_AFmax", "B_loss_AFmax",
                "B_ins_AFmax", "B_inv_AFmax")
  add_optional_columns(x, optional, "AnnotSV table")
}

run_checked <- function(command, args = character(), stdin = "") {
  executable <- Sys.which(command)
  if (!nzchar(executable)) stop(sprintf("Required executable not found: %s", command), call. = FALSE)
  output <- suppressWarnings(system2(executable, args = vapply(args, shQuote, character(1)),
                                     stdin = stdin, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("%s failed (exit %s):\n%s", command, status, paste(output, collapse = "\n")), call. = FALSE)
  }
  invisible(output)
}

filter_vcf_checked <- function(input, output, ids_file) {
  if (!file.exists(input)) stop(sprintf("Input VCF does not exist: %s", input), call. = FALSE)
  if (!file.exists(ids_file)) stop(sprintf("Variant ID list does not exist: %s", ids_file), call. = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  staged <- tempfile(pattern = ".filtered-", tmpdir = dirname(output), fileext = ".vcf")
  on.exit(unlink(staged), add = TRUE)
  # An empty selection must produce a valid header-only VCF, not all records.
  if (file.info(ids_file)$size == 0) {
    args <- c("view", "--no-version", "-h", "-Ov", "-o", staged, input)
  } else {
    # BCFtools' @filename expression parser cannot handle spaces in filenames.
    # Feed the ID list on stdin instead of embedding a path in its expression.
    args <- c("filter", "--no-version", "-i", "ID=@-", "-Ov", "-o", staged, input)
  }
  run_checked("bcftools", args, stdin = if (file.info(ids_file)$size) ids_file else "")
  if (!file.exists(staged) || !file.info(staged)$size) {
    stop("bcftools did not produce a nonempty VCF", call. = FALSE)
  }
  if (!file.rename(staged, output)) stop(sprintf("Could not install filtered VCF: %s", output), call. = FALSE)
  invisible(output)
}
