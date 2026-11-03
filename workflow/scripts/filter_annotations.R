# Use the active R environment's library paths (including its installed packages).
snakemake@source("annotation_io.R")
snakemake@source("utils.R")

file_annotsv <- snakemake@input[["annotsv"]]
file_maf <- snakemake@input[["maf"]]
file_tsv <- snakemake@input[["tsv"]]
file_tab <- snakemake@input[["tab"]]
vep <- snakemake@input[["vep"]]
snpeff <- snakemake@input[["snpeff"]]
vcf <- snakemake@input[["vcf"]]
vep_final <- snakemake@output[["vep"]]
snpeff_final <- snakemake@output[["snpeff"]]
vcf_final <- snakemake@output[["vcf"]]
path_ids <- snakemake@output[["ids"]]
table <- snakemake@output[["table"]]
callers <- sort(snakemake@params[["callers"]])
terms_relative <- snakemake@params[["terms_relative"]]
missing_population_af <- snakemake@params[["missing_population_af"]]
threads <- snakemake@threads
s <- snakemake@wildcards[["sample"]]
t <- snakemake@wildcards[["caller"]]
v <- snakemake@wildcards[["type_sv"]]


annotsv <- file_annotsv |>
  my_vroom(quote = "") |>
  prepare_annotsv_annotations(sv_type = v)
id_info <- file_tab |>
  my_vroom(col_names = c("chr", "pos", "ref", "alt", callers))
if (nrow(annotsv) != 0) {
  tsv_annotsv <- annotsv |>
    dplyr::filter(ID %in% na.omit(id_info[[t]])) |>
    add_overlap_flag_annotsv(sv_type = v, workers = threads)
} else {
  tsv_annotsv <- annotsv
}
maf_vep <- file_maf |>
  my_vroom(na_append = ".") |>
  prepare_vep_annotations(missing_population_af = missing_population_af) |>
  tidyr::drop_na(Start_Position)
tsv_snpeff <- file_tsv |>
  my_vroom() |>
  prepare_snpeff_annotations()

message(glue("[INFO] Applying filters on {s} | {t} | {v} ..."))
ids <- filter_svid(
  caller = t, sv_type = v, annotsv_result = tsv_annotsv, vep_result = maf_vep,
  snpeff_result = tsv_snpeff, terms_relative = terms_relative,
  missing_population_af = missing_population_af, stringent = FALSE
)

list_annotation <- list(
  VEP = maf_vep, SnpEff = tsv_snpeff, AnnotSV = tsv_annotsv
)

annotation_modified <- lapply(seq_along(list_annotation), function(i) {
  x <- list_annotation[[i]]
  if ("vcf_id" %in% colnames(x)) {
    x <- dplyr::mutate(x, ID = vcf_id)
  }
  if (length(ids) != 0 && nrow(x) != 0) {
    tmp <- x[x$ID %in% ids, ]
    colnames(tmp) <- glue::glue("{names(list_annotation)[i]}.{colnames(tmp)}")
    tmp
  } else {
    NULL
  }
}) |>
  dplyr::bind_rows()

write.table(
  annotation_modified, table, quote = FALSE, col.names = TRUE,
  row.names = FALSE, sep = "\t", na = ""
)

write.table(ids, path_ids, quote = FALSE, col.names = FALSE, row.names = FALSE)
filter_vcf_checked(vep, vep_final, path_ids)
filter_vcf_checked(snpeff, snpeff_final, path_ids)
filter_vcf_checked(vcf, vcf_final, path_ids)

save.image(snakemake@output[["rdata"]])
