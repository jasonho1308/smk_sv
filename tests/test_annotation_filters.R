args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
source(file.path(root, "workflow/scripts/annotation_io.R"))
source(file.path(root, "workflow/scripts/utils.R"))

expect_error <- function(expression, pattern) {
  message <- tryCatch({force(expression); NULL}, error = function(e) conditionMessage(e))
  stopifnot(!is.null(message), grepl(pattern, message))
}

work <- tempfile("filter-tests-")
dir.create(work)

# AnnotSV widths are header-driven, extra fields and apostrophes are preserved.
wide <- as.data.frame(setNames(rep(list("value"), 131), paste0("field", seq_len(131))))
wide[[1]] <- "O'Brien"
path <- file.path(work, "wide.tsv")
write.table(wide, path, sep = "\t", quote = FALSE, row.names = FALSE)
parsed <- my_vroom(path, quote = "")
stopifnot(ncol(parsed) == 131, parsed[[1]][[1]] == "O'Brien")
# Use the workflow's AnnotSV reader: # and quotes are literal annotation text,
# not comment/quote delimiters. Later columns and subsequent rows must survive.
annotsv_path <- file.path(work, "annotsv.tsv")
writeLines(c("ID\tSV_chrom\tSV_start\tSV_end\tOMIM_phenotype\tACMG_class",
             "del1\tchr1\t10\t20\tO'Brien disease #123456\t3",
             "del2\tchr2\t30\t40\t\"Disease\" #654321\t4"), annotsv_path)
for (sv_type in c("DEL", "INV")) {
  parsed_annotsv <- suppressWarnings(read_annotsv_annotations(annotsv_path, sv_type))
  stopifnot(nrow(parsed_annotsv) == 2L,
            identical(parsed_annotsv$ID, c("del1", "del2")),
            identical(parsed_annotsv$OMIM_phenotype, c("O'Brien disease #123456", "\"Disease\" #654321")),
            identical(parsed_annotsv$ACMG_class, c("3", "4")))
}
# INV benign-overlap evidence comes after free-text annotations. A literal hash
# must not truncate that evidence or silently change the filtering decision.
writeLines(c("ID\tSV_chrom\tSV_start\tSV_end\tOMIM_phenotype\tB_inv_coord\tB_inv_AFmax",
             "inv1\tchr1\t10\t20\tDisease #123456\tchr1:10-20\t0.1"), annotsv_path)
parsed_inv <- suppressWarnings(read_annotsv_annotations(annotsv_path, "INV"))
parsed_inv <- add_overlap_flag_annotsv(parsed_inv, "INV", workers = 1)
stopifnot(identical(parsed_inv$OMIM_phenotype, "Disease #123456"),
          identical(parsed_inv$B_inv_AFmax, "0.1"), parsed_inv$flag_inv,
          identical(filter_annotsv("INV", parsed_inv, "blood")$id_exclude, "inv1"))
writeLines(c("a\tb", "one\ttwo\tthree"), file.path(work, "bad.tsv"))
expect_error(suppressWarnings(my_vroom(file.path(work, "bad.tsv"))), "Malformed annotation")
expect_error(suppressWarnings(my_vroom(file.path(work, "bad.tsv"))),
             "row 2, column 3: expected 2 columns; actual 3 columns")
writeLines(c("ID\tSV_chrom\tSV_start\tSV_end", "del1\tchr1\t10\t20\textra"), annotsv_path)
for (sv_type in c("DEL", "INV")) {
  expect_error(suppressWarnings(read_annotsv_annotations(annotsv_path, sv_type)), "Malformed annotation")
}
# Other readers still accept metadata comments (e.g. VEP MAF headers).
commented_path <- file.path(work, "commented.tsv")
writeLines(c("# metadata", "a\tb", "one\ttwo"), commented_path)
stopifnot(identical(my_vroom(commented_path)$a, "one"))
expect_error(prepare_annotsv_annotations(data.frame(foo = "bar"), "DEL"), "required columns")

# Absent frequencies are unknown, AF=0 is measured, and comma-separated AF works.
maf <- data.frame(vcf_id = c("unknown", "common", "rare", "zero", "multi"),
                  Start_Position = as.character(1:5), vcf_qual = "30",
                  AF = c(NA, "0.02", "0.001", "0", "0.001,0.2"), stringsAsFactors = FALSE)
prepared <- suppressWarnings(prepare_vep_annotations(maf))
stopifnot(identical(prepared$population_AF_status, c("unknown", rep("observed", 4))))
stopifnot(is.na(prepared$population_AF_max[[1]]), prepared$population_AF_max[[4]] == 0)
stopifnot(setequal(filter_vep("DEL", prepared, "cutesv")$id_exclude, c("common", "multi")))
expect_error(prepare_vep_annotations(maf, "error"), "frequency unknown")
invalid_af <- maf[2, ]
invalid_af$AF <- "1.2"
expect_error(prepare_vep_annotations(invalid_af), "outside")
no_af <- maf[, c("vcf_id", "Start_Position", "vcf_qual")]
stopifnot(all(suppressWarnings(prepare_vep_annotations(no_af))$population_AF_status == "unknown"))
expect_error(prepare_vep_annotations(data.frame(bad = "schema")), "required columns")

# Missing optional AnnotSV evidence must not invent benign overlap.
annotsv <- suppressWarnings(prepare_annotsv_annotations(
  data.frame(ID = "unknown", SV_chrom = "chr1", SV_start = "1", SV_end = "100"), "DEL"))
annotsv <- add_overlap_flag_annotsv(annotsv, "DEL", workers = 1)
stopifnot(!annotsv$flag_loss, !annotsv$flag_po_loss)
annotsv$B_loss_coord <- "chr1:1-100"
annotsv$B_loss_AFmax <- "0.1"
annotsv <- add_overlap_flag_annotsv(annotsv, "DEL", workers = 1)
stopifnot(annotsv$flag_loss, identical(filter_annotsv("DEL", annotsv, "blood")$id_exclude, "unknown"))
ins <- suppressWarnings(prepare_annotsv_annotations(data.frame(ID = "x", SV_chrom = "1", SV_start = "2"), "INS"))
stopifnot(!suppressWarnings(add_overlap_flag_annotsv(ins, "INS", workers = 1))$flag_ins)

# Strict missing-AF policy cannot be bypassed by SnpEff's evidence branch.
empty_annotsv <- prepare_annotsv_annotations(data.frame(), "DEL")
snpeff <- prepare_snpeff_annotations(data.frame(ID = maf$vcf_id, QUAL = "30"))
ids <- filter_svid("cutesv", "DEL", annotsv_result = empty_annotsv,
                  vep_result = prepared, snpeff_result = snpeff, terms_relative = "blood",
                  missing_population_af = "exclude", stringent = FALSE)
stopifnot(!"unknown" %in% ids, "zero" %in% ids)
no_vep <- prepare_vep_annotations(data.frame())
ids <- filter_svid("cutesv", "DEL", annotsv_result = empty_annotsv,
                  vep_result = no_vep, snpeff_result = snpeff, terms_relative = "blood",
                  missing_population_af = "exclude", stringent = FALSE)
stopifnot(length(ids) == 0)
expect_error(filter_svid("cutesv", "DEL", annotsv_result = empty_annotsv,
                        vep_result = no_vep, snpeff_result = snpeff, terms_relative = "blood",
                        missing_population_af = "error", stringent = FALSE), "frequency unknown")

# Command failures are fatal and diagnostics retain the exit status.
bad <- file.path(work, "bad command")
writeLines(c("#!/bin/sh", "echo deliberate-failure >&2", "exit 7"), bad)
Sys.chmod(bad, "0755")
expect_error(run_checked(bad, "argument with spaces"), "exit 7.*|deliberate-failure")
expect_error(run_checked("nonexistent-smk-sv-test-executable"), "not found")

# Actual bcftools checks cover selected IDs, empty IDs, and paths containing spaces.
if (nzchar(Sys.which("bcftools"))) {
  input <- file.path(work, "input variants.vcf")
  output <- file.path(work, "output variants.vcf")
  ids_file <- file.path(work, "selected ids.txt")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=100>",
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
               "chr1\t1\tid1\tA\tG\t30\tPASS\t.",
               "chr1\t2\tid2\tA\tT\t30\tPASS\t."), input)
  writeLines("id1", ids_file)
  filter_vcf_checked(input, output, ids_file)
  records <- readLines(output)
  stopifnot(any(grepl("id1", records)), !any(grepl("id2", records)))
  file.create(ids_file)
  filter_vcf_checked(input, output, ids_file)
  stopifnot(all(grepl("^#", readLines(output))))
  original <- readLines(output)
  tools <- file.path(work, "tools")
  dir.create(tools)
  file.copy(bad, file.path(tools, "bcftools"))
  Sys.chmod(file.path(tools, "bcftools"), "0755")
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste(tools, old_path, sep = .Platform$path.sep))
  expect_error(filter_vcf_checked(input, output, ids_file), "failed")
  stopifnot(identical(readLines(output), original), !length(list.files(work, pattern = "^\\.filtered-")))
  Sys.setenv(PATH = old_path)
}
unlink(work, recursive = TRUE)
cat("Annotation schema, missing-frequency policies, and checked-command tests passed\n")
