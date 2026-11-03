# AnnotSV / cuteSV insertion-length mismatch

## Observed failure

For `cutesv/h1/cutesv.INS.vcf`, cuteSV records contain `INFO/SVLEN`, but many rows in `cutesv/h1/cutesv.INS.annotsv.tsv` have an empty `SV_length` (`NA` when read by R). In one run, **9,292 of 13,857 AnnotSV rows** had an unavailable `SV_length` (counts include repeated annotation rows). The downstream `filter_annotations` rule used `SV_start + SV_length - 1` for insertions, causing `IRanges` to fail with `'start' or 'end' cannot contain NAs`.

Example from that run:

| VCF ID | VCF POS | INFO/SVLEN | INFO/CIPOS | AnnotSV SV_start | AnnotSV SV_end | AnnotSV SV_length |
| --- | ---: | ---: | --- | ---: | ---: | --- |
| `cuteSV.INS.0` | 90320 | 115 | `-34,34` | 90286 | 90320 | blank |
| `cuteSV.INS.1` | 136954 | 276 | `-19,19` | 136935 | 136954 | blank |

`CIPOS` explains the start shift: 90320 − 34 = 90286. It represents breakpoint uncertainty, **not** insertion length. This is not evidence of a GRCh19/GRCh38 mismatch: `config/config.yaml` selects GRCh38, and the [AnnotSV rule](../workflow/rules/annotator/annotsv.smk) passes that value via `-genomeBuild`.

## Cause and scope

In the [AnnotSV v3.5.10 VCF reader](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-vcf.tcl), `SVLEN` is read from INFO and stored under an internal variant ID built from the original VCF `posVCF`. The [v3.5.10 writer](https://github.com/lgmgeo/AnnotSV/blob/v3.5.10/share/tcl/AnnotSV/AnnotSV-write.tcl) reconstructs that ID using the `CIPOS`-adjusted `SVstart_VCF` to retrieve `SV_length`. When those starts differ, the lookup misses; for INS, the writer cannot calculate length from the span, so it emits an empty value. This is an inference from the versioned source and the observed VCF/TSV values, not an upstream-confirmed bug report. [AnnotSV issue #37](https://github.com/lgmgeo/AnnotSV/issues/37) confirms the intended `includeCI` coordinate expansion but does not report this specific lost-length issue.

The workflow [pins AnnotSV 3.5.10](../workflow/envs/annotsv.yaml); **upgrading to 3.5.10 alone does not repair the mismatch**. `annotation.annotsv.version: v3.5` in config refers to annotation resources, not the genome build or the installed AnnotSV patch version.

A separate problem surfaced while restoring lengths: VCF `ID` is not necessarily unique (for example `chr10_3314715_N_<INS>` appeared with different `SVLEN` values). Never infer the length of an ambiguous ID from the first matching VCF record.

## Workflow mitigation

The [restore rule](../workflow/rules/filtration/restore_annotsv_length.smk) writes `{caller}.{type_sv}.annotsv.length.tsv` between the original AnnotSV output and `filter_annotations`. Its [script](../workflow/scripts/restore_annotsv_length.py) fills only missing INS `SV_length` cells, primarily using that **row's own INFO/SVLEN**; it falls back to a VCF lookup only for an unambiguous ID and fails instead of guessing. It leaves the original AnnotSV TSV and CI-adjusted `SV_start`/`SV_end` unchanged. This is a workflow-level repair, not a fix in AnnotSV. Do not disable `-includeCI` merely to hide the issue: doing so changes genomic overlap annotations and potentially filtering decisions.

Validation: `python3 -m unittest discover -s tests -p 'test_restore_annotsv_length.py'`. On real data, compare the raw and restored AnnotSV TSVs: only previously missing INS `SV_length` cells should differ. If a row lacks both a usable `INFO/SVLEN` and an unambiguous VCF `SVLEN`, the restore rule intentionally stops and requires investigation. An empty AnnotSV TSV (observed for `svision/*.BND.annotsv.tsv`) is passed through as an empty file; there are no rows to restore.
