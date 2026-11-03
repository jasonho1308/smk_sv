# Workflow package update review

## Applied update — 2026-11-03

The environment pins and `pixi.lock` have now been upgraded against the live
Bioconda/conda-forge and PyPI registries. All 21 workflow environments resolve for
`linux-64`, including NanoVar's pip dependencies (resolution only, not execution).

| Component | New version |
| --- | --- |
| SAMtools / BCFtools / HTSlib | 1.25 |
| minimap2 | 2.31 |
| Sniffles | 2.8.1 |
| Clair3 | 2.0.3 (latest Bioconda); pinned container 2.0.5 |
| DELLY / Severus / NanoVar | 2.7.0 / 1.7 / 1.8.3 |
| VEP / cache release | 116.3 / 116 |
| SnpEff / SnpSift | 5.4.0c (requires Java >=21) |
| awscli / curl / gdown | 2.38.3 / 8.23.0 / 6.5.1 |
| perl-compress-raw-zlib | 2.214 |

Already-current caller pins and AnnotSV 3.5.10 remain unchanged. SVision only
publishes a mutable `latest` image; it is now pinned to its current digest.
Python 3.11 is retained where required by Clair3 and NanoVar/TensorFlow.
Unpinned support/R dependencies resolve to the latest compatible versions.
The Pixi Snakemake/PBS forks are preserved and pinned to their current commits,
rather than silently replaced with upstream packages.

Defaults now use the T2T-CHM13v2.0 analysis set (including chrY), its matching
Severus repeat BED, and a checksummed PyTorch model archive. SnpEff now builds its
custom CHM13 database from RefSeq/Liftoff v5.3 automatically; VEP now uses local
CHM13 GFF/GTF models without an Ensembl cache: see [README](../README.md#reference-resources).
Do not reuse GRCh38 alignments/results, and benchmark caller changes before
interpreting biological differences.

Validation: all 21 Linux environments solved, `pixi lock --check --no-install`
passed, and all 9 tests in `tests/test_t2t_configuration.py` passed using the
pinned Snakemake fork and deployment plugins. Tests cover configuration/schema
consistency, PyTorch archive layout and legacy-model rejection, complete
analysis/resource DAG dry-runs, missing-cache guidance, and existing-cache reuse.
No whole-genome analysis, Linux tool execution, or Docker image rebuild was run.

Subsequent SnpEff build integration: **23 tests pass**, including real SnpEff
5.4c toy-genome builds, CDS/protein verification, rejection of mismatched CDS,
and build-to-annotation execution through Snakemake. The full CHM13 database
has not been built here. The implementation follows the
[official database-building guide](https://pcingola.github.io/SnpEff/snpeff/build_db/).

Offline VEP/filter integration: **35 tests pass** and all **22 Linux environment
specifications resolve**. Tests include HTSlib sorting/indexing, a real VEP 116.3
GFF and GTF annotation in a network-disabled Linux container, conversion with vcf2maf 1.6.22,
and real R/BCFtools filtering under keep/exclude/error missing-frequency policies.
Annotated TSV widths are header-driven and command failures are fatal. A temporary
native macOS VEP install failed at startup; the real VEP execution test therefore
uses Linux, the production target. No whole-genome analysis was run.

Version sources: [Bioconda registry](https://api.anaconda.org/package/bioconda/samtools),
[conda-forge registry](https://api.anaconda.org/package/conda-forge/awscli),
[NanoVar PyPI](https://pypi.org/project/nanovar/),
[Clair3 README](https://github.com/HKU-BAL/Clair3#latest-updates),
[Clair3 model archive](https://github.com/HKU-BAL/Clair3/releases/tag/v2.0.2),
[T2T analysis sets](https://github.com/marbl/CHM13#genome-assembly).

## Historical review — 2026-09-24

Reviewed 2026-09-24 against the `linux-64` packages available from Bioconda/conda-forge and the upstream release notes linked below. These are recommendations, **not** changes to environment files. Update and validate packages in small groups; newer variant callers can change biological results.

## Recommended first: HTS utilities

| Package | Current workflow pins | Available version | Assessment |
| --- | --- | --- | --- |
| bcftools | 1.22 (`bcftools.yaml`, `duphold.yaml`) | 1.24 | **Update.** [1.23.1](https://github.com/samtools/bcftools/releases/tag/1.23.1) fixed silent truncation when input reads fail in commands including `filter` and `query`; [1.24](https://github.com/samtools/bcftools/releases/tag/1.24) fixes additional VCF handling. Test sorting, filtering, reheadering and querying. |
| samtools | 1.22 (`samtools.yaml`, `vcf2maf.yaml`), 1.20 (`nanovar.yaml`), 1.24 (`minimap2.yaml`); container Dockerfile builds 1.20 | 1.24 | **Align the older pins with 1.24**, subject to env solves. [1.23](https://github.com/samtools/samtools/releases/tag/1.23) and [1.24](https://github.com/samtools/samtools/releases/tag/1.24) include bug fixes. The most notable HTSlib 1.23 long-alignment region-iterator fix concerns CRAM; this workflow normally uses BAM, so it is *not* an explanation for the duphold crash. Check whether the container Dockerfile is used before changing its independent build. |
| htslib | 1.22 explicitly in `vcf2maf.yaml`; otherwise resolved as a dependency | 1.24 | **Update the explicit pin alongside the tools.** [1.23](https://github.com/samtools/htslib/releases/tag/1.23) fixed CRAM region iteration and other issues; [1.24](https://github.com/samtools/htslib/releases/tag/1.24) adds BCF validation and security fixes. Avoid changing the underlying library without testing consumers. |

Run a small smoke test of BAM sort/index and the workflow's `bcftools sort`, `filter`, `query`, and `reheader` steps; compare record counts and verify outputs can be read to EOF. Update all relevant environment definitions deliberately; `pixi.toml` is the workflow driver, not a replacement for Snakemake's per-rule Conda environments.

## Test separately: updates that can change calls

| Package | Current pin → available | Assessment |
| --- | --- | --- |
| minimap2 | 2.30 → 2.31 | **Worth testing on a new analysis run.** [2.31](https://github.com/lh3/minimap2/releases/tag/v2.31) fixes occasionally wrong supplementary/secondary flags and rare inversion-alignment out-of-bounds access. Upstream says alignments may change; the mapper's BAM feeds every caller. Validate the reference index and compare alignment/call metrics before adoption. |
| Sniffles | 2.2 → 2.8.1 | **Worth benchmarking.** [2.8](https://github.com/fritzsedlazeck/Sniffles/releases/tag/v2.8.0) reworks BND calling and improves DUP genotyping. This workflow explicitly passes `--minsvlen 100`; [2.8.1](https://github.com/fritzsedlazeck/Sniffles/releases/tag/v2.8.1) treats explicit values as hard limits, while `~100` preserves the former leeway. Review this option and downstream filters when comparing callsets. |
| cuteSV | 2.1.1 committed → 2.1.4 | **Test on `h8` first.** The working tree currently has an uncommitted change to `workflow/envs/cutesv.yaml` specifying 2.1.4. The [latest release](https://github.com/tjiangHIT/cuteSV/releases/tag/cuteSV-v2.1.4) does not document a fix for the observed POS=0 *insertion*. The [POS=0 issue fixed upstream](https://github.com/tjiangHIT/cuteSV/issues/178) concerns a BND coordinate instead. Count POS=0 records in the new VCF and compare calls. |
| Clair3 | 1.1.2 → 2.0.3 | **Treat as a migration.** [v2](https://github.com/HKU-BAL/Clair3#v200--feb-9-2026-major-release) switched from TensorFlow to PyTorch; v1 models are incompatible. Verify the configured model at `config/config.yaml`, select a compatible model using the [official model catalogs](https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/) and [converted Rerio catalog](https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/), and pin the mutable `docker://hkubal/clair3:latest` in `workflow/rules/caller/clair3.smk` before changing versions. The currently selected caller list in `config/config.yaml` does not include Clair3. |
| Severus | 1.6 → 1.7 | Optional single-caller benchmark; no workflow-specific bug fix verified. |
| DELLY | 1.2.6 → 2.6.0 | Major caller change; validate `delly lr`, the output VCF and downstream filters before adopting. [Release](https://github.com/dellytools/delly/releases/tag/v2.6.0). |
| NanoVar (PyPI) | 1.7.0 → 1.8.3 | Test independently, including its fixed Python and tool dependencies. It is not in the default caller list. |

## Defer unless addressing a specific need

- **VEP:** 114.2 → 116.2. `config/config.yaml` specifies VEP cache version **114**, and `awscli.yaml` also pins VEP 114.2. Update executable, cache paths/version and annotation tests together; do not bump the executable alone.
- **SnpEff and SnpSift:** 5.2 → 5.4.0c. Update together and verify database compatibility and extracted annotation fields. The annotation config specifies its own database version (`87`).
- **AnnotSV:** 3.5 → 3.5.10. Validate tabular output and the workflow's R-side annotation parsing.
- **awscli / curl / gdown:** 2.27.53 → 2.36.47 / 8.14.1 → 8.22.0 / 5.2.0 → 6.4.0. These can be tested in the resource-download steps independently; they are unrelated to the SV crash.
- **Already current in checked channels:** SVIM 2.0.0, deBreak 1.3, NanoSV 1.2.4, SURVIVOR 1.0.7, vcf2maf 1.6.22 and bedtools 2.31.1. The separately downloaded [duphold v0.2.3](https://github.com/brentp/duphold/releases/tag/v0.2.3) is also the latest GitHub release.

## Specific `h8` failure

The observed cuteSV VCF record has `POS=0`, `END=0`, `SVTYPE=INS` and a 9,225-base ALT. It is invalid VCF and precedes duphold's crash; a package upgrade is **not** proven to repair it. At the time of this review, `workflow/rules/filtration/duphold.smk` has no POS=0 guard. Validate or filter malformed records before `bcftools sort`/duphold, or demonstrate on the same BAM and reference that the upgraded caller no longer emits them. Keep raw caller output for comparison and record the number of dropped records.
