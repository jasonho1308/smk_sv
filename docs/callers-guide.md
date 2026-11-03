# Variant callers: a simple guide

A **variant caller** finds differences between a sample's DNA and a reference
genome. This pipeline mainly looks for **structural variants (SVs)**: relatively
large pieces of DNA that are missing, added, copied, reversed, or rearranged.

The links below are the tools' original project repositories. Descriptions were
checked against their official READMEs; enabled/disabled status comes from this
repository's `config/config.yaml`.

### Reading the T2T column

Here, **T2T** means the human **T2T-CHM13v2.0** reference, not every T2T assembly.

- **Documented:** upstream explicitly provides a CHM13 example or resource.
- **Expected:** the tool accepts a supplied reference or reference-aligned BAM,
  but its checked README does not explicitly establish T2T support. This is an
  inference, **not** a validated T2T benchmark or a guarantee.
- **Conditional:** extra resources or workflow changes are needed.

None of these callers has been tested on T2T in this working session. Tool-level
compatibility does **not** mean this entire pipeline is ready to run unchanged;
see [T2T requirements below](#5-t2t-requirements-for-this-pipeline).

## 1. SV callers enabled in this pipeline

These five callers examine the same sample using different approaches. Their
results are filtered and combined rather than treating one caller as definitive.

| Caller | What it does, in simple terms | T2T-CHM13 reference? | Original URL |
| --- | --- | --- | --- |
| **cuteSV** | Finds SV evidence in long-read alignments, groups similar evidence, and refines it into variant calls. A general-purpose long-read SV caller. | **Expected.** Supply the T2T FASTA and BAM aligned to it. | <https://github.com/tjiangHIT/cuteSV> |
| **Sniffles2** (`sniffles`) | Finds and genotypes SVs from long-read alignments. Also offers population and mosaic calling modes, although this pipeline does not enable those modes. | **Expected.** Use T2T-aligned reads and a matching repeat BED when supplying `--tandem-repeats`; do not reuse GRCh38 intervals. | <https://github.com/fritzsedlazeck/Sniffles> |
| **SVIM** | Finds and classifies several SV types, including distinguishing copies placed next to the original DNA from copies inserted elsewhere. | **Expected.** Its alignment interface accepts a supplied genome FASTA; use the matching T2T-aligned BAM. | <https://github.com/eldariont/svim> |
| **SVision** | Uses deep learning to recognize simple and complex rearrangements from image-like representations of aligned reads. Here, it runs separately on 24 chromosomes and then gathers the results. | **Expected; workflow conditional.** Supply T2T FASTA + index and matching BAM; chromosome selection must match the reference. T2T-specific model accuracy is not established here. | <https://github.com/xjtu-omics/SVision> |
| **Severus** | Uses long-read phasing and a graph of DNA breakpoints to find SVs, including complex cancer-related rearrangements. Here, it receives phased data from Clair3. | **Documented.** Upstream supplies CHM13 VNTR annotations and discusses CHM13 panels of normals. Use T2T-aligned/phased inputs and matching tracks. | <https://github.com/KolmogorovLab/Severus> |

**Important for Severus:** the original tool supports tumor–normal analysis, but
this pipeline supplies only the tumor/sample BAM, with no `--control-bam`. It uses
Severus's `all_SVs/severus_all.vcf`, not a matched-normal somatic call set. The
upstream README notes that without a control BAM, only germline-mode calls are
made; this does **not** mean every variant in a tumor sample is biologically
inherited. These calls are not confirmed somatic variants.

## 2. Clair3: a helper caller, not an SV caller

| Tool | What it does here | T2T-CHM13 reference? | Original URL |
| --- | --- | --- | --- |
| **Clair3** | Calls small DNA changes—single-base changes and small insertions/deletions. This workflow also uses WhatsHap through Clair3 to phase variants and tag reads by chromosome copy, providing inputs for Severus. | **Expected.** Accepts an indexed FASTA and matching BAM. Check contig selection and a sequencing-model/software-version match; changing reference alone does not establish T2T accuracy. | <https://github.com/HKU-BAL/Clair3> |

**Phasing** groups variants belonging to the same chromosome copy.
**Haplotagging** labels reads with that chromosome-copy assignment.
Clair3's small-variant calls are **not** included as a sixth SV caller in the
SURVIVOR consensus.

## 3. Additional SV caller modules—not currently enabled

These modules exist under `workflow/rules/caller/`, but are not selected in the
current configuration. Their presence does not mean their integration has been
validated; see the [detailed DAG notes](workflow-dag.md#4-other-available-caller-branches)
before enabling them.

| Caller | What it does, in simple terms | T2T-CHM13 reference? | Original URL |
| --- | --- | --- | --- |
| **DeBreak** | Finds SVs from long reads and can refine breakpoint positions using read consensus; it can also rescue very large insertions using local assembly. | **Expected.** Accepts a supplied reference, which must match the input BAM. | <https://github.com/Maggi-Chen/DeBreak> |
| **Delly** | Finds SVs using short-read or long-read data. This repository's rule uses its long-read mode, `delly lr`, for ONT reads. | **Documented CHM13 use.** Upstream shows `chm13.fa` in its alternate-alignment workflow; this does not validate this repository's inactive Delly module. | <https://github.com/dellytools/delly> |
| **NanoSV** | Finds rearrangement junctions by grouping reads whose aligned pieces have similar positions and orientations. | **Conditional / unverified.** Use T2T-aligned reads and a T2T sampling BED if depth-support filtering is enabled; the bundled default BED is hg19. Its local formatter also needs repair. | <https://github.com/mroosmalen/nanosv> |
| **NanoVar** | Finds SVs from long-read data, with a focus on low-depth sequencing. | **Expected at tool level; not supported unchanged by this workflow.** The tool accepts a supplied FASTA and optional custom filter BED, but this rule only maps GRCh37/GRCh38 to built-in filters. It needs changes for `CHM13`. | <https://github.com/cytham/nanovar> |

## 4. How the callers work together

```text
ONT reads → Minimap2 alignment → cuteSV / Sniffles2 / SVIM / SVision
                             ↘ Clair3 phasing → Severus

Each SV caller's VCF
  → depth annotation and hard filtering
  → split into SV types
  → first SURVIVOR consensus + caller-specific annotation
  → annotation filtering
  → final SURVIVOR consensus VCFs
```

Using multiple callers provides complementary evidence, but agreement does not
guarantee a true variant: callers share the same reads and can share errors.
The final output is a set of **candidate SVs**, not proof that they are
cancer-specific. Tumor-only analysis cannot reliably separate all inherited
variants from somatic variants without additional evidence.

### SV abbreviations in the output

- **DEL:** deletion—a piece of DNA is missing.
- **INS:** insertion—extra DNA is added.
- **DUP:** duplication—a piece of DNA is copied.
- **INV:** inversion—a piece of DNA is reversed.
- **BND:** breakend—a record describing one side of a rearrangement junction,
  often used for translocations or complex rearrangements.

## 5. T2T requirements for this pipeline

A T2T configuration exists at [`config/config_T2T.yaml`](../config/config_T2T.yaml),
but it is **not evidence of a successful end-to-end T2T run**.

1. **Use one exact assembly throughout.** Align reads to the selected T2T FASTA
   and rebuild its Minimap2 and FASTA indexes. All BAMs, VCFs, BED tracks, and
   annotations must use matching coordinates and chromosome names. GRCh38 files
   cannot be reused just by renaming them.
2. **Use T2T repeat tracks.** The T2T configuration points Sniffles and Severus to
   Severus's CHM13 BED. This is an assembly-matched resource, not proof that the
   same track is optimal for both callers.
3. **Check chromosome coverage.** `workflow/rules/constants.smk` fixes SVision's
   list to `chr1`–`chr22`, `chrX`, and `chrY`; other contigs are omitted. The
   linked T2T analysis-set `chm13v2.0.fa.gz` includes HG002-derived Y, whereas
   `chm13v2.0_noY.fa.gz` does not. Adjust the list for the actual FASTA; switching
   to a no-Y or accession-named assembly without changes is unsafe.
4. **Prepare T2T annotation databases.** `snakemake build_snpeff_cache` builds
   SnpEff `CHM13.v2.0` from the matching FASTA and RefSeq/Liftoff v5.3 annotation;
   `snakemake resources` also includes this build. Review its log and configure
   matching CDS/protein FASTAs to enable sequence cross-checks. VEP now uses
   local, sorted/indexed CHM13 GFF/GTF models and the reference FASTA, without an
   Ensembl cache. Missing population frequencies are explicitly reported as
   unknown, not zero. AnnotSV
   explicitly accepts `-genomeBuild CHM13`, but still needs the corresponding
   annotation data. Do not substitute GRCh38 caches.
5. **Keep Clair3's software and model compatible.** The T2T config names an older
   R9.4.1 PyTorch model; the rule pins container v2.0.5 (Bioconda currently
   provides v2.0.3). Verify the installed
   release and model format as well as chemistry/basecaller compatibility;
   T2T reference compatibility does not resolve a model-version mismatch.

### T2T evidence and further reading

- [Official T2T reference assemblies and analysis sets](https://github.com/marbl/CHM13#genome-assembly).
- [Severus CHM13 repeat-track directory](https://github.com/KolmogorovLab/Severus/tree/main/vntrs)
  and its [README](https://github.com/KolmogorovLab/Severus).
- [Delly README](https://github.com/dellytools/delly): alternate-alignment example
  using `chm13.fa`.
- [SVIM reference argument](https://github.com/eldariont/svim/blob/master/src/svim/SVIM_input_parsing.py).
- [AnnotSV command-line options](https://github.com/lgmgeo/AnnotSV/blob/master/commandLineOptions.txt):
  accepted genome builds include `CHM13`.
- For other **Expected** entries, the original repository links in the tables
  document the custom-reference or reference-aligned-input interface, not an
  explicit T2T validation claim.

For the exact rule dependencies and output paths, see the
[detailed workflow DAG](workflow-dag.md).
