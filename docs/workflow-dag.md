# Detailed Snakemake DAG

This is a **source-derived draft**, not a successfully executed `snakemake --dag`
export. It follows `workflow/Snakefile`, the included rule files, and the five
callers and three annotators selected in `config/config.yaml`.

The editable Graphviz equivalent is [workflow-dag.dot](workflow-dag.dot).

## Scope and notation

- Arrows point from prerequisites to consumers. Boxes name actual Snakemake rules;
  labels identify their principal products. Input/output nodes are artifacts, not
  additional jobs.
- `{sample}` expands over PEP `sample_name` values.
- `{caller}` expands over `cutesv`, `severus`, `sniffles`, `svim`, and `svision`.
- `{type_sv}` expands over `DEL`, `DUP`, `INV`, `BND`, and `INS`.
- SVision additionally expands over `chr1`–`chr22`, `chrX`, and `chrY`.
- The common caller-processing block represents **five parallel copies per
  sample**; its type-specific stages represent **25 copies per sample**. Each
  SURVIVOR rule represents five jobs per sample, gathering all five callers for
  one SV type. Shared reference/model/cache preparation is not repeated per sample.
- The diagrams show declared file dependencies, not every file a tool may open
  internally. Temporary files, logs, and most reports are omitted for readability.

## 1. Main analysis DAG

Reference, model, and annotation-cache inputs are separated into the resource
view below to keep the main graph readable.

```mermaid
flowchart TD
    FASTQ[/"{dir_data}/{sample}{suffix_fastq}"/]
    MMI[/"Minimap2 reference index (.mmi)"/]

    subgraph ALIGN["Alignment — per sample"]
        minimap2["minimap2<br/>ONT alignment and coordinate sorting<br/>{sample}.sorted.bam + .bam.csi"]
    end

    subgraph CALL["Configured caller branches — per sample"]
        cutesv["cutesv<br/>cutesv.vcf"]
        sniffles["sniffles<br/>sniffles.vcf"]
        svim["svim<br/>variants.vcf + svim.vcf symlink"]
        clair3["clair3<br/>phased_merge_output.vcf.gz<br/>haplotagged BAM + BAI"]
        severus["severus<br/>all_SVs/severus_all.vcf"]
        format_severus["format_severus<br/>severus.vcf"]
        svision["svision × 24 chromosomes<br/>chroms/{chrom}/*.graph.vcf"]
        format_svision["format_svision<br/>gather chromosomes and normalize VCF<br/>svision.vcf"]
    end

    subgraph PREP["Common processing — per caller × sample"]
        duphold["duphold<br/>remove invalid POS; sort; annotate depth<br/>{caller}.duphold.vcf"]
        hard_filter["hard_filter<br/>caller-specific support, coverage,<br/>size and depth-fold-change filters<br/>{caller}.duphold.filtered.vcf"]
        separate_types["separate_types × 5 SV types<br/>{caller}.{type_sv}.vcf"]
    end

    subgraph CONSENSUS["First consensus — per sample × SV type"]
        survivor["survivor<br/>gather type-specific VCFs from all callers<br/>minimum support: n_callers[type_sv]<br/>{sample}.{type_sv}.merged.vcf"]
    end

    subgraph ANNOTATE["Parallel annotation — per caller × sample × SV type"]
        vep["vep<br/>{caller}.{type_sv}.vep.vcf + HTML"]
        snpeff["snpeff<br/>{caller}.{type_sv}.snpeff.vcf + HTML"]
        annotsv["annotsv<br/>{caller}.{type_sv}.annotsv.tsv"]
        restore_annotsv_length["restore_annotsv_length<br/>restore SV length using original VCF<br/>{caller}.{type_sv}.annotsv.length.tsv"]
    end

    subgraph REFINE["Consensus selection and annotation filtering — per caller × sample × SV type"]
        extract_annotations["extract_annotations<br/>map first-consensus records back to caller IDs<br/>subset original, VEP and SnpEff VCFs<br/>write cross-caller ID table"]
        convert_vep["convert_vep<br/>normalize breakpoints; vcf2maf<br/>merged/{caller}.{type_sv}.vep.maf"]
        convert_snpeff["convert_snpeff<br/>SnpSift extractFields<br/>merged/{caller}.{type_sv}.snpeff.tsv"]
        filter_annotations["filter_annotations<br/>combine AnnotSV, VEP, SnpEff and caller IDs<br/>select IDs; filter three VCFs<br/>write TSV + RData"]
    end

    survivor_final["survivor_final — per sample × SV type<br/>gather annotation-filtered caller VCFs<br/>minimum support: n_callers[type_sv] - 1"]
    FINAL[/"survivor/{sample}/final/{sample}.{type_sv}.merged.vcf<br/>DEL · DUP · INV · BND · INS"/]
    all["all<br/>requires all five final VCFs for every sample"]

    FASTQ --> minimap2
    MMI --> minimap2
    minimap2 --> cutesv
    minimap2 --> sniffles
    minimap2 --> svim
    minimap2 --> clair3
    minimap2 --> svision
    clair3 -->|"phased VCF + haplotagged BAM"| severus
    severus --> format_severus
    svision --> format_svision

    cutesv --> duphold
    sniffles --> duphold
    svim --> duphold
    format_severus --> duphold
    format_svision --> duphold
    minimap2 -->|"original sorted BAM"| duphold
    duphold --> hard_filter --> separate_types

    separate_types -->|"all callers, same sample and SV type"| survivor
    separate_types --> vep
    separate_types --> snpeff
    separate_types --> annotsv
    separate_types -->|"original type-specific VCF"| restore_annotsv_length
    annotsv --> restore_annotsv_length

    survivor -->|"consensus membership"| extract_annotations
    separate_types -->|"original type-specific VCF"| extract_annotations
    vep --> extract_annotations
    snpeff --> extract_annotations
    extract_annotations --> convert_vep
    extract_annotations --> convert_snpeff
    extract_annotations -->|"three selected VCFs + ID table"| filter_annotations
    restore_annotsv_length --> filter_annotations
    convert_vep --> filter_annotations
    convert_snpeff --> filter_annotations

    filter_annotations -->|"all callers, same sample and SV type"| survivor_final
    survivor_final --> FINAL --> all

    classDef input fill:#e2e8f0,stroke:#64748b,color:#0f172a;
    classDef caller fill:#dbeafe,stroke:#2563eb,color:#0f172a;
    classDef process fill:#ffedd5,stroke:#ea580c,color:#0f172a;
    classDef annotation fill:#ede9fe,stroke:#7c3aed,color:#0f172a;
    classDef consensus fill:#dcfce7,stroke:#16a34a,color:#0f172a;
    class FASTQ,MMI input;
    class cutesv,sniffles,svim,clair3,severus,format_severus,svision,format_svision caller;
    class minimap2,duphold,hard_filter,separate_types process;
    class vep,snpeff,annotsv,restore_annotsv_length,extract_annotations,convert_vep,convert_snpeff,filter_annotations annotation;
    class survivor,survivor_final,FINAL,all consensus;
```

### Important dependency details

1. **Annotation is caller-specific, not annotation of the merged VCF.** `vep`,
   `snpeff`, and `annotsv` all consume `{caller}/{sample}/{caller}.{type_sv}.vcf`.
   They can run concurrently with the first `survivor` job. `extract_annotations`
   waits for that consensus and the VEP/SnpEff outputs, then subsets caller VCFs
   by their consensus-associated IDs.
2. **Clair3 is a prerequisite of Severus, not an SV consensus contributor.** It
   supplies phasing and haplotagging. It is included in the Snakefile but is not
   listed in the configured `CALLERS` collection.
3. **duphold always uses the original mapper BAM**, including for Severus-derived
   VCFs; only Severus itself uses the Clair3 haplotagged BAM.
4. **AnnotSV length restoration uses the original type-specific VCF**, not a
   consensus-subset VCF. `filter_annotations` subsequently restricts AnnotSV rows
   using the consensus ID table.
5. **Two consensus thresholds are used.** In the current configuration the first
   merge requires three callers and the final merge requires two, independently
   for each SV type. Both use `distance_sv[type_sv] = 10`, ignore type and strand
   in SURVIVOR's matching flags, and enable size-based distance estimation via
   `estimate_distance`. The distance setting should therefore not be described
   simply as a fixed 10-bp window.
6. **There is no separate somatic/germline classification output in this DAG.**
   The final targets are annotation-filtered consensus SV VCFs. Tumor-only input
   and annotation filtering do not establish somatic status on their own.

## 2. Shared resource DAG

These rules provide resources for the main graph. Existing resources can be
reused; missing ones are prepared through the corresponding producer rules.
`resources` is an alternative aggregate target, **not** a prerequisite of `all`.
The normal analysis depends directly on the resource files it needs.

```mermaid
flowchart LR
    download_reference["download_reference<br/>reference FASTA"]
    minimap2_index["minimap2_index<br/>reference .mmi"]
    faidx["faidx<br/>reference .fai"]
    download_tandem_repeats["download_tandem_repeats<br/>shared tandem-repeat / VNTR BED"]
    download_clair3_model["download_clair3_model<br/>Clair3 model directory"]
    download_svision_model["download_svision_model<br/>checkpoint .data / .index / .meta"]
    download_vep_cache["download_vep_cache"]
    download_snpeff_cache["download_snpeff_cache"]
    download_annotsv_human["download_annotsv_human<br/>Annotations_Human"]
    download_annotsv_exomiser["download_annotsv_exomiser<br/>versioned Exomiser directory"]

    minimap2["minimap2"]
    reference_users["FASTA consumers:<br/>cutesv, sniffles, svim, svision,<br/>clair3, duphold, vep, snpeff, convert_vep"]
    clair3["clair3"]
    sniffles["sniffles"]
    severus["severus"]
    svision["svision"]
    vep["vep"]
    snpeff["snpeff"]
    annotsv["annotsv"]
    resources["resources<br/>aggregate preparation target"]

    download_reference --> minimap2_index --> minimap2
    download_reference --> faidx --> clair3
    download_reference --> reference_users
    download_clair3_model --> clair3
    download_svision_model --> svision
    download_tandem_repeats --> sniffles
    download_tandem_repeats --> severus
    download_vep_cache --> vep
    download_snpeff_cache --> snpeff
    download_annotsv_human --> annotsv
    download_annotsv_exomiser --> annotsv

    download_reference --> resources
    minimap2_index --> resources
    faidx --> resources
    download_tandem_repeats --> resources
    download_clair3_model --> resources
    download_svision_model --> resources
    download_vep_cache --> resources
    download_snpeff_cache --> resources
    download_annotsv_human --> resources
    download_annotsv_exomiser --> resources
```

`caller_settings.nvtr` and `caller_settings.tandem_repeats` currently point to the
same BED, so `download_tandem_repeats` serves both Severus and Sniffles. If those
paths differ, a separate `download_nvtr` rule supplies Severus and the `resources`
target. NanoSV-only resources are outside the configured analysis.

## 3. Key artifacts

Paths below are relative to `run.dir`, except configured input/cache paths.

| Stage | Rule | Principal output |
| --- | --- | --- |
| Alignment | `minimap2` | `minimap2/{sample}/{sample}.sorted.bam` and `.bam.csi` |
| Phasing | `clair3` | `clair3/{sample}/phased_merge_output.vcf.gz`, phased BAM/BAI, and `{sample}.sorted.haplotagged.bam`/BAI symlinks |
| Caller normalization | caller rule or formatter | `{caller}/{sample}/{caller}.vcf` |
| Depth annotation | `duphold` | `{caller}/{sample}/{caller}.duphold.vcf` |
| Hard filtering | `hard_filter` | `{caller}/{sample}/{caller}.duphold.filtered.vcf` |
| Type splitting | `separate_types` | `{caller}/{sample}/{caller}.{type_sv}.vcf` |
| First consensus | `survivor` | `survivor/{sample}/{sample}.{type_sv}.merged.vcf` |
| Annotation | `vep`, `snpeff`, `annotsv` | `{caller}/{sample}/{caller}.{type_sv}.{vep.vcf,snpeff.vcf,annotsv.tsv}` |
| Length correction | `restore_annotsv_length` | `{caller}/{sample}/{caller}.{type_sv}.annotsv.length.tsv` |
| Consensus-associated subsets | `extract_annotations` | `{caller}/{sample}/merged/{caller}.{type_sv}.{vcf,vep.vcf,snpeff.vcf}` and `survivor/{sample}/{sample}.{caller}.{type_sv}.tab` |
| Annotation conversion | `convert_vep`, `convert_snpeff` | `{caller}/{sample}/merged/{caller}.{type_sv}.{vep.maf,snpeff.tsv}` |
| Annotation filtering | `filter_annotations` | `{caller}/{sample}/merged/filtered/{caller}.{type_sv}.{vcf,vep.vcf,snpeff.vcf,tsv,RData}` |
| Final consensus | `survivor_final` | `survivor/{sample}/final/{sample}.{type_sv}.merged.vcf` |

Brace suffix lists in this table are shorthand, not literal filenames.

## 4. Other available caller branches

The following modules exist but are **not selected** in the current configuration.
They are deliberately excluded from the main DAG rather than presented as active
or validated branches.

| Caller | Source-declared chain before `duphold` | Caveat |
| --- | --- | --- |
| DeBreak | mapper BAM + FASTA → `debreak` → `format_debreak` → `debreak.vcf` | Formatter rename-table path overlaps other formatter outputs if enabled together. |
| Delly | mapper BAM + FASTA → `delly` → `format_delly` → `delly.vcf` | Formatter writes `debreak/{sample}/rename.tab`, also used by `format_severus`. |
| NanoSV | mapper BAM + NanoSV config → `nanosv` → **broken formatter dependency** | `nanosv` produces `tmp.vcf`, but `format_nanosv` declares its own final VCF as both input and output. Sampling BED is passed via `params`, not declared as an analysis input. |
| NanoVar | mapper BAM + FASTA → `nanovar` → formatter → `nanovar.vcf` | Its formatter is named `format_nanosv`, colliding with NanoSV's formatter if both modules are included. |

Although the schema permits reduced annotator lists, the downstream
`extract_annotations` and `filter_annotations` rules unconditionally require all
three annotation tools' products. This draft therefore assumes all three remain
enabled. Similarly, adding `clair3` to `workflow.callers` is not represented as a
valid SV branch: Clair3 is already included separately and does not produce the
`{caller}/{sample}/{caller}.vcf` expected by `duphold`.

## 5. Render or regenerate

### Render the source-derived draft

Open this Markdown file in a Mermaid-capable viewer, or use Graphviz for the
combined analysis/resource graph:

```bash
# From the repository root; requires Graphviz's `dot` executable.
dot -Tsvg docs/workflow-dag.dot -o docs/workflow-dag.svg
dot -Tpdf docs/workflow-dag.dot -o docs/workflow-dag.pdf
```

### Export Snakemake's actual instantiated DAG

After configuring the environment, real samples, and resources:

```bash
# Explicitly select `all`: resource rules are included before rule all.
snakemake --snakefile workflow/Snakefile all --dry-run

# One node per concrete job, expanded across samples/callers/types/chromosomes.
# Separate generation and rendering so a failed Snakemake command does not
# silently leave an empty PDF through a shell pipeline.
snakemake --snakefile workflow/Snakefile all --dag > dag.dot &&
    dot -Tpdf dag.dot -o dag.pdf

# Compact rule-level view.
snakemake --snakefile workflow/Snakefile all --rulegraph > rulegraph.dot &&
    dot -Tsvg rulegraph.dot -o rulegraph.svg

# Resource-preparation graph only.
snakemake --snakefile workflow/Snakefile resources --dag > resources.dot &&
    dot -Tsvg resources.dot -o resources.svg
```

`--rulegraph` collapses wildcard instances; use `--dag` when the actual job fan-out
is important. Depending on existing outputs, an instantiated DAG may differ from
the full source-derived view.

### Local export limitations observed while drafting

- The installed `snakemake` entry point fails during import with
  `ModuleNotFoundError: No module named 'tenacity'`.
- Graphviz's `dot` executable was not found on `PATH`.
- `config/pep/config.yaml` and its sample table are absent; only the PEP template
  `config/pep/.config.yaml` is present.
- `config/config.yaml` has `resources.downloads: {}`. The AnnotSV resource output
  helper requires `annotsv_exomiser.version` even when caches already exist;
  missing downloadable inputs also need appropriate URLs. The configuration
  template `config/.config.yaml` contains example resource definitions.
- The configured input and run paths refer to an external analysis environment.

No pipeline rules, configuration, or existing `dag.pdf` were changed to produce
this documentation.
