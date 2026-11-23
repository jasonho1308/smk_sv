# Targeted DUX4 rearrangement screen (T2T-CHM13v2.0)

Use this opt-in screen for ONT DNA reads after minimap2 alignment. It reports junction candidates at the 21 requested partner loci and separately searches soft clips across the BAM. It does not run the genome-wide SV callers. An existing pipeline BAM and CSI are reused.

```bash
# All configured samples (add -n to preview)
snakemake -s workflow/Snakefile --use-conda --cores 16 dux4

# Just one sample (replace SAMPLE with a PEP sample_name)
snakemake -s workflow/Snakefile --use-conda --cores 16 dux4/SAMPLE/SAMPLE.junctions.tsv
```

`scan_dux4` requests up to 16 threads per sample and passes the allocation to `blastn`. `--cores` sets the total CPU budget for a local run; Snakemake scales jobs down if the budget is smaller. To benchmark a different per-sample allocation, use `--set-threads scan_dux4=N` with `--cores` at least `N`.

Files are relative to the configured run directory. The indexed partner windows are the [21 BED intervals](../workflow/resources/dux4/chm13v2.0.partners.bed) plus 10 kb of padding on each side. The table uses the pipeline's pinned [RefSeq/Liftoff v5.3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz), including the full annotated IGH locus. LINC02193 is absent from that file; its chr16 primary interval comes from the pinned [GENCODEv35 CAT/Liftoff annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz). Its CAT alternate contig placement on chr2 is not included. The [D4Z4 target](../workflow/resources/dux4/chm13v2.0.sources.bed) is **chr4q35 only**, from 2 kb upstream of DUX4L8 through the annotated terminal DUX4 gene. The chr10q26 D4Z4 array is not searched as a source, but all chr10 partner loci remain in the scan.

## Evidence and interpretation

- The scan starts with segments with ≥200 matched/aligning CIGAR bases and MAPQ ≥20 near a partner. It checks for junction-adjacent D4Z4 alignments in minimap2 `SA` tags, including low-MAPQ repeat segments. Unresolved soft clips (≥100 bp, up to 3 kb from the clip boundary) and internal insertions with ≥200 bp of aligned sequence on each side are searched by local `blastn` against the chr4q35 interval. BLAST requires ≥100 bp and ≥80% identity with the hit within 100 bp of a junction edge.
- `*.junctions.tsv` lists clusters by partner and breakend side, within 100 bp. `review` means ≥2 distinct supporting read names, not a validated call. `nearby_opposite_ends` links L/R clusters within 100 bp supported by the same read(s) (`junction_id:read_count`). `insertion_like` indicates DUX4-region similarity within a CIGAR insertion; paired ends, read inspection, or local assembly can help resolve the structure. One junction alone does not distinguish an insertion from a translocation or a complex event.
- `*.reads.tsv` gives read names, evidence methods, candidate source placements and breakend sides for review in IGV. L/R is the reference edge of the aligned segment facing the junction; I is an internal insertion. Cuts are **0-based between-base positions**: cut `550` is between 1-based bases `550` and `551`. A `flank_placement` requires an `SA` alignment with MAPQ ≥20 extending at least 250 aligned bp beyond the source interval on one side. It is a provisional genome-wide reference placement, not proof of donor origin. All other hits, including BLAST-only hits, are `ambiguous`; a match to the chr4-only target does not establish chr4 as the donor, because homologous repeats on chr10 can also match.
- `*.loci.tsv` lists each examined locus, the scan window, distinct high-quality anchor read names, and candidate counts, including zeros. Anchor counts are not local breakpoint coverage. At 15×, no candidates does not rule out a rearrangement, especially with low tumor purity, a breakpoint outside the windows, or poor mappability. A DNA junction does not by itself establish an expressed fusion.

## Soft-clip BLAST comparison

- `dux4/SAMPLE/SAMPLE.blast.tsv` reports the best significant BLAST hit per soft-clipped alignment end across the BAM. It searches the full soft-clipped sequence at either end of each mapped alignment that stores that sequence, including alignments outside the partner windows and low-MAPQ alignments. The partner anchor, clip-length, identity, and junction-edge filters are not applied. The `partner` column is `.` because this search is genome-wide.
- `dux4/SAMPLE/SAMPLE.blast_filtered.tsv` reports only soft-clip BLAST hits from the original targeted screen, using the criteria above. Split-read and insertion evidence remain in `*.reads.tsv` but are not included in this soft-clip comparison.

Both files use the same chr4q35 target and BLAST E-value cutoff (1e-8). Clips with no qualifying hit are not listed. The genome-wide search does not change the junction, read, or locus reports. `logs/SAMPLE/dux4.log` records split evidence and query counts for both searches. Benchmark one sample on an existing indexed BAM before planning a larger run.

## Optional DRAGEN v4.3 attempt

The `dux4_dragen` target runs Illumina's [DUX4 Rearrangement Caller](https://help.dragen.illumina.com/dragen-v4.3/product-guide/dragen-v4.3/dragen-dna-pipeline/dux4-rearrangement-caller) separately. It does not change the `dux4` ONT screen, the default `all` target, or the SV consensus. Run it where licensed DRAGEN is available as `dragen` on `PATH`.

For the default CHM13v2.0 pipeline BAM, no extra configuration is needed. The first `dux4_dragen` run downloads and extracts the [Illumina v4.3 CHM13_v2 v4 linear hash table](https://support.illumina.com/sequencing/sequencing_software/dragen-bio-it-platform/product_files.html) into `references/dragen/chm13_v2-linear-r4.0/` under the run directory. Later samples reuse it. The archive is about 14 GB; allow space for both the archive and the extracted files during installation. The default `all`, `resources`, and `dux4` targets do not download it. To prefetch without running DRAGEN:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 1 download_dragen_dux4_hash_table
```

To use an already installed table instead, set `caller_settings.dragen_dux4.hash_table` in `config/config.yaml`. For an external BAM, set both `hash_table` and `bam_pattern`:

```yaml
  dragen_dux4:
    hash_table: /path/to/dragen/hg38_hash_table
    bam_pattern: /path/to/hg38/{sample}.sorted.markdup.bam
    skip_sanity_check: true
```

The hash table must match the BAM's reference exactly. The caller primarily supports hg38. If your pipeline BAM was aligned to a custom reference, use a table built from that exact FASTA rather than the automatic download. Never pair a CHM13 BAM with an hg38 table. The rule passes an already aligned and sorted BAM to DRAGEN without remapping or marking duplicates. Use a duplicate-marked BAM when available, as in Illumina's example.

`skip_sanity_check` defaults to `true`, using Illumina's documented `--dux4-skip-santiy-check` option, so this branch does not enforce a 25× cutoff. Illumina validated tumor-only, paired-end WGS at 25–120× mean coverage, with 100–151 bp reads and 300–500 bp fragments. At ~15×, ONT reads, lack of paired-end evidence, and CHM13 all fall outside the documented setup. The bypass does not make these inputs supported; DRAGEN may still reject them. A VCF with no calls does not rule out a rearrangement. Set the option to `false` to keep DRAGEN's checks enabled.

```bash
# One sample first (add -n to preview; --use-conda is for any upstream ONT alignment).
snakemake -s workflow/Snakefile --use-conda --cores 4 dux4_dragen/SAMPLE/attempt.json
# All configured samples, one DRAGEN job at a time on a single DRAGEN host.
snakemake -s workflow/Snakefile --use-conda --cores 4 dux4_dragen --resources dragen_slots=1
```

Check `dux4_dragen/SAMPLE/attempt.json` for `vcf_written`, `no_vcf`, `failed`, or `not_started`. `vcf_written` means a VCF file exists, not that a rearrangement was found. DRAGEN outputs are under `dux4_dragen/SAMPLE/dragen/` and its log is `logs/SAMPLE/dux4_dragen.log`. Failed invocations record status and a warning without publishing a partial or placeholder VCF. After fixing a recorded failure, rerun the sample with `--force`.
