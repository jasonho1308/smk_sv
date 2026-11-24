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

## Run the screen and assembly together

Request one `dux4_files.txt` per suspected sample (use the PEP `sample_name`). This runs the CHM13 read-level screen, both soft-clip BLAST reports, and local assembly, then lists the report and contig file paths relative to `run.dir`. No config sample list is needed, and the `all` and `dux4` targets are unchanged.

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16 \
  dux4/SAMPLE/dux4_files.txt
```

Pass multiple `dux4_files.txt` paths to run several samples. For GRCh38 inputs, request local assembly alone as below; the read-level screen requires CHM13v2.0.

## Optional local assembly of DUX4 partners

Local assembly is a follow-up for suspected DUX4 rearrangements. It accepts HiFi or ONT WGS aligned to a **coordinate-sorted, indexed BAM**. The default input is the workflow's minimap2 BAM/CSI; CHM13v2.0 defaults include the [chr4q35 and chr10q26 D4Z4 arrays](../workflow/resources/dux4/chm13v2.0.assembly_sources.bed) and **all 21 [partner loci](../workflow/resources/dux4/chm13v2.0.partners.bed)**, not just IGH. The read-level `scan_dux4` still uses a chr4-only source as described above.

Run assembly manually by requesting the breakpoint report for a suspected sample. Replace `SAMPLE` with its PEP `sample_name`, or pass several report paths to follow up on several samples. Neither `dux4` nor `all` runs local assembly, and no sample list is needed in `config/config.yaml`.

```bash
snakemake -s workflow/Snakefile --software-deployment-method conda --cores 12 \
  dux4/local_asm/SAMPLE/breakpoints.tsv
```

No additional assembly settings are needed for the default CHM13v2.0 ONT pipeline input. For HiFi or an external BAM/reference, use the same per-sample target with these optional assembly settings:

```yaml
caller_settings:  # merge into existing settings
  dux4_local_assembly:
    platform: hifi      # or ont (default)
    ont_mode: nano-hq   # or nano-raw for older ONT data
    bam_pattern: 'results/aligned/{sample}.bam'
    bam_index_pattern: 'results/aligned/{sample}.bam.bai'  # optional; can also point to CSI
    ref: /reference/GRCh38.fa  # omit to reuse inputs.fasta
    regions:
      sources: /reference/GRCh38.d4z4.bed
      partners: /reference/GRCh38.partners.bed
    partner_padding: 10000  # bp; set to 0 for exact BED windows
```

Both regions files must be BED4 (`chrom`, `start`, `end`, `locus_name`; 0-based, end-exclusive); the names are printed in the reports. Supply both files for GRCh38 or other custom FASTAs. Use exact reference coordinates and contig names matching the BAM, include both D4Z4 arrays where represented, and list any partners of interest (IGH is not required). Choose these windows independently for each assembly rather than reusing the CHM13 coordinates. Paths are relative to the configured run directory unless absolute. An alternate FASTA is indexed by the branch; use an uncompressed FASTA. This branch uses Conda environments and does not require Docker or Apptainer.

- **Extraction:** `samtools view` recruits read names from primary or supplementary alignments with a mapped footprint overlapping a target, including records with soft clips. A second pass over the whole BAM retrieves complete records even if the primary is off-target. Exactly one full sequence per read name is retained in `target_reads.fastq.gz`, with the original strand and quality order restored. Hard-clipped records are never exported; if the BAM has no full record for a selected read, it is omitted and counted in `extraction_stats.tsv`. Missing qualities use `!` and are counted. Breakpoint reads with no mapped segment overlapping any window will not be recruited.
- **Assembly:** HiFi uses hifiasm primary contigs without duplicate purging; ONT uses Flye `--meta --keep-haplotypes` for uneven coverage. Neither assembler guarantees recovery of low-VAF haplotypes. By default, fewer than `min_reads: 3` reads of `min_read_length: 1000` bp skips assembly; this is a minimum input check, not a coverage or VAF estimate. `assembly_status.tsv` distinguishes `no_reads`, `insufficient_reads`, `no_contigs`, and `assembled`. No contigs produces an empty `contigs.fa`, a valid indexed header-only BAM, and zero breakpoint rows. Unexpected tool errors still fail the job.
- **Evidence:** Load `contigs_aligned.bam` and `.bai` in IGV or Ribbon, and optionally load `breakpoints.bed` as a review track. `breakpoints.tsv` lists paired cut positions (0-based between-base), contig query spans in original orientation, strands, and MAPQs. A candidate must have neighboring primary/supplementary segments at a D4Z4 source and a configured partner, with at least `min_anchor_bases: 100` aligned bp inside each window and the junction-facing edge inside the window (or on its boundary). The partner needs `min_partner_mapq: 20`; D4Z4 may have low MAPQ. `max_query_gap: 1000` and `max_query_overlap: 100` define the allowed separation on the contig. The BED marks one aligned base beside each cut, grouped by candidate ID. Inspect the underlying reads to verify support: a split contig may be chimeric or ambiguous in repeats, and a D4Z4 placement does not prove chr4 vs. chr10 donor origin. No candidate does not rule out a rearrangement.

Sample outputs are under `dux4/local_asm/SAMPLE/`, relative to `run.dir`. Rule logs are in `logs/SAMPLE/`, and per-rule benchmarks are in `benchmarks/SAMPLE/`.
