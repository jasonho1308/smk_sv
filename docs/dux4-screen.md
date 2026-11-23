# Targeted DUX4 rearrangement screen (T2T-CHM13v2.0)

Use this opt-in screen for ONT DNA reads after minimap2 alignment. It scans the 21 requested partner loci first; it does not run the genome-wide SV callers or search every soft clip. An existing pipeline BAM and CSI are reused.

```bash
# All configured samples (add -n to preview)
snakemake -s workflow/Snakefile --use-conda --cores 4 dux4

# Just one sample (replace SAMPLE with a PEP sample_name)
snakemake -s workflow/Snakefile --use-conda --cores 4 dux4/SAMPLE/SAMPLE.junctions.tsv
```

Files are relative to the configured run directory. The indexed partner windows are the [21 BED intervals](../workflow/resources/dux4/chm13v2.0.partners.bed) plus 10 kb of padding on each side. The table uses the pipeline's pinned [RefSeq/Liftoff v5.3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RefSeq_Liftoff_v5.3.gff.gz), including the full annotated IGH locus. LINC02193 is absent from that file; its chr16 primary interval comes from the pinned [GENCODEv35 CAT/Liftoff annotation](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_GENCODEv35_CAT_Liftoff.vep.gff3.gz). Its CAT alternate contig placement on chr2 is not included. The [D4Z4 target](../workflow/resources/dux4/chm13v2.0.sources.bed) is **chr4q35 only**, from 2 kb upstream of DUX4L8 through the annotated terminal DUX4 gene. The chr10q26 D4Z4 array is not searched as a source, but all chr10 partner loci remain in the scan.

## Evidence and interpretation

- The scan starts with segments with ≥200 matched/aligning CIGAR bases and MAPQ ≥20 near a partner. It checks for junction-adjacent D4Z4 alignments in minimap2 `SA` tags, including low-MAPQ repeat segments. Unresolved soft clips (≥100 bp, up to 3 kb from the clip boundary) and internal insertions with ≥200 bp of aligned sequence on each side are searched by local `blastn` against the chr4q35 interval. BLAST requires ≥100 bp and ≥80% identity with the hit within 100 bp of a junction edge.
- `*.junctions.tsv` lists clusters by partner and breakend side, within 100 bp. `review` means ≥2 distinct supporting read names, not a validated call. `nearby_opposite_ends` links L/R clusters within 100 bp supported by the same read(s) (`junction_id:read_count`). `insertion_like` indicates DUX4-region similarity within a CIGAR insertion; paired ends, read inspection, or local assembly can help resolve the structure. One junction alone does not distinguish an insertion from a translocation or a complex event.
- `*.reads.tsv` gives read names, evidence methods, candidate source placements and breakend sides for review in IGV. L/R is the reference edge of the aligned segment facing the junction; I is an internal insertion. Cuts are **0-based between-base positions**: cut `550` is between 1-based bases `550` and `551`. A `flank_placement` requires an `SA` alignment with MAPQ ≥20 extending at least 250 aligned bp beyond the source interval on one side. It is a provisional genome-wide reference placement, not proof of donor origin. All other hits, including BLAST-only hits, are `ambiguous`; a match to the chr4-only target does not establish chr4 as the donor, because homologous repeats on chr10 can also match.
- `*.loci.tsv` lists each examined locus, the scan window, distinct high-quality anchor read names, and candidate counts, including zeros. Anchor counts are not local breakpoint coverage. At 15×, no candidates does not rule out a rearrangement, especially with low tumor purity, a breakpoint outside the windows, or poor mappability. A DNA junction does not by itself establish an expressed fusion.

`logs/SAMPLE/dux4.log` records the split-evidence and BLAST query counts. Benchmark one sample on an existing indexed BAM before planning a larger run.
