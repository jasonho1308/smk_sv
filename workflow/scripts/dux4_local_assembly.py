#!/usr/bin/env python3
"""Read recovery, local assembly, and split-contig evidence for the DUX4 workflow."""

import argparse
from collections import defaultdict
from dataclasses import dataclass
import gzip
from itertools import groupby
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import TemporaryDirectory


COMPLEMENT = str.maketrans("ACMGRSVTWYHKDBNacmgrsvtwyhkdbn", "TGKCYSBAWRDMHVNtgkcysbawrdmhvn")
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str


def reference_lengths(fai):
    with open(fai) as handle:
        return {fields[0]: int(fields[1]) for fields in map(str.split, handle)}


def read_regions(path, lengths, padding=0):
    regions = []
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                chrom, start, end, name, *_ = line.split()
                start, end = int(start), int(end)
                if not 0 <= start < end <= lengths[chrom]:
                    raise ValueError("interval outside reference bounds")
            except (ValueError, KeyError) as exc:
                raise ValueError(f"Invalid BED4 interval at {path}:{number}: {line.strip()}") from exc
            regions.append(Region(chrom, max(0, start - padding), min(lengths[chrom], end + padding), name))
    if not regions:
        raise ValueError(f"No intervals in {path}")
    return regions


def targets(args):
    lengths = reference_lengths(args.fai)
    return (lengths, read_regions(args.sources, lengths),
            read_regions(args.partners, lengths, args.partner_padding))


def sam_records(lines):
    return (line.rstrip("\n").split("\t") for line in lines if not line.startswith("@"))


def write_metrics(path, **metrics):
    text = "metric\tvalue\n" + "".join(f"{key}\t{value}\n" for key, value in metrics.items())
    path.write_text(text)
    print(text, end="", flush=True)


def recover_fastq(lines, fastq):
    recovered = missing_qual = 0
    with gzip.open(fastq, "wt") as handle:
        for name, group in groupby(sam_records(lines), key=lambda row: row[0]):
            # Choose one sequence per unpaired molecule. Soft clips retain bases; hard clips do not.
            eligible = (row for row in group if "H" not in row[5] and row[9] != "*")
            row = max(eligible, key=lambda r: (len(r[9]), r[10] != "*", not (int(r[1]) & 0x800)), default=None)
            if row is None:
                continue
            seq, qual = row[9:11]
            if qual == "*":
                qual = "!" * len(seq)
                missing_qual += 1
            if int(row[1]) & 0x10:
                seq, qual = seq.translate(COMPLEMENT)[::-1], qual[::-1]
            handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
            recovered += 1
    return recovered, missing_qual


def extract(args):
    lengths, sources, partners = targets(args)
    header = subprocess.check_output(["samtools", "view", "-H", args.bam], text=True)
    bam_lengths = {}
    for line in header.splitlines():
        if line.startswith("@SQ\t"):
            tags = dict(field.split(":", 1) for field in line.split("\t")[1:])
            bam_lengths[tags["SN"]] = int(tags["LN"])
    for chrom in {region.chrom for region in sources + partners}:
        if bam_lengths.get(chrom) != lengths[chrom]:
            raise ValueError(f"BAM and FASTA target contig names/lengths differ for {chrom}")
    with TemporaryDirectory(prefix="dux4-extract-") as temp:
        tmp = Path(temp)
        bed = tmp / "targets.bed"
        bed.write_text("".join(f"{r.chrom}\t{r.start}\t{r.end}\n" for r in sources + partners))
        # Indexed first pass: mapped primary or supplementary alignments, no MAPQ filter.
        with subprocess.Popen(["samtools", "view", "-F", "0x104", "-M", "-L", bed,
                               "-X", args.bam, args.index], stdout=subprocess.PIPE, text=True) as proc:
            names = {row[0] for row in sam_records(proc.stdout)}
            if proc.wait():
                raise subprocess.CalledProcessError(proc.returncode, proc.args)
        recovered = missing_qual = 0
        if not names:
            with gzip.open(args.fastq, "wb"):
                pass
        else:
            # Full BAM pass: the complete sequence may be on a distant primary or unmapped record.
            names_file, sorted_bam = tmp / "names.txt", tmp / "recovered.bam"
            names_file.write_text("\n".join(sorted(names)) + "\n")
            with subprocess.Popen(["samtools", "view", "-u", "-F", "0x100", "-N", names_file,
                                   args.bam], stdout=subprocess.PIPE) as proc:
                subprocess.run(["samtools", "sort", "-N", "-@", str(max(0, args.threads - 2)),
                                "-m", "256M", "-T", tmp / "sort", "-o", sorted_bam, "-"],
                               stdin=proc.stdout, check=True)
                proc.stdout.close()
                if proc.wait():
                    raise subprocess.CalledProcessError(proc.returncode, proc.args)
            with subprocess.Popen(["samtools", "view", sorted_bam], stdout=subprocess.PIPE, text=True) as proc:
                recovered, missing_qual = recover_fastq(proc.stdout, args.fastq)
                if proc.wait():
                    raise subprocess.CalledProcessError(proc.returncode, proc.args)
    write_metrics(args.stats, selected_read_names=len(names), recovered_reads=recovered,
                  missing_full_sequence=len(names) - recovered, missing_qualities=missing_qual)


def fastq_counts(path, min_length):
    reads = long_reads = 0
    with gzip.open(path, "rt") as handle:
        for header in handle:
            seq = handle.readline().rstrip("\n")
            plus = handle.readline()
            qual = handle.readline().rstrip("\n")
            if not header.startswith("@") or not seq or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record {reads + 1} in {path}")
            reads += 1
            long_reads += len(seq) >= min_length
    return reads, long_reads


def assemble(args):
    reads, long_reads = fastq_counts(args.fastq, args.min_read_length)
    status = "assembled"
    args.contigs.write_text("")
    if not reads:
        status = "no_reads"
    elif long_reads < args.min_reads:
        status = "insufficient_reads"
    else:
        with TemporaryDirectory(prefix="dux4-assembly-") as work:
            tmp = Path(work)
            if args.platform == "hifi":
                prefix = tmp / "hifiasm"
                subprocess.run(["hifiasm", "-t", str(args.threads), "-f0", "-l0", "-o", prefix,
                                args.fastq], check=True)
                with args.contigs.open("w") as output:
                    subprocess.run(["awk", "-F", "\t", '$1 == "S" && $3 != "*" && $3 != "" {print ">"$2; print $3}',
                                    str(prefix) + ".bp.p_ctg.gfa"], stdout=output, check=True)
            else:
                out = tmp / "flye"
                cmd = ["flye", "--" + args.ont_mode, args.fastq, "--meta", "--keep-haplotypes",
                       "--threads", str(args.threads), "--out-dir", out]
                # Preserve the command output in the rule log, even on failure.
                run = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                print(run.stdout, end="", flush=True)
                if run.returncode:
                    if not re.search(r"No (?:disjointigs|contigs) were assembled", run.stdout):
                        run.check_returncode()
                else:
                    shutil.copyfile(out / "assembly.fasta", args.contigs)
        if not args.contigs.stat().st_size:
            status = "no_contigs"
    write_metrics(args.status, status=status, platform=args.platform, reads=reads,
                  reads_at_least_min_length=long_reads, min_read_length=args.min_read_length,
                  min_reads=args.min_reads)


@dataclass(frozen=True)
class Segment:
    chrom: str
    strand: str
    mapq: int
    qstart: int
    qend: int
    blocks: tuple

    def facing_cut(self, before):
        # If before the junction in query coordinates, use the high query edge.
        return self.blocks[-1][1] if before == (self.strand == "+") else self.blocks[0][0]


def segment(row):
    flag = int(row[1])
    if flag & 0x104:
        return None
    ref_pos, qpos = int(row[3]) - 1, 0
    blocks, query = [], []
    for count, op in CIGAR.findall(row[5]):
        n = int(count)
        if op in "M=X":
            blocks.append((ref_pos, ref_pos + n))
            query.append((qpos, qpos + n))
        if op in "MDN=X":
            ref_pos += n
        if op in "MISH=X":
            qpos += n
    if not blocks:
        return None
    lo, hi = query[0][0], query[-1][1]
    if flag & 0x10:
        lo, hi = qpos - hi, qpos - lo
    return Segment(row[2], "-" if flag & 0x10 else "+", int(row[4]), lo, hi, tuple(blocks))


def anchors(seg, regions, min_bases):
    return [r for r in regions if seg.chrom == r.chrom and min_bases <= sum(
        max(0, min(end, r.end) - max(start, r.start)) for start, end in seg.blocks)]


def split_pairs(records, sources, partners, min_bases, min_mapq, max_gap, max_overlap):
    by_query = defaultdict(list)
    for row in records:
        seg = segment(row)
        if seg is None:
            continue
        source_hits = anchors(seg, sources, min_bases)
        partner_hits = anchors(seg, partners, min_bases) if seg.mapq >= min_mapq else []
        if source_hits or partner_hits:
            by_query[row[0]].append((seg, source_hits, partner_hits))
    for name, hits in sorted(by_query.items()):
        for src, source_hits, _ in hits:
            for dst, _, partner_hits in hits:
                if not source_hits or not partner_hits:
                    continue
                first, second = sorted((src, dst), key=lambda s: s.qstart)
                gap = second.qstart - first.qend
                if first.qstart >= second.qstart or first.qend >= second.qend or not -max_overlap <= gap <= max_gap:
                    continue
                src_before = src is first
                src_cut, dst_cut = src.facing_cut(src_before), dst.facing_cut(not src_before)
                if src.chrom == dst.chrom and src.strand == dst.strand:
                    ref_gap = second.facing_cut(False) - first.facing_cut(True)
                    if src.strand == "-":
                        ref_gap = -ref_gap
                    if abs(ref_gap - gap) <= max(max_gap, max_overlap):
                        continue  # Nearby collinear segments are not rearrangement evidence.
                for a in source_hits:
                    for b in partner_hits:
                        if (a.chrom == b.chrom and a.start < b.end and b.start < a.end
                                or not a.start <= src_cut <= a.end or not b.start <= dst_cut <= b.end):
                            continue
                        yield name, src, a, src_cut, dst, b, dst_cut, gap, src_before


def breakpoint_reports(sam, args, sources, partners):
    header = ("candidate\tcontig\tsource_locus\tsource_chrom\tsource_cut\tsource_strand\tsource_mapq\t"
              "source_qstart\tsource_qend\tpartner_locus\tpartner_chrom\tpartner_cut\tpartner_strand\t"
              "partner_mapq\tpartner_qstart\tpartner_qend\tquery_gap\n")
    bed_rows = []
    with sam.open() as handle, args.summary.open("w") as tsv:
        tsv.write(header)
        pairs = split_pairs(sam_records(handle), sources, partners, args.min_anchor_bases,
                            args.min_partner_mapq, args.max_query_gap, args.max_query_overlap)
        for i, (name, src, source, src_cut, dst, partner, dst_cut, gap, src_before) in enumerate(pairs, 1):
            candidate = f"split_{i}"
            row = [candidate, name, source.name, src.chrom, src_cut, src.strand, src.mapq, src.qstart, src.qend,
                   partner.name, dst.chrom, dst_cut, dst.strand, dst.mapq, dst.qstart, dst.qend, gap]
            tsv.write("\t".join(map(str, row)) + "\n")
            for seg, region, cut, before in ((src, source, src_cut, src_before),
                                            (dst, partner, dst_cut, not src_before)):
                start = cut - int(before == (seg.strand == "+"))
                bed_rows.append((seg.chrom, start, start + 1, f"{candidate}|{name}|{region.name}", 0, seg.strand))
    with args.bed.open("w") as handle:
        for row in sorted(bed_rows):
            handle.write("\t".join(map(str, row)) + "\n")
    print(f"Split-contig candidate pairs: {len(bed_rows) // 2}", flush=True)


def align(args):
    lengths, sources, partners = targets(args)
    with TemporaryDirectory(prefix="dux4-align-") as work:
        sam = Path(work) / "contigs.sam"
        with sam.open("w") as output:
            if args.contigs.stat().st_size:
                subprocess.run(["minimap2", "-ax", args.preset, "--cs", "-Y", "-t", str(args.threads),
                                args.ref, args.contigs], stdout=output, check=True)
            else:
                print("No assembled contigs; producing a header-only BAM and empty evidence reports.")
                output.write("@HD\tVN:1.6\tSO:coordinate\n" + "".join(
                    f"@SQ\tSN:{chrom}\tLN:{length}\n" for chrom, length in lengths.items()))
        workers = str(max(0, args.threads - 1))
        subprocess.run(["samtools", "sort", "-@", workers, "-m", "256M", "-T", Path(work) / "sort",
                        "-o", args.bam, sam], check=True)
        subprocess.run(["samtools", "index", "-@", workers, "-b", args.bam, args.index], check=True)
        breakpoint_reports(sam, args, sources, partners)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for func, paths in ((extract, "bam index fastq stats"), (assemble, "fastq contigs status"),
                        (align, "contigs ref bam index summary bed")):
        cmd = commands.add_parser(func.__name__)
        cmd.set_defaults(func=func)
        cmd.add_argument("--threads", type=int, required=True)
        for path in paths.split():
            cmd.add_argument("--" + path, type=Path, required=True)
        if func in (extract, align):
            for path in ("fai", "sources", "partners"):
                cmd.add_argument("--" + path, type=Path, required=True)
            cmd.add_argument("--partner-padding", type=int, required=True)
        if func is assemble:
            cmd.add_argument("--platform", choices=("hifi", "ont"), required=True)
            cmd.add_argument("--ont-mode", choices=("nano-hq", "nano-raw"), required=True)
            cmd.add_argument("--min-reads", type=int, required=True)
            cmd.add_argument("--min-read-length", type=int, required=True)
        if func is align:
            cmd.add_argument("--preset", choices=("asm5", "asm20"), required=True)
            for threshold in ("min-anchor-bases", "min-partner-mapq", "max-query-gap", "max-query-overlap"):
                cmd.add_argument("--" + threshold, type=int, required=True)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
