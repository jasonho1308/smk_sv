#!/usr/bin/env python3
"""Screen gene-anchored ONT reads for DUX4/D4Z4 junction evidence."""

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass, fields
from pathlib import Path
import re
from statistics import median_low
import subprocess
from tempfile import TemporaryDirectory

import pysam


MATCH = set("M=X")
CLIP = set("SH")
REF = MATCH | set("DN")
CHM13_LENGTHS = {"chr4": 193574945, "chr10": 134758134, "chr14": 101161492}


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int
    name: str


def read_bed(path):
    regions = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name = line.split()[:4]
            if not 0 <= int(start) < int(end):
                raise ValueError(f"Invalid BED interval: {line.strip()}")
            regions.append(Region(chrom, int(start), int(end), name))
    if not regions or len({r.name for r in regions}) != len(regions):
        raise ValueError(f"Expected nonempty BED with unique names: {path}")
    return regions


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    strand: str
    mapq: int
    cigar: str

    @property
    def ops(self):
        return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", self.cigar)]

    def layout(self):
        ops = self.ops
        left = right = 0
        while ops and ops[0][1] in CLIP:
            left += ops.pop(0)[0]
        while ops and ops[-1][1] in CLIP:
            right += ops.pop()[0]
        aligned = sum(n for n, op in ops if op in MATCH)
        span = aligned + sum(n for n, op in ops if op == "I")
        offset = left if self.strand == "+" else right
        return offset, offset + span, left + span + right, aligned

    def blocks(self):
        ref = self.start
        for n, op in self.ops:
            if op in MATCH:
                yield ref, ref + n
            if op in REF:
                ref += n

    def cut(self, side):
        return self.start if side == "L" else self.start + sum(n for n, op in self.ops if op in REF)


@dataclass(frozen=True)
class Anchor:
    partner: str
    read: str
    partner_chrom: str
    partner_cut: int
    partner_side: str
    partner_mapq: int


@dataclass(frozen=True)
class Evidence:
    anchor: Anchor
    method: str
    source_status: str
    source_target: str
    source_chrom: str
    source_cut: int | None
    source_side: str
    source_mapq: int | None = None
    query_gap: int | None = None
    blast_identity: float | None = None
    blast_length: int | None = None


def split_evidence(anchor, segment, sa_segments, sources, args):
    qstart, qend, qlen, _ = segment.layout()
    after = (anchor.partner_side == "R") == (segment.strand == "+")
    best = None
    for sa in sa_segments:
        start, end, length, aligned = sa.layout()
        if length != qlen or aligned < args.min_clip:
            continue
        # Both alignments must extend away from the junction in original-read coordinates.
        if (after and end <= qend) or (not after and start >= qstart):
            continue
        gap = start - qend if after else qstart - end
        if abs(gap) > args.max_query_gap:
            continue
        side = "L" if after == (sa.strand == "+") else "R"
        cut = sa.cut(side)
        for source in sources:
            if sa.chrom != source.chrom or not source.start - args.source_flank <= cut <= source.end + args.source_flank:
                continue
            blocks = list(sa.blocks())
            core = sum(max(0, min(e, source.end) - max(s, source.start)) for s, e in blocks)
            if core < args.min_clip:
                continue
            west = sum(max(0, min(e, source.start) - s) for s, e in blocks)
            east = sum(max(0, e - max(s, source.end)) for s, e in blocks)
            placed = sa.mapq >= args.min_mapq and sa.mapq != 255 and max(west, east) >= args.min_flank
            rank = (placed, -abs(gap), sa.mapq if sa.mapq != 255 else -1, core)
            if best is None or rank > best[0]:
                evidence = Evidence(anchor, "split_read", "flank_placement" if placed else "ambiguous",
                                    source.name, sa.chrom, cut, side, sa.mapq, gap)
                best = (rank, evidence)
    return best[1] if best else None


def insertion_sequences(read, min_anchor, min_clip):
    ops = Segment(read.reference_name, read.reference_start, "+", read.mapping_quality, read.cigarstring).ops
    ref, query, before = read.reference_start, 0, 0
    after = sum(n for n, op in ops if op in MATCH)
    for n, op in ops:
        if op == "I" and n >= min_clip and min(before, after) >= min_anchor:
            yield ref, read.query_sequence[query:query + n]
        if op in MATCH:
            before += n
            after -= n
        if op in REF:
            ref += n
        if op in MATCH or op in "IS":
            query += n


def scan_anchors(bam, partners, sources, args, query_path):
    evidence, queries, counts = [], [], {}
    with query_path.open("w") as fasta:
        def queue(anchor, seq):
            if seq and len(seq) >= args.min_clip:
                fasta.write(f">{len(queries)}\n{seq}\n")
                queries.append((anchor, len(seq)))

        for region in partners:
            start = max(0, region.start - args.padding)
            end = min(bam.get_reference_length(region.chrom), region.end + args.padding)
            names = set()
            for read in bam.fetch(region.chrom, start, end):
                if (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail
                        or not args.min_mapq <= read.mapping_quality < 255 or not read.cigarstring):
                    continue
                strand = "-" if read.is_reverse else "+"
                segment = Segment(region.chrom, read.reference_start, strand, read.mapping_quality, read.cigarstring)
                if segment.layout()[3] < args.min_anchor:
                    continue
                names.add(read.query_name)
                def anchor(cut, side):
                    return Anchor(region.name, read.query_name, region.chrom, cut, side, read.mapping_quality)

                sa = []
                if read.has_tag("SA"):
                    for entry in read.get_tag("SA").rstrip(";").split(";"):
                        chrom, pos, orientation, cigar, mapq, _ = entry.split(",")
                        sa.append(Segment(chrom, int(pos) - 1, orientation, int(mapq), cigar))
                for side, ops in (("L", segment.ops), ("R", segment.ops[::-1])):
                    cut = segment.cut(side)
                    if not start <= cut <= end:
                        continue
                    clip = next((n for n, op in ops[:2] if op == "S"), 0) if ops[0][1] in CLIP else 0
                    if clip < args.min_clip and not sa:
                        continue
                    support = split_evidence(anchor(cut, side), segment, sa, sources, args)
                    if support:
                        evidence.append(support)
                    elif clip >= args.min_clip and read.query_sequence:
                        n = min(clip, args.max_clip)
                        seq = read.query_sequence[:clip][-n:] if side == "L" else read.query_sequence[-clip:][:n]
                        queue(anchor(cut, side), seq)
                if read.query_sequence:
                    for cut, seq in insertion_sequences(read, args.min_anchor, args.min_clip):
                        if start <= cut <= end:
                            queue(anchor(cut, "I"), seq)
            counts[region.name] = (start, end, len(names))
    return evidence, queries, counts


def blast_evidence(reference, sources, queries, tmp, args, *, filtered):
    if not queries:
        return []
    subject = tmp / "dux4.fa"
    db = tmp / "dux4_db"
    if not subject.exists():
        with subject.open("w") as handle:
            for region in sources:
                handle.write(f">{region.name}\n{reference.fetch(region.chrom, region.start, region.end)}\n")
        subprocess.run(["makeblastdb", "-dbtype", "nucl", "-in", str(subject), "-out", str(db)], check=True)
    hits = tmp / "hits.tsv"
    subprocess.run([
        "blastn", "-task", "blastn", "-db", str(db), "-query", str(tmp / "queries.fa"),
        "-strand", "plus", "-dust", "no", "-evalue", "1e-8",
        "-max_hsps", "100", "-num_threads", str(args.threads), "-out", str(hits),
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore",
        *(["-perc_identity", str(args.min_identity)] if filtered else []),
    ], check=True)
    by_name = {r.name: r for r in sources}
    best = {}
    with hits.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            qid, source_name, ident, length, qs, qe, ss, se, score = row
            idx = int(qid)
            anchor, qlen = queries[idx]
            length, qs, qe, ss, se = map(int, (length, qs, qe, ss, se))
            if filtered and min(length, qe - qs + 1, abs(se - ss) + 1) < args.min_clip:
                continue
            if anchor.partner_side == "R":
                gap, edge, subject_side = qs - 1, ss, "L" if ss <= se else "R"
            elif anchor.partner_side == "L":
                gap, edge, subject_side = qlen - qe, se, "R" if ss <= se else "L"
            else:
                gap, edge, subject_side = min(qs - 1, qlen - qe), None, "."
            if filtered and gap > args.max_query_gap:
                continue
            rank = (float(score), -gap, length)
            if idx not in best or rank > best[idx][0]:
                region = by_name[source_name]
                cut = region.start + edge - (subject_side == "L") if edge is not None else None
                evidence = Evidence(anchor, "cigar_insertion_blast" if edge is None else "softclip_blast",
                                    "ambiguous", source_name, region.chrom, cut, subject_side,
                                    query_gap=gap, blast_identity=float(ident), blast_length=length)
                best[idx] = (rank, evidence)
    return [hit for _, hit in best.values()]


def all_softclip_hits(bam, reference, sources, tmp, args, *, batch_bases=2_000_000):
    """Search full soft-clipped ends genome-wide without the partner screen's gates."""
    queries = []
    bases = total = 0
    with (tmp / "queries.fa").open("w") as fasta:
        for read in bam.fetch(until_eof=True):
            ops = read.cigartuples
            if read.is_unmapped or not ops:
                continue
            # CIGAR codes: S=4, H=5. A terminal S can sit inside a hard clip.
            left = ops[1] if ops[0][0] == 5 and len(ops) > 1 else ops[0]
            right = ops[-2] if ops[-1][0] == 5 and len(ops) > 1 else ops[-1]
            if left[0] != 4 and right[0] != 4:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            for side, (op, n) in (("L", left), ("R", right)):
                if op != 4:
                    continue
                clip = seq[:n] if side == "L" else seq[-n:]
                cut = read.reference_start if side == "L" else read.reference_end
                anchor = Anchor(".", read.query_name, read.reference_name, cut, side, read.mapping_quality)
                fasta.write(f">{len(queries)}\n{clip}\n")
                queries.append((anchor, len(clip)))
                bases += len(clip)
                total += 1
                if bases >= batch_bases:
                    fasta.flush()
                    yield from blast_evidence(reference, sources, queries, tmp, args, filtered=False)
                    fasta.seek(0)
                    fasta.truncate()
                    queries.clear()
                    bases = 0
        fasta.flush()
        yield from blast_evidence(reference, sources, queries, tmp, args, filtered=False)
    print(f"Genome-wide soft-clip BLAST queries: {total}", flush=True)


def write_blast_report(path, hits):
    anchor_fields = [f.name for f in fields(Anchor)]
    evidence_fields = [f.name for f in fields(Evidence) if f.name != "anchor"]
    write_table(path, (*anchor_fields, *evidence_fields),
                ((*(getattr(h.anchor, name) for name in anchor_fields),
                  *(getattr(h, name) for name in evidence_fields)) for h in hits))


def write_table(path, header, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_reports(partners, evidence, counts, args):
    unique = set(evidence)
    grouped = defaultdict(list)
    for hit in unique:
        a = hit.anchor
        grouped[a.partner, a.partner_chrom, a.partner_side].append(hit)
    junction_rows, read_rows = [], []
    ends_by_read = defaultdict(set)
    num_clusters, num_review = defaultdict(int), defaultdict(int)
    for (partner, chrom, side), hits in sorted(grouped.items()):
        clusters = []
        for hit in sorted(hits, key=lambda h: (h.anchor.partner_cut, h.anchor.read, h.method, repr(h))):
            if not clusters or hit.anchor.partner_cut - clusters[-1][0].anchor.partner_cut > args.cluster_bp:
                clusters.append([])
            clusters[-1].append(hit)
        for cluster in clusters:
            jid = f"J{len(junction_rows) + 1:04d}"
            reads = {h.anchor.read for h in cluster}
            cuts = sorted({(h.anchor.read, h.anchor.partner_cut) for h in cluster})
            per_read = defaultdict(list)
            for name, cut in cuts:
                per_read[name].append(cut)
            representative = median_low([median_low(values) for values in per_read.values()])
            placed = sorted({h.source_chrom for h in cluster if h.source_status == "flank_placement"})
            if side == "I":
                geometry = "insertion_like"
            elif len(placed) > 1:
                geometry = "conflicting_flank_placements"
            elif not placed:
                geometry = "repeat_associated_adjacency"
            elif placed[0] == chrom:
                geometry = "intrachromosomal_adjacency"
            else:
                geometry = "interchromosomal_adjacency"
            num_clusters[partner] += 1
            review = len(reads) >= args.min_reads
            num_review[partner] += review
            junction_rows.append((jid, partner, chrom, representative, side, geometry, ",".join(placed) or ".",
                                  len(reads), len({h.anchor.read for h in cluster if h.method == "split_read"}),
                                  len({h.anchor.read for h in cluster if h.method != "split_read"}),
                                  "review" if review else "low_support"))
            for h in cluster:
                if side in {"L", "R"}:
                    ends_by_read[partner, h.anchor.read].add((side, jid, h.anchor.partner_cut))
                read_rows.append((jid, *(getattr(h.anchor, f.name) for f in fields(Anchor)),
                                  *(getattr(h, f.name) for f in fields(Evidence) if f.name != "anchor")))
    links = defaultdict(lambda: defaultdict(set))
    for (_, name), ends in ends_by_read.items():
        for side, left, cut in ends:
            if side != "L":
                continue
            for other_side, right, other_cut in ends:
                if other_side == "R" and abs(cut - other_cut) <= args.cluster_bp:
                    links[left][right].add(name)
                    links[right][left].add(name)
    write_table(args.junctions, ("junction_id", "partner", "chrom", "cut0", "side", "geometry", "flank_placements",
                                "reads", "split_reads", "blast_reads", "status", "nearby_opposite_ends"),
                ((*row, ",".join(f"{jid}:{len(names)}" for jid, names in sorted(links[row[0]].items())) or ".")
                 for row in junction_rows))
    write_table(args.reads, ("junction_id", *(f.name for f in fields(Anchor)),
                             *(f.name for f in fields(Evidence) if f.name != "anchor")), read_rows)
    write_table(args.loci, ("partner", "chrom", "start0", "end0", "scan_start0", "scan_end0", "anchor_reads",
                            "candidate_clusters", "review_clusters"),
                ((r.name, r.chrom, r.start, r.end, *counts[r.name], num_clusters[r.name], num_review[r.name])
                 for r in partners))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("bam", "bam-index", "reference", "partners", "sources", "junctions", "reads", "loci",
                 "blast", "filtered-blast"):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--require-chm13", action="store_true")
    for name, default in (("threads", 4), ("padding", 10000), ("min-anchor", 200), ("min-mapq", 20),
                          ("min-clip", 100), ("max-clip", 3000), ("max-query-gap", 100),
                          ("source-flank", 5000), ("min-flank", 250), ("cluster-bp", 100), ("min-reads", 2)):
        p.add_argument(f"--{name}", type=int, default=default)
    p.add_argument("--min-identity", type=float, default=80.0)
    args = p.parse_args()
    if (args.max_clip < args.min_clip or args.min_clip < 1 or args.min_anchor < 1 or args.threads < 1
            or args.min_reads < 1 or not 0 <= args.min_mapq <= 254 or not 0 < args.min_identity <= 100
            or min(args.padding, args.max_query_gap, args.source_flank, args.cluster_bp) < 0 or args.min_flank < 1):
        p.error("Invalid scan thresholds")
    partners, sources = read_bed(args.partners), read_bed(args.sources)
    with pysam.AlignmentFile(args.bam, "rb", index_filename=args.bam_index, require_index=True) as bam, \
            pysam.FastaFile(args.reference) as reference, TemporaryDirectory(prefix="dux4-") as tmp:
        ref_lengths = dict(zip(reference.references, reference.lengths))
        bam_lengths = dict(zip(bam.references, bam.lengths))
        if args.require_chm13 and any(ref_lengths.get(c) != length for c, length in CHM13_LENGTHS.items()):
            p.error("Expected the chr-prefixed T2T-CHM13v2.0 reference FASTA")
        for r in partners + sources:
            if ref_lengths.get(r.chrom, 0) < r.end or ref_lengths.get(r.chrom) != bam_lengths.get(r.chrom):
                p.error(f"Reference/BAM mismatch or out-of-bounds interval: {r.name} ({r.chrom})")
        tmp = Path(tmp)
        evidence, queries, counts = scan_anchors(bam, partners, sources, args, tmp / "queries.fa")
        print(f"Scanned {len(partners)} loci; {len(evidence)} split observations; {len(queries)} filtered BLAST queries "
              "(unresolved soft clips and insertions)", flush=True)
        filtered_hits = blast_evidence(reference, sources, queries, tmp, args, filtered=True)
        evidence.extend(filtered_hits)
        write_blast_report(args.filtered_blast,
                           (hit for hit in dict.fromkeys(filtered_hits) if hit.method == "softclip_blast"))
        with pysam.AlignmentFile(args.bam, "rb") as all_reads:
            write_blast_report(args.blast, all_softclip_hits(all_reads, reference, sources, tmp, args))
    write_reports(partners, evidence, counts, args)


if __name__ == "__main__":
    main()
