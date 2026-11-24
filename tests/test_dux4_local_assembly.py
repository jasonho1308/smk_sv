"""Synthetic DUX4 local assembly cases; tool-dependent checks skip if not on PATH."""

import csv
import gzip
from pathlib import Path
import random
import shutil
import subprocess
import sys

import pytest

HELPER = Path(__file__).resolve().parents[1] / "workflow/scripts/dux4_local_assembly.py"


def run(*args, **kwargs):
    result = subprocess.run(args, capture_output=True, text=True, timeout=60, **kwargs)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


def metrics(path):
    return dict(row.split("\t") for row in path.read_text().splitlines()[1:])


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def require(*tools):
    for tool in tools:
        if not shutil.which(tool):
            pytest.skip(f"requires {tool}")


def intervals(tmp, *, source="chr10\t80\t400\tD4Z4_10q26", partner="chr12\t80\t400\tETV6"):
    sources, partners = tmp / "sources.bed", tmp / "partners.bed"
    sources.write_text(source + "\n")
    partners.write_text(partner + "\n")
    return sources, partners


def sam_row(name, flag, chrom, pos, cigar, seq="*", qual="*", mapq=60):
    return f"{name}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n"


def sam_to_bam(tmp, body):
    sam, bam, index = tmp / "in.sam", tmp / "in.bam", tmp / "custom.csi"
    header = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(f"@SQ\tSN:{chrom}\tLN:1000\n" for chrom in ("chr10", "chr12"))
    sam.write_text(header + "".join(body))
    run("samtools", "sort", "-o", bam, sam)
    run("samtools", "index", "-c", bam, index)
    return bam, index


def test_extract_recruits_non_igh_supplementary_and_recovers_full_reads(tmp_path):
    require("samtools")
    seq, qual = "AACCCGGT" * 15, "!" * 40 + "I" * 80
    bam, index = sam_to_bam(tmp_path, [
        sam_row("rescue", 0x810, "chr10", 100, "60H60M", "A" * 60),
        sam_row("rescue", 0x10, "chr12", 700, "12S108M", reverse_complement(seq), qual[::-1]),
        sam_row("duplicate", 0, "chr12", 100, "60M60S", seq),
        sam_row("duplicate", 0x800, "chr10", 101, "60H60M", "A" * 60),
        sam_row("missing", 0, "chr10", 200, "60M60H", "C" * 60),
        sam_row("not_seeded", 0x100, "chr12", 200, "120M", seq),
        sam_row("not_seeded", 0, "chr12", 500, "120M", seq),
    ])
    fai = tmp_path / "ref.fai"
    fai.write_text("chr10\t1000\nchr12\t1000\n")
    sources, partners = intervals(tmp_path)
    fastq, stats = tmp_path / "target.fq.gz", tmp_path / "stats.tsv"
    run(sys.executable, HELPER, "extract", "--bam", bam, "--index", index, "--threads", "2",
        "--sources", sources, "--partners", partners, "--fai", fai, "--partner-padding", "0",
        "--fastq", fastq, "--stats", stats)
    with gzip.open(fastq, "rt") as handle:
        assert handle.read().splitlines() == ["@duplicate", seq, "+", "!" * 120, "@rescue", seq, "+", qual]
    assert metrics(stats) == {"selected_read_names": "3", "recovered_reads": "2",
                              "missing_full_sequence": "1", "missing_qualities": "1"}


@pytest.mark.parametrize("n, status", [(0, "no_reads"), (2, "insufficient_reads")])
def test_no_assembly_writes_empty_fasta(n, status, tmp_path):
    fastq, contigs, status_file = tmp_path / "test.fastq.gz", tmp_path / "contigs.fa", tmp_path / "status.tsv"
    with gzip.open(fastq, "wt") as handle:
        for i in range(n):
            handle.write(f"@read{i}\n" + "A" * 1000 + "\n+\n" + "I" * 1000 + "\n")
    run(sys.executable, HELPER, "assemble", "--fastq", fastq, "--contigs", contigs, "--status", status_file,
        "--platform", "hifi", "--ont-mode", "nano-hq", "--threads", "1", "--min-reads", "3", "--min-read-length", "1000")
    assert contigs.read_text() == ""
    assert metrics(status_file)["status"] == status


@pytest.fixture
def toy_ref(tmp_path):
    require("samtools")
    rng = random.Random(400)
    sequences = {chrom: "".join(rng.choices("ACGT", k=30000)) for chrom in ("chr10", "chr12")}
    ref = tmp_path / "ref.fa"
    ref.write_text("".join(f">{chrom}\n{seq}\n" for chrom, seq in sequences.items()))
    run("samtools", "faidx", ref)
    sources, partners = intervals(tmp_path, source="chr10\t5000\t18500\tD4Z4_10q26", partner="chr12\t16000\t30000\tETV6")
    return sequences, ref, sources, partners


def align(tmp, ref, sources, partners, contigs, preset="asm5"):
    bam, bai, summary, bed = (tmp / name for name in ("contigs.bam", "contigs.bam.bai", "breakpoints.tsv", "breakpoints.bed"))
    run(sys.executable, HELPER, "align", "--ref", ref, "--fai", str(ref) + ".fai", "--contigs", contigs,
        "--bam", bam, "--index", bai, "--summary", summary, "--bed", bed, "--min-anchor-bases", "100",
        "--min-partner-mapq", "20", "--max-query-gap", "1000", "--max-query-overlap", "100",
        "--threads", "2", "--preset", preset, "--sources", sources, "--partners", partners, "--partner-padding", "0")
    run("samtools", "quickcheck", bam)
    with summary.open() as handle:
        return bam, bai, list(csv.DictReader(handle, delimiter="\t")), bed.read_text()


def test_empty_fasta_yields_indexed_header_only_bam(tmp_path, toy_ref):
    _, ref, sources, partners = toy_ref
    contigs = tmp_path / "empty.fa"
    contigs.touch()
    bam, bai, rows, bed = align(tmp_path, ref, sources, partners, contigs)
    assert bai.stat().st_size > 0
    assert run("samtools", "view", "-c", bam, "chr10:1-30000").strip() == "0"
    assert rows == [] and bed == ""


def test_split_alignment_reverse_anchor_and_cut_coordinates(tmp_path, toy_ref):
    require("minimap2")
    seq, ref, sources, partners = toy_ref
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(f">joined\n{seq['chr10'][5000:18000]}{reverse_complement(seq['chr12'][17000:30000])}\n"
                       f">native\n{seq['chr12']}\n")
    bam, _, rows, bed = align(tmp_path, ref, sources, partners, contigs)
    assert len(rows) == 1
    row = rows[0]
    assert (row["contig"], row["source_locus"], row["partner_locus"], row["partner_strand"]) == (
        "joined", "D4Z4_10q26", "ETV6", "-")
    assert abs(int(row["source_cut"]) - 18000) < 10
    assert abs(int(row["partner_cut"]) - 30000) < 10
    assert abs(int(row["partner_qstart"]) - 13000) < 10
    assert "split_1|joined|ETV6" in bed
    assert "cs:Z:" in run("samtools", "view", bam)


@pytest.mark.parametrize("platform, tool, preset", [("hifi", "hifiasm", "asm5"), ("ont", "flye", "asm20")])
def test_assemble_contig_with_non_igh_partner(platform, tool, preset, tmp_path, toy_ref):
    require("minimap2", tool)
    seq, ref, sources, partners = toy_ref
    chimera = seq["chr10"][5000:18000] + seq["chr12"][17000:30000]
    fastq, contigs, status_file = tmp_path / "test.fastq.gz", tmp_path / "contigs.fa", tmp_path / "status.tsv"
    with gzip.open(fastq, "wt") as handle:
        for i in range(100):
            offset = (i * 337) % 14600
            read = chimera[offset:offset + 11400]
            handle.write(f"@read{i}\n{read}\n+\n" + "I" * len(read) + "\n")
    run(sys.executable, HELPER, "assemble", "--fastq", fastq, "--contigs", contigs, "--status", status_file,
        "--platform", platform, "--ont-mode", "nano-hq", "--threads", "2", "--min-reads", "3", "--min-read-length", "1000")
    assert metrics(status_file)["status"] == "assembled"
    _, _, rows, _ = align(tmp_path, ref, sources, partners, contigs, preset)
    assert any(row["source_locus"] == "D4Z4_10q26" and row["partner_locus"] == "ETV6" for row in rows)
