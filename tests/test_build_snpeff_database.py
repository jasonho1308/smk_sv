"""SnpEff build checks. Set SNPEFF_EXECUTABLE for a real build/annotation smoke test."""

import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("build_snpeff", ROOT / "workflow/scripts/build_snpeff_database.py")
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


@pytest.fixture
def inputs(tmp_path):
    sequence = "ATG" + "AAA" * 9 + "TAA"
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chr1\n{sequence}\n")
    fai = tmp_path / "genome.fa.fai"
    fai.write_text(f"chr1\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n")
    genes = tmp_path / "genes.gff3"
    genes.write_text(
        "##gff-version 3\n"
        "chr1\tfixture\tgene\t1\t33\t.\t+\t.\tID=g1;gene_name=GENE1;gene=GENE1;gene_biotype=protein_coding\n"
        "chr1\tfixture\ttranscript\t1\t33\t.\t+\t.\tID=tx1;Parent=g1;gbkey=mRNA;transcript_biotype=mRNA\n"
        "chr1\tfixture\texon\t1\t33\t.\t+\t.\tID=ex1;Parent=tx1\n"
        "chr1\tfixture\tCDS\t1\t33\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )
    return dict(fasta=fasta, fai=fai, genes=genes,
                destination=tmp_path / "cache/Toy.1", database="Toy.1", mem_mb=512)


def write_predictor(command, **kwargs):
    data = Path(command[command.index("-dataDir") + 1])
    (data / command[-1] / "snpEffectPredictor.bin").write_bytes(b"fixture database")
    return subprocess.CompletedProcess(command, 0)


def test_build_layout_registration_and_explicit_skipped_checks(inputs):
    with patch.object(builder.subprocess, "run", side_effect=write_predictor) as run:
        builder.build_database(**inputs)
    destination = inputs["destination"]
    assert (destination / "sequences.fa").resolve() == inputs["fasta"]
    assert (destination / "genes.gff").resolve() == inputs["genes"]
    config = (destination / "snpEff.config").read_text()
    assert "Toy.1.genome" in config
    assert "Toy.1.M.codonTable : Vertebrate_Mitochondrial" in config
    assert "-gff3" in run.call_args.args[0]
    assert "-noCheckCds" in run.call_args.args[0]
    assert "-noCheckProtein" in run.call_args.args[0]
    manifest = json.loads((destination / "build_manifest.json").read_text())
    assert manifest["features"] == 4
    assert manifest["skipped_sequence_checks"] == ["CDS", "protein"]
    assert len(manifest["genes_sha256"]) == 64


def test_supplied_sequences_enable_verification(inputs, tmp_path):
    cds = tmp_path / "cds.fa"
    protein = tmp_path / "protein.fa"
    cds.write_text(">tx1\n" + "ATG" + "AAA" * 9 + "TAA\n")
    protein.write_text(">tx1\nMKKKKKKKKK*\n")
    with patch.object(builder.subprocess, "run", side_effect=write_predictor) as run:
        builder.build_database(**inputs, cds=cds, protein=protein)
    command = run.call_args.args[0]
    assert "-noCheckCds" not in command and "-noCheckProtein" not in command
    manifest = json.loads((inputs["destination"] / "build_manifest.json").read_text())
    assert manifest["skipped_sequence_checks"] == []
    assert (inputs["destination"] / "cds.fa").resolve() == cds


@pytest.mark.parametrize("bad_row", [
    "1\tfixture\tgene\t1\t33\t.\t+\t.\tID=g1\n",
    "chr1\tfixture\tgene\t1\t34\t.\t+\t.\tID=g1\n",
    "chr1\tfixture\tgene\t0\t33\t.\t+\t.\tID=g1\n",
    "##FASTA\n>chr1\nACGT\n",
    "# No annotations\n",
])
def test_mismatched_or_invalid_annotations_do_not_publish_database(inputs, bad_row):
    inputs["genes"].write_text(bad_row)
    with patch.object(builder.subprocess, "run") as run:
        with pytest.raises(ValueError):
            builder.build_database(**inputs)
    run.assert_not_called()
    assert not inputs["destination"].exists()


def test_compressed_gtf_layout_is_supported(inputs, tmp_path):
    genes = tmp_path / "genes.gtf.gz"
    with gzip.open(genes, "wt") as handle:
        handle.write('chr1\ttest\texon\t1\t33\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n')
    inputs["genes"] = genes
    with patch.object(builder.subprocess, "run", side_effect=write_predictor) as run:
        builder.build_database(**inputs, format="gtf22")
    assert "-gtf22" in run.call_args.args[0]
    assert (inputs["destination"] / "genes.gtf.gz").resolve() == genes


def test_failed_build_does_not_publish_partial_database(inputs):
    with patch.object(builder.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "snpEff")):
        with pytest.raises(subprocess.CalledProcessError):
            builder.build_database(**inputs)
    assert not inputs["destination"].exists()
    assert not list(inputs["destination"].parent.glob(".Toy.1.*"))


def test_success_without_predictor_is_rejected(inputs):
    with patch.object(builder.subprocess, "run"):
        with pytest.raises(ValueError, match="snpEffectPredictor.bin"):
            builder.build_database(**inputs)
    assert not inputs["destination"].exists()


def test_existing_database_is_not_overwritten(inputs):
    inputs["destination"].mkdir(parents=True)
    predictor = inputs["destination"] / "snpEffectPredictor.bin"
    predictor.write_bytes(b"existing")
    with pytest.raises(FileExistsError):
        builder.build_database(**inputs)
    assert predictor.read_bytes() == b"existing"


@pytest.mark.skipif(not os.environ.get("SNPEFF_EXECUTABLE"), reason="Set SNPEFF_EXECUTABLE for real SnpEff smoke test")
def test_real_snpeff_rejects_mismatched_cds(inputs, tmp_path):
    cds = tmp_path / "incorrect_cds.fa"
    cds.write_text(">tx1\n" + "CCC" * 11 + "\n")
    with pytest.raises((subprocess.CalledProcessError, ValueError)):
        builder.build_database(**inputs, cds=cds, executable=os.environ["SNPEFF_EXECUTABLE"])
    assert not inputs["destination"].exists()


@pytest.mark.skipif(not os.environ.get("SNPEFF_EXECUTABLE"), reason="Set SNPEFF_EXECUTABLE for real SnpEff smoke test")
def test_real_snpeff_build_and_annotation(inputs, tmp_path):
    executable = os.environ["SNPEFF_EXECUTABLE"]
    cds = tmp_path / "cds.fa"
    protein = tmp_path / "protein.fa"
    cds.write_text(">tx1\n" + "ATG" + "AAA" * 9 + "TAA\n")
    protein.write_text(">tx1\nMKKKKKKKKK*\n")
    builder.build_database(**inputs, cds=cds, protein=protein, executable=executable)
    build_log = (inputs["destination"].parent / "Toy.1.build.log").read_text()
    assert "CDS check:" in build_log
    assert "Protein check:" in build_log
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=33>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t4\t.\tA\tG\t.\tPASS\t.\n"
    )
    destination = inputs["destination"]
    result = subprocess.run([
        executable, "-Xmx512M", "-c", str(destination / "snpEff.config"),
        "-dataDir", str(destination.parent), "-nodownload", "-noStats",
        "Toy.1", str(vcf),
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "missense_variant" in result.stdout
    assert "GENE1" in result.stdout
