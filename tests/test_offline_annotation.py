"""Offline annotation preprocessing and R filtering regression tests."""

import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("prepare_vep", ROOT / "workflow/scripts/prepare_vep_annotation.py")
preparer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparer)


@pytest.fixture
def annotation(tmp_path):
    genes = tmp_path / "source genes.gff3"
    genes.write_text(
        "##gff-version 3\n"
        "chr2\tfixture\tgene\t20\t30\t.\t+\t.\tID=g2\n"
        "chr1\tfixture\tgene\t10\t20\t.\t+\t.\tID=g1\n"
        "chr1\tfixture\texon\t2\t5\t.\t+\t.\tParent=t1\n"
    )
    fai = tmp_path / "genome.fa.fai"
    fai.write_text("chr1\t100\t6\t100\t101\nchr2\t100\t113\t100\t101\n")
    return genes, fai, tmp_path / "prepared annotation.gff3.gz"


@pytest.mark.skipif(not (shutil.which("bgzip") and shutil.which("tabix")), reason="HTSlib tools required")
def test_sort_bgzip_tabix_and_query(annotation):
    genes, fai, destination = annotation
    preparer.prepare_annotation(genes, fai, destination)
    with gzip.open(destination, "rt") as handle:
        lines = handle.readlines()
    assert [line.split("\t")[0:4:3] for line in lines] == [["chr1", "2"], ["chr1", "10"], ["chr2", "20"]]
    assert Path(f"{destination}.tbi").stat().st_size > 0
    query = subprocess.run(["tabix", str(destination), "chr2:20-30"], capture_output=True, text=True, check=True)
    assert "ID=g2" in query.stdout


@pytest.mark.parametrize("bad", [
    "1\ttest\tgene\t1\t10\t.\t+\t.\tID=g1\n",
    "chr1\ttest\tgene\t1\t101\t.\t+\t.\tID=g1\n",
    "##FASTA\n>chr1\nACGT\n",
    "# empty annotation\n",
])
def test_invalid_annotations_do_not_publish(annotation, bad):
    genes, fai, destination = annotation
    genes.write_text(bad)
    with pytest.raises(ValueError):
        preparer.prepare_annotation(genes, fai, destination)
    assert not destination.exists()
    assert not Path(f"{destination}.tbi").exists()


def test_compression_failure_preserves_existing_files(annotation):
    genes, fai, destination = annotation
    destination.write_bytes(b"existing annotation")
    index = Path(f"{destination}.tbi")
    index.write_bytes(b"existing index")
    with patch.object(preparer.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "bgzip")):
        with pytest.raises(subprocess.CalledProcessError):
            preparer.prepare_annotation(genes, fai, destination)
    assert destination.read_bytes() == b"existing annotation"
    assert index.read_bytes() == b"existing index"


@pytest.mark.skipif(not shutil.which("Rscript"), reason="Rscript required")
def test_r_filtering_and_command_failures():
    result = subprocess.run(["Rscript", str(ROOT / "tests/test_annotation_filters.R"), str(ROOT)],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "tests passed" in result.stdout
