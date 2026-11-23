"""Configuration and dry-run regression checks; no reference downloads required.

Run with: python -m pytest tests/test_t2t_configuration.py
Requires pytest, pyyaml, jsonschema, snakemake, peppy, and eido.
The configured Snakemake fork also requires its conda/container deployment
plugins, as declared in pixi.toml.
"""

from copy import deepcopy
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from unittest.mock import patch

import jsonschema
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


def template():
    return yaml.safe_load((ROOT / "config/.config.yaml").read_text())


def test_t2t_template_matches_schema_and_resources():
    config = template()
    schema = json.loads((ROOT / "workflow/schemas/config.schema.json").read_text())
    jsonschema.validate(config, schema)
    assert config["annotation"]["genome"] == "CHM13"
    assert "T2T-CHM13v2.0" in config["inputs"]["fasta"]
    assert "T2T-CHM13v2.0" in config["inputs"]["index_minimap2"]
    assert config["caller_settings"]["tandem_repeats"] == config["caller_settings"]["nvtr"]
    assert config["caller_settings"]["nanosv"]["bed"] is None
    downloads = config["resources"]["downloads"]
    assert set(downloads["reference"]["assemblies"]) == {"CHM13"}
    assert downloads["reference"]["assemblies"]["CHM13"]["checksum"] == (
        "md5:06caded0055ec647a69acb708e13beff"
    )
    assert "chm13.bed" in downloads["tandem_repeats"]["assemblies"]["CHM13"]["url"]
    assert downloads["clair3_model"]["archive_root"] == "models/r941_prom_sup_g5014"
    assert "pytorch" in config["caller_settings"]["clair3"]["model"]


def test_annotation_versions_match_executable_environments():
    config = template()
    for name in ("vep", "awscli"):
        deps = yaml.safe_load((ROOT / f"workflow/envs/{name}.yaml").read_text())["dependencies"]
        release = next(dep.split("=")[1].split(".")[0] for dep in deps if dep.startswith("ensembl-vep="))
        assert int(release) == config["annotation"]["vep"]["version"]
    assert config["annotation"]["snpeff"]["version"] == "v2.0"
    deps = yaml.safe_load((ROOT / "workflow/envs/snpeff.yaml").read_text())["dependencies"]
    assert "openjdk>=21" in deps
    assert "snpeff=5.4.0c" in deps and "snpsift=5.4.0c" in deps


def test_pytorch_archive_installs_selected_model(tmp_path):
    """The bundled archive contains models/, not a model at its root."""
    import tarfile

    spec = importlib.util.spec_from_file_location("install_archive", ROOT / "workflow/scripts/install_archive.py")
    installer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(installer)
    source = tmp_path / "source/models/r941_prom_sup_g5014"
    source.mkdir(parents=True)
    for name in ("pileup.pt", "full_alignment.pt"):
        (source / name).write_bytes(b"model fixture")
    archive = tmp_path / "models.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(source.parent, arcname="models")
    destination = tmp_path / "installed"
    with patch.object(installer, "download", return_value=archive):
        installer.install_archive("https://example.org/models.tar.gz", destination, "tar", "models/r941_prom_sup_g5014")
    assert sorted(path.name for path in destination.iterdir()) == ["full_alignment.pt", "pileup.pt"]


@pytest.fixture
def project(tmp_path):
    shutil.copytree(ROOT / "workflow", tmp_path / "workflow", ignore=shutil.ignore_patterns("__pycache__"))
    config = deepcopy(template())
    config["run"]["dir"] = str(tmp_path / "results_t2t")
    config["inputs"]["fasta"] = str(tmp_path / "resources/chm13v2.0.fa")
    config["inputs"]["index_minimap2"] = str(tmp_path / "resources/chm13v2.0.mmi")
    config["inputs"]["dir_data"] = str(tmp_path / "data")
    config["annotation"]["snpeff"]["cache"] = str(tmp_path / "resources/snpeff")
    config["annotation"]["snpeff"]["build"]["genes"] = str(tmp_path / "resources/genes.gff3")
    config["annotation"]["vep"]["annotation"] = str(tmp_path / "resources/vep.gff3")
    config["annotation"]["annotsv"]["cache"] = str(tmp_path / "resources/annotsv")
    settings = config["caller_settings"]
    settings["tandem_repeats"] = settings["nvtr"] = str(tmp_path / "resources/chm13.bed")
    settings["clair3"]["model"] = str(tmp_path / "resources/pytorch/r941_prom_sup_g5014")
    settings["svision"]["model"] = str(tmp_path / "resources/svision/model.ckpt")
    settings["nanosv"]["config"] = str(tmp_path / "resources/nanosv/config.ini")
    (tmp_path / "config/pep").mkdir(parents=True)
    (tmp_path / "config/config.yaml").write_text(yaml.safe_dump(config))
    (tmp_path / "config/pep/config.yaml").write_text("pep_version: 2.1.0\nsample_table: samples.csv\n")
    (tmp_path / "config/pep/samples.csv").write_text("sample_name\nexample\n")
    (tmp_path / "data").mkdir()
    (tmp_path / "data/example.fq.gz").touch()
    return tmp_path, config


def run_snakemake(project, *args, cores=1):
    path, _ = project
    executable = shutil.which("snakemake")
    assert executable, "Install Snakemake to run the workflow regression checks"
    return subprocess.run(
        [executable, "--snakefile", "workflow/Snakefile", "--workflow-profile", "none", "--cores", str(cores), *args],
        cwd=path, capture_output=True, text=True, timeout=90,
    )


def test_resources_dry_run_uses_t2t_and_pytorch(project):
    proc = run_snakemake(project, "resources", "--dry-run", "--printshellcmds")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "download_reference" in output
    assert "chm13v2.0.fa" in output
    assert "chm13.bed" in output
    assert "prepare_vep_annotation" in output
    assert "vep.gff3.sorted.gz" in output
    assert "download_vep_cache" not in output
    assert "CHM13.v2.0" in output
    assert "build_snpeff_cache" in output
    assert "download_snpeff_genes" in output
    assert "pytorch/r941_prom_sup_g5014" in output
    assert "GRCh38.primary_assembly" not in output
    assert "download_dragen_dux4_hash_table" not in output


def test_full_analysis_dry_run_builds_t2t_dag(project):
    proc = run_snakemake(project, "all", "--dry-run", "--printshellcmds")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "rule minimap2:" in output
    assert "rule clair3:" in output
    assert "rule severus:" in output
    assert "rule survivor_final:" in output
    assert "run_dragen_dux4" not in output
    assert "download_dragen_dux4_hash_table" not in output
    assert "CHM13" in output
    assert "--gff" in output
    # VEP's literal --offline flag implies a cache; custom GFF mode uses neither.
    assert "--offline" not in output
    assert "--cache_version" not in output
    assert "--filter_common" not in output
    assert "GRCh38.primary_assembly" not in output


def test_dux4_target_uses_indexed_t2t_alignment_without_consensus(project):
    proc = run_snakemake(project, "dux4", "--dry-run", "--printshellcmds", cores=16)
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "rule scan_dux4:" in output
    assert "chm13v2.0.partners.bed" in output and "chm13v2.0.sources.bed" in output
    assert "example.sorted.bam.csi" in output
    assert "example.junctions.tsv" in output
    assert "--blast dux4/example/example.blast.tsv" in output
    assert "--filtered-blast dux4/example/example.blast_filtered.tsv" in output
    assert "--threads 16 > logs/example/dux4.log" in output
    assert "rule survivor_final:" not in output
    assert "run_dragen_dux4" not in output
    assert "download_dragen_dux4_hash_table" not in output


def configure_dragen(project, **settings):
    path, config = project
    hash_table = path / "DRAGEN hash table"
    hash_table.mkdir()
    config["caller_settings"]["dragen_dux4"] = {"hash_table": str(hash_table), **settings}
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))


def test_dragen_dux4_external_bam_requires_matching_hash_table():
    config = template()
    config["caller_settings"]["dragen_dux4"] = {"bam_pattern": "/data/{sample}.bam"}
    schema = json.loads((ROOT / "workflow/schemas/config.schema.json").read_text())
    with pytest.raises(jsonschema.ValidationError, match="hash_table.*dependency"):
        jsonschema.validate(config, schema)


def test_dragen_dux4_uses_pipeline_bam_and_downloads_hash_table_by_default(project):
    proc = run_snakemake(project, "dux4_dragen", "--dry-run", "--printshellcmds")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "rule minimap2:" in output
    assert "rule download_dragen_dux4_hash_table:" in output
    assert "--hash-table references/dragen/chm13_v2-linear-r4.0" in output
    assert "--bam minimap2/example/example.sorted.bam" in output
    assert "--skip-sanity-check" in output
    assert "rule survivor_final:" not in output and "rule scan_dux4:" not in output


@pytest.mark.skipif(not shutil.which("curl"), reason="requires curl")
@pytest.mark.parametrize("valid", [True, False])
def test_dragen_hash_table_download_and_cache(project, valid):
    import tarfile

    # The standalone target remains available even when no caller would use it.
    configure_dragen(project)
    path, config = project
    config["annotation"]["genome"] = "hg38"
    source = path / "source"
    source.mkdir()
    filenames = ["hash_table.cfg.bin", "reference.bin"] if valid else ["hash_table.cfg.bin"]
    for name in filenames:
        (source / name).write_bytes(b"fixture")
    archive = path / "dragen.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        handle.add(source, arcname=".")
    config["resources"]["downloads"]["dragen_dux4_hash_table"] = {"url": archive.as_uri()}
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))

    target = "download_dragen_dux4_hash_table"
    proc = run_snakemake(project, target)
    output = proc.stdout + proc.stderr
    installed = Path(config["run"]["dir"]) / "references/dragen/chm13_v2-linear-r4.0"
    if not valid:
        assert proc.returncode != 0
        assert "Required file(s) missing from archive: reference.bin" in output
        assert not installed.exists()
        return
    assert proc.returncode == 0, output
    assert all((installed / name).read_bytes() == b"fixture" for name in filenames)
    archive.unlink()
    cached = run_snakemake(project, target)
    assert cached.returncode == 0, cached.stdout + cached.stderr
    assert "Nothing to be done" in cached.stdout + cached.stderr


@pytest.mark.parametrize("mode, outcome, exit_code", [
    ("success", "vcf_written", 0),
    ("empty", "no_vcf", 0),
    ("reject", "failed", 12),
])
def test_dragen_dux4_external_bam_attempt(project, monkeypatch, mode, outcome, exit_code):
    path, config = project
    bam = path / "data/example tumor.bam"
    bam.touch()
    bypass = mode == "success"
    configure_dragen(project, bam_pattern=str(path / "data/{sample} tumor.bam"), skip_sanity_check=bypass)
    tools = path / "tools"
    tools.mkdir()
    stub = tools / "dragen"
    stub.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "import os, sys\n"
        "out = next(arg.split('=', 1)[1] for arg in sys.argv if arg.startswith('--output-dir='))\n"
        "mode = os.environ['DRAGEN_TEST_MODE']\n"
        "if mode != 'empty':\n"
        "    (Path(out) / 'example.dux4.vcf').write_text('##fileformat=VCFv4.2\\n')\n"
        "if mode == 'reject':\n"
        "    print('fixture rejects input')\n"
        "    sys.exit(12)\n"
    )
    stub.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("DRAGEN_TEST_MODE", mode)
    proc = run_snakemake(project, "dux4_dragen/example/attempt.json")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "rule minimap2:" not in output
    assert "download_dragen_dux4_hash_table" not in output
    run_dir = Path(config["run"]["dir"])
    attempt = run_dir / "dux4_dragen/example"
    vcfs = ["dragen/example.dux4.vcf"] if mode == "success" else []
    assert json.loads((attempt / "attempt.json").read_text()) == {
        "outcome": outcome, "exit_code": exit_code, "skip_sanity_check": bypass, "vcfs": vcfs,
    }
    assert sorted(str(vcf.relative_to(attempt)) for vcf in (attempt / "dragen").glob("*.vcf")) == vcfs
    log = (run_dir / "logs/example/dux4_dragen.log").read_text()
    assert "--enable-dux4-caller=true" in log
    assert f"--tumor-bam-input '{bam}'" in log
    assert ("--dux4-skip-santiy-check true" in log) == bypass
    if mode == "reject":
        assert "fixture rejects input" in log


def test_dragen_dux4_missing_binary_records_not_started(tmp_path):
    status = tmp_path / "attempt.json"
    log = tmp_path / "dux4_dragen.log"
    results = tmp_path / "dragen"
    proc = subprocess.run(
        [sys.executable, str(ROOT / "workflow/scripts/run_dragen_dux4.py"),
         "--sample", "example", "--bam", str(tmp_path / "fixture.bam"),
         "--hash-table", str(tmp_path), "--status", str(status),
         "--log", str(log), "--results", str(results)],
        env={**os.environ, "PATH": str(tmp_path)}, capture_output=True, text=True, timeout=10,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert json.loads(status.read_text())["outcome"] == "not_started"
    assert "Could not start DRAGEN" in log.read_text()
    assert results.is_dir() and not list(results.iterdir())


def test_clair3_rejects_existing_tensorflow_model(project):
    path, config = project
    fasta = Path(config["inputs"]["fasta"])
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\nACGT\n")
    Path(f"{fasta}.fai").write_text("chr1\t4\t6\t4\t5\n")
    bam = Path(config["run"]["dir"]) / "minimap2/example/example.sorted.bam"
    bam.parent.mkdir(parents=True)
    bam.touch()
    model = Path(config["caller_settings"]["clair3"]["model"])
    model.mkdir(parents=True)
    (model / "pileup.index").touch()
    proc = run_snakemake(project, "clair3/example/phased_merge_output.vcf.gz", "--allowed-rules", "clair3")
    output = proc.stdout + proc.stderr
    assert proc.returncode != 0
    log = path / "results_t2t/logs/example/clair3.log"
    assert log.exists(), output
    assert "requires PyTorch" in log.read_text()
    assert "run_clair3.sh: command not found" not in log.read_text()


def test_legacy_chm13_cache_mode_fails_with_offline_guidance(project):
    path, config = project
    config["annotation"]["vep"]["mode"] = "cache"
    config["annotation"]["vep"]["cache"] = str(path / "resources/vep")
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))
    proc = run_snakemake(project, "download_vep_cache")
    output = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "No official CHM13" in output
    assert "README.md" in output


def test_prepared_custom_caches_are_not_downloaded(project):
    _, config = project
    annotation = config["annotation"]
    snpeff = Path(annotation["snpeff"]["cache"]) / "CHM13.v2.0"
    vep = Path(f"{annotation['vep']['annotation']}.sorted.gz")
    snpeff.mkdir(parents=True)
    vep.touch()
    Path(f"{vep}.tbi").touch()
    proc = run_snakemake(project, "build_snpeff_cache", "prepare_vep_annotation", "--dry-run")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "Nothing to be done" in output


@pytest.mark.skipif(not os.environ.get("SNPEFF_EXECUTABLE"), reason="Set SNPEFF_EXECUTABLE for real SnpEff workflow test")
def test_real_snpeff_build_target(project, monkeypatch):
    path, config = project
    tools = path / "tools"
    tools.mkdir()
    (tools / "snpEff").symlink_to(os.environ["SNPEFF_EXECUTABLE"])
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ['PATH']}")
    fasta = Path(config["inputs"]["fasta"])
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\n" + "ATG" + "AAA" * 9 + "TAA\n")
    Path(f"{fasta}.fai").write_text("chr1\t33\t6\t33\t34\n")
    genes = Path(config["annotation"]["snpeff"]["build"]["genes"])
    genes.write_text(
        "##gff-version 3\n"
        "chr1\ttest\tgene\t1\t33\t.\t+\t.\tID=g1;Name=GENE1\n"
        "chr1\ttest\tmRNA\t1\t33\t.\t+\t.\tID=tx1;Parent=g1\n"
        "chr1\ttest\texon\t1\t33\t.\t+\t.\tParent=tx1\n"
        "chr1\ttest\tCDS\t1\t33\t.\t+\t0\tParent=tx1\n"
    )
    proc = run_snakemake(project, "build_snpeff_cache", "--set-resources", "build_snpeff_cache:mem_mb=512")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    database = Path(config["annotation"]["snpeff"]["cache"]) / "CHM13.v2.0"
    assert (database / "snpEffectPredictor.bin").stat().st_size > 0
    assert (database / "snpEff.config").exists()
    manifest = json.loads((database / "build_manifest.json").read_text())
    assert manifest["skipped_sequence_checks"] == ["CDS", "protein"]
    second = run_snakemake(project, "build_snpeff_cache", "--dry-run")
    assert second.returncode == 0, second.stdout + second.stderr
    assert "Nothing to be done" in second.stdout + second.stderr
    # Exercise the annotation rule too: it must pass the custom -c registration.
    vcf = path / "results_t2t/cutesv/example/cutesv.DEL.vcf"
    vcf.parent.mkdir(parents=True)
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=33>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t4\t.\tA\tG\t.\tPASS\t.\n"
    )
    annotation = run_snakemake(
        project, "cutesv/example/cutesv.DEL.snpeff.vcf",
        "--allowed-rules", "snpeff", "--set-resources", "snpeff:mem_mb=512",
    )
    assert annotation.returncode == 0, annotation.stdout + annotation.stderr
    result = vcf.with_name("cutesv.DEL.snpeff.vcf").read_text()
    assert "missense_variant" in result and "GENE1" in result


@pytest.mark.skipif(not os.environ.get("VEP_EXECUTABLE"), reason="Set VEP_EXECUTABLE for real offline VEP test")
@pytest.mark.parametrize("mode", ["gff", "gtf"])
def test_real_offline_vep_and_maf_conversion(project, monkeypatch, mode):
    path, config = project
    config["annotation"]["vep"]["mode"] = mode
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))
    tools = path / "tools"
    tools.mkdir()
    (tools / "vep").symlink_to(os.environ["VEP_EXECUTABLE"])
    if os.environ.get("VCF2MAF_EXECUTABLE"):
        (tools / "vcf2maf.pl").symlink_to(os.environ["VCF2MAF_EXECUTABLE"])
    monkeypatch.setenv("PATH", f"{tools}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("VEP_TEST_ROOT", str(path))
    fasta = Path(config["inputs"]["fasta"])
    fasta.parent.mkdir(parents=True)
    fasta.write_text(">chr1\n" + "ATG" + "AAA" * 9 + "TAA\n")
    Path(f"{fasta}.fai").write_text("chr1\t33\t6\t33\t34\n")
    genes = Path(config["annotation"]["vep"]["annotation"])
    if mode == "gff":
        genes.write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t33\t.\t+\t.\tID=g1;Name=GENE1;biotype=protein_coding\n"
            "chr1\ttest\ttranscript\t1\t33\t.\t+\t.\tID=tx1;Parent=g1;biotype=protein_coding\n"
            "chr1\ttest\texon\t1\t33\t.\t+\t.\tID=ex1;Parent=tx1\n"
            "chr1\ttest\tCDS\t1\t33\t.\t+\t0\tID=cds1;Parent=tx1\n"
        )
    else:
        attrs = 'gene_id "g1"; transcript_id "tx1"; gene_name "GENE1"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
        genes.write_text("".join(
            f"chr1\ttest\t{feature}\t1\t33\t.\t+\t{phase}\t{attrs}\n"
            for feature, phase in [("gene", "."), ("transcript", "."), ("exon", "."), ("CDS", "0")]
        ))
    vcf = path / "results_t2t/cutesv/example/cutesv.DEL.vcf"
    vcf.parent.mkdir(parents=True)
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=33>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t4\tsv1\tA\tG\t30\tPASS\t.\n"
    )
    proc = run_snakemake(project, "cutesv/example/cutesv.DEL.vep.vcf", "--allowed-rules", "vep", "prepare_vep_annotation")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    annotated = vcf.with_name("cutesv.DEL.vep.vcf")
    text = annotated.read_text()
    assert "missense_variant" in text and "GENE1" in text
    assert not (path / "resources/vep/homo_sapiens/116_CHM13").exists()
    if os.environ.get("VCF2MAF_EXECUTABLE"):
        merged = vcf.parent / "merged/cutesv.DEL.vep.vcf"
        merged.parent.mkdir()
        shutil.copyfile(annotated, merged)
        conversion = run_snakemake(project, "cutesv/example/merged/cutesv.DEL.vep.maf", "--allowed-rules", "convert_vep")
        assert conversion.returncode == 0, conversion.stdout + conversion.stderr
        maf = merged.with_suffix(".maf").read_text()
        assert "GENE1" in maf and "sv1" in maf


@pytest.mark.skipif(not (shutil.which("Rscript") and shutil.which("bcftools")), reason="R and bcftools required")
@pytest.mark.parametrize("policy", ["keep", "exclude", "error"])
def test_filter_annotations_missing_frequency_policy(project, policy):
    path, config = project
    config["filtering"]["missing_population_af"] = policy
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))
    caller = path / "results_t2t/cutesv/example"
    merged = caller / "merged"
    merged.mkdir(parents=True)
    (caller / "cutesv.DEL.annotsv.length.tsv").write_text("ID\tSV_chrom\tSV_start\tSV_end\nsv1\t1\t1\t10\n")
    (merged / "cutesv.DEL.vep.maf").write_text("vcf_id\tStart_Position\tvcf_qual\nsv1\t1\t30\n")
    (merged / "cutesv.DEL.snpeff.tsv").write_text("ID\tQUAL\nsv1\t30\n")
    tab = path / "results_t2t/survivor/example/example.cutesv.DEL.tab"
    tab.parent.mkdir(parents=True)
    tab.write_text("chr1\t1\tA\tG\t" + "\t".join("sv1" for _ in config["workflow"]["callers"]) + "\n")
    vcf_text = ("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t1\tsv1\tA\tG\t30\tPASS\t.\n")
    for suffix in ("vcf", "vep.vcf", "snpeff.vcf"):
        (merged / f"cutesv.DEL.{suffix}").write_text(vcf_text)
    proc = run_snakemake(project, "cutesv/example/merged/filtered/cutesv.DEL.vcf", "--allowed-rules", "filter_annotations")
    output = proc.stdout + proc.stderr
    if policy == "error":
        assert proc.returncode != 0
        assert "frequency unknown" in output
        return
    assert proc.returncode == 0, output
    final = merged / "filtered/cutesv.DEL.vcf"
    records = [line for line in final.read_text().splitlines() if not line.startswith("#")]
    assert len(records) == (1 if policy == "keep" else 0)
    if policy == "keep":
        table = (merged / "filtered/cutesv.DEL.tsv").read_text()
        assert "VEP.population_AF_status" in table and "unknown" in table


def test_prebuilt_snpeff_mode_still_uses_download_rule(project):
    path, config = project
    del config["annotation"]["snpeff"]["build"]
    config["annotation"]["genome"] = "GRCh38"
    config["annotation"]["snpeff"]["version"] = "87"
    (path / "config/config.yaml").write_text(yaml.safe_dump(config))
    proc = run_snakemake(project, "download_snpeff_cache", "--dry-run", "--printshellcmds")
    output = proc.stdout + proc.stderr
    assert proc.returncode == 0, output
    assert "snpEff download GRCh38.87" in output
    assert "build_snpeff_cache" not in output
