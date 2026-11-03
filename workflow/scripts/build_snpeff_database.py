#!/usr/bin/env python3
"""Validate and atomically build a custom SnpEff GFF3/GTF database.

Only the caller's data directory is written; installed SnpEff config is untouched.
"""

import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


# SnpEff's Vertebrate_Mitochondrial table (NCBI translation table 2).
# https://github.com/pcingola/SnpEff/blob/master/config/snpEff.config
MITOCHONDRIAL_TABLE = (
    "TTT/F,TTC/F,TTA/L,TTG/L,TCT/S,TCC/S,TCA/S,TCG/S,TAT/Y,TAC/Y,TAA/*,TAG/*,"
    "TGT/C,TGC/C,TGA/W,TGG/W,CTT/L,CTC/L,CTA/L,CTG/L,CCT/P,CCC/P,CCA/P,CCG/P,"
    "CAT/H,CAC/H,CAA/Q,CAG/Q,CGT/R,CGC/R,CGA/R,CGG/R,ATT/I+,ATC/I+,ATA/M+,ATG/M+,"
    "ACT/T,ACC/T,ACA/T,ACG/T,AAT/N,AAC/N,AAA/K,AAG/K,AGT/S,AGC/S,AGA/*,AGG/*,"
    "GTT/V,GTC/V,GTA/V,GTG/V+,GCT/A,GCC/A,GCA/A,GCG/A,GAT/D,GAC/D,GAA/E,GAG/E,"
    "GGT/G,GGC/G,GGA/G,GGG/G"
)


def validate_annotations(genes, fai):
    """Reject mismatched contig names, invalid coordinates, or empty annotations."""
    lengths = {}
    with Path(fai).open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError("Reference FASTA index is empty")

    count = 0
    opener = gzip.open if str(genes).endswith(".gz") else open
    with opener(genes, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                raise ValueError("Use a GFF without embedded FASTA; the configured reference must be used")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Annotation line {number}: expected 9 GFF/GTF columns")
            chrom, start, end = fields[0], int(fields[3]), int(fields[4])
            if chrom not in lengths:
                raise ValueError(f"Annotation contig {chrom!r} is absent from the reference FASTA index")
            if not 1 <= start <= end <= lengths[chrom]:
                raise ValueError(f"Annotation line {number}: coordinates outside reference contig {chrom}")
            count += 1
    if not count:
        raise ValueError("No annotation features found")
    return count


def build_database(fasta, fai, genes, destination, database, format="gff3",
                   mem_mb=32000, cds=None, protein=None, log=None, executable="snpEff"):
    if format not in ("gff3", "gtf22"):
        raise ValueError(f"Unsupported annotation format: {format}")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", database):
        raise ValueError("Invalid SnpEff database identifier")
    if mem_mb < 1:
        raise ValueError("mem_mb must be positive")
    destination = Path(destination).resolve()
    if destination.name != database:
        raise ValueError("Destination directory must be named after the database")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite an existing database: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    log = Path(log) if log else destination.parent / f"{database}.build.log"
    log.parent.mkdir(parents=True, exist_ok=True)

    with log.open("w") as handle:
        count = validate_annotations(genes, fai)
        print(f"Validated {count} annotation features against {fai}", file=handle, flush=True)
        with tempfile.TemporaryDirectory(prefix=f".{database}.", dir=destination.parent) as temporary:
            data_dir = Path(temporary)
            staged = data_dir / database
            staged.mkdir()
            suffix = ".gz" if str(genes).endswith(".gz") else ""
            annotation_name = ("genes.gff" if format == "gff3" else "genes.gtf") + suffix
            sources = {"sequences.fa": fasta, annotation_name: genes}
            for name, source in (("cds.fa", cds), ("protein.fa", protein)):
                if source:
                    sources[name + (".gz" if str(source).endswith(".gz") else "")] = source
            for name, source in sources.items():
                source = Path(source).resolve(strict=True)
                (staged / name).symlink_to(source)

            # Standard code is built into SnpEff. Register the vertebrate mitochondrial
            # code as well; SnpEff normalizes chrM to M internally.
            config = staged / "snpEff.config"
            config.write_text(
                f"{database}.genome : Custom reference {database}\n"
                f"codon.Vertebrate_Mitochondrial : {MITOCHONDRIAL_TABLE}\n"
                f"{database}.M.codonTable : Vertebrate_Mitochondrial\n"
                f"{database}.MT.codonTable : Vertebrate_Mitochondrial\n"
            )
            # Leave room for the JVM's native memory within the scheduler request.
            heap_mb = max(1, int(mem_mb * 0.9))
            command = [executable, f"-Xmx{heap_mb}M", "build", f"-{format}",
                       "-v", "-c", str(config), "-dataDir", str(data_dir)]
            skipped = []
            for source, flag, name in ((cds, "-noCheckCds", "CDS"),
                                       (protein, "-noCheckProtein", "protein")):
                if not source:
                    command.append(flag)
                    skipped.append(name)
                    print(f"WARNING: {name} sequence cross-check skipped: no matching FASTA supplied",
                          file=handle, flush=True)
            command.append(database)
            print(shlex.join(command), file=handle, flush=True)
            subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
            predictor = staged / "snpEffectPredictor.bin"
            if not predictor.is_file() or predictor.stat().st_size == 0:
                raise ValueError("SnpEff exited without producing a nonempty snpEffectPredictor.bin")
            with Path(genes).open("rb") as annotation_handle:
                digest = hashlib.file_digest(annotation_handle, "sha256").hexdigest()
            (staged / "build_manifest.json").write_text(json.dumps({
                "database": database, "format": format,
                "fasta": str(Path(fasta).resolve()), "fai": str(Path(fai).resolve()),
                "genes": str(Path(genes).resolve()), "genes_sha256": digest,
                "features": count, "skipped_sequence_checks": skipped,
                "cds": str(Path(cds).resolve()) if cds else None,
                "protein": str(Path(protein).resolve()) if protein else None,
                "command": command,
            }, indent=2) + "\n")
            os.replace(staged, destination)


if __name__ == "__main__":
    build_database(
        fasta=snakemake.input.fasta,
        fai=snakemake.input.fai,
        genes=snakemake.input.genes,
        destination=snakemake.output.cache,
        database=snakemake.params.database,
        format=snakemake.params.format,
        mem_mb=snakemake.resources.mem_mb,
        cds=snakemake.input.cds[0] if snakemake.input.cds else None,
        protein=snakemake.input.protein[0] if snakemake.input.protein else None,
        log=snakemake.log[0],
    )
