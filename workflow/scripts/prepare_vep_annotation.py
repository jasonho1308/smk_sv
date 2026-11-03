#!/usr/bin/env python3
"""Validate, coordinate-sort, bgzip and tabix-index an offline VEP GFF/GTF."""

import gzip
import os
from pathlib import Path
import subprocess
import tempfile


def prepare_annotation(genes, fai, destination, threads=1, log=None, mem_mb=4000):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    log = Path(log) if log else destination.with_name(f"{destination.name}.log")
    log.parent.mkdir(parents=True, exist_ok=True)
    with Path(fai).open() as handle:
        lengths = {fields[0]: int(fields[1]) for fields in (line.split("\t") for line in handle)}
    if not lengths:
        raise ValueError("Reference FASTA index is empty")

    with log.open("w") as errors, tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.", dir=destination.parent
    ) as temporary:
        work = Path(temporary)
        # GNU/BSD sort spill to disk rather than retaining the entire human GFF.
        raw = work / "records"
        opener = gzip.open if str(genes).endswith(".gz") else open
        count = 0
        with opener(genes, "rt") as source, raw.open("w") as target:
            for number, line in enumerate(source, 1):
                if line.startswith("##FASTA"):
                    raise ValueError("Embedded FASTA is not supported; use the configured reference FASTA")
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Annotation line {number}: expected 9 GFF/GTF columns")
                chrom, start, end = fields[0], int(fields[3]), int(fields[4])
                if chrom not in lengths:
                    raise ValueError(f"Annotation contig {chrom!r} is absent from the reference FASTA index")
                if not 1 <= start <= end <= lengths[chrom]:
                    raise ValueError(f"Annotation line {number}: coordinates outside reference contig {chrom}")
                target.write(line if line.endswith("\n") else line + "\n")
                count += 1
        if not count:
            raise ValueError("No GFF/GTF annotation features found")
        sorted_path = work / "sorted"
        with sorted_path.open("wb") as target:
            subprocess.run(
                ["sort", "-S", f"{max(1, int(mem_mb * 0.5))}M", "-s", "-k1,1", "-k4,4n", str(raw)], check=True,
                stdout=target, stderr=errors,
                env={**os.environ, "LC_ALL": "C", "TMPDIR": str(work)},
            )
        raw.unlink()
        compressed = work / "annotation.gz"
        with compressed.open("wb") as target:
            subprocess.run(["bgzip", "-c", "-@", str(threads), str(sorted_path)],
                           check=True, stdout=target, stderr=errors)
        subprocess.run(["tabix", "-p", "gff", str(compressed)], check=True, stderr=errors)
        index = Path(f"{compressed}.tbi")
        if not index.is_file() or not index.stat().st_size:
            raise ValueError("tabix did not produce a nonempty index")
        os.replace(compressed, destination)
        os.replace(index, f"{destination}.tbi")
        print(f"Validated and indexed {count} annotation features", file=errors)


if __name__ == "__main__":
    prepare_annotation(
        genes=snakemake.input.genes,
        fai=snakemake.input.fai,
        destination=snakemake.output.annotation,
        threads=snakemake.threads,
        log=snakemake.log[0],
        mem_mb=snakemake.resources.mem_mb,
    )
