#!/usr/bin/env python3
"""Download one configured resource atomically, with optional verification."""

from __future__ import annotations

import gzip
import hashlib
import os
from pathlib import Path
import shutil
import subprocess


def verify_checksum(path: Path, checksum: str | None) -> None:
    if not checksum:
        return

    try:
        algorithm, expected = checksum.split(":", 1)
        digest = hashlib.new(algorithm)
    except (ValueError, TypeError) as error:
        raise ValueError(
            "checksum must have the form 'md5:<digest>' or 'sha256:<digest>'"
        ) from error

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    observed = digest.hexdigest()
    if observed.lower() != expected.lower():
        raise ValueError(
            f"Checksum mismatch for {path}: expected {expected}, observed {observed}"
        )


def download(url: str, destination: Path) -> Path:
    if not url:
        raise ValueError("No download URL configured for this resource")

    artifact = destination.with_name(f".{destination.name}.download")
    subprocess.run(
        [
            "curl",
            "--fail",
            "--location",
            "--retry",
            "5",
            "--retry-delay",
            "5",
            "--continue-at",
            "-",
            "--output",
            str(artifact),
            url,
        ],
        check=True,
    )
    return artifact


def postprocess_file(path: Path, operation: str | None) -> None:
    if operation in (None, "", "none"):
        return
    if operation != "ucsc_canonical":
        raise ValueError(f"Unsupported postprocessing operation: {operation}")

    canonical = {str(chromosome) for chromosome in range(1, 23)} | {"X", "Y"}
    processed = path.with_name(f"{path.name}.processed")
    try:
        with path.open() as source, processed.open("w") as target:
            for line in source:
                fields = line.split("\t", 1)
                if fields[0] in canonical:
                    target.write(f"chr{line}")
                elif fields[0] in {"M", "MT"} and len(fields) == 2:
                    target.write(f"chrM\t{fields[1]}")
                else:
                    target.write(line)
        os.replace(processed, path)
    finally:
        processed.unlink(missing_ok=True)


def install_file(
    url: str,
    destination: str | Path,
    checksum: str | None = None,
    compression: str | None = None,
    postprocess: str | None = None,
) -> None:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    artifact = download(url, destination)
    verify_checksum(artifact, checksum)

    staged = destination.with_name(f".{destination.name}.tmp")
    staged.unlink(missing_ok=True)
    try:
        if compression in (None, "", "none"):
            os.replace(artifact, staged)
        elif compression == "gzip":
            with gzip.open(artifact, "rb") as source, staged.open("wb") as target:
                shutil.copyfileobj(source, target)
            artifact.unlink()
        else:
            raise ValueError(f"Unsupported compression type: {compression}")
        postprocess_file(staged, postprocess)
        os.replace(staged, destination)
    finally:
        staged.unlink(missing_ok=True)


if __name__ == "__main__":
    install_file(
        url=snakemake.params.url,  # type: ignore[name-defined]
        destination=snakemake.output[0],  # type: ignore[name-defined]
        checksum=snakemake.params.checksum,  # type: ignore[name-defined]
        compression=snakemake.params.compression,  # type: ignore[name-defined]
        postprocess=snakemake.params.postprocess,  # type: ignore[name-defined]
    )
