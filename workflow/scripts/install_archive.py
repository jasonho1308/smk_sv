#!/usr/bin/env python3
"""Download and atomically install a tar or zip resource directory."""

from __future__ import annotations

import os
from pathlib import Path
import hashlib
import shutil
import subprocess
import tarfile
import tempfile
import zipfile


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
    if digest.hexdigest().lower() != expected.lower():
        raise ValueError(f"Checksum mismatch for {path}")


def download(url: str, destination: Path) -> Path:
    if not url:
        raise ValueError("No download URL configured for this resource")

    artifact = destination.parent / f".{destination.name}.download"
    artifact.parent.mkdir(parents=True, exist_ok=True)
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


def safe_extract_tar(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive) as handle:
        for member in handle.getmembers():
            member_path = (destination / member.name).resolve()
            if destination not in member_path.parents and member_path != destination:
                raise ValueError(f"Unsafe archive member: {member.name}")
        handle.extractall(destination, filter="data")


def safe_extract_zip(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            member_path = (destination / member.filename).resolve()
            if destination not in member_path.parents and member_path != destination:
                raise ValueError(f"Unsafe archive member: {member.filename}")
        handle.extractall(destination)


def install_archive(
    url: str,
    destination: str | Path,
    archive_type: str,
    archive_root: str | None = None,
    checksum: str | None = None,
) -> None:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)

    artifact = download(url, destination.parent / f"{destination.name}.archive")
    verify_checksum(artifact, checksum)

    with tempfile.TemporaryDirectory(
        prefix=f".{destination.name}.", dir=destination.parent
    ) as temporary_directory:
        unpacked = Path(temporary_directory) / "unpacked"
        unpacked.mkdir()

        if archive_type == "tar":
            safe_extract_tar(artifact, unpacked)
        elif archive_type == "zip":
            safe_extract_zip(artifact, unpacked)
        else:
            raise ValueError(f"Unsupported archive type: {archive_type}")

        source = unpacked / archive_root if archive_root else unpacked
        if not source.exists():
            raise ValueError(f"Expected archive path was not found: {archive_root}")

        staged = Path(temporary_directory) / "installed"
        os.replace(source, staged)
        if destination.exists():
            shutil.rmtree(destination)
        os.replace(staged, destination)

    artifact.unlink(missing_ok=True)


if __name__ == "__main__":
    install_archive(
        url=snakemake.params.url,  # type: ignore[name-defined]
        destination=snakemake.params.destination,  # type: ignore[name-defined]
        archive_type=snakemake.params.archive_type,  # type: ignore[name-defined]
        archive_root=snakemake.params.archive_root,  # type: ignore[name-defined]
        checksum=snakemake.params.checksum,  # type: ignore[name-defined]
    )
