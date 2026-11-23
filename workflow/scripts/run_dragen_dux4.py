#!/usr/bin/env python3
"""Run DRAGEN DUX4 as an exploratory branch, not as part of SV consensus."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
from tempfile import TemporaryDirectory


def run(args):
    args.status.parent.mkdir(parents=True, exist_ok=True)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    exit_code = None
    vcfs = []
    with TemporaryDirectory(prefix=".dux4-", dir=args.status.parent) as tmp:
        stage = Path(tmp) / "dragen"
        stage.mkdir()
        command = [
            "dragen", "-r", str(args.hash_table),
            "--enable-map-align=false", "--enable-map-align-output=false",
            "--enable-sort=false", "--tumor-bam-input", str(args.bam),
            "--enable-dux4-caller=true", f"--output-dir={stage}",
            f"--output-file-prefix={args.sample}",
        ]
        if args.skip_sanity_check:
            # Illumina's v4.3 option is spelled "santiy" in the product guide.
            command.extend(["--dux4-skip-santiy-check", "true"])

        with args.log.open("w", encoding="utf-8") as log:
            print(shlex.join(command), file=log, flush=True)
            if not args.hash_table.is_dir():
                outcome = "not_started"
                print("The DRAGEN hash table must be a directory.", file=log)
            else:
                try:
                    exit_code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
                except OSError as error:
                    outcome = "not_started"
                    print(f"Could not start DRAGEN: {error}", file=log)
                else:
                    if exit_code == 0:
                        vcfs = sorted(
                            path.name for path in stage.iterdir()
                            if path.name.startswith(f"{args.sample}.")
                            and path.name.endswith((".vcf", ".vcf.gz"))
                            and path.is_file() and path.stat().st_size > 0
                        )
                        outcome = "vcf_written" if vcfs else "no_vcf"
                    else:
                        outcome = "failed"
            print(f"DRAGEN DUX4 outcome: {outcome} (exit code: {exit_code})", file=log)

        if exit_code == 0:
            stage.rename(args.results)
        else:
            # Do not expose partial VCFs from a failed run as caller results.
            args.results.mkdir()

    args.status.write_text(json.dumps({
        "outcome": outcome,
        "exit_code": exit_code,
        "skip_sanity_check": args.skip_sanity_check,
        "vcfs": [f"dragen/{name}" for name in vcfs],
    }, indent=2) + "\n", encoding="utf-8")
    if outcome != "vcf_written":
        print(f"WARNING: DRAGEN DUX4 for {args.sample}: {outcome}. See {args.log}.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    for name in ("bam", "hash-table", "status", "results", "log"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--skip-sanity-check", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
