"""Fill missing AnnotSV insertion lengths from the VCF given to AnnotSV.

AnnotSV 3.5.10 can lose SVLEN when CIPOS expands the reported start.
Do not replace SV_start/SV_end: they are the CI-adjusted annotation interval.
"""

import sys
from pathlib import Path


def restore_lengths(annotsv_path, vcf_path, output_path):
    lengths = {}
    with open(vcf_path, encoding="utf-8") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line[:80]!r}")
            variant_id = fields[2]
            svlens = [item[6:] for item in fields[7].split(";") if item.startswith("SVLEN=")]
            if not svlens or svlens[0] in ("", "."):
                continue
            if variant_id in lengths and lengths[variant_id] != svlens[0]:
                raise ValueError(f"Ambiguous SVLEN for VCF ID {variant_id}")
            lengths[variant_id] = svlens[0]

    restored = 0
    with open(annotsv_path, encoding="utf-8", newline="") as source, open(
        output_path, "w", encoding="utf-8", newline=""
    ) as output:
        header = source.readline()
        columns = header.rstrip("\r\n").split("\t")
        required = ("ID", "SV_type", "SV_length")
        if not all(column in columns for column in required):
            raise ValueError(f"AnnotSV header missing one of {required}")
        id_col, type_col, length_col = (columns.index(column) for column in required)
        output.write(header)
        for line in source:
            if not line.strip() or line.startswith("#"):
                output.write(line)
                continue
            newline = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError(f"AnnotSV row has {len(fields)} fields; expected {len(columns)}")
            if fields[type_col] == "INS" and fields[length_col] in ("", "NA", "NaN"):
                variant_id = fields[id_col]
                svlen = lengths.get(variant_id)
                if svlen is None:
                    raise ValueError(f"No VCF SVLEN for AnnotSV INS ID {variant_id}")
                try:
                    if int(svlen) <= 0:
                        raise ValueError
                except ValueError:
                    raise ValueError(f"Invalid VCF SVLEN {svlen!r} for {variant_id}") from None
                fields[length_col] = svlen
                restored += 1
                output.write("\t".join(fields) + newline)
            else:
                output.write(line)
    print(f"Restored SV_length in {restored} AnnotSV INS rows from VCF SVLEN")


if __name__ == "__main__":
    if "snakemake" in globals():
        restore_lengths(snakemake.input.annotsv, snakemake.input.vcf, snakemake.output.tsv)
    elif len(sys.argv) == 4:
        restore_lengths(*map(Path, sys.argv[1:]))
    else:
        raise SystemExit("Usage: restore_annotsv_length.py annotsv.tsv input.vcf output.tsv")
