"""Fill missing AnnotSV insertion lengths from their retained INFO/SVLEN.

AnnotSV 3.5.10 can lose SVLEN when CIPOS expands the reported start.
Do not replace SV_start/SV_end: they are the CI-adjusted annotation interval.
"""

import sys
from pathlib import Path


def svlen_from_info(info):
    values = [item[6:] for item in info.split(";") if item.startswith("SVLEN=")]
    if len(values) > 1 and len(set(values)) > 1:
        raise ValueError(f"Conflicting SVLEN values in INFO: {info!r}")
    return values[0] if values and values[0] not in ("", ".") else None


def restore_lengths(annotsv_path, vcf_path, output_path):
    if Path(annotsv_path).stat().st_size == 0:
        Path(output_path).write_bytes(b"")
        print("AnnotSV input is empty; no insertion lengths to restore")
        return

    lengths = {}
    with open(vcf_path, encoding="utf-8") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line[:80]!r}")
            variant_id = fields[2]
            svlen = svlen_from_info(fields[7])
            if variant_id not in lengths:
                lengths[variant_id] = svlen
            elif lengths[variant_id] != svlen:
                # Duplicate VCF IDs cannot be used as a lookup key. Prefer
                # the INFO column belonging to each AnnotSV row below.
                lengths[variant_id] = None

    restored = 0
    with open(annotsv_path, encoding="utf-8", newline="") as source, open(
        output_path, "w", encoding="utf-8", newline=""
    ) as output:
        header = source.readline()
        columns = header.rstrip("\r\n").split("\t")
        required = ("ID", "SV_type", "SV_length")
        missing = [column for column in required if column not in columns]
        if missing:
            raise ValueError(
                f"AnnotSV input {annotsv_path}: missing header fields {missing}; "
                f"first line: {header[:200]!r} (file size: {Path(annotsv_path).stat().st_size} bytes)"
            )
        id_col, type_col, length_col = (columns.index(column) for column in required)
        info_col = columns.index("INFO") if "INFO" in columns else None
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
                svlen = svlen_from_info(fields[info_col]) if info_col is not None else None
                if svlen is None:
                    svlen = lengths.get(variant_id)
                if svlen is None:
                    raise ValueError(f"No unambiguous SVLEN for AnnotSV INS ID {variant_id}")
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
    print(f"Restored SV_length in {restored} AnnotSV INS rows from INFO/SVLEN")


if __name__ == "__main__":
    if "snakemake" in globals():
        restore_lengths(snakemake.input.annotsv, snakemake.input.vcf, snakemake.output.tsv)
    elif len(sys.argv) == 4:
        restore_lengths(*map(Path, sys.argv[1:]))
    else:
        raise SystemExit("Usage: restore_annotsv_length.py annotsv.tsv input.vcf output.tsv")
