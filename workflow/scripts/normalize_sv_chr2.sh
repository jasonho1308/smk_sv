#!/usr/bin/env bash
# CHR2 means the chromosome of the second breakpoint, NOT chromosome 2.
# For same-chromosome DEL/DUP/INV calls, Sniffles may omit it even though
# vcf2maf requires it; use the record's CHROM when CHR2 is missing.
set -euo pipefail

if (( $# != 2 )); then
    echo "Usage: $0 input.vcf output.vcf" >&2
    exit 2
fi

awk -F '\t' 'BEGIN { OFS = FS }
    /^##INFO=<ID=CHR2,/ { has_chr2_header = 1 }
    /^#CHROM/ {
        if (!has_chr2_header)
            print "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of the second SV breakpoint\">"
        print
        next
    }
    /^#/ { print; next }
    {
        if (NF < 8 || $1 == "") {
            print "Invalid VCF record at line " NR > "/dev/stderr"
            exit 1
        }
        n = split($8, tags, ";")
        type = ""
        chr2 = ""
        chr2_index = 0
        for (i = 1; i <= n; i++) {
            if (tags[i] ~ /^SVTYPE=/) type = substr(tags[i], 8)
            if (tags[i] ~ /^CHR2=/) {
                chr2 = substr(tags[i], 6)
                chr2_index = i
            }
        }
        if ((type == "DEL" || type == "DUP" || type == "INV") && chr2 == "") {
            if (chr2_index) {
                tags[chr2_index] = "CHR2=" $1
                $8 = tags[1]
                for (i = 2; i <= n; i++) $8 = $8 ";" tags[i]
            } else {
                $8 = ($8 == "." ? "" : $8 ";") "CHR2=" $1
            }
        }
        print
    }
' "$1" > "$2"
