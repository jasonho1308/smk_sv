#!/usr/bin/env bash
# CHR2 means the chromosome of the second breakpoint, NOT chromosome 2.
# vcf2maf expects CHR2 and END for breakpoint events. Same-chromosome
# DEL/DUP/INV can use CHROM; BND/TRA must obtain their mate from ALT instead.
set -euo pipefail

if (( $# != 2 )); then
    echo "Usage: $0 input.vcf output.vcf" >&2
    exit 2
fi

awk -F '\t' 'BEGIN { OFS = FS }
    /^##INFO=<ID=CHR2,/ { has_chr2_header = 1 }
    /^##INFO=<ID=END,/ { has_end_header = 1 }
    /^#CHROM/ {
        if (!has_chr2_header)
            print "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome of the second SV breakpoint\">"
        if (!has_end_header)
            print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Position of the second SV breakpoint\">"
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
        end = ""
        chr2_index = end_index = 0
        for (i = 1; i <= n; i++) {
            if (tags[i] ~ /^SVTYPE=/) type = substr(tags[i], 8)
            if (tags[i] ~ /^CHR2=/) { chr2 = substr(tags[i], 6); chr2_index = i }
            if (tags[i] ~ /^END=/) { end = substr(tags[i], 5); end_index = i }
        }
        if (type == "DEL" || type == "DUP" || type == "INV") {
            if (end !~ /^[0-9]+$/ || end + 0 < 1) {
                print "Missing/invalid END at VCF line " NR > "/dev/stderr"
                exit 1
            }
            if (chr2 == "") {
                if (chr2_index) tags[chr2_index] = "CHR2=" $1
                else tags[++n] = "CHR2=" $1
            }
        } else if (type == "BND" || type == "TRA") {
            # Breakend ALT: N[chr:pos[, N]chr:pos], [chr:pos[N or ]chr:pos]N.
            bracket = index($5, "[") ? "[" : (index($5, "]") ? "]" : "")
            first = bracket == "" ? 0 : index($5, bracket)
            tail = first ? substr($5, first + 1) : ""
            second = bracket == "" ? 0 : index(tail, bracket)
            mate = second ? substr(tail, 1, second - 1) : ""
            if (mate ~ /^[^:]+:[0-9]+$/) {
                split(mate, parts, ":")
                if ((chr2 != "" && chr2 != parts[1]) ||
                    (end != "" && end != parts[2])) {
                    print "BND mate ALT disagrees with CHR2/END at VCF line " NR > "/dev/stderr"
                    exit 1
                }
                if (chr2 == "") {
                    if (chr2_index) tags[chr2_index] = "CHR2=" parts[1]
                    else tags[++n] = "CHR2=" parts[1]
                }
                if (end == "") {
                    if (end_index) tags[end_index] = "END=" parts[2]
                    else tags[++n] = "END=" parts[2]
                }
            } else if (chr2 == "" || end !~ /^[0-9]+$/ || end + 0 < 1) {
                print "BND has no usable mate in ALT or CHR2/END at VCF line " NR > "/dev/stderr"
                exit 1
            }
        }
        if (type == "DEL" || type == "DUP" || type == "INV" || type == "BND" || type == "TRA") {
            $8 = tags[1]
            for (i = 2; i <= n; i++) $8 = $8 ";" tags[i]
        }
        print
    }
' "$1" > "$2"
