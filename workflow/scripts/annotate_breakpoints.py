"""Report genes overlapping each breakpoint in a consensus SV VCF."""

import gzip
import re
from collections import defaultdict
from urllib.parse import unquote


MATE = re.compile(r"([\[\]])([^:\[\]]+):(\d+)\1")


def open_text(path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def attributes(text):
    if re.match(r"[^=;\s]+=", text):
        # GFF3 delimiters are split before decoding escaped values.
        return {key: unquote(value) for key, value in
                (item.split("=", 1) for item in text.split(";") if "=" in item)}
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def load_genes(path):
    genes = defaultdict(list)
    with open_text(path) as handle:
        for line in handle:
            if line.rstrip() == "##FASTA":
                break
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "pseudogene"}:
                continue
            info = attributes(fields[8])
            name = next((info[key] for key in ("gene_name", "gene", "Name", "gene_id", "ID")
                         if info.get(key)), "unknown")
            genes[fields[0]].append((int(fields[3]), int(fields[4]), name))
    if not genes:
        raise ValueError(f"No gene or pseudogene features found in {path}")
    for intervals in genes.values():
        intervals.sort()
    return genes


def info_map(text):
    return dict(item.split("=", 1) if "=" in item else (item, "")
                for item in text.split(";") if item)


def breakpoint_records(vcf_path):
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, record_id, alt = fields[0], int(fields[1]), fields[2], fields[4]
            info = info_map(fields[7])
            yield record_id, chrom, pos, "left"

            # BND mates can be encoded only in ALT; take both coordinates together.
            mate = MATE.search(alt)
            if mate:
                yield record_id, mate[2], int(mate[3]), "right"
                continue

            end = info.get("END", "")
            chrom2 = info.get("CHR2") or chrom
            if chrom2 == ".":
                chrom2 = chrom
            if end.isdigit() and int(end) > 0 and (chrom2, int(end)) != (chrom, pos):
                yield record_id, chrom2, int(end), "right"


def genes_at(genes, chrom, position, flank):
    start = max(1, position - flank)
    end = position + flank
    names = set()
    for gene_start, gene_end, name in genes.get(chrom, []):
        if gene_start > end:
            break
        if gene_end >= start:
            names.add(name)
    return sorted(names)


def annotate_breakpoints(vcf, gene_models, output, flank=0):
    genes = load_genes(gene_models)
    with open(output, "w", encoding="utf-8") as handle:
        handle.write("record_id\tchrom\tposition\tbreakend\tgene_names\n")
        for record_id, chrom, position, breakend in breakpoint_records(vcf):
            names = genes_at(genes, chrom, position, flank)
            handle.write(f"{record_id}\t{chrom}\t{position}\t{breakend}\t"
                         f"{','.join(names) if names else 'intergenic'}\n")


if "snakemake" in globals():
    annotate_breakpoints(
        snakemake.input.vcf, snakemake.input.genes, snakemake.output[0],
    )
