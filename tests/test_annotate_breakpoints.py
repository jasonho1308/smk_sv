import gzip
import runpy
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/annotate_breakpoints.py"


def test_breakpoint_tsv(tmp_path):
    gtf = tmp_path / "genes.gtf.gz"
    with gzip.open(gtf, "wt") as handle:
        handle.write(
            'chr1\tref\tgene\t8\t10\t.\t+\t.\tgene_name "ZETA";\n'
            'chr1\tref\tgene\t9\t12\t.\t+\t.\tgene_name "ALPHA";\n'
            'chr1\tref\tgene\t10\t20\t.\t+\t.\tgene_name "ALPHA";\n'
            'chr1\tref\texon\t98\t100\t.\t+\t.\tgene_name "EXON";\n'
            'chr2\tref\tgene\t42\t45\t.\t+\t.\tgene_id "gene2";\n'
        )
    vcf = tmp_path / "consensus.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t12\tdel\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=50\n"
        "chr1\t100\tbnd-alt\tN\tN]chr2:40]\t.\tPASS\tSVTYPE=BND\n"
        "chr1\t40\tbnd-info\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr2;END=40\n"
        "chr1\t12\tsingle\tN\tN.\t.\tPASS\tSVTYPE=BND\n"
        "chr1\t99\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=99\n"
    )
    out = tmp_path / "breakpoints.tsv"
    annotate = runpy.run_path(str(SCRIPT))["annotate_breakpoints"]
    annotate(vcf, gtf, out, flank=2)

    assert out.read_text().splitlines() == [
        "record_id\tchrom\tposition\tbreakend\tgene_names",
        "del\tchr1\t12\tleft\tALPHA,ZETA",
        "del\tchr1\t50\tright\tintergenic",
        "bnd-alt\tchr1\t100\tleft\tintergenic",
        "bnd-alt\tchr2\t40\tright\tgene2",
        "bnd-info\tchr1\t40\tleft\tintergenic",
        "bnd-info\tchr2\t40\tright\tgene2",
        "single\tchr1\t12\tleft\tALPHA,ZETA",
        "ins\tchr1\t99\tleft\tintergenic",
    ]
