"""Synthetic indexed BAM + real BLAST; run in workflow/envs/dux4.yaml."""

import csv
from pathlib import Path
from random import Random
from shutil import which
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest

try:
    import pysam
except ImportError:
    pysam = None

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/scan_dux4.py"


@unittest.skipUnless(pysam and which("blastn") and which("makeblastdb"), "requires pysam and BLAST+")
class TestDux4Scan(unittest.TestCase):
    def test_partner_first_scan_and_repeat_ambiguity(self):
        with TemporaryDirectory() as path:
            folder = Path(path)
            rng = Random(23)

            def synth(n):
                return "".join(rng.choices("ACGT", k=n))

            repeat = synth(400)
            partner = synth(2000)
            reference = folder / "test.fa"
            chr4_flank = synth(1100)
            reference.write_text(f">chr4\n{synth(500)}{repeat}{chr4_flank}\n"
                                 f">chr10\n{synth(500)}{repeat}{synth(1100)}\n"
                                 f">chr14\n{partner}\n")
            pysam.faidx(str(reference))
            panel = folder / "partners.bed"
            panel.write_text("chr14\t300\t1550\tIGH\nchr14\t1900\t1970\tEmpty\n")
            sources = folder / "sources.bed"
            sources.write_text("chr4\t500\t900\tD4Z4_4q35\n")
            bam = folder / "sample.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": c, "LN": 2000} for c in ("chr4", "chr10", "chr14")]}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as output:
                def add(name, pos, cigar, seq, sa=None, flag=0, mapq=60):
                    record = pysam.AlignedSegment(output.header)
                    record.query_name = name
                    record.flag = flag
                    record.reference_id = 2
                    record.reference_start = pos
                    record.mapping_quality = mapq
                    record.cigarstring = cigar
                    record.query_sequence = seq
                    if sa:
                        record.set_tag("SA", sa)
                    output.write(record)

                add("split", 300, "250M400S", partner[300:550] + repeat,
                    "chr4,501,+,250S400M,2,0;")
                add("reverse", 305, "250M400S", partner[305:555] + repeat,
                    "chr4,501,+,400M250H,0,0;", flag=16)
                add("left_clip", 550, "400S250M", repeat + partner[550:800])
                add("paired", 550, "400S250M400S", repeat + partner[550:800] + repeat)
                add("paired", 550, "400S250M400S", repeat + partner[550:800] + repeat, flag=2048)
                add("insertion", 610, "250M400I250M", partner[610:860] + repeat + partner[860:1110])
                add("chr10_sa_only", 650, "250M400H", partner[650:900],
                    "chr10,501,+,250S400M,60,0;", flag=2048)
                linked = partner[850:1100] + repeat + partner[1102:1352]
                add("linked", 850, "250M650S", linked)
                add("flank", 975, "250M650S", partner[975:1225] + repeat + chr4_flank[:250],
                    "chr4,501,+,250S650M,60,0;")
                add("linked", 1102, "650S250M", linked, flag=2048)
                add("distant_clip", 1125, "250M650S", partner[1125:1375] + synth(250) + repeat)
                add("ordinary", 1150, "300M", partner[1150:1450])
                add("short_clip", 1200, "10H90S250M", repeat[:90] + partner[1200:1450])
                add("tiny_clip", 1220, "250M5S", partner[1220:1470] + repeat[:5])
                add("short_anchor", 1250, "75M400S", partner[1250:1325] + repeat)
                add("low_mapq", 1300, "220M400S", partner[1300:1520] + repeat, mapq=0)
                add("long_clip", 1310, "200M3600S", partner[1310:1510] + synth(3200) + repeat)
                add("off_target", 1500, "250M400S", partner[1500:1750] + repeat)
                add("outside_windows", 1640, "250M400S", partner[1640:1890] + repeat, flag=256)
            pysam.index("-c", str(bam))
            reports = {n: folder / f"{n}.tsv" for n in ("junctions", "reads", "loci", "blast")}
            reports["filtered-blast"] = folder / "blast_filtered.tsv"
            command = [sys.executable, str(SCRIPT), "--bam", str(bam), "--bam-index", f"{bam}.csi",
                       "--reference", str(reference), "--partners", str(panel), "--sources", str(sources),
                       "--padding", "0", "--threads", "1"]
            for name, report in reports.items():
                command.extend((f"--{name}", str(report)))
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

            def rows(name):
                with reports[name].open() as handle:
                    return list(csv.DictReader(handle, delimiter="\t"))

            loci = {r["partner"]: r for r in rows("loci")}
            self.assertEqual(loci["Empty"]["candidate_clusters"], "0")
            by_cut = {(r["cut0"], r["side"]): r for r in rows("junctions")}
            split = by_cut["550", "R"]
            self.assertEqual((split["reads"], split["split_reads"], split["geometry"]),
                             ("2", "2", "repeat_associated_adjacency"))
            self.assertEqual(split["flank_placements"], ".")
            self.assertEqual(split["status"], "review")
            left = by_cut["550", "L"]
            self.assertEqual((left["reads"], left["blast_reads"], left["status"]), ("2", "2", "review"))
            self.assertEqual(left["nearby_opposite_ends"], ".")
            self.assertEqual(by_cut["1102", "L"]["nearby_opposite_ends"], f"{by_cut['1100', 'R']['junction_id']}:1")
            self.assertEqual(by_cut["860", "I"]["geometry"], "insertion_like")
            self.assertEqual(by_cut["1225", "R"]["flank_placements"], "chr4")
            self.assertEqual(by_cut["1225", "R"]["geometry"], "interchromosomal_adjacency")
            self.assertNotIn(("1375", "R"), by_cut)
            self.assertNotIn(("1750", "R"), by_cut)
            read_rows = rows("reads")
            self.assertFalse(any(r["read"] == "chr10_sa_only" for r in read_rows))
            clip = next(r for r in read_rows if r["read"] == "left_clip")
            self.assertEqual((clip["source_chrom"], clip["source_status"]), ("chr4", "ambiguous"))
            reverse = next(r for r in read_rows if r["read"] == "reverse")
            self.assertEqual((reverse["source_side"], reverse["source_cut"], reverse["query_gap"]), ("R", "900", "0"))
            broad = {(r["read"], r["partner_side"]): r for r in rows("blast")}
            filtered = rows("filtered-blast")
            self.assertCountEqual(filtered, [{k: v for k, v in r.items() if k != "junction_id"}
                                             for r in read_rows if r["method"] == "softclip_blast"])
            self.assertEqual(set(broad), {
                ("split", "R"), ("reverse", "R"), ("left_clip", "L"), ("paired", "L"), ("paired", "R"),
                ("linked", "L"), ("linked", "R"), ("flank", "R"), ("distant_clip", "R"), ("short_clip", "L"),
                ("short_anchor", "R"), ("low_mapq", "R"), ("long_clip", "R"), ("off_target", "R"),
                ("outside_windows", "R"),
            })
            self.assertEqual(broad["short_clip", "L"]["blast_length"], "90")
            self.assertGreaterEqual(int(broad["long_clip", "R"]["query_gap"]), 3000)
            self.assertEqual(broad["outside_windows", "R"]["partner_cut"], "1890")
            self.assertTrue(all(r["partner"] == "." for r in broad.values()))
            self.assertIn("Genome-wide soft-clip BLAST queries: 18", run.stdout)

            guarded = subprocess.run([*command, "--require-chm13"], capture_output=True, text=True)
            self.assertNotEqual(guarded.returncode, 0)
            self.assertIn("Expected the chr-prefixed T2T-CHM13v2.0", guarded.stderr)


if __name__ == "__main__":
    unittest.main()
