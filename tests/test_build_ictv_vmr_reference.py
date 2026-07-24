from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_ictv_vmr_reference.py"
SPEC = importlib.util.spec_from_file_location("build_ictv_vmr_reference", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class IctvVmrReferenceBuilderTest(unittest.TestCase):
    def test_accession_extraction_handles_segment_labels_and_versions(self) -> None:
        self.assertEqual(MODULE.accessions_from_cell("Seg1: OP436269; Seg2: OP436270"), ["OP436269", "OP436270"])
        self.assertEqual(MODULE.accessions_from_cell("NC_001802.1; PP467602"), ["NC_001802", "PP467602"])

    def test_baltimore_mapping_covers_vmr_genome_notation(self) -> None:
        self.assertEqual(MODULE.baltimore_group("dsDNA"), "I")
        self.assertEqual(MODULE.baltimore_group("ssDNA(+/-)"), "II")
        self.assertEqual(MODULE.baltimore_group("ssRNA(+/-)"), "V")
        self.assertEqual(MODULE.baltimore_group("ssRNA-RT"), "VI")
        self.assertEqual(MODULE.baltimore_group("unknown"), "unclassified")

    def test_ncbi_cds_header_is_linked_to_source_accession_and_protein(self) -> None:
        header = "lcl|PP467602.1_prot_WYC14516.1_1 [protein=DNA-binding protein] [protein_id=WYC14516.1]"
        self.assertEqual(MODULE.accession_from_header(header), "PP467602")
        self.assertEqual(MODULE.protein_from_header(header, 1), "WYC14516.1")
