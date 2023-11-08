import unittest
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, ".")

import unittest
from CRISPR_TAPE.general_function import context, find_motifs_in_amino_acid_sequence, general_function, make_distance_matrix, find_sign_change
from CRISPR_TAPE.shared_functions import PAMposition

class TestScriptFunctions(unittest.TestCase):

    def test_context(self):
        # execution
        actual_context_1 = context(2, ["M", "V", "K", "P", "L", "F", "C"], "K")
        actual_context_2 = context(1, ["M", "V", "K", "P", "L", "F", "C"], "V")
        actual_context_3 = context(6, ["M", "V", "K", "P", "L", "F", "C"], "C")
        actual_context_4 = context(2, ["M", "V", "K", "P", "L", "F", "C"], "KPL")
        actual_context_5 = context(0, ["M", "V", "K", "P", "L", "F", "C"], "MVKP")
        # assertion
        self.assertEqual(actual_context_1, "MV(K)PL")
        self.assertEqual(actual_context_2, "M(V)KP")
        self.assertEqual(actual_context_3, "LF(C)")
        self.assertEqual(actual_context_4, "MV(KPL)FC")
        self.assertEqual(actual_context_5, "(MVKP)LF")
        self.assertRaises(AssertionError, context, 0, ["M", "V", "K", "P", "L", "F", "C"], "PLF")

    def test_find_motifs_in_amino_acid_sequence(self):
        # setup
        protein_dict = {0: {'Amino Acid': 'M', 'base_1': (1, 'A'), 'base_2': (2, 'C'), 'base_3': (3, 'A')},
                        1: {'Amino Acid': 'V', 'base_1': (4, 'G'), 'base_2': (5, 'T'), 'base_3': (6, 'C')},
                        2: {'Amino Acid': 'A', 'base_1': (7, 'G'), 'base_2': (8, 'C'), 'base_3': (9, 'C')},
                        3: {'Amino Acid': 'D', 'base_1': (10, 'G'), 'base_2': (11, 'A'), 'base_3': (12, 'T')},
                        4: {'Amino Acid': 'L', 'base_1': (13, 'C'), 'base_2': (14, 'T'), 'base_3': (15, 'T')},
                        5: {'Amino Acid': 'A', 'base_1': (16, 'G'), 'base_2': (17, 'C'), 'base_3': (18, 'C')},
                        6: {'Amino Acid': 'S', 'base_1': (19, 'T'), 'base_2': (20, 'C'), 'base_3': (21, 'C')},
                        7: {'Amino Acid': 'D', 'base_1': (22, 'G'), 'base_2': (23, 'A'), 'base_3': (24, 'T')},
                        8: {'Amino Acid': 'E', 'base_1': (25, 'G'), 'base_2': (26, 'A'), 'base_3': (27, 'G')},
                        9: {'Amino Acid': 'W', 'base_1': (28, 'T'), 'base_2': (29, 'G'), 'base_3': (30, 'G')}}
        # execution
        actual_amino_position_1, actual_amino_acids_1 = find_motifs_in_amino_acid_sequence("M", protein_dict)
        actual_amino_position_2, actual_amino_acids_2 = find_motifs_in_amino_acid_sequence("W", protein_dict)
        actual_amino_position_3, actual_amino_acids_3 = find_motifs_in_amino_acid_sequence("A", protein_dict)
        actual_amino_position_4, actual_amino_acids_4 = find_motifs_in_amino_acid_sequence("DLA", protein_dict)
        # assertion
        self.assertTrue(all(p in actual_amino_position_1 for p in [0]))
        self.assertEqual(actual_amino_acids_1, ["M", "V", "A", "D", "L", "A", "S", "D", "E", "W"])
        self.assertTrue(all(p in actual_amino_position_2 for p in [9]))
        self.assertEqual(actual_amino_acids_2, ["M", "V", "A", "D", "L", "A", "S", "D", "E", "W"])
        self.assertTrue(all(p in actual_amino_position_3 for p in [2, 5]))
        self.assertEqual(actual_amino_acids_3, ["M", "V", "A", "D", "L", "A", "S", "D", "E", "W"])
        self.assertTrue(all(p in actual_amino_position_4 for p in [3]))
        self.assertEqual(actual_amino_acids_4, ["M", "V", "A", "D", "L", "A", "S", "D", "E", "W"])

    def test_make_distance_matrix(self):
        # setup
        base_list = [10, 20, 30]
        guide_array = np.array([5, 15, 25, 35])
        # execution
        actual_output = make_distance_matrix(base_list, guide_array)
        # assertion
        expected_output = [
            [-5, 5, 15, 25],
            [-15, -5, 5, 15],
            [-25, -15, -5, 5]
        ]
        self.assertEqual(list(actual_output[0]), expected_output[0])
        self.assertEqual(list(actual_output[1]), expected_output[1])
        self.assertEqual(list(actual_output[2]), expected_output[2])

    def test_no_sign_change(self):
        self.assertIsNone(find_sign_change([1, 2, 3, 4, 5]))

    def test_sign_change_positive_to_negative(self):
        self.assertEqual(find_sign_change([1, 2, -3, 4]), 1)

    def test_sign_change_negative_to_positive(self):
        self.assertEqual(find_sign_change([-1, -2, 3, 4]), 1)

    def test_sign_change_at_end(self):
        self.assertEqual(find_sign_change([1, 2, 3, -4]), 2)

    def test_sign_change_at_beginning(self):
        self.assertEqual(find_sign_change([-1, 2, 3, 4]), 0)

    def test_empty_list(self):
        self.assertIsNone(find_sign_change([]))

    def test_single_element(self):
        self.assertIsNone(find_sign_change([1]))

    def test_zero_in_list(self):
        self.assertEqual(find_sign_change([-1, 0, 1]), 0)

    def test_general_function_single_aa_NGG(self):
        # setup
        motif = "L"
        PAM = "NGG"
        DNA = "aaaatgttggATGGTCTCCGAGCTGCAGCGCCAGCTGGtttggggccctttCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTTCGTGCTCAaaatttagg"
        CDS = "ATGGTCTCCGAGCTGCAGCGCCAGCTGGCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTTCGTGCTCA"
        reference_genome = "AAAATGTTGGATGGTCTCCGAGCTGCAGCGCCAGCTGGTTTGGGGCCCTTTCGCTGCATCGGCAGACCCGCGGTGTAGGGTCTTCGTCGACTGCTTCGTGCTCAAAATTTAGGAACCAGCTGGCGCTGCAGCTCGGCCGATGCAGCGAAAGGGCCCCAA"
        # execution
        actual_guide_positions, actual_gRNA_list, actual_guide_strands = PAMposition(DNA.upper(), PAM)
        actual_guides = general_function(motif,
                                    PAM,
                                    DNA,
                                    CDS,
                                    "",
                                    "",
                                    reference_genome)
        # assertion
        expected_gRNA_list = ["CTCCGAGCTGCAGCGCCAGCTGG",
                        "AGCTGCAGCGCCAGCTGGTTTGG",
                        "GCTGCAGCGCCAGCTGGTTTGGG",
                        "CTGCAGCGCCAGCTGGTTTGGGG",
                        "TTGGGGCCCTTTCGCTGCATCGG",
                        "TCGCTGCATCGGCAGACCCGCGG",
                        "CATCGGCAGACCCGCGGTGTAGG",
                        "ATCGGCAGACCCGCGGTGTAGGG",
                        "CTGCTTCGTGCTCAAAATTTAGG",
                        "AACCAGCTGGCGCTGCAGCTCGG",
                        "CGAAAGGGCCCCAAACCAGCTGG",
                        "GGTCTGCCGATGCAGCGAAAGGG",
                        "GGGTCTGCCGATGCAGCGAAAGG",
                        "CGACGAAGACCCTACACCGCGGG",
                        "TCGACGAAGACCCTACACCGCGG"]
        self.assertEqual(len(actual_gRNA_list), len(expected_gRNA_list))
        self.assertTrue(all(g in expected_gRNA_list for g in actual_gRNA_list) and all(g in actual_gRNA_list for g in expected_gRNA_list))
        self.assertEqual(len(actual_guides), 3)
        self.assertEqual(list(actual_guides["Amino Acid Position"]), [5, 9, 11])
        self.assertEqual(list(actual_guides["Context"]), ["SE(L)QR", "RQ(L)AL", "LA(L)HR"])
        actual_guides.to_csv("unittest_guide.csv")
        # test the 5' guides
        self.assertEqual(list(actual_guides["5' gRNA Sequence"]), ["AACCAGCTGGCGCTGCAGCT",
                                                                "CTCCGAGCTGCAGCGCCAGC",
                                                                "GGGTCTGCCGATGCAGCGAA"])
        self.assertEqual(list(actual_guides["5' PAM"]), ["CGG", "TGG", "AGG"])
        self.assertEqual(list(actual_guides["5' gRNA Strand"]), ["reverse", "forward", "reverse"])
        self.assertEqual(list(actual_guides["Distance of 5' Cut Site from Amino Acid (bp)"]), [-1, 2, 1])
        self.assertEqual(list(actual_guides["5' gRNA G/C Content (%)"]), [69.57, 73.91, 65.22])
        self.assertEqual(list(actual_guides["5' Notes"]), ["No leading G.", "No leading G.", ""])
        self.assertEqual(list(actual_guides["5' gRNA Off Target Count"]), [1, 0, 0])
        # test the 3' guides"
        self.assertEqual(list(actual_guides["3' gRNA Sequence"]), ["CTCCGAGCTGCAGCGCCAGC",
                                                                "AGCTGCAGCGCCAGCTGGTT",
                                                                "TTGGGGCCCTTTCGCTGCAT"])
        self.assertEqual(list(actual_guides["3' PAM"]), ["TGG", "TGG", "CGG"])
        self.assertEqual(list(actual_guides["3' gRNA Strand"]), ["forward", "forward", "forward"])
        self.assertEqual(list(actual_guides["Distance of 3' Cut Site from Amino Acid (bp)"]), [7, 0, 0])
        self.assertEqual(list(actual_guides["3' gRNA G/C Content (%)"]), [73.91, 65.22, 65.22])
        self.assertEqual(list(actual_guides["3' Notes"]), ["No leading G.", "No leading G.", "No leading G."])
        self.assertEqual(list(actual_guides["3' gRNA Off Target Count"]), [1, 0, 1])

if __name__ == "__main__":
    unittest.main()