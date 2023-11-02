import unittest
import sys
sys.path.insert(0, ".")

from CRISPR_TAPE.shared_functions import translate, clean_inputs, PAMposition, get_codon_index, analyse_text, notes, pamcolumn, remove_pam, correct_distance, list_search, get_count

class TestScript(unittest.TestCase):

    def test_translate(self):
        # assertion
        self.assertEqual(translate("ATGCGA"), "MR")
        self.assertEqual(translate("ATGGGGTTT"), "MGF")
        self.assertEqual(translate(""), "")

    def test_clean_inputs(self):
        # execution
        loci_exon, loci, loci_edited, coding_sequence, organism_genome = clean_inputs("ATGC acgtTGCA",
                                                                                    "ATGC TGCA",
                                                                                    "ATGNTTGCCAAGGNG",
                                                                                    "aaa",
                                                                                    "ttt")
        # assertion
        self.assertEqual(loci_exon, "ATGCTGCA")
        self.assertEqual(loci, "aaaATGCacgtTGCAttt")
        self.assertEqual(loci_edited, "AAAATGCACGTTGCATTT")
        self.assertEqual(coding_sequence, "ATGCTGCA")
        self.assertEqual(organism_genome, "ATGTTGCCAAGGG")

    def test_PAMposition_NGG(self):
        # setup
        example_loci = "GGCCATGGATGCTAGCTACGGCTAACTGGGCTTAACGGATGCTAGCTGG"
        PAM = "NGG"
        # execution
        positions, entries, strands = PAMposition(example_loci, PAM)
        print(entries)
        # assertion
        self.assertTrue(5 == len(entries) == len(positions) == len(strands))
        self.assertEqual([8, 22, 23, 31, 42], positions)
        expected_entries = ['CCATGGATGCTAGCTACGGCTAA', 'GGATGCTAGCTACGGCTAACTGG', 'GATGCTAGCTACGGCTAACTGGG', 'CTACGGCTAACTGGGCTTAACGG', 'TGGGCTTAACGGATGCTAGCTGG']
        self.assertEqual(expected_entries, entries)
        expected_strands = ["reverse", "forward", "forward", "forward", "forward"]
        self.assertEqual(expected_strands, strands)

    def test_get_codon_index(self):
        dna = "ATGCGA"
        cds = "ATGCGA"
        codon_index = get_codon_index(dna, cds)
        self.assertIsInstance(codon_index, dict)
        self.assertIn("Amino Acid", codon_index[0])
        self.assertEqual(codon_index[0]["Amino Acid"], "M")

    def test_analyse_text(self):
        self.assertEqual(analyse_text("CGCG"), 100.00)
        self.assertEqual(analyse_text("ATAT"), 0.00)

    def test_notes(self):
        # execution
        actual_note_1 = notes("TTTTGCGC", 80)
        actual_note_2 = notes("GCGC", 50)
        actual_note_3 = notes("GTTTTGCGC", 40)
        actual_note_4 = notes("CGCGC", 40)
        actual_note_5 = notes("GCGCGC", 100)
        # assertion
        self.assertEqual(actual_note_1, 'PolyT present. No leading G. G/C content over 75%.')
        self.assertEqual(actual_note_2, '')
        self.assertEqual(actual_note_3, 'PolyT present.')
        self.assertEqual(actual_note_4, 'No leading G.')
        self.assertEqual(actual_note_5, 'G/C content over 75%.')

    def test_pamcolumn(self):
        # execution
        actual_pam_one = pamcolumn("CTGACTGGATTGCTGGTCAATGG", "forward", "CCATTGACCAGCAATCCAGTCAG", "NGG")
        actual_pam_two = pamcolumn("CCATTGACCAGCAATCCAGTCAG", "reverse", "CTGACTGGATTGCTGGTCAATGG", "NGG")
        actual_pam_three = pamcolumn("CTGACTGGATTGCTGGTCAATAG", "forward", "CTATTGACCAGCAATCCAGTCAG", "YG")
        actual_pam_four = pamcolumn("CTATTGACCAGCAATCCAGTCAG", "reverse", "CTGACTGGATTGCTGGTCAATAG", "YG")
        actual_pam_five = pamcolumn("TTTACTATTGACCAGCAATCCAG", "forward", "CTGGATTGCTGGTCAATAGTAAA", "TTTN")
        actual_pam_six = pamcolumn("CTGGATTGCTGGTCAATAGTAAA", "reverse", "TTTACTATTGACCAGCAATCCAG", "TTTN")
        # assertion
        self.assertEqual(actual_pam_one, "TGG")
        self.assertEqual(actual_pam_two, "TGG")
        self.assertEqual(actual_pam_three, "AG")
        self.assertEqual(actual_pam_four, "AG")
        self.assertEqual(actual_pam_five, "TTTA")
        self.assertEqual(actual_pam_six, "TTTA")
        self.assertRaises(AssertionError, pamcolumn, "CTGGATTGCTGGTCAATAGTGA", "reverse", "TCACTATTGACCAGCAATCCAG", "NGG")

    def test_remove_pam(self):
        # setup
        actual_guide_one = remove_pam("CTGACTGGATTGCTGGTCAATGG", "forward", "CCATTGACCAGCAATCCAGTCAG", "NGG")
        actual_guide_two = remove_pam("CCATTGACCAGCAATCCAGTCAG", "reverse", "CTGACTGGATTGCTGGTCAATGG", "NGG")
        actual_guide_three = remove_pam("CTGACTGGATTGCTGGTCAATAG", "forward", "CTATTGACCAGCAATCCAGTCAG", "YG")
        actual_guide_four = remove_pam("CTATTGACCAGCAATCCAGTCAG", "reverse", "CTGACTGGATTGCTGGTCAATAG", "YG")
        actual_guide_five = remove_pam("TTTACTATTGACCAGCAATCCAG", "forward", "CTGGATTGCTGGTCAATAGTAAA", "TTTN")
        actual_guide_six = remove_pam("CTGGATTGCTGGTCAATAGTAAA", "reverse", "TTTACTATTGACCAGCAATCCAG", "TTTN")
        # assertion
        self.assertEqual(actual_guide_one, "CTGACTGGATTGCTGGTCAA")
        self.assertEqual(actual_guide_two, "CTGACTGGATTGCTGGTCAA")
        self.assertEqual(actual_guide_three, "CTGACTGGATTGCTGGTCAAT")
        self.assertEqual(actual_guide_four, "CTGACTGGATTGCTGGTCAAT")
        self.assertEqual(actual_guide_five, "CTATTGACCAGCAATCCAG")
        self.assertEqual(actual_guide_six, "CTATTGACCAGCAATCCAG")
        self.assertRaises(AssertionError, pamcolumn, "CTGGATTGCTGGTCAATAGTGA", "reverse", "TCACTATTGACCAGCAATCCAG", "NGG")

    def test_correct_distance(self):
        self.assertEqual(correct_distance(10, "reverse"), 9)
        self.assertEqual(correct_distance(-5, "forward"), -4)

    def test_list_search(self):
        counts = list_search(["ATG"], "ATGCATGCATG")
        self.assertIn("ATG", counts)

    def test_get_count(self):
        self.assertEqual(get_count("ATG", "CAT", {"ATG": 2, "CAT": 1}), 2)

if __name__ == '__main__':
    unittest.main()
