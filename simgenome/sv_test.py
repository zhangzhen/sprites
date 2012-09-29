import unittest
from bx.seq.fasta import FastaFile
import sv

validseq = "CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGTTAGGGAGCTGTGGACCCTGCAGC" \
         + "CTGGCTGTGGGGGCCGCAGTGGCTGAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCAC" \
         + "AGGGGCTTAACCTCTGGTGACTGCCAGAGCTGCTGGCAAGCTAGAGTCCCATTTGGAGCC" \
         + "CCTCTAAGCCGTTCTATTTGTAATGAAAACTATATTTATGCTATTCAGTTCTAAATATAG" \
         + "AAATTGAAACAGCTGTGTTTAGTGCCTTTGTTCAACCCCCTTGCAACAACCTTGAGAACC" \
         + "CCAGGGAATTTGTCAATGTCAGGGAAGGAGCATTTTGTCAGTTACCAAATGTGTTTATTA" \
         + "CCAGAGGGATGGAGGGAAGAGGGACGCTGAAGAACTTTGATGCCCTCTTCTTCCAAAGAT" \
         + "GAAACGCGTAACTGCGCTCTCATTCACTCCAGCTCCCTGTCACCCAATGGACCTGAATTT" \
         + "CCCAGAATCCAGATATCACTTCATCCTGGACCCTGAGAGATTCTGCAGCCCAGCTCCAGA" \
         + "TTGCTTGTGGTCTGACAGGCTGCAACTGTGAGCCATCACAATGAACAACAGGAAGAAAAG" \
         + "GTCTTTCAAAAGGTGATGTGTGTTCTCATCAACCTCATACACACACATGGTTTAGGGGTA" \
         + "TAATACCTCTACATGGCTGATTATGCTCTTCTGGCTTGTAAGGTTTCTGCTGAGAAGCCC" \
         + "ATTGTTAAAACAATGTTCCCCAGATAC"


class SvTest(unittest.TestCase):
    def setUp(self):
        self.ref = FastaFile(file("data/ex1.fa", "r"))
        self.reader = sv.GffReader(file("data/svlist.gff", "r"))
    def test_reversecomplement(self):
        self.assertEqual(sv.reversecomplement("GTGATATCTGGATTCTGGGAAATTC"), "GAATTTCCCAGAATCCAGATATCAC")
    def test_loadevents(self):
        events = sv.loadevents(self.reader, "seq1")
        self.assertEqual(len(events), 3)
	self.assertTrue(isinstance(events[0], sv.Deletion))
	self.assertTrue(isinstance(events[1], sv.Inversion))
	self.assertTrue(isinstance(events[2], sv.Insertion))
    def test_buildchrom(self):
        events = sv.loadevents(self.reader, "seq1")
        seq = sv.buildchrom(self.ref, events)
        self.assertEqual(seq.text[0:747], validseq)
    def test_isgffsorted(self):
        self.assertTrue(sv.isgffsorted(sv.GffReader(file("data/svlist.gff", "r"))))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(SvTest)
    unittest.TextTestRunner(verbosity=2).run(suite)


