"""
07 Nov 2012


"""

"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from melting_nt.nearest_neighbor import oligo_Tm


SEQUENCES = ['ATATATATATATATATATA',
             'ATATATATATATATATAT',
             'GCATATATGCA',
             'GCGCGCGCGCGCGCGCGCG']

class test_oligo_Tm(unittest.TestCase):
    """
    test main tadbit funtion
    """
    def test_SAL1996(self):
        tms = []
        for seq in SEQUENCES:
            tms.append(oligo_Tm(seq, method='SAL1996'))
        self.assertEqual(tms, [19.830917544829276, 17.45583181343744, 18.270824015439985, 81.40827106738902])


    def test_SAL2004(self):
        tms = []
        for seq in SEQUENCES:
            tms.append(oligo_Tm(seq, method='SAL2004'))
        self.assertEqual(tms, [28.81964469222504, 25.852532749514864, 22.010722145333943, 81.42812254178148])

        
if __name__ == "__main__":
    unittest.main()
