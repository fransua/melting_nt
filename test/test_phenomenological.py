"""
08 Nov 2012


"""

"""
07 Nov 2012


"""

"""
30 Oct 2012

unittest for pytadbit functions
"""

import unittest
from melting_nt.phenomenological import oligo_Tm


SEQUENCES = ['GACGACAAGACCGCG',
             'ATGCCATTGCAATGCCATTGCAATGCCATTGCAATGCCATTGCAATGCCATTGCAATGCCATTGCAATGC']

class test_oligo_Tm(unittest.TestCase):
    """
    test main tadbit funtion
    """
    def test_calculator(self):
        tms = []
        for seq in SEQUENCES:
            tms.append(oligo_Tm(seq))
        self.assertEqual(tms, [65.06861352737879, 87.96493053740234])

        
if __name__ == "__main__":
    unittest.main()
