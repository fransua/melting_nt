"""
07 Nov 2012

computes DNA melting point according to Phenomenological method [Khandelwal2010]_


.. [Khandelwal2010] Khandelwal, G., & Bhyravabhotla, J. (2010). A phenomenological model for predicting melting temperatures of DNA sequences. PloS one, 5(8), e12433. doi:10.1371/journal.pone.0012433

"""

from math import log

ENERGY = {'AA':  5, 'AT':  7, 'AG':  8, 'AC': 10,
          'TA':  4, 'TT':  5, 'TG':  7, 'TC':  8,
          'GA':  8, 'GT': 10, 'GG': 11, 'GC': 13,
          'CA':  7, 'CT':  8, 'CG': 10, 'CC': 11}


def oligo_Tm(seq, salt_conc=0.22, dna_conc=2.0):
    """
    :argument seq: a 5'->3' DNA sequence. It assumes that we work with double strand DNA
    :argument 0.22 salt_conc: concentration in Na+ in M
    :argument 2.0 dna_conc: concentration in Na+ in nM
    """
    strength = sum([ENERGY[seq[i:i+2]] for i in xrange(len(seq)-1)])
    return 7.35 * float(strength)/len(seq) + 17.34 * log(len(seq)) \
           + 4.96 * log(salt_conc) + 0.89 * log(dna_conc*1000000) - 25.42
