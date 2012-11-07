"""
07 Nov 2012

computes DNA melting point according to Nearest-Neighbor method [SanataLucia1998]_


.. [SantaLucia1998] SantaLucia, J., Allawi, H. T., & Seneviratne, P. a. (1996). Improved nearest-neighbor parameters for predicting DNA duplex stability. Biochemistry, 35(11), 3555–62. doi:10.1021/bi951907q

Note:
* conversion from paper's tables:
  * TT -> AA
  * AG -> CT
  * AC -> GT
  * TG -> CA
  * TC -> GA
  * CC -> GG

"""

from operator import div
from math import log
from string import maketrans

REV = maketrans('ATGC','TACG')

SAL1996 = {
    'AA':   {'H':  -8400.0, 'S': -23.6},
    'AT':   {'H':  -6500.0, 'S': -18.8},
    'AG':   {'H':  -6100.0, 'S': -16.1},
    'AC':   {'H':  -8600.0, 'S': -23.0},
    'TA':   {'H':  -6300.0, 'S': -18.5},
    'TT':   {'H':  -8400.0, 'S': -23.6},
    'TG':   {'H':  -7400.0, 'S': -19.3},
    'TC':   {'H':  -7700.0, 'S': -20.3},
    'GA':   {'H':  -7700.0, 'S': -20.3},
    'GT':   {'H':  -8600.0, 'S': -23.0},
    'GG':   {'H':  -6700.0, 'S': -15.6},
    'GC':   {'H': -11100.0, 'S': -28.4},
    'CA':   {'H':  -7400.0, 'S': -19.3},
    'CT':   {'H':  -6100.0, 'S': -16.1},
    'CG':   {'H': -10100.0, 'S': -25.5},
    'CC':   {'H':  -6700.0, 'S': -15.6},
    'stAT': {'H':      0.0, 'S':  -9.0},
    'stGC': {'H':      0.0, 'S':  -5.9},
    'ndAT': {'H':    400.0, 'S':   0.0},
    'symm': {'H':      0.0, 'S': -1.40}
    }

SAL2004 = {
    'AA':   {'H':  -7600.0, 'S': -21.3},
    'AT':   {'H':  -7200.0, 'S': -20.4},
    'AG':   {'H':  -7800.0, 'S': -21.0},
    'AC':   {'H':  -8400.0, 'S': -22.4},
    'TA':   {'H':  -7200.0, 'S': -21.3},
    'TT':   {'H':  -7600.0, 'S': -21.3},
    'TG':   {'H':  -8500.0, 'S': -22.7},
    'TC':   {'H':  -8200.0, 'S': -22.2},
    'GA':   {'H':  -8200.0, 'S': -22.2},
    'GT':   {'H':  -8400.0, 'S': -22.4},
    'GG':   {'H':  -8000.0, 'S': -19.9},
    'GC':   {'H':  -9800.0, 'S': -24.4},
    'CA':   {'H':  -8500.0, 'S': -22.7},
    'CT':   {'H':  -7800.0, 'S': -21.0},
    'CG':   {'H': -10600.0, 'S': -27.2},
    'CC':   {'H':  -8000.0, 'S': -19.9},    
    'stAT': {'H':    200.0, 'S': -5.70},
    'stGC': {'H':    200.0, 'S': -5.70},
    'ndAT': {'H':   2200.0, 'S':  6.90},
    'symm': {'H':      0.0, 'S': -1.40}
    }

def delta_G(delta_H, delta_S, temp=37):
    return (delta_H - delta_S * (temp + 273.15)) / 1000

def reverse(seq):
    return seq.translate(REV)

def is_symmetric(seq, len_seq):
    half = int(len_seq/2)
    return seq[:half] == reverse(seq[-half:])[::-1]

def get_nn_params(method):
    if method=='SAL1996':
        return SAL1996
    elif method=='SAL2004':
        return SAL2004
    Exception('Method %s not implemented yet\n' % method)

def oligo_Tm(seq, method='SAL1996', dna_conc=.2, corrector=4, verbose=False):
    """
    computes melting temperature of a given sequence.
    TODO: take into account Na+ concentration?

    :argument seq: a 5'->3' DNA sequence. I it assumed that we work with double strand DNA
    :argument SAL1996 method: parameters for the Nearest Neighbor algorithm
    """
    nn_params = get_nn_params(method)
    len_seq = len(seq)
    # R * ln(concentration of DNA/1 or 4 if complementary strands or not)
    const = 1.9872 * log(float(str(dna_conc)+'e-9')/corrector)
    # deltaHd and deltaSd
    delta_Hd = 0
    delta_Sd = 0
    for pos in xrange(len_seq-1):
        bp = seq[pos:pos+2]
        delta_Hd += nn_params[bp]['H']
        delta_Sd += nn_params[bp]['S']
        print bp, nn_params[bp]['S']
    # delta Hi and delta Si
    delta_Hi = delta_Si = 0
    if not len(seq.translate(None, 'AT')):
        delta_Hi = nn_params['stAT']['H']
        delta_Si = nn_params['stAT']['S']
    else:
        delta_Hi = nn_params['stGC']['H']
        delta_Si = nn_params['stGC']['S']
    # delta self if self complementary sequence
    delta_self = nn_params['symm']['S'] if is_symmetric(seq, len_seq) else 0
    # terminal penality
    delta_Ht = delta_St = 0
    if method is 'SAL1996':
        if seq[-1] is 'A':
            delta_Ht = nn_params['ndAT']['H']
            delta_St = nn_params['ndAT']['S']
    elif method is 'SAL2004':
        if seq[0] in 'AT':
            delta_Ht = nn_params['ndAT']['H']
            delta_St = nn_params['ndAT']['S']
        if seq[-1] in 'AT':
            delta_Ht = nn_params['ndAT']['H']
            delta_St = nn_params['ndAT']['S']

    if verbose:
        print 'Enthalpy: {:,}'.format(delta_Hd + delta_Hi + delta_Ht)
        print 'Entropy: {:,}'.format(delta_Sd + delta_Si + delta_St + delta_self)
        
    return div(delta_Hd + delta_Hi + delta_Ht,
               delta_Sd + delta_Si + delta_St + delta_self + const) - 273.15

