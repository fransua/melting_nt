"""
14 Nov 2012


.. [Khandelwal2010] Khandelwal, G., & Bhyravabhotla, J. (2010). A phenomenological model for predicting melting temperatures of DNA sequences. PloS one, 5(8), e12433. doi:10.1371/journal.pone.0012433
.. [SantaLucia1998] SantaLucia, J., Allawi, H. T., & Seneviratne, P. a. (1996). Improved nearest-neighbor parameters for predicting DNA duplex stability. Biochemistry, 35(11), 3555-62. doi:10.1021/bi951907q
.. [SantaLucia2004] SantaLucia, J., & Hicks, D. (2004). The thermodynamics of DNA structural motifs. Annual review of biophysics and biomolecular structure, 33, 415–40. doi:10.1146/annurev.biophys.32.110601.141800

"""

from melting_nt.nearest_neighbor import oligo_Tm as oligo_Tm_nn
from melting_nt.phenomenological import oligo_Tm as oligo_Tm_ph

def oligo_Tm(seq, method='SAL1996', dna_conc=.2, salt_conc=0.22, corrector=4,
             verbose=False):
    """
    :argument seq: a 5'->3' DNA sequence. It assumes that we work with double strand DNA
    :argument SAL1996 method: can be either SAL1996, SAL2004 or PHENOME for using
    respectively methods from [SantaLucia1998]_, [SantaLucia2004]_ or [Khandelwal2010]_
    :argument 0.22 salt_conc: concentration in Na+ in M (only PHENOME method)
    :argument 2.0 dna_conc: concentration in Na+ in nM
    :argument 4 corrector: only for SAL methods
    :argument False verbose: only for SAL methods
    """
    if 'SAL' in method:
        return oligo_Tm_nn(seq, method=method, dna_conc=dna_conc,
                           corrector=corrector, verbose=verbose)
    elif method is 'PHENOME':
        return oligo_Tm_ph(seq, dna_conc=dna_conc, salt_conc=salt_conc)
