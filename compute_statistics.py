from scipy.stats import ks_2samp
from random import sample
import numpy as np


SAMPLESIZE = 1000
MINSAMPLESIZE = 100


LOG10_PVAL_TRH = 10
EFFSIZE_TRH = 0.25


def cohend(d1, d2):
    n1, n2 = len(d1), len(d2)
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.mean(d1), np.mean(d2)
    return (u1 - u2) / s


def get_statistics(
    native_motifs, 
    wga_motifs, 
    maxsamplesize=SAMPLESIZE,
    minsamplesize=MINSAMPLESIZE,

):

    print('Getting difsignals...')
    ks_stats = []
    effsizes = []

    motifs_line = []

    cnt = 1

    motifs = list(native_motifs.keys())

    for MOTIF in motifs:
        print('{} from {}'.format(cnt, len(motifs)), end='')
        
        try:
            s1 = sample(native_motifs[MOTIF], maxsamplesize)
        except ValueError:
            s1 = native_motifs[MOTIF]
            
        try:
            s2 = sample(wga_motifs[MOTIF], maxsamplesize)
        except ValueError:
            s2 = wga_motifs[MOTIF]
        
        if len(s1) < minsamplesize or len(s2) < minsamplesize:
            continue
        
        ks_stats.append(ks_2samp(s1,s2, mode='asymp')[1])
        effsizes.append(cohend(s1, s2)) 
        
        motifs_line.append(MOTIF)
        cnt += 1
        print('\r', end ='')


def get_difsignals(
    motifs_line, 
    ks_stats, 
    effsizes,
    log10_pval_thr = LOG10_PVAL_TRH,
    effsize_thr = EFFSIZE_TRH

    ):
    

    passed_motifs = []

    for i in range(len(ks_stats)):
        if -np.log10(ks_stats[i]) >= log10_pval_thr and np.abs(effsizes[i]) >= effsize_thr:
            passed_motifs.append(motifs_line[i])
    return passed_motifs




def save_results (motifs, out_fasta):

    with open(out_fasta) as f:
        
        cnt = 1

        for m in motifs:
            out_fasta.write('>MOTIF_{}\n{}\n'.format(cnt, m))



