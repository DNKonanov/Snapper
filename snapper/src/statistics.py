from scipy.stats import ks_2samp
from random import sample
import numpy as np


SAMPLESIZE = 200

# just a patch
MINSAMPLESIZE = 10

LOG10_PVAL_TRH = 5



def get_statistics(
    native_motifs, 
    wga_motifs, 
    maxsamplesize=SAMPLESIZE,
    minsamplesize=MINSAMPLESIZE,

):

    print('Getting difsignals...')
    

    contigs = list(native_motifs.keys())

    motifs_lines = {}
    ks_stat_lines = {}

    for contig in contigs:
        ks_stat_line = []
        motifs_line = []
    
        cnt = 1


        for MOTIF in native_motifs[contig]:
            print('{}: {} out of {}'.format(contig, cnt, len(native_motifs[contig])), end='')
            
            if MOTIF not in wga_motifs[contig]:
                cnt += 1
                print('\r', end ='')
                continue

            
            try:
                s1 = sample(native_motifs[contig][MOTIF], maxsamplesize)
            except ValueError:
                s1 = native_motifs[contig][MOTIF]
                
            try:
                s2 = sample(wga_motifs[contig][MOTIF], maxsamplesize)
            except ValueError:
                s2 = wga_motifs[contig][MOTIF]
            
            if len(s1) < minsamplesize or len(s2) < minsamplesize:
                print('\r', end ='')
                cnt += 1
                continue
            
            ks_stat_line.append(ks_2samp(s1,s2, mode='asymp')[1])
            
            motifs_line.append(MOTIF)
            cnt += 1
            print('\r', end ='')
        print()

        ks_stat_lines[contig] = ks_stat_line
        motifs_lines[contig] = motifs_line

    return motifs_lines, ks_stat_lines


def get_difsignals(
    motifs_line, 
    ks_stats, 
    log10_pval_thr = LOG10_PVAL_TRH,

    ):
    

    passed_motifs = []

    for i in range(len(ks_stats)):
        if -np.log10(ks_stats[i]) >= log10_pval_thr:
            passed_motifs.append(motifs_line[i])
    return passed_motifs




