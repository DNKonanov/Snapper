from scipy.stats import ks_2samp
from random import sample
import numpy as np
import os
from multiprocessing import Process, Manager
from tqdm import tqdm

SAMPLESIZE = 200

# just a patch
MINSAMPLESIZE = 10

LOG10_PVAL_TRH = 5

def get_ks_test(
    motif_subset, 
    ks_stat_line,
    motifs_line, 
    native_motifs, 
    wga_motifs, 
    contig, 
    minsamplesize,
    maxsamplesize):

    ks_test_batch = []
    motif_batch = []

    for MOTIF in motif_subset:
        if MOTIF not in wga_motifs[contig].keys():
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
            continue

        ks_test_batch.append(ks_2samp(s1,s2, mode='asymp')[1])
        motif_batch.append(MOTIF)
    #return_data[os.getpid()] = {}
    ks_stat_line += ks_test_batch
    motifs_line += motif_batch

def get_statistics(
    native_motifs, 
    wga_motifs, 
    maxsamplesize=SAMPLESIZE,
    minsamplesize=MINSAMPLESIZE,
    threads=1,
    ):

    print('Getting difsignals...')
    contigs = list(native_motifs.keys())
    motifs_lines = {}
    ks_stat_lines = {}

    for contig in contigs:

      #  ks_stat_line = []
       # motifs_line = []
        procs = []
        KS_manager = Manager()
        ks_stat_line = KS_manager.list()

        motif_manager = Manager()
        motifs_line = motif_manager.list()
        
        interval_coordinates = np.linspace(0, len(native_motifs[contig].keys()), threads+1)
        intervals = [(interval_coordinates[idx], interval_coordinates[idx+1]) for idx,_ in  list(enumerate(interval_coordinates))[: -1]]
        motifs = list(native_motifs[contig].keys())

        # Thread filtering
        threads_limit = len(motifs)
        if threads > threads_limit:
            threads = threads_limit

        for thread_process in tqdm(intervals):
            
            motif_subset = motifs[int(thread_process[0]): int(thread_process[1])]
            #print(motif_subset)
            proc = Process(target=get_ks_test, args=(motif_subset, ks_stat_line, motifs_line, native_motifs, wga_motifs, contig, minsamplesize, maxsamplesize))
            procs.append(proc)
            proc.start()
        
        for proc in procs:
            proc.join()
        
        print()
        ks_stat_lines[contig] = list(ks_stat_line)
        motifs_lines[contig] = list(motifs_line)

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