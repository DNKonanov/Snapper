from multiprocessing import Pool
from unittest import result
from webbrowser import get
import numpy as np
import h5py
from itertools import product
from Bio.SeqIO import parse
from Bio.Seq import reverse_complement
import os
from tqdm import tqdm

letters = ['A','G','C','T']


def get_reference(reference_file, target_chr='all'):

    ref_file = parse(reference_file, format='fasta')

    
    if target_chr == 'all':

        refs = {}
        reverse_refs = {}        
        for rec in ref_file:
            
            seq = str(rec.seq)
            contig = str(rec.description).split(' ')[0]
            refs[contig] = seq
            reverse_refs[contig] = reverse_complement(seq)
        return refs, reverse_refs
    
    if target_chr == 'longest':
        length = 0

        refs = {}
        reverse_refs = {}

        for rec in ref_file:
            
            contig = str(rec.description).split(' ')[0]
            if len(rec.seq) > length:
                seq = str(rec.seq)
                
                refs = {contig : seq}
                reverse_refs = {contig : reverse_complement(seq)}
                length = str(rec.seq)

        return refs, reverse_refs

    else:
        refs = {}
        reverse_refs = {}
        for rec in ref_file:
            contig = str(rec.description).split(' ')[0]
            if contig == target_chr:
                seq = str(rec.seq)

                refs[target_chr] = seq
                reverse_refs[target_chr] = reverse_complement(seq)
        
        if len(refs) == 0:
            raise KeyError('{} contig does not exist!'.format(target_chr))

        return refs, reverse_refs



def get_max_replicon(refs):

    length = 0
    max_chrom = 0
    for chrom in refs:
        if len(refs[chrom]) > length:
            length = len(refs[chrom])
            max_chrom = chrom

    return max_chrom
        


def _get_shifts(MOTIF_LEN):

    left_shift = int(np.floor(MOTIF_LEN/2))
    right_shift = int(np.ceil(MOTIF_LEN/2))

    return left_shift, right_shift


def parse_data(fast5dir, reference_file, target_chr='all', required_coverage=30, MOTIF_LEN=11, LONG_MOTIF_LEN=29):

    refs, reverse_refs = get_reference(reference_file, target_chr)

    coverages = {
        ref: np.zeros(len(refs[ref])) for ref in refs
    }
    
    rev_coverages = {
        ref: np.zeros(len(refs[ref])) for ref in reverse_refs
    }

    l, r = _get_shifts(MOTIF_LEN)

    long_l, long_r = _get_shifts(LONG_MOTIF_LEN)

    motifs = {}
    reverse_motifs = {}

    long_motifs = {}
    long_reverse_motifs = {}




    for ref in refs:
        motifs[ref] = {}
        reverse_motifs[ref] = {}
        
        long_motifs[ref] = {}
        long_reverse_motifs[ref] = {}


    files = [file for file in os.listdir(fast5dir) if '.fast5' in file]

    batch = 1

    max_chrom = get_max_replicon(refs)


    for f in files:
        print('Batch {} out of {}...'.format(batch, len(files)))
        
        batch += 1
        

        try:

            with h5py.File('{}/{}'.format(fast5dir, f), 'r', rdcc_nbytes=1024**3) as file:

                for i in tqdm(list(file.items()), leave=False, ncols=75):

                    
                    readname = i[0]
                    try:
                        trace = file['/{}/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'.format(readname)][:]

                    except KeyError:
                        continue
                        
                    chrom = file['/{}/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'.format(readname)].attrs['mapped_chrom']  
                    
                    if chrom not in motifs:
                        continue
                    seq = [t[4].decode() for t in trace]

                    str_seq = ''.join(seq).upper()


                    f = refs[chrom].find(str_seq)


                    if f != -1:

                        for i in range(l, len(seq)-r):
                            context = str_seq[i-l:i+r]
                            
                            if context not in motifs[chrom]:
                                motifs[chrom][context] = []
                                
                            motifs[chrom][context].append(trace[i][0])

                        
                        for i in range(long_l, len(seq) - long_r):
                            long_context = str_seq[i-long_l:i+long_r]

                            if long_context not in long_motifs[chrom]:
                                long_motifs[chrom][long_context] = []

                            long_motifs[chrom][long_context].append(trace[i][0])

                        if chrom == max_chrom:
                            coverages[chrom][f:f+len(seq)] += 1

                        continue


                    
                
                    f_reverse = reverse_refs[chrom].find(str_seq)
                    if f_reverse != -1:
                        

                        for i in range(l, len(seq)-r):
                            context = str_seq[i-l:i+r]

                            if context not in reverse_motifs[chrom]:
                                reverse_motifs[chrom][context] = []

                            reverse_motifs[chrom][context].append(trace[i][0])
                        

                        for i in range(long_l, len(seq) - long_r):
                            long_context = str_seq[i-long_l:i+long_r]

                            if long_context not in long_reverse_motifs[chrom]:
                                long_reverse_motifs[chrom][long_context] = []

                            long_reverse_motifs[chrom][long_context].append(trace[i][0])


                        
                        if chrom == max_chrom:
                            rev_coverages[chrom][f_reverse:f_reverse+len(seq)] += 1
                        continue
        except KeyboardInterrupt:
            import sys
            sys.exit()
            pass
        except:
            print('Invalid batch!')
            continue    
        
        
        current_forward_coverage = np.round(np.mean(coverages[max_chrom]), 2)
        current_reverse_coverage = np.round(np.mean(rev_coverages[max_chrom]), 2)

        
        if min(current_forward_coverage, current_reverse_coverage) > required_coverage:
            break

        print(f'Current forward coverage {current_forward_coverage}X ; reverse coverage {current_reverse_coverage}X')

    print(f'Final coverage depth: forward {current_forward_coverage}X ; reverse {current_reverse_coverage}X (with {required_coverage}X threshold)')
    return motifs, reverse_motifs, long_motifs, long_reverse_motifs, coverages, rev_coverages



    

