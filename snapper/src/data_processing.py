from multiprocessing import Pool
from unittest import result
from webbrowser import get
import numpy as np
import h5py
from itertools import product
from Bio.SeqIO import parse
from Bio.Seq import reverse_complement
import os

letters = ['A','G','C','T']


TARGET_CHRS = ['all', 'divided', 'longest', 'specified']


def _get_motifs(k = 6):

    return [''.join(s) for s in list(product(letters, repeat=k))]




def get_reference(reference_file, target_chr):

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



def parse_data(fast5dir, reference_file, target_chr='all', n_batches=np.inf, MOTIF_LEN=6):

    refs, reverse_refs = get_reference(reference_file, target_chr)
    
    motifs = {}
    shifts = {}

    reverse_motifs = {}
    reverse_shifts = {}

    for ref in refs:
        motifs[ref] = {
            motif: [] for motif in _get_motifs(k = MOTIF_LEN)
        }
        reverse_motifs[ref] = {
            motif: [] for motif in _get_motifs(k = MOTIF_LEN)
        }

        shifts[ref] = {
            i: [] for i in range(len(refs[ref]))
        }
        reverse_shifts[ref] = {
            i: [] for i in range(len(reverse_refs[ref]))
        }

    files = [file for file in os.listdir(fast5dir) if '.fast5' in file]

    batch = 1
    for f in files:

        if batch > n_batches:
            break
        print('Batch {} out of {}...'.format(batch, min(n_batches, len(files))))
        batch += 1
        
        with h5py.File('{}/{}'.format(fast5dir, f), 'r', rdcc_nbytes=1024**3) as file:

            for i in list(file.items()):
                
                readname = i[0]
                try:
                    trace = file['/{}/Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'.format(readname)][:]

                except KeyError:
                    continue
                    
                chrom = file['/{}/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment'.format(readname)].attrs['mapped_chrom']  
                 
                    
                seq = [t[4].decode() for t in trace]

                str_seq = ''.join(seq).upper()

                
            
                f = refs[chrom].find(str_seq)
                
                if f != -1:

                    for i in range(3, len(seq) - 3):
                        context = str_seq[i-3:i+3]

                        shifts[chrom][f + i].append(trace[i][0])

                        motifs[chrom][context].append(trace[i][0])

                    continue
                
            
                f_reverse = reverse_refs[chrom].find(str_seq)
                
                if f_reverse != -1:
                    

                    for i in range(3, len(seq) - 3):
                        context = str_seq[i-3:i+3]

                        reverse_shifts[chrom][f + i].append(trace[i][0])

                        reverse_motifs[chrom][context].append(trace[i][0])


    return shifts, reverse_shifts, motifs, reverse_motifs



    

