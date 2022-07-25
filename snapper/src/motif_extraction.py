from audioop import reverse
from email import iterators
from functools import reduce
from xxlimited import new
from Bio.SeqIO import parse
import numpy as np


letter_codes = {
    'A': set(['A']),
    'C': set(['C']),
    'G': set(['G']),
    'T': set(['T']),
    'M': set(['A','C']),
    'R': set(['A','G']),
    'W': set(['A','T']),
    'S': set(['C','G']),
    'Y': set(['C','T']),
    'K': set(['G','T']),
    'V': set(['A','C','G']),
    'H': set(['A','C','T']),
    'D': set(['A','G','T']),
    'B': set(['C','G','T']),
    'N': set(['A','C','G','T'])
}

letter_anticodes = {
    'A': set(['C', 'G', 'T']),
    'C': set(['A', 'G', 'T']),
    'G': set(['A', 'C', 'T']),
    'T': set(['A', 'C', 'G']),
    'M': set(['G','T']),
    'R': set(['C','T']),
    'W': set(['C','G']),
    'S': set(['A','T']),
    'Y': set(['C','T']),
    'K': set(['A','C']),
    'V': set(['T']),
    'H': set(['G']),
    'D': set(['C']),
    'B': set(['A']),
    'N': set([])
}


letter_codes_rev = {
    ('A',): 'A',
    ('C',): 'C',
    ('G',): 'G',
    ('T',): 'T',
    ('A', 'C'): 'M',
    ('A', 'G'): 'R',
    ('A', 'T'): 'W',
    ('C', 'G'): 'S',
    ('C', 'T'): 'Y',
    ('G', 'T'): 'K',
    ('A', 'C', 'G'): 'V',
    ('A', 'C', 'T'): 'H',
    ('A', 'G', 'T'): 'D',
    ('C', 'G', 'T'): 'B',
    ('A', 'C', 'G', 'T'): 'N'
}

regular_letters = set(['A','G','C','T'])


def _get_dist(seq1, seq2):
    
    dist = 0
    
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
            
    return dist


def get_seq_clusters(seqs, allowed_dist=2):

    if len(seqs) == 0:
        return []
    
    clusters = [[]]

    for seq in seqs:
    
        cluster_hit = False
        
        for cluster in clusters:
            
            hit = True
            for seq2 in cluster:
                if _get_dist(seq, seq2) > allowed_dist:
                    
                    hit = False
                    
                    break
                    
            if hit == True:
                cluster_hit = True
                cluster.append(seq)
        if cluster_hit == False:
            clusters.append([])
            clusters[-1].append(seq)

    for cluster in clusters:

        for seq in seqs:
            if seq in cluster:
                continue
            hit = True
            for seq2 in cluster:
                if _get_dist(seq, seq2) > allowed_dist:
                    hit = False
                    break

            if hit:
                cluster.append(seq)


    clusters.sort(key=len, reverse=True)
    return clusters




def get_seq_clusters(seqs, allowed_dist=2):

    if len(seqs) == 0:
        return []
    
    clusters = [[]]

    for seq in seqs:
    
        cluster_hit = False
        
        for cluster in clusters:
            
            hit = True
            for seq2 in cluster:
                if _get_dist(seq, seq2) > allowed_dist:
                    
                    hit = False
                    
                    break
                    
            if hit == True:
                cluster_hit = True

                current_cluster = list(cluster)

                cluster.append(seq)

        if cluster_hit:
            clusters.append(current_cluster)

        if cluster_hit == False:
            clusters.append([])
            clusters[-1].append(seq)

    for cluster in clusters:

        for seq in seqs:
            if seq in cluster:
                continue
            hit = True
            for seq2 in cluster:
                if _get_dist(seq, seq2) > allowed_dist:
                    hit = False
                    break

            if hit:
                cluster.append(seq)


    clusters = _remove_duplicated_clusters(clusters)

    clusters.sort(key=len, reverse=True)
    return clusters


def _remove_duplicated_clusters(clusters):

    new_clusters = []

    for c in clusters:
        vec_c = tuple(sorted(c))

        if vec_c not in new_clusters:
            new_clusters.append(vec_c)

    for i in range(len(new_clusters)):
        new_clusters[i] = list(new_clusters[i])

    return new_clusters




def get_consensus(cluster, min_lngth=4):

    consensus_var = [tuple(sorted(list(set([seq[i] for seq in cluster])))) for i in range(len(cluster[0]))]

    consensus = ''.join([letter_codes_rev[consensus_var[i]] for i in range(len(consensus_var))])

    cm = True
    while cm:
        
        cm = False
        
        if consensus[0] == 'N':
            consensus = consensus[1:]
            cm = True
        
        if consensus[-1] == 'N':
            consensus = consensus[:-1]
            cm = True
    
    
    if len(consensus) <= 3:
        return None

    return consensus

def update_clusters(clusters, update_seq):

    if update_seq is None:
        return clusters[1:]

    new_seqs = []
    for c in clusters:
        for seq in c:
            if _in(update_seq, seq) == False:
                new_seqs.append(seq)

    new_seqs = list(set(new_seqs))
    new_clusters = get_seq_clusters(new_seqs)

    return new_clusters


def _gen_variants(seq):
    variants = ['']
    
    for i in range(len(seq)):
        new_variants = []
        for l in letter_codes[seq[i]]:
            
            for v in variants:
                new_variants.append(v + l)
                
        variants = new_variants
            
    return variants


def _in(query, template):
    
    for v in _gen_variants(query):
        if v in template:
            return True
    return False


MOTIF_S_LETTERS = 4
def check_reduce(clusters, update_seqs):

    reduce_result = []

    for i in range(len(update_seqs)):
        update_seq = update_seqs[i]

        if update_seq is None:
            reduce_result.append((0, update_seq, i))
            continue

        _pass = 0

        for letter in update_seq:
            if letter in regular_letters:
                _pass += 1

        if _pass < MOTIF_S_LETTERS:
            continue


        cnt = 0
        for c in clusters:
            for seq in c:
                if _in(update_seq, seq):
                    cnt += 1

        reduce_result.append((cnt, update_seq, i))

    reduce_result.sort(reverse=True)

    if len(reduce_result) == 0:
        return (None, None, None)

    return reduce_result[0]


TOP_CLUSTERS = 5

def extract_motifs(seqs):

   
    motifs = []



    seqs.sort()
    clusters = get_seq_clusters(seqs)

    iteration = 1

    while len(clusters) > 0:
        print('\n-------Iteration {}-------'.format(iteration))



        iteration += 1
        

        putative_consensuses = []

        for i in range(min(len(clusters), TOP_CLUSTERS)):


            putative_consensuses.append(get_consensus(clusters[i]))

        _, top_consensus, idx = check_reduce(clusters, putative_consensuses)

        print('Motifs cluster: {}'.format(clusters[idx]))
        print('Extracted motif: {}'.format(top_consensus))
        print()

        if top_consensus is None:
            break

        clusters = update_clusters(clusters, top_consensus)
        motifs.append(top_consensus)

    return motifs


def get_seqs(file):

    fasta_m = parse(file, format='fasta')

    seqs = [str(rec.seq) for rec in fasta_m]

    return seqs


def merge_motifs(motifs):

    for i in range(len(motifs)):

        m1 = motifs[i]

        for j in range(i + 1, len(motifs)):
            m2 = motifs[j]



            if len(m1) != len(m2):
                continue
            if _get_dist(m1, m2) < 2:
                print(m1, m2)




def crop_variant(variant):

    new_variant = variant
    for i in range(len(variant)):
        if variant[i] not in regular_letters:
            new_variant = new_variant[1:]
        else:
            break

    for i in range(len(variant)-1, -1, -1):
        if variant[i] not in regular_letters:
            new_variant = new_variant[:-1]
        else:
            break
    
    
    return new_variant



def collapse_motifs(motifs):

    print('Post-alignment...')

    ref_motifs = motifs
    new_motifs = []

    motifs_to_drop = []

    for motif in motifs:
        for letter in motif:
            if letter in regular_letters:
                continue
            for var in letter_anticodes[letter]:
                
                variant = motif.replace(letter, var)
                
                for ref_motif in ref_motifs:
                    if variant in ref_motif or ref_motif in variant:
                        
                        cropped_motif = crop_variant(motif)
                        
                        new_motifs.append(cropped_motif)
                        print('{} was merged to {}'.format(ref_motif, cropped_motif))

                        motifs_to_drop += [motif, ref_motif]

    
    filtered_motifs = []

    for m in ref_motifs:
        if m in motifs_to_drop:
            continue
        filtered_motifs.append(m)

    filtered_motifs += new_motifs

    return filtered_motifs                     
 
