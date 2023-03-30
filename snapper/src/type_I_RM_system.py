import numpy as np
from itertools import product
from snapper.src.seq_processing import letter_codes, gen_variants
from scipy.stats import ks_2samp
from tqdm import tqdm
from random import sample
from scipy.stats import chi2_contingency
from tqdm import tqdm

bases = ['A', 'G', 'T', 'C', 'N']


MAXSAMPLESIZE = 200


def gen_template(motif_variant, pos_variant, lenmotif):
    
    template = ['N',]*lenmotif
    
    for i, pos in enumerate(pos_variant):
        template[pos] = motif_variant[i]
        
    template =  ''.join(template)

    return template


def get_delta(k_size, long_k_size):
    return int(long_k_size/2) - int(k_size/2)





#### THE PROBLEM IS HERERERERERE!!!!

def filter_long_kmers(target, long_kmers, lenmotif, long_k_size):
    

    
    delta = get_delta(lenmotif, long_k_size)
    
    
    
    if target[2][-1] + delta >= long_k_size:
        return []

    long_kmers = np.array([list(l) for l in long_kmers])


    target_seq = ''.join(target[1])

    target_variants = gen_variants(target_seq)




    filtered_long_kmers = []
    
    for target_variant in target_variants:

        current_long_kmers = long_kmers.copy()
    
        for i in range(len(target[1])):
            try:
                current_long_kmers = current_long_kmers[current_long_kmers[:, target[2][i] + delta] == target_variant[i]]
            except IndexError:
                continue

        filtered_long_kmers += [''.join(l) for l in current_long_kmers]

    return filtered_long_kmers






def check_for_completeness(
    motif, 
    sample_motifs, 
    control_motifs, 
    long_sample_motifs, 
    long_control_motifs,
    lenmotif, 
    long_k_size, 
    reference, 
    outputdir,
    log_threshold=2,
    long_motif_confidence=100_000,
    ):

    
    #return motif

    print(f'Checking for {"".join(motif[1])} motif completeness...')

    _sample = []
    _control = []


    ancMOTIF = gen_template(motif[1], motif[2], lenmotif)


    for MOTIF in gen_variants(ancMOTIF):
            
        
        if MOTIF not in sample_motifs or MOTIF not in control_motifs:
            continue
            
    
        _sample += sample_motifs[MOTIF]
        _control += control_motifs[MOTIF]

    # THINKKKKK!!!!

    if len(_sample) > MAXSAMPLESIZE:
        _sample = sample(_sample, MAXSAMPLESIZE)
    
    if len(_control) > MAXSAMPLESIZE:
        _control = sample(_control, MAXSAMPLESIZE)

    global_ks_results = ks_2samp(_sample, _control)[1]

    if -np.log10(global_ks_results) > log_threshold:
        print(f'Motif {"".join(motif[1])} seems complete (p-val = {global_ks_results})')
        print()
        return motif


    print(f'Motif {"".join(motif[1])} is probably incomplete (p-val = {global_ks_results}). Extending enrichment heuristics...')
    
    _sample = []
    _control = []

    significant_contexts = []

    filtered_long_motifs = filter_long_kmers(motif, list(long_control_motifs.keys()), lenmotif, long_k_size)

    print('Collecting extended contexts...')
    for MOTIF in tqdm(filtered_long_motifs):
        
        
        if MOTIF not in long_sample_motifs or MOTIF not in long_control_motifs:
            continue   
        _sample += long_sample_motifs[MOTIF]
        _control += long_control_motifs[MOTIF]

        if -np.log10(ks_2samp(long_sample_motifs[MOTIF], long_control_motifs[MOTIF])[1]) > log_threshold:
            significant_contexts.append(MOTIF)

    import os

    os.makedirs(outputdir)

    context_cnt = 1
    with open(outputdir + '/long_contexts.fasta', 'w') as f_contexts:
        for context in significant_contexts:
            f_contexts.write(f'>long_context_{context_cnt}\n{context}\n')
            context_cnt += 1


    delta = get_delta(lenmotif, long_k_size)
    trd1_pos = motif[2][0] + delta
    
    
    long_motif_veriants = find_possible_trd2(motif, significant_contexts, reference, trd1_pos, long_k_size)

    with open(outputdir + '/long_motif_variants.tsv', 'w') as f_variants:
        for motif in long_motif_veriants:
            confidence, possible_long_motif = motif

            f_variants.write(f'{possible_long_motif}\t{confidence}\n')


    print(f'{possible_long_motif} has shown the best statistics. All data have been saved to {outputdir}.')
    print()

    return motif
    
        
    


def create_long_motif_template(long_motif, trd1_pos, confidence):

    motif_template = []
    pos_template = []
    for i in range(len(long_motif)):

        pos_template.append(i + trd1_pos)
        motif_template.append(long_motif[i])


    while motif_template[0] == 'N':
        motif_template = motif_template[1:]
        pos_template = pos_template[1:]
    
    while motif_template[-1] == 'N':
        motif_template = motif_template[:-1]
        pos_template = pos_template[:-1]
    

    return (
        confidence, tuple(motif_template), tuple(pos_template),
    )




def filter_trd2(trd2_variants, min_nondegenerate_letters=2, max_N_letters=2):
    
    filtered_trd2_variants = []
    
    for v in trd2_variants:
        
        Ncnt = v.count('N')
        
        if v[0] == 'N':
            continue
            
        if v[-1] == 'N':
            continue
        
        if Ncnt > max_N_letters:
            continue
            
        if len(v) - Ncnt < min_nondegenerate_letters:
            continue
        
        filtered_trd2_variants.append(v)
    return filtered_trd2_variants
    



def generate_RM_type_I_templates(trd1, N_lens=(5,6,7,8), trd2_lens=(2,3,4,5,6)):
        
        
    trd2_variants = []

    for trd2_len in trd2_lens:
        trd2_variants += [''.join(v) for v in list(product(bases, repeat=trd2_len))]

    trd2_variants = filter_trd2(trd2_variants)
    
    
    N_variants = ['N'*N_len for N_len in N_lens]

    while trd1[0] == 'N':
        trd1 = trd1[1:]

    while trd1[-1] == 'N':
        trd1 = trd1[:-1]
    
    templates = list(product(gen_variants(trd1), N_variants,trd2_variants))
    
    templates = [''.join(t) for t in templates]
    return templates
    
def find_possible_trd2(trd1, significant_contexts, reference, trd1_pos, lenmotif, N_lens=(5,6,7,8), trd2_lens=(2,3,4,5,6)):



    
    trd1 = ''.join(trd1[1])

    seq_array = np.array([list(l) for l in significant_contexts])
    ref_array = list(set([reference[i:i + lenmotif]  for i in range(len(reference)-lenmotif)]))

    ref_array = np.array([list(l) for l in ref_array])

    N_ref = len(ref_array)
    filtered_ref_array = None
    for short_motif in gen_variants(trd1):
        tmp_ref_array = ref_array.copy()

        for pos in range(trd1_pos, trd1_pos + len(trd1)):
            #print(ref_array, pos, trd1, trd1_pos)
            tmp_ref_array = tmp_ref_array[tmp_ref_array[:,pos] == short_motif[pos-trd1_pos]]
        if filtered_ref_array is None:

            filtered_ref_array = tmp_ref_array
        else:
            filtered_ref_array = np.concatenate((filtered_ref_array, tmp_ref_array))
        
    ref_array = filtered_ref_array
    
    N_seqset = len(seq_array)

    templates = generate_RM_type_I_templates(trd1, N_lens=N_lens, trd2_lens=trd2_lens)

    trd2_testing_results = []
    
    # THINKKKKK!!!
    print('TRD2 sequence optimization...')

    for template in tqdm(templates):
        
        if trd1_pos + len(template) >= lenmotif:
            continue

        subseq = seq_array.copy()
        ref_subseq = ref_array.copy()
        
        for pos in range(trd1_pos, trd1_pos + len(template)):
            if template[pos - trd1_pos] == 'N':
                continue
            
            if len(subseq) != 0:
                subseq = subseq[subseq[:, pos+1] == template[pos - trd1_pos]]
            
            if len(ref_subseq) != 0:
                ref_subseq = ref_subseq[ref_subseq[:, pos+1] == template[pos - trd1_pos]]

        subseq_N = len(subseq)
        ref_subseq_N = len(ref_subseq)

        try:
            chi_res = chi2_contingency(
                [
                    [subseq_N, N_seqset-subseq_N],
                    [ref_subseq_N, N_ref-ref_subseq_N],
                ]
            )[0]
        except ValueError:
            chi_res = 0

        trd2_testing_results.append(
            (
                chi_res, template
            )
        )
        

    trd2_testing_results.sort(reverse=True)

    print('Best results:')
    for t in trd2_testing_results[:10]:
        print(t[0], t[1])


    return trd2_testing_results[:20]



                
            
