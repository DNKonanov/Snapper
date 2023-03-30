import numpy as np
from itertools import combinations, product
import re
from snapper.src.seq_processing import gen_variants, letter_codes_rev, letter_anticodes
from scipy.stats import chi2_contingency, mode
from tqdm import tqdm


from multiprocessing import Process, Manager



regular_letters = ['A','G','C','T']
non_regular_letters = ['M', 'R', 'W', 'S', 'Y','K', 'V', 'H', 'D','B']


def filter_pos_variants_l3(pos_variants):

    filtered_pos_variants = []

    for pos_variant in pos_variants:
        _p = sorted(pos_variant)
        
        if _p[-1] - _p[0] >= 6:
            continue
        
        if tuple(pos_variant) not in filtered_pos_variants:
            filtered_pos_variants.append(tuple(_p))
    
    return filtered_pos_variants



def filter_pos_variants(pos_variants):


    # custom filtering for pos_variants with length of 3
    if len(pos_variants[0]) == 3:
        return filter_pos_variants_l3(pos_variants)



    filtered_pos_variants = []
    for pos_variant in pos_variants:
        
        _p = sorted(pos_variant)
        
        if _p[-1] - _p[0] >= 6:
            continue
            
        filtered_pos_variants.append(_p)
    
    _2_filtered_pos_variants = []

    
    for pos_variant in filtered_pos_variants:
        
        #for i in range(1, len(pos_variant) - 1):
        #    if (pos_variant[i] - pos_variant[i-1] > 1) and (pos_variant[i+1] - pos_variant[i] > 1):
        #        continue
        #    
        if tuple(pos_variant) in _2_filtered_pos_variants:
            continue
        
        #if pos_variant[1] - pos_variant[0] > 1 or  pos_variant[-1] - pos_variant[-2] > 1:
        #    continue
        
        _2_filtered_pos_variants.append(tuple(pos_variant))
    
    
    return _2_filtered_pos_variants
    

    
def filter_motifs(motif_variants):
    filtered_motifs  = []
    
    for motif in motif_variants:
        if 'C' not in motif and 'A' not in motif:
            continue
        filtered_motifs.append(motif)
        
    return filtered_motifs




def extract_template_subset(pos_variant, motif_variant, seq_array):
    subseq = seq_array
    for i in range(len(pos_variant)):
        if motif_variant[i] == '.':
            continue

        subseq = subseq[subseq[:,pos_variant[i]] == motif_variant[i]]

    return subseq


def extract_template_count(pos_variant, motif_variant, seq_array):
    
    subseq = seq_array
    for i in range(len(pos_variant)):
        subseq = subseq[subseq[:,pos_variant[i]] == motif_variant[i]]
        
    return len(subseq)


def gen_regexp_template(motif_variant, pos_variant, length=6):
    
    template = ['.',]*length
    
    base_pos = pos_variant[0]
    for i, pos in enumerate(pos_variant):
        template[pos-base_pos] = motif_variant[i]
        
    return ''.join(template)

def normalized_variation(array):
    return np.std(array)/np.mean(array)


def local_filter_seqs(seqs, pos_variant, motif_variant):
    
    new_seqs = []
    template = ''.join(motif_variant)
    template = template.replace('.', 'N')

    template_subvariants = gen_variants(template)
    
    for s in seqs:
        
        str_vec = ''.join([s[i] for i in pos_variant])
        if str_vec in template_subvariants:
            continue
    
        new_seqs.append(s)
            
    return new_seqs


def modify_seq(seq, pos, target_letter):

    newseq = list(seq)
    newseq[pos] = target_letter

    return ''.join(newseq)


def generate_reference_freqs_parallel(seq_array, batch, dict_per_length):

    for pos_variant, motif_variant in batch:
        variant_count = extract_template_count(pos_variant, motif_variant, seq_array)
        dict_per_length[(motif_variant, pos_variant)] = variant_count
     
    
def generate_reference_freqs(reference, length, threads, lengths=(4,5,6)):

    variants_counter = {}

    seqs = list(set([
        reference[i:i+length] for i in range(len(reference) - length)
    ]))

    seq_array = np.array([list(s) for s in seqs])

    print(len(seq_array))
    for LENGTH in lengths: 

        print('Reference indexing with length of {}...'.format(LENGTH))

        manager = Manager()
        dict_per_length = manager.dict()

        pos_variants = list(combinations(range(0,length), r=LENGTH))
        pos_variants = filter_pos_variants(pos_variants)

        motif_variants = list(product(regular_letters, repeat=LENGTH))
        motif_variants = filter_motifs(motif_variants)
        
        batch_len = len(pos_variants)*len(motif_variants)//threads
        
        processes = [] #all processes
        for i in range(threads+1):
            try:
                batch = list(product(pos_variants, motif_variants))[(i)*batch_len:(i+1)*batch_len]
            except IndexError: 
                batch = list(product(pos_variants, motif_variants))[(i)*batch_len:]
            p = Process(target=generate_reference_freqs_parallel, 
                                 args = (seq_array, batch, dict_per_length,))
            
            processes.append(p)
            p.start()

        #join processes    
        [p.join() for p in processes]
        
        variants_counter[LENGTH] = dict(dict_per_length)
              
    return variants_counter, len(seq_array)




def add_N(motif):
    
    if motif[0] != 'N':
        motif = 'N' + motif
    
    if motif[-1] != 'N':
        motif += 'N'     
    
    return motif
    

def is_superset(motif1, motif2, edgelength=2):
    
    motif1 = add_N(motif1)
    motif2 = add_N(motif2)
    
    if len(motif2) <= len(motif1):
        extended_motif1 = motif1

    else:
        extended_motif1 = 'N' * edgelength + motif1 + 'N' * edgelength
    
    motif1_variants = gen_variants(extended_motif1)
    motif2_variatns = gen_variants(motif2)
    
    global_match = True
    for variant2 in motif2_variatns:
        match = False
        for variant1 in motif1_variants:
            if variant2 in variant1:
                match = True
                break
        if match == False:
            global_match = False
            break

    
    
    return global_match


def get_alternate_variants(motif_variant, lenmotif=11, range_of_filtering=5):

    seq_variant, pos_variant  = motif_variant[1], motif_variant[2]

    while seq_variant[0] == 'N':
        seq_variant = seq_variant[1:]
        pos_variant = pos_variant[1:]

    while seq_variant[-1] == 'N':
        seq_variant = seq_variant[:-1]
        pos_variant = pos_variant[:-1]
        
    
    alternate_variants = []


    for i in range(
        max(0,        pos_variant[0]  - range_of_filtering), 
        min(lenmotif, pos_variant[-1] + range_of_filtering)
    ):
        
        shift = i - pos_variant[0]
        
        pos_alternate = tuple(j+shift for j in pos_variant)
        if pos_alternate[-1] >= lenmotif:
            break
            
        alternate_variants.append((motif_variant[0], seq_variant, pos_alternate))


    return alternate_variants

            



def is_subset(motif1, motif2, edgelength=2):
    return is_superset(motif2, motif1, edgelength=edgelength)


def variant_counts_parallel(seq_array, ref_motifs_counter, N_REF, batch, LENGTH, total_variants_counter_list):
    variants_counter_list = [] 
    N_VARIANT = len(seq_array)
    for pos_variant, motif_variant in batch:
        try:
            reference_count = ref_motifs_counter[LENGTH][(motif_variant, pos_variant)]
                
        except KeyError:
            variants_counter_list.append((0, motif_variant, pos_variant))
            
        else:
            variant_count = extract_template_count(pos_variant, motif_variant, seq_array)

            
            if variant_count == 0 and reference_count == 0:
                    variants_counter_list.append((0, motif_variant, pos_variant))

            else:
                chi2_result = chi2_contingency(
                    [
                        [variant_count, N_VARIANT-variant_count],
                        [reference_count, N_REF-reference_count],
                    ]
                )

                # chi2_log_pval = -np.log10(chi2_result[1])
                chi2_statistic = chi2_result[0]

                variants_counter_list.append((chi2_statistic, motif_variant, pos_variant))

    total_variants_counter_list+=variants_counter_list

               
def collect_variant_counts(seq_array, ref_motifs_counter, N_REF, threads, lengths=(4,5,6), lenmotif=11):
    merged_variants_counter_list = []
    
    for LENGTH in lengths:

        print('\tOBSERVING ANCHOR MOTIFS WITH LENGTH OF', LENGTH)

        pos_variants = list(combinations(range(0,lenmotif), r=LENGTH))
        pos_variants = filter_pos_variants(pos_variants)

        motif_variants = list(product(regular_letters, repeat=LENGTH))
        motif_variants = filter_motifs(motif_variants)
        
        #create batch
        batch_len = len(pos_variants)*len(motif_variants)//threads
        total_variants_counter_list = Manager().list() #for all outputs
        
        processes = [] #all processes
        args_list = list(product(pos_variants, motif_variants))

        
        for i in range(threads+1):
            try:
                batch = args_list[i*batch_len:(i+1)*batch_len]
            except IndexError: 
                batch = args_list[i*batch_len:]

            p = Process(target=variant_counts_parallel, 
                                 args = (seq_array, ref_motifs_counter, N_REF, batch, LENGTH, total_variants_counter_list))
            
            processes.append(p)
            p.start()
        
        [p.join() for p in processes]
        
        merged_variants_counter_list+=list(total_variants_counter_list) # add to 

    merged_variants_counter_list.sort(reverse=True)

        
    return merged_variants_counter_list


def get_significant_letters(sub_seq_array, top_variant, pos, reference, threshold_ratio):

    print('\tLocal motif adjustment...')

    reference_letter_freqs = {'A':0, 'G':0, 'T':0, 'C':0}
    variant_subset_letter_freqs = {'A':0, 'G':0, 'T':0, 'C':0}
    ref_vs_variant_ratios = {'A':0, 'G':0, 'T':0, 'C':0}

    variant_length = (top_variant[2][-1] - top_variant[2][0] + 1)
    re_variant = gen_regexp_template(top_variant[1], top_variant[2], length=variant_length)

    pos_letters = sub_seq_array[:,pos]

    for letter in reference_letter_freqs:
        re_variant_mod = modify_seq(re_variant, pos-top_variant[2][0], letter)
        ref_letter_count = len(re.findall(re_variant_mod, reference))
        
        variant_subset_letter_count = len(pos_letters[pos_letters == letter])

        reference_letter_freqs[letter] += ref_letter_count
        variant_subset_letter_freqs[letter] +=  variant_subset_letter_count

    
    list_variant_letter_freqs = [
        (variant_subset_letter_freqs[k], k) for k in variant_subset_letter_freqs
    ]

    list_variant_letter_freqs.sort(reverse=True)

    # consider the first letter to be presented apriori
    the_first_letter = list_variant_letter_freqs[0][1]

    ref_vs_variant_ratios[the_first_letter] = 1

    significant_letters = set([the_first_letter])

    for record in list_variant_letter_freqs[1:]:

        try:
            ref_letter_ratio = reference_letter_freqs[the_first_letter]/reference_letter_freqs[record[1]]
        except ZeroDivisionError:
            ref_letter_ratio = np.inf

        try:
            variant_subset_letter_ratio = variant_subset_letter_freqs[the_first_letter]/variant_subset_letter_freqs[record[1]]
        except ZeroDivisionError:
            variant_subset_letter_ratio = np.inf

        ref_vs_variant_ratio = variant_subset_letter_ratio/ref_letter_ratio

        ref_vs_variant_ratios[record[1]] = round(ref_vs_variant_ratio, 4)

        if ref_vs_variant_ratio > threshold_ratio:
            break

        significant_letters.add(record[1])


    return tuple(sorted(list(significant_letters)))


def adjust_letter(seq_array, top_variant, pos, reference, threshold_ratio=5):

    sub_seq_array = extract_template_subset(top_variant[2], top_variant[1], seq_array)

    pos_letters = get_significant_letters(sub_seq_array, top_variant, pos, reference, threshold_ratio=threshold_ratio) 

    return letter_codes_rev[pos_letters]


def change_subset_motif(supermotif, submotif, edgelength=2):
    
    extended_supermotif = 'N'*edgelength + ''.join(supermotif[1]) + 'N'*edgelength
    
    super_variants = gen_variants(extended_supermotif)
    sub_variants = gen_variants(''.join(submotif[1]))
    
    shifts = []  
    for subvariant in sub_variants:
        for supervariant in super_variants:
            if subvariant in supervariant:
                shift = edgelength - supervariant.find(subvariant)
                shifts.append(shift)
    
    shift = mode(shifts).mode[0]

    left_pos = max(0, submotif[2][0] + shift)
    right_pos = min(11, submotif[2][0] + shift + len(supermotif[2]))


    # check left edge case 
    if shift < 0:
        adjusted_subvariant = (
            submotif[0],
            supermotif[1][-shift:],
            tuple(range(submotif[2][0], submotif[2][0] + len(supermotif[1][-shift:])))
        )

    
    # check rigth edge case
    elif submotif[1][-1] in regular_letters and submotif[2][-1] == 10 and supermotif[1][-1] == 'N':
        adjusted_subvariant = (
            submotif[0],
            supermotif[1][:-1],
            tuple(range(left_pos, 11))
        )

    # common case
    else:
        adjusted_subvariant = (
            submotif[0],
            supermotif[1],
            tuple(range(left_pos, right_pos))
        )

    
    # just a patch, must be formalized!!
    if len(adjusted_subvariant[1]) != len(adjusted_subvariant[2]):
        adjusted_subvariant = [
            submotif[0],
            supermotif[1],
            tuple(range(left_pos, left_pos + len(supermotif[1])))
        ]  
        while adjusted_subvariant[2][-1] > 10:
            adjusted_subvariant[1] = adjusted_subvariant[1][:-1]
            adjusted_subvariant[2] = adjusted_subvariant[2][:-1]
    
    return tuple(adjusted_subvariant)



def extend_template(top_variant, maxlength=11):

    extended_top_variant = [ top_variant[0], list(top_variant[1]), list(top_variant[2])]

    if top_variant[2][0] != 0:
        extended_top_variant[2] = [extended_top_variant[2][0] - 1] + extended_top_variant[2]
        extended_top_variant[1] = ['.'] + extended_top_variant[1]

    if top_variant[2][-1] != maxlength-1:
        extended_top_variant[2] = extended_top_variant[2] + [extended_top_variant[2][-1] + 1]
        extended_top_variant[1] = extended_top_variant[1] + ['.']

    
    variant_length = (extended_top_variant[2][-1] - extended_top_variant[2][0] + 1)
    re_variant = gen_regexp_template(extended_top_variant[1], extended_top_variant[2], length=variant_length)

    extended_top_variant = (
        top_variant[0],
        tuple(re_variant), 
        list(range(extended_top_variant[2][0], extended_top_variant[2][-1] + 1))
    )

    return extended_top_variant



def save_results (motifs, out_fasta):

    with open(out_fasta, 'w') as f:
        
        cnt = 1

        for m in motifs:
            f.write('>MOTIF_{} conflevel={}\n{}\n'.format(cnt, m[0], ''.join(m[1])))
            cnt += 1


def save_k_mers (motifs, out_fasta):
    with open(out_fasta, 'w') as f:
        
        cnt = 1

        for m in motifs:
            f.write('>MOTIF_{}\n{}\n'.format(cnt, m))
            cnt += 1
