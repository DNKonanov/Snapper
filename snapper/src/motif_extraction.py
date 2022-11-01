import numpy as np
import os
from Bio.SeqIO import parse
from pickle import dump, load
from snapper.src.methods import collect_variant_counts, is_superset, is_subset, local_filter_seqs, adjust_letter, extend_template, generate_reference_freqs, change_subset_motif
from snapper.src.methods import get_alternate_variants


def extract_motifs(
    seqs, 
    reference, 
    savepath, 
    max_motifs,
    min_conf, 
    contig_name,
    threads=10):


    print()
    print('Motif enrichment')
    print()


    N_REF = len(set(
        [reference[i:i+11] for i in range(len(reference) - 11)]
    ))

    lengths = [3,4,5,6]

    print('Reference indexing...')
    ref_motifs_counter, N_REF = generate_reference_freqs(reference, 11, threads, lengths=lengths)





    ITERATION = 1


    new_seqs = seqs.copy()
    try:
        os.mkdir(savepath + '/seq_iter')

    except:
        pass


    os.mkdir(savepath + '/seq_iter/{}/'.format(contig_name))

    with open(savepath + '/seq_iter/{}/seqs_iter_{}.fasta'.format(contig_name, ITERATION), 'w') as fseqiter:

        for seq in new_seqs:
            fseqiter.write('>')
            fseqiter.write(seq)
            fseqiter.write('\n')
            fseqiter.write(seq)
            fseqiter.write('\n')

    seq_array = np.array([list(s) for s in new_seqs])

    initial_seq_array = seq_array.copy()


    MOTIFS_SET = []
    DETAILED_MOTIF_SET = []

    print('ITERATION 1 ({} unexplained 11mers):'.format(len(seq_array)))




    variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, threads=threads, lengths=lengths)


    ITERATION = 2
    while variants_counter_list[0][0] > min_conf:

        
        
        top_variant = variants_counter_list[0]

        extended_top_variant = extend_template(top_variant, maxlength=11)

        positions_to_adjust = []

        for i, pos in enumerate(extended_top_variant[2]):
            if extended_top_variant[1][i] == '.':
                positions_to_adjust.append((pos, i))
        
        modifiable_extended_top_variant = [
            extended_top_variant[0],
            list(extended_top_variant[1]),
            list(extended_top_variant[2])
        ]

        
        for pos in positions_to_adjust:

            adjusted_pos_letter = adjust_letter(initial_seq_array, extended_top_variant, pos[0], reference)
            modifiable_extended_top_variant[1][pos[1]] = adjusted_pos_letter

        extended_top_variant = (
            extended_top_variant[0],
            tuple(modifiable_extended_top_variant[1]),
            tuple(modifiable_extended_top_variant[2]),
        )

        print(extended_top_variant)

        is_superset_check = False
        is_subset_check = False

        for i, motif in enumerate(MOTIFS_SET):
            is_superset_check = is_superset(motif, ''.join(extended_top_variant[1]))
            is_subset_check = is_subset(motif, ''.join(extended_top_variant[1]))

            if is_subset_check:
                break

            if is_superset_check:
                break

        
        if is_superset_check == False:
            MOTIFS_SET.append(''.join(extended_top_variant[1]))
            DETAILED_MOTIF_SET.append(extended_top_variant)
        
        
        else:
            print('{} already has a supermotif!'.format(extended_top_variant))
            
            extended_top_variant = change_subset_motif(
                DETAILED_MOTIF_SET[i],
                extended_top_variant,
                edgelength=2
            )
            print('Changed to {}'.format(extended_top_variant))
            

        if len(MOTIFS_SET) == max_motifs:
            break


        print(MOTIFS_SET)

        alternate_variants = get_alternate_variants(extended_top_variant)

        print('Filtering seq_set...')

        n_seqs = len(new_seqs)
        
        for variant in alternate_variants:
            if variant[0] > min_conf:

                new_seqs = local_filter_seqs(new_seqs, variant[2], variant[1])
        

        # filter seq_set by top_variant to prevent infinite loop
        if len(new_seqs) == n_seqs:
            alternate_variants = get_alternate_variants(top_variant)    
            for variant in alternate_variants:
                if variant[0] > min_conf:

                    new_seqs = local_filter_seqs(new_seqs, variant[2], variant[1])
            
        
        
        with open(savepath + '/seq_iter/{}/seqs_iter_{}.fasta'.format(contig_name, ITERATION), 'w') as fseqiter:

            for seq in new_seqs:
                fseqiter.write('>')
                fseqiter.write(seq)
                fseqiter.write('\n')
                fseqiter.write(seq)
                fseqiter.write('\n')


        seq_array = np.array([list(s) for s in new_seqs])
        
        print('ITERATION {} ({} unexplained 11mers):'.format(ITERATION, len(seq_array)))
        ITERATION += 1
        
        variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, threads=threads, lengths=lengths)

    
    return DETAILED_MOTIF_SET