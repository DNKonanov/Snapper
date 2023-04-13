import numpy as np
import os
from Bio.SeqIO import parse
from pickle import dump, load
from snapper.src.methods import collect_variant_counts, is_superset, is_subset, local_filter_seqs, adjust_letter, extend_template, generate_reference_freqs, change_subset_motif
from snapper.src.methods import get_alternate_variants
from snapper.src.type_I_RM_system import check_for_completeness


def extract_motifs(
    seqs, 
    reference, 
    savepath, 
    max_motifs,
    min_conf, 
    contig_name,

    
    sample_motifs, 
    control_motifs, 
    sample_long_motifs, 
    control_long_motifs, 
    k_size, 
    long_k_size,  
    ks_t,

    threads=10,
    lenmotif=11,

    ):


    print()
    print('Motif enrichment')
    print()


    N_REF = len(set(
        [reference[i:i+lenmotif] for i in range(len(reference) - lenmotif)]
    ))

    lengths = [3,4,5,6]

    print('Reference indexing...')
    ref_motifs_counter, N_REF = generate_reference_freqs(reference, lenmotif, threads, lengths=lengths)





    ITERATION = 1


    new_seqs = seqs.copy()
    try:
        os.mkdir(savepath + '/seq_iter')
        os.mkdir(savepath + '/motif_refine')

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

    print(f'ITERATION 1 ({len(seq_array)} unexplained {lenmotif}-mers):')




    variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, threads=threads, lengths=lengths, lenmotif=lenmotif)


    ITERATION = 2
    while variants_counter_list[0][0] > min_conf and len(seq_array) > 0:
        
        for v in variants_counter_list[:15]:
            print('\t', v)

        top_variant = variants_counter_list[0]

        extended_top_variant = extend_template(top_variant, maxlength=lenmotif)

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

        refine_outdir = f'{savepath}/motif_refine/{contig_name}/{"".join(extended_top_variant[1])}' 
        
        complete_motif = check_for_completeness(
                extended_top_variant, 
                sample_motifs, 
                control_motifs, 
                sample_long_motifs, 
                control_long_motifs, 
                k_size, 
                long_k_size, 
                reference, 
                outputdir=refine_outdir,
                log_threshold=ks_t
                )
        

        alternate_variants = get_alternate_variants(extended_top_variant, lenmotif=lenmotif)

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

        
        MOTIFS_SET.append(''.join(extended_top_variant[1]))
        DETAILED_MOTIF_SET.append(extended_top_variant)
            
        
        
        print(MOTIFS_SET)

        with open(savepath + '/seq_iter/{}/seqs_iter_{}.fasta'.format(contig_name, ITERATION), 'w') as fseqiter:

            for seq in new_seqs:
                fseqiter.write('>')
                fseqiter.write(seq)
                fseqiter.write('\n')
                fseqiter.write(seq)
                fseqiter.write('\n')

        
        if len(MOTIFS_SET) == max_motifs:
            break

        seq_array = np.array([list(s) for s in new_seqs])
        
        print(f'ITERATION {ITERATION} ({len(seq_array)} unexplained {lenmotif}-mers):')
        ITERATION += 1
        
        variants_counter_list = collect_variant_counts(seq_array, ref_motifs_counter, N_REF, threads=threads, lengths=lengths)

    
    return DETAILED_MOTIF_SET