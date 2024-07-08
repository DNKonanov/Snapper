from argparse import ArgumentParser
import os
import sys

import warnings
warnings.filterwarnings("ignore")




def main():

    parser = ArgumentParser()

    parser.add_argument('-sample_fast5dir', type=str, help='sample multi fast5 dir', required=True)
    parser.add_argument('-control_fast5dir', type=str, help='control multi fast5 dir', required=True)
    parser.add_argument('-reference', type=str, help='reference genome in the fasta format', required=True)
    parser.add_argument('-ks_t', type=int, default=3, help='-log ks_test p-value (default 3).')
    parser.add_argument('-outdir', type=str, default='default', help='output directory name')
    parser.add_argument('-coverage', type=float, help='minimal genome coverage depth (default 40)', default=40)
    parser.add_argument('-threads', type=int, default=8, help='number of threads used (default 8)')
    parser.add_argument('-k_size', type=int, default=15, help='k-mer size, must be odd, should not be less than 11 (default 15)')
    parser.add_argument('-long_k_size', type=int, default=29, help='long k-mer size, must be odd, should not be less than 21 (default 29)')
    parser.add_argument('-max_motifs', help='the maximum expected number of motifs extracted (default 20)', default=20, type=int)
    parser.add_argument('-min_conf', help='the minimal confidence value (default is 100)', type=float, default=100)
    parser.add_argument('-target_chr', help='target chromosome name (by default all contigs/replicons are considered)', type=str, default='all')
    


    from snapper.src.motif_extraction import extract_motifs
    from snapper.src.plotting import plot_motif, plot_coverage, plot_dist
    from snapper.src.data_processing import get_reference, parse_data
    from snapper.src.statistics_methods import get_difsignals, get_statistics
    from snapper.src.methods import save_results, save_k_mers
    from snapper.src.statistics_methods import SAMPLESIZE, MINSAMPLESIZE

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)


    

    args = parser.parse_args()

    if args.k_size%2 == 0 or args.long_k_size%2 == 0:
        raise ValueError('Both -k_size and -long_k_size must be odd numbers')
    
    if args.k_size < 11:
        raise ValueError('-k_size parameter should not be less than 11')

    if args.long_k_size < 21:
        raise ValueError('-long_k_size parameter should not be less than 21')
    
    if args.k_size >= args.long_k_size:
        raise ValueError('K_SIZE should be less than LONG_K_SIZE')
    



    if args.outdir == 'default':
        import datetime
        
        sp = str(datetime.datetime.now()
            ).replace(' ', '_').replace(':', '').replace('-', '_').split('.')[0]
        outdir = 'Results_' + sp

    else:
        outdir  = args.outdir

    try:
        os.mkdir(outdir)
    except:
        raise FileExistsError('The specified output dir already exists!')

    print('\nSample data collecting...')

    sample_motifs, sample_reverse_motifs, sample_long_motifs, sample_long_reverse_motifs, sample_coverages, sample_rev_coverages = parse_data(
        args.sample_fast5dir, 
        args.reference, 
        target_chr=args.target_chr, 
        required_coverage=args.coverage,
        MOTIF_LEN=args.k_size,
        LONG_MOTIF_LEN=args.long_k_size,
    )



    print('\nControl data collecting...')
    control_motifs, control_reverse_motifs, control_long_motifs, control_long_reverse_motifs, control_coverages, control_rev_coverages = parse_data(
        args.control_fast5dir, 
        args.reference, 
        target_chr=args.target_chr, 
        required_coverage=args.coverage,
        MOTIF_LEN=args.k_size,
        LONG_MOTIF_LEN=args.long_k_size,
    )


    refs, reverse_refs = get_reference(
        args.reference,
        target_chr=args.target_chr
    )

    for contig in refs:
        print(contig, len(refs[contig]))


    print('\nForward strand signals processing...')
    motifs_lines, ks_stat_lines = get_statistics(
        sample_motifs, 
        control_motifs, 
        maxsamplesize=SAMPLESIZE,
        minsamplesize=MINSAMPLESIZE,
        threads=args.threads
    )


    print('\nReverse strand signals processing...')
    reverse_motifs_lines, reverse_ks_stat_lines = get_statistics(
        sample_reverse_motifs, 
        control_reverse_motifs, 
        maxsamplesize=SAMPLESIZE,
        minsamplesize=MINSAMPLESIZE,
        threads=args.threads
    )

    


    # MOTIFS EXTRACTION

    for contig in motifs_lines:


        print('Processing forward motifs {}...'.format(contig))

            
        
        contig_passed_motifs = get_difsignals(
            motifs_lines[contig], 
            ks_stat_lines[contig], 
            log10_pval_thr = args.ks_t,
        )

        if len(contig_passed_motifs) < 100:
            print('---The number of k-mers is insufficient for the enrichment process. {} is skipped.---'.format(contig))
            continue


        



        plotdir = outdir + '/plots_forward_{}'.format(contig) 
        os.mkdir(plotdir)

        save_k_mers(contig_passed_motifs, outdir + '/passed_motifs_forward_{}.fasta'.format(contig))
        motifs = extract_motifs(contig_passed_motifs, 
                                refs[contig], 
                                outdir, 
                                args.max_motifs,
                                args.min_conf, 
                                'forward_' + contig,
                                
                                sample_motifs[contig], 
                                control_motifs[contig], 
                                sample_long_motifs[contig], 
                                control_long_motifs[contig], 
                                args.k_size, 
                                args.long_k_size,  
                                args.ks_t,

                                threads=args.threads,
                                lenmotif=args.k_size
                                )


        for motif in motifs:
            plot_dist(motif, sample_motifs[contig], control_motifs[contig], plotdir, lenmotif=args.k_size)



        save_results(motifs, outdir + '/final_motifs_forward_{}.fasta'.format(contig))
        plot_coverage(sample_coverages[contig], control_coverages[contig], contig, f'{outdir}/coverage_forward_{contig}.pdf')


    for contig in reverse_motifs_lines:

        print(contig, len(reverse_refs[contig]))
        print('Processing reversed motifs {}...'.format(contig))
            

        contig_passed_motifs = get_difsignals(
            reverse_motifs_lines[contig], 
            reverse_ks_stat_lines[contig], 
            log10_pval_thr = args.ks_t,
        )

        
        if len(contig_passed_motifs) < 100:
            print('---The number of k-mers is insufficient for the enrichment process. {}(reverse) is skipped.---'.format(contig))
            continue
        

        plotdir = outdir + '/plots_reverse_{}'.format(contig) 
        os.mkdir(plotdir)

        save_k_mers(contig_passed_motifs, outdir + '/passed_motifs_reverse_{}.fasta'.format(contig))
        motifs = extract_motifs(contig_passed_motifs, 
                                reverse_refs[contig], 
                                outdir, 
                                args.max_motifs,
                                args.min_conf, 
                                'reverse_' + contig,

                                sample_reverse_motifs[contig], 
                                control_reverse_motifs[contig], 
                                sample_long_reverse_motifs[contig], 
                                control_long_reverse_motifs[contig], 
                                args.k_size, 
                                args.long_k_size,  
                                args.ks_t,

                                threads=args.threads,
                                lenmotif=args.k_size
                                )

        


        for motif in motifs:
            plot_dist(motif, sample_reverse_motifs[contig], control_reverse_motifs[contig], plotdir, lenmotif=args.k_size)



        save_results(motifs, outdir + '/final_motifs_reverse_{}.fasta'.format(contig))
        plot_coverage(sample_rev_coverages[contig], control_rev_coverages[contig], contig, f'{outdir}/coverage_reverse_{contig}.pdf')


    print('Done!')


if __name__ == '__main__':
    main()