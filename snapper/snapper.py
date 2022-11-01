from argparse import ArgumentParser


import warnings
warnings.filterwarnings("ignore")

import os
import sys




def main():

    parser = ArgumentParser()

    parser.add_argument('-sample_fast5dir', type=str, help='sample multi fast5 dir')
    parser.add_argument('-control_fast5dir', type=str, help='control multi fast54 dir')
    parser.add_argument('-reference', type=str, help='reference fasta')
    parser.add_argument('-ks_t', type=int, default=5, help='-log ks_test p-value (default 5)')
    parser.add_argument('-outdir', type=str, default='default', help='output directory name')
    parser.add_argument('-n_batches', type=int, default=100, help='number of parsed fast5 batches')
    parser.add_argument('-threads', type=int, default=8, help='number of threads used (derfault is 8)')
    parser.add_argument('-max_motifs', help='the maximum expected number of motifs extracted', default=20, type=int)
    parser.add_argument('-min_conf', help='the minimal confidence value. Default is 1000', type=float, default=1000)
    parser.add_argument('-target_chr', help='target chromosome name. By default is "all"', type=str, default='all')
    


    from snapper.src.motif_extraction import extract_motifs
    from snapper.src.plotting import plot_motif
    from snapper.src.data_processing import get_reference, parse_data
    from snapper.src.statistics import get_difsignals, get_statistics
    from snapper.src.methods import save_results, save_k_mers
    from snapper.src.statistics import SAMPLESIZE, MINSAMPLESIZE

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()


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

    sample_motifs, sample_reverse_motifs = parse_data(
        args.sample_fast5dir, 
        args.reference, 
        target_chr=args.target_chr, 
        n_batches=args.n_batches, 
    )


    print('\nControl data collecting...')
    control_motifs, control_reverse_motifs = parse_data(
        args.control_fast5dir, 
        args.reference, 
        target_chr=args.target_chr, 
        n_batches=args.n_batches, 
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
    )


    print('\nReverse strand signals processing...')
    reverse_motifs_lines, reverse_ks_stat_lines = get_statistics(
        sample_reverse_motifs, 
        control_reverse_motifs, 
        maxsamplesize=SAMPLESIZE,
        minsamplesize=MINSAMPLESIZE,
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
                                threads=args.threads
                                )


        for motif in motifs:
            plot_motif(motif, sample_motifs[contig], control_motifs[contig], plotdir)



        save_results(motifs, outdir + '/final_motifs_forward_{}.fasta'.format(contig))


    for contig in reverse_motifs_lines:

        print(contig, len(reverse_refs[contig]))
        print('Processing reversed motifs {}...'.format(contig))
            

        contig_passed_motifs = get_difsignals(
            reverse_motifs_lines[contig], 
            reverse_ks_stat_lines[contig], 
            log10_pval_thr = args.ks_t,
        )
        

        plotdir = outdir + '/plots_reverse_{}'.format(contig) 
        os.mkdir(plotdir)

        save_k_mers(contig_passed_motifs, outdir + '/passed_motifs_reverse_{}.fasta'.format(contig))
        motifs = extract_motifs(contig_passed_motifs, 
                                reverse_refs[contig], 
                                outdir, 
                                args.max_motifs,
                                args.min_conf, 
                                'reverse_' + contig,
                                threads=args.threads
                                )


        for motif in motifs:
            plot_motif(motif, sample_motifs[contig], control_motifs[contig], plotdir)



        save_results(motifs, outdir + '/final_motifs_reverse_{}.fasta'.format(contig))



    print('Done!')


if __name__ == '__main__':
    main()