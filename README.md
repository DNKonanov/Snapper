# Snapper: nanopore-based modification motifs caller

This tool is designed to perform a high-sensitive search for modification sites using ONT sequencing data.
This tool provide a new two-sided and balanced approch to compute statistic for each motif which is likely to be modified.

## Dependencies
- ont-tombo
- h5py
- biopython


## Usage

Firstly, fast5 files should be resquiggled using [Tombo](https://github.com/nanoporetech/tombo) software. 
After resquiggling, fast5 files should be converted to multi-fast5 format.

```
usage: snapper [-h] [-sample_fast5dir SAMPLE_FAST5DIR] [-control_fast5dir CONTROL_FAST5DIR] [-reference REFERENCE] [-ks_t KS_T] [-eff_size EFF_SIZE] [-outdir OUTDIR] [-n_batches N_BATCHES]
                      [-n_threads N_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -sample_fast5dir SAMPLE_FAST5DIR
                        sample multi fast5 dir
  -control_fast5dir CONTROL_FAST5DIR
                        control multi fast54 dir
  -reference REFERENCE  reference fasta
  -ks_t KS_T            -log ks_test p-value (default 50)
  -eff_size EFF_SIZE    Cohen d-effect size (default 0.25)
  -outdir OUTDIR        output directory name (default Results_%yyyy_%mm_%dd_%hhmmss)
  -n_batches N_BATCHES  number of parsed fast5 batches (default all)
  -n_threads N_THREADS  number of threads used (default 8)

```


NB! Depending on the data size and the number of threds used, the process can consume up to 1Tb of RAM.
The authors recommend to specify the same values for `-n_batches` and `-n_threads`.

## Output explanation

The output directory contains the following files:
- `passed_motifs_[contig_name].fasta` - all k-mers defined as modifed
- `final_motifs_[contig_name].fasta` - optimal set of motifs obtained from the passed motifs with greedy algorithm
- `plots_[contig_name]` - signal distributions for all motifs detected 
