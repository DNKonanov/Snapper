# Snapper: nanopore-based modification motifs caller

This tool is designed to perform a high-sensitive search for modification sites using ONT sequencing data.
This tool provide a new two-sided and balanced approch to compute statistic for each motif which is likely to be modified.

## Dependencies
- python 3.7 (later versions might be incompatible because of inner biopython dependencies)
- ont-tombo
- h5py
- biopython
- matplotlib
- scipy

## Installation

```pip install snapper-ont```

## Usage

Firstly, fast5 files should be resquiggled using [Tombo](https://github.com/nanoporetech/tombo) software. 
After resquiggling, fast5 files should be converted to multi-fast5 format using [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api).

```
usage: snapper [-h] [-sample_fast5dir SAMPLE_FAST5DIR] [-control_fast5dir CONTROL_FAST5DIR] [-reference REFERENCE] [-ks_t KS_T]
               [-outdir OUTDIR] [-n_batches N_BATCHES] [-threads THREADS] [-max_motifs MAX_MOTIFS] [-min_conf MIN_CONF]

optional arguments:
  -h, --help            show this help message and exit
  -sample_fast5dir SAMPLE_FAST5DIR
                        sample multi fast5 dir
  -control_fast5dir CONTROL_FAST5DIR
                        control multi fast54 dir
  -reference REFERENCE  reference fasta
  -ks_t KS_T            -log ks_test p-value (default 5)
  -outdir OUTDIR        output directory name
  -n_batches N_BATCHES  number of parsed fast5 batches
  -threads THREADS      number of threads used (derfault is 8)
  -max_motifs MAX_MOTIFS
                        the maximum expected number of motifs extracted
  -min_conf MIN_CONF    the minimal confidence value. Default is 1000

```


Example run command:
```
snapper -sample_fast5dir ../HelicobacterMod/fast5/J99_multi/ -control_fast5dir ../HelicobacterMod/fast5/J99_wga_multi/ -reference ../HelicobacterMod/genome/J99.fasta -n_batches 5
```

## Output explanation

The output directory contains the following files:
- `passed_motifs_[contig_name].fasta` - all k-mers defined as modifed
- `final_motifs_[contig_name].fasta` - optimal set of motifs obtained from the passed motifs with greedy algorithm
- `plots_[contig_name]` - signal distributions for all motifs detected 
