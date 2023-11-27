<img src="logo.png" align="left" width=150> 

# Snapper: nanopore-based modification motifs caller

This tool is designed to efficiently detect methylation sites using ONT sequencing data.
Snapper uses balanced approach to compute statistics for each k-mer which is likely to be modified.
The core feature of Snapper in comparison with other tools is a new high-sensitive greedy algorithm that is used 
for position-specific motif enrichment. This repository contains not the Snapper tool itself but its pip distribution.

## Dependencies
- python 3.7 (later versions might be incompatible because of inner biopython dependencies)
- ont-tombo
- h5py
- biopython
- matplotlib
- scipy
- seaborn

## Installation

```
(base) $ conda create -n snapper python=3.7
(base) $ conda activate snapper
(snapper) $ conda install -c bioconda ont-fast5-api ont-tombo
(snapper) $ pip install snapper-ont
```

## Usage

Firstly, fast5 files should be resquiggled using [Tombo](https://github.com/nanoporetech/tombo) software. 
After resquiggling, fast5 files should be converted to the multi-fast5 format using [ont_fast5_api](https://github.com/nanoporetech/ont_fast5_api).

A more detailed usage guideline and few usercases are available in [Snapper's documentation](https://snapper-tutorial.readthedocs.io/en/latest/index.html)

```
usage: snapper [-h] -sample_fast5dir SAMPLE_FAST5DIR -control_fast5dir
               CONTROL_FAST5DIR -reference REFERENCE [-ks_t KS_T]
               [-outdir OUTDIR] [-coverage COVERAGE] [-threads THREADS]
               [-k_size K_SIZE] [-long_k_size LONG_K_SIZE]
               [-max_motifs MAX_MOTIFS] [-min_conf MIN_CONF]
               [-target_chr TARGET_CHR]

optional arguments:
  -h, --help            show this help message and exit
  -sample_fast5dir SAMPLE_FAST5DIR
                        sample multi fast5 dir
  -control_fast5dir CONTROL_FAST5DIR
                        control multi fast5 dir
  -reference REFERENCE  reference genome in the fasta format
  -ks_t KS_T            -log ks_test p-value (default 3).
  -outdir OUTDIR        output directory name
  -coverage COVERAGE    minimal genome coverage depth (default 40)
  -threads THREADS      number of threads used (default 8)
  -k_size K_SIZE        k-mer size, must be odd, 
                        should not be less than 11 (default 15)
  -long_k_size LONG_K_SIZE
                        k-mer size, must be odd, 
                        should not be less than 21 (default 29)
  -max_motifs MAX_MOTIFS
                        the maximum expected number of motifs extracted
  -min_conf MIN_CONF    the minimal confidence value (default is 100)
  -target_chr TARGET_CHR
                        target chromosome name (by default all
                        contigs/replicons are considered)


```


Typical run command:
```
snapper -sample_fast5dir ../HelicobacterMod/fast5/J99_multi/ -control_fast5dir ../HelicobacterMod/fast5/J99_wga_multi/ -reference ../HelicobacterMod/genome/J99.fasta
```

## Output explanation

The output directory contains the following files:
- `passed_motifs_[strand]_[contig_name].fasta` - all k-mers that most likely bring a modified base
- `final_motifs_[strand]_[contig_name].fasta` - optimal set of motifs generated from the passed motifs by the Snapper greedy algorithm
- `plots_[strand]_[contig_name]` - signal distribution plots for each extracted motif  

## Citation

Dmitry N Konanov, Vladislav V Babenko, Aleksandra M Belova, Arina G Madan, Daria I Boldyreva, Oksana E Glushenko, Ivan O Butenko, Dmitry E Fedorov, Alexander I Manolov, Danil V Krivonos, Vassilii N Lazarev, Vadim M Govorun, Elena N Ilina, [Snapper: high-sensitive detection of methylation motifs based on Oxford Nanopore reads](https://doi.org/10.1093/bioinformatics/btad702), Bioinformatics, 2023
