# NanoCaller: nanopore-based modification motifs caller

## Dependencies
- ont-tombo
- h5py
- biopython


## Usage
### Motifs extraction

```
usage: cli_main.py [-h] [-sample_fast5dir SAMPLE_FAST5DIR]
                   [-control_fast5dir CONTROL_FAST5DIR] [-reference REFERENCE]
                   [-ks_t KS_T] [-eff_size EFF_SIZE] [-outdir OUTDIR]
                   [-n_batches N_BATCHES]

optional arguments:
  -h, --help            show this help message and exit
  -sample_fast5dir SAMPLE_FAST5DIR
                        sample multi fast5 dir
  -control_fast5dir CONTROL_FAST5DIR
                        control multi fast54 dir
  -reference REFERENCE  reference fasta
  -ks_t KS_T            -log ks_test p-value (default 10)
  -eff_size EFF_SIZE    Cohen d-effect size (default 0.25)
  -outdir OUTDIR        output directory name
  -n_batches N_BATCHES  number of parsed fast5 batche
```


### Difsignals extraction

```
usage: cli.py [-h] [-motifs MOTIFS]

optional arguments:
  -h, --help      show this help message and exit
  -motifs MOTIFS  fasta file returned by NanoCaller
```