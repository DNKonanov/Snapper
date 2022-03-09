from argparse import ArgumentParser
from extract_motifs import get_seqs, extract_motifs, merge_motifs, collapse_motifs
from plot_motifs import get_variants

parser = ArgumentParser()

parser.add_argument('-motifs', type=str, help='fasta file returned by NanoCaller')


args = parser.parse_args()

seqs = get_seqs(args.motifs)

motifs = collapse_motifs(extract_motifs(seqs))


print()
print('--MOTIFS--')
for m in motifs:
    print(m)

print(get_variants('TGRCAR'))

#print()
#print('--PAIRS---')
#merge_motifs(motifs)