import matplotlib.pyplot as plt
from snapper.src.seq_processing import letter_codes, gen_variants
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import numpy as np


def cohend(d1, d2):
    n1, n2 = len(d1), len(d2)
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.mean(d1), np.mean(d2)
    return (u1 - u2) / s


regular_letters = ['A','G','C','T']




def gen_template(motif_variant, pos_variant):
    
    template = ['N',]*11
    
    for i, pos in enumerate(pos_variant):
        template[pos] = motif_variant[i]
        
    return ''.join(template)



def plot_motif(motif, sample_motifs, control_motifs, savepath):

    _sample = []
    _control = []


    ancMOTIF = gen_template(motif[1], motif[2])

    effect_size_dist = []

    for MOTIF in gen_variants(ancMOTIF):
            
        
        if MOTIF not in sample_motifs or MOTIF not in control_motifs:
            continue
            
        
        
        _sample += sample_motifs[MOTIF]
        _control += control_motifs[MOTIF]

        effect_size_dist.append(np.abs(cohend(sample_motifs[MOTIF], control_motifs[MOTIF])))
                


    plt.figure(figsize=(8,5))

    plt.grid()

    sns.distplot(x = _control, hist=False, label='Control', color='red')
    sns.distplot(x = _sample, hist=False, label='Sample', color='green')
    #plt.savefig('tnp/check3.png', dpi=400)

    plt.title('{}, confidence = {}\nmed effsize = {}'.format(ancMOTIF, motif[0], np.median(effect_size_dist)))

    plt.xlim(-5,5)
    plt.xlabel('Normalized signal')
    plt.legend()

    plt.tight_layout()
    plt.savefig(savepath + '/{}.png'.format(ancMOTIF), format='png', dpi=400)

    plt.close()