import matplotlib.pyplot as plt
from snapper.src.seq_processing import letter_codes, gen_variants
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns

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


    for MOTIF in gen_variants(ancMOTIF):
            
        
        if MOTIF not in sample_motifs:
            continue
            
        
        
        _sample += sample_motifs[MOTIF]
        _control += control_motifs[MOTIF]
                


    plt.figure(figsize=(8,5))

    plt.grid()

    sns.distplot(x = _control, hist=False, label='Control', color='red')
    sns.distplot(x = _sample, hist=False, label='Sample', color='green')
    #plt.savefig('tnp/check3.png', dpi=400)

    plt.title(ancMOTIF)

    plt.xlim(-5,5)
    plt.xlabel('Normalized signal')
    plt.legend()

    plt.tight_layout()
    plt.savefig(savepath + '/{}.png'.format(ancMOTIF), format='png', dpi=400)

    plt.close()