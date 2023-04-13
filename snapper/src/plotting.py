import matplotlib.pyplot as plt
from snapper.src.seq_processing import letter_codes, gen_variants
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import numpy as np

from random import sample


def cohend(d1, d2):
    n1, n2 = len(d1), len(d2)
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.mean(d1), np.mean(d2)
    return (u1 - u2) / s


regular_letters = ['A','G','C','T']



def gen_template(motif_variant, pos_variant, lenmotif):
    
    template = ['N',]*lenmotif
    
    for i, pos in enumerate(pos_variant):
        template[pos] = motif_variant[i]
        
    return ''.join(template)



def get_anc_variants(ancMOTIF, lenmotif=15):
    
    ancMOTIF_c = ancMOTIF
    while ancMOTIF_c[0] == 'N':
        ancMOTIF_c = ancMOTIF_c[1:]
        
    while ancMOTIF_c[-1] == 'N':
        ancMOTIF_c = ancMOTIF_c[:-1]
        
    ext_length = lenmotif - len(ancMOTIF_c)
    anc_variants = []
    print(ancMOTIF_c)
    for i in range(ext_length + 1):
        anc_variants.append('N'*i + ancMOTIF_c + 'N'*(ext_length - i))
        
    return anc_variants

def plot_dist(motif, native_motifs, wga_motifs, savepath, lenmotif=15, MAXSAMPLESIZE = 2000):
    
    print(f'Rendering {"".join(motif[1])}...')
    
    ancMOTIF_init = gen_template(motif[1], motif[2], lenmotif)

    anc_variants = get_anc_variants(ancMOTIF_init, lenmotif=lenmotif)
    
    N = len(anc_variants)
    fig, axs = plt.subplots(N, 2, figsize=(14, 4*N))
    
    motif_cnt = 0
    for ancMOTIF in anc_variants:
        
        _wga = []
        _native = []
        
        lens = []
        effect_size_dist = []

        cnt = 0
        for MOTIF in gen_variants(ancMOTIF):


            if MOTIF not in wga_motifs or MOTIF not in native_motifs:
                continue

            _wga += wga_motifs[MOTIF]
            _native += native_motifs[MOTIF]

            effect_size_dist.append(np.abs(cohend(native_motifs[MOTIF], wga_motifs[MOTIF])))

            #print(len(_wga), len(_native))
        
        if len(effect_size_dist) == 0:
            continue
        if len(_native) > MAXSAMPLESIZE:
            _native = sample(_native, MAXSAMPLESIZE)

        if len(_wga) > MAXSAMPLESIZE:
            _wga = sample(_wga, MAXSAMPLESIZE)
        
        
        sns.distplot(x = _wga, hist = False, label='WGA', color='red', ax=axs[motif_cnt][0])
        sns.distplot(x = _native, hist = False, label='native', color='green', ax=axs[motif_cnt][0])
        axs[motif_cnt][0].grid()

        
        axs[motif_cnt][0].set_title('{}, confidence = {}\nmed effsize = {}'.format(ancMOTIF, motif[0], np.median(effect_size_dist)))
        axs[motif_cnt][0].set_xlabel('Normalized signal')

        axs[motif_cnt][1].hist(effect_size_dist, bins=50, density=True, rwidth=0.8, color='black')
        axs[motif_cnt][1].set_xlabel('eff size')
        
        title = np.round([
                np.percentile(effect_size_dist, 10),
                np.percentile(effect_size_dist, 25),
                np.percentile(effect_size_dist, 50),
                np.percentile(effect_size_dist, 75),
                np.percentile(effect_size_dist, 90)
                ], 2)
        
        title = list(map(str, title))
        axs[motif_cnt][1].set_title(
            f'{" ".join(title)}'
        )

        axs[motif_cnt][0].legend()
        
        motif_cnt += 1

    plt.tight_layout()

    plt.savefig(savepath + '/{}.png'.format("".join(motif[1])), format='png', dpi=400)
    plt.show()




def plot_motif(motif, sample_motifs, control_motifs, savepath, lenmotif=11):

    _sample = []
    _control = []


    ancMOTIF = gen_template(motif[1], motif[2], lenmotif)


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



def plot_coverage(sample_coverage, control_coverage, chrom, output):

    plt.figure(figsize=(25,7))

    plt.xlabel('chrom position')
    plt.ylabel('depth')

    sample_mean_cov = np.round(np.mean(sample_coverage), 0)
    control_mean_cov = np.round(np.mean(control_coverage), 0)


    plt.title(f'{chrom}\nNative mean cov = {sample_mean_cov}X\nControl mean cov = {control_mean_cov}X')

    plt.plot(sample_coverage, label='Native')
    plt.plot(control_coverage, label='Control')

    plt.legend()

    plt.tight_layout()

    plt.savefig(output, format='pdf')