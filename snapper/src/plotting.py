import matplotlib.pyplot as plt


regular_letters = ['A','G','C','T']



def get_variants(seq):

    if len(seq) == 4:
        seq = 'N' + seq + 'N'

    elif len(seq) == 5:
        seq = seq + 'N'


    variants = ['']
    for letter in seq:

        if letter in regular_letters:
            for i in range(len(variants)):
                variants[i] += letter

        if letter not in regular_letters:
            
            new_variants = []
            for i in range(len(variants)):
                for l in regular_letters:
                    new_variants.append(variants[i] + l)
            variants = new_variants


    return variants


def plot_motif(native_signals: dict, wga_signals: dict, motif: str, outputdir):
        

    motif_variants = get_variants(motif)

    print('Creating pdf for {}'.format(motif))
    fig, axs = plt.subplots(
        1, len(motif_variants), 
        figsize=(5*len(motif_variants),5)
    )


    for i in range(len(motif_variants)):
        MOTIF = motif_variants[i]

        try:
            axs[i].hist(
                    native_signals[MOTIF], alpha=0.5, bins=100, density=True, label='native'
                )
            axs[i].hist(
                    wga_signals[MOTIF], alpha=0.5, bins=100, density=True, label='wga'
                )
            
            axs[i].set_title(MOTIF)
        except TypeError:
            
            axs.hist(
                    native_signals[MOTIF], alpha=0.5, bins=100, density=True, label='native'
                )
            axs.hist(
                    wga_signals[MOTIF], alpha=0.5, bins=100, density=True, label='wga'
                )
            
            axs.set_title(MOTIF)

    plt.savefig('{}/{}.pdf'.format(outputdir, motif))
    plt.close()