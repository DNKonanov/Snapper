letter_codes = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'M': ['A','C'],
    'R': ['A','G'],
    'W': ['A','T'],
    'S': ['C','G'],
    'Y': ['C','T'],
    'K': ['G','T'],
    'V': ['A','C','G'],
    'H': ['A','C','T'],
    'D': ['A','G','T'],
    'B': ['C','G','T'],
    'N': ['A','C','G','T']
}


letter_anticodes = {
    'A': set(['C', 'G', 'T']),
    'C': set(['A', 'G', 'T']),
    'G': set(['A', 'C', 'T']),
    'T': set(['A', 'C', 'G']),
    'M': set(['G','T']),
    'R': set(['C','T']),
    'W': set(['C','G']),
    'S': set(['A','T']),
    'Y': set(['C','T']),
    'K': set(['A','C']),
    'V': set(['T']),
    'H': set(['G']),
    'D': set(['C']),
    'B': set(['A']),
    'N': set([])
}

letter_codes_rev = {
    ('A',): 'A',
    ('C',): 'C',
    ('G',): 'G',
    ('T',): 'T',
    ('A', 'C'): 'M',
    ('A', 'G'): 'R',
    ('A', 'T'): 'W',
    ('C', 'G'): 'S',
    ('C', 'T'): 'Y',
    ('G', 'T'): 'K',
    ('A', 'C', 'G'): 'V',
    ('A', 'C', 'T'): 'H',
    ('A', 'G', 'T'): 'D',
    ('C', 'G', 'T'): 'B',
    ('A', 'C', 'G', 'T'): 'N'
}





def gen_variants(seq):
    variants = ['']
    
    for i in range(len(seq)):
        new_variants = []
        for l in letter_codes[seq[i]]:
            
            for v in variants:
                new_variants.append(v + l)
                
        variants = new_variants
            
    return variants