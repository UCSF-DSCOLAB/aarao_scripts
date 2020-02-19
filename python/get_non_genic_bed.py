chr_lengths = {
    '1': 248956422,
    '2': 242193529,
    '3': 198295559,
    '4': 190214555,
    '5': 181538259,
    '6': 170805979,
    '7': 159345973,
    '8': 145138636,
    '9': 138394717,
    '10': 133797422,
    '11': 135086622,
    '12': 133275309,
    '13': 114364328,
    '14': 107043718,
    '15': 101991189,
    '16': 90338345,
    '17': 83257441,
    '18': 80373285,
    '19': 58617616,
    '20': 64444167,
    '21': 46709983,
    '22': 50818468,
    'X': 156040895,
    'Y': 57227415,
    'MT': 16569,
}
genic_intervals = {str(x):[] for x in list(range(1, 23)) + ['X', 'Y', 'MT']}

with open('/cbc/resources/human/hg38/Homo_sapiens.GRCh38.85.gtf') as iff:
    for line in iff:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[0] not in genic_intervals or line[2] != 'gene':
            continue
        genic_intervals[line[0]].append((int(line[3]), int(line[4])))

genic_intervals = {x: sorted(genic_intervals[x], key=lambda y: y[0]) for x in genic_intervals}

with open('genic_intervals.bed', 'w') as off:
    for chrom in genic_intervals:
        for istart, istop in genic_intervals[chrom]:
                print(chrom, istart, istop, sep='\t', end='\n', file=off)


with open('non_genic_intervals.bed', 'w') as off:
    for chrom in genic_intervals:
        start = 1
        for istart, istop in genic_intervals[chrom]:
            if start < istart:
                print(chrom, start, istart-1, sep='\t', end='\n', file=off)
            start = istop+1
        print(chrom, start, chr_lengths[chrom], sep='\t', end='\n', file=off)


