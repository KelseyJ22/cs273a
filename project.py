import sys
import re
from itertools import *
import pandas as pd
import numpy as np
import math

def reverse_complement(sequence):
    rc = ''
    for base in sequence:
        if base == 'A':
            rc = 'T' + rc
        elif base == 'T':
            rc = 'A' + rc
        elif base == 'C':
            rc = 'G' + rc
        else:
            rc = 'C' + rc
    return rc

def build_map_of_indices(k):
    list_of_permutations = [''.join(permutation) for permutation in product('ATCG', repeat=k)]
    map_of_indices = {}
    index = 0
    num_permutations = 0
    for permutation in list_of_permutations:
        if reverse_complement(permutation) not in map_of_indices:
            num_permutations += 1
            map_of_indices[permutation] = index
            index += 1
    return map_of_indices

def window_to_string(window):
    to_return = ''
    for base in window:
        to_return += base
    return to_return

def get_averages(lines):
    results = list()
    for line in lines:
        pieces = line.split('\t')
        if len(pieces) >= 4:
            avg = float(pieces[4].rstrip())
            if not math.isnan(avg):
                if avg != float('Inf'):
                    results.append(avg)
                else:
                    results.append(0.0)
            else:
                results.append(0.0)
        else:
            results.append(0.0)
    return results


def main():
    filename = 'vistaEnhancerBrowser.txt'

    print 'Parsing...'
    k = 6

    print 'Building map of indices...'
    map_of_indices = build_map_of_indices(k)

    # counts number of sequences so DataFrame can be preallocated, since
    # adding rows to DF one by one is very expensive
    print 'Preallocating DataFrame...'
    f = open(filename)
    num_seqs = 0
    for line in f.readlines():
        if len(line) > 0 and line[0] == '>':
            num_seqs += 1
    f.close()

    data = pd.DataFrame(index=np.arange(0, num_seqs), columns=[x for x in xrange(len(map_of_indices) + 13)])
    f = open('vistaEnhancerBrowser.txt')
    t1 = open('H3K27me3.tab', 'r')
    t2 = open('H3K36me3.tab', 'r')
    t3 = open('H3K4me3.tab', 'r')
    t4 = open('H3K4me1.tab', 'r')
    t5 = open('H3K9ac.tab', 'r')

    l1 = t1.readline()
    l2 = t2.readline()
    l3 = t3.readline()
    l4 = t4.readline()
    l5 = t5.readline()
    line = f.readline()
    labels = []
    window = []
    features = [0 for x in xrange(len(map_of_indices))]
    sequence = 0
    print 'Computing features...'
    while(line):
        if '>' in line:
            if labels != []:
                # normalizes frequencies
                divisor = float(sum(features))
                features = [feature / divisor for feature in features]

                averages = get_averages([l1, l2, l3, l4, l5])
                
                full_vector = features + averages + labels

                data.loc[sequence] = full_vector
                features = [0 for x in xrange(len(map_of_indices))]
                sequence += 1
                window = []
                l1 = t1.readline()
                l2 = t2.readline()
                l3 = t3.readline()
                l4 = t4.readline()
                l5 = t5.readline()

            # labels are [enhancer, brain, forebrain, midbrain, hindbrain, limb, neural tube, heart]
            labels = [0, 0, 0, 0, 0, 0, 0, 0]
            if 'positive' in line:
                labels[0] = 1
                pieces = line.strip().split('|')
                for x in xrange(4, len(pieces)):
                    label = re.sub('\[.*\]', '', re.sub('\(.*\)', '', pieces[x])).strip()
                    if label == 'forebrain':
                        labels[2] = 1
                        labels[1] = 1
                    elif label == 'midbrain':
                        labels[3] = 1
                        labels[1] = 1
                    elif label == 'hindbrain':
                        labels[4] = 1
                        labels[1] = 1
                    elif label == 'limb':
                        labels[5] = 1
                    elif label == 'neural tube':
                        labels[6] = 1
                    elif label == 'heart':
                        labels[7] = 1
        else:
            for base in line.strip():
                if len(window) < k:
                    window.append(base.upper())
                else:
                    window = window[1:] + [base.upper()]
                if (len(window) == k):
                    window_string = window_to_string(window)
                    if (window_to_string(window)) in map_of_indices:
                        features[map_of_indices[window_string]] += 1
                    else:
                        features[map_of_indices[reverse_complement(window_string)]] += 1
        line = f.readline()

    # normalizes frequencies
    divisor = float(sum(features))
    features = [feature / divisor for feature in features]
    averages = get_averages([l1, l2, l3, l4, l5])
    full_vector = features + averages + labels
    data.loc[sequence] = full_vector

    f.close()
    print 'Done parsing! Data ready to use.'
    data.reindex(np.random.permutation(data.index)).to_pickle('pickled_vista_data_k_' + str(k) + '.pkl')
    # To load the data, do 'data = pd.read_pickle(filename)'

if __name__ == "__main__":
    main()