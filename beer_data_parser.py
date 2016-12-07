import sys
import re
from itertools import *
import pandas as pd
import numpy as np

# Removed element 1517 because it contained 'n' as a base
# start w 6mers
# uppercase everything
# feature vectors are counts of kmer appearance in sequence
# then have variable for enhancer or not, then labels
# fold in reverse complements
# normalize
# group mouse and human
# ignore #s after labels
# remember class weight (try without too)
# use forebrain, midbrain, hindbrain, limb, neural tube, heart (count 243) | branchial arch (count 155)
# make map of permutations to index, use that to build feature vectors for each sequence (remove complements here)
# then try different classifiers

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

def main():
    pos_file_1 = 'sequence_files/enh_fb.fa'
    pos_file_2 = ''
    pse_file_3 = ''
    neg_file = 'sequence_files/random4000.fa'
    out_file = 'enh_fb_random4000_k_3.pkl'

    print 'Parsing...'
    k = 3

    print 'Building map of indices...'
    map_of_indices = build_map_of_indices(k)

    # counts number of sequences so DataFrame can be preallocated, since
    # adding rows to DF one by one is very expensive
    print 'Preallocating DataFrame...'
    f = open(pos_file_1)
    num_seqs = 0
    for line in f.readlines():
        if len(line) > 0 and line[0] == '>':
            num_seqs += 1
    f.close()
    f = open(neg_file)
    for line in f.readlines():
        if len(line) > 0 and line[0] == '>':
            num_seqs += 1
    f.close()

    data = pd.DataFrame(index=np.arange(0, num_seqs), columns=[x for x in xrange(len(map_of_indices) + 1)])

    f = open(pos_file_1)
    line = f.readline()
    labels = []
    window = []
    features = [0 for x in xrange(len(map_of_indices))]
    sequence = 0
    print 'Computing positive features...'

    # POSITIVE SAMPLES
    while(line):
        if '>' in line:
            if labels != []:
                # normalizes frequencies
                divisor = float(sum(features))
                features = [feature / divisor for feature in features]
                full_vector = features + labels
                data.loc[sequence] = full_vector
                features = [0 for x in xrange(len(map_of_indices))]
                sequence += 1
                window = []
            labels = [1]
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
    full_vector = features + labels
    data.loc[sequence] = full_vector
    sequence += 1

    f.close()

    labels = []
    window = []
    features = [0 for x in xrange(len(map_of_indices))]
    f = open(neg_file)
    line = f.readline()
    print 'Computing negative features...'

    # NEGATIVE SAMPLES
    while(line):
        if '>' in line:
            if labels != []:
                # normalizes frequencies
                divisor = float(sum(features))
                features = [feature / divisor for feature in features]
                full_vector = features + labels
                data.loc[sequence] = full_vector
                features = [0 for x in xrange(len(map_of_indices))]
                sequence += 1
                window = []
            labels = [-1]
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
    full_vector = features + labels
    data.loc[sequence] = full_vector

    f.close()

    data.reindex(np.random.permutation(data.index)).to_pickle(out_file)
    # To load the data, do 'data = pd.read_pickle(filename)'
    print 'Done parsing! Data ready to use.'

if __name__ == "__main__":
    main()