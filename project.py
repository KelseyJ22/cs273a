import sys
import re
from itertools import *


# start w 6mers
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


def build_map_of_indices(k):
    #BUILDING THE MAP
    list_of_permutations = [''.join(permutation) for permutation in product('ATCG', repeat=k)]
    print list_of_permutations
    return 'hi'


def main():

    labels = {}

    total_samples = 0

    k = 4

    map_of_indices = build_map_of_indices(k)

    # THIS STARTS PARSING
    # f = open('vistaEnhancerBrowser.txt')
    # line = f.readline()
    # while(line):
    #     if '>' in line:
    #         # labels are [enhancer, brain, forebrain, midbrain, hindbrain, limb, neural tube, heart]
    #         labels = [0, 0, 0, 0, 0, 0, 0, 0]
    #         if 'positive' in line:
    #             labels[0] = 1
    #             pieces = line.strip().split('|')
    #             for x in xrange(4, len(pieces)):
    #                 label = re.sub('\[.*\]', '', re.sub('\(.*\)', '', pieces[x])).strip()
    #                 if label == 'forebrain':
    #                     labels[2] = 1
    #                     labels[1] = 1
    #                 elif label == 'midbrain':
    #                     labels[3] = 1
    #                     labels[1] = 1
    #                 elif label == 'hindbrain':
    #                     labels[4] = 1
    #                     labels[1] = 1
    #                 elif label == 'limb':
    #                     labels[5] = 1
    #                 elif label == 'neural tube':
    #                     labels[6] = 1
    #                 elif label == 'heart':
    #                     labels[7] = 1
    #         else:
                #do things
    # for key, value in sorted(labels.iteritems(), key=lambda (k,v): (v,k)):
    #     print "%s: %s" % (key, value)
    # print total_samples
    # f.close()

if __name__ == "__main__":
    main()