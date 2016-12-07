#!/usr/local/bin/python

"""
	svm_enh_train.py; train a support vector machine using shogun toolbox
	Copyright (C) 2011 Dongwon Lee

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import numpy
import optparse

from shogun.PreProc import SortWordString, SortUlongString
from shogun.Kernel import CommWordStringKernel, CommUlongStringKernel

from shogun.Features import StringWordFeatures, StringUlongFeatures, \
		StringCharFeatures, DNA, Labels

from shogun.Classifier import SVMLight, LibSVM, MSG_INFO

"""
global variables
"""
g_kmers_list = []
g_revcomp_mapping = []

def generate_kmers_list(kmerlen):
	nts = ['A', 'C', 'G', 'T']
	kmers = [] 
	kmers.append('')
	l = 0
	while l < kmerlen: 
		imers = [] 
		for imer in kmers:
			for nt in nts: 
				imers.append(imer+nt)
		kmers = imers
		l += 1 
	
	return kmers


def generate_revcomp_mapping_table(kmerlen):
	kmers = g_kmers_list

	rc = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
	revcomp = lambda s: ''.join([rc[s[i]] for i in xrange(len(s)-1, -1, -1)])

	kmer_id_dict = {}
	for i in xrange(len(kmers)): 
		kmer_id_dict[kmers[i]] = i

	revcomp_mapping_table = []
	for kmerid in xrange(len(kmers)): 
		rc_id = kmer_id_dict[revcomp(kmers[kmerid])]
		if rc_id < kmerid:
			revcomp_mapping_table.append(rc_id)
		else:
			revcomp_mapping_table.append(kmerid)
	
	return revcomp_mapping_table


def read_fastafile(filename):
	"""Read a file in FASTA format

	Arguments:
	filename -- string, the name of the sequence file in FASTA format

	Return: 
	list of sequences, list of sequence ids
	"""
	id = '' 
	ids = []
	seqs = []

	try:
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)

	seq = [] 
	for line in lines:
		if line[0] == '>':
			ids.append(line[1:].rstrip('\n'))
			if seq != []: seqs.append("".join(seq))
			seq = []
		else:
			seq.append(line.rstrip('\n').upper())

	if seq != []:
		seqs.append("".join(seq))

	return seqs, ids


def non_redundant_word_features(feats, kmerlen):
	revcomp_mapping = g_revcomp_mapping

	for idx in xrange(feats.get_num_vectors()):
		fvec = feats.get_feature_vector(idx)

		fvec2 = []
		for kmerid in fvec:
			fvec2.append(revcomp_mapping[int(kmerid)])

		feats.set_feature_vector(numpy.array(fvec2, numpy.dtype('u2')), idx)

	preproc = SortWordString()
	preproc.init(feats)
	feats.add_preproc(preproc)
	feats.apply_preproc()

	return feats


def non_redundant_ulong_features(feats, kmerlen):
	revcomp_mapping = g_revcomp_mapping

	for idx in xrange(feats.get_num_vectors()):
		fvec = feats.get_feature_vector(idx)

		fvec2 = []
		for kmerid in fvec:
			fvec2.append(revcomp_mapping[int(kmerid)])

		feats.set_feature_vector(numpy.array(fvec2, numpy.dtype('u8')), idx)

	preproc = SortUlongString()
	preproc.init(feats)
	feats.add_preproc(preproc)
	feats.apply_preproc()

	return feats


def svm_learn(kernel, labels, svmC, epsilon, weight):
	try: 
		svm=SVMLight(svmC, kernel, Labels(labels))
	except NameError:
		print 'No support for SVMLight available.'
		return

	svm.io.set_loglevel(MSG_INFO)
	svm.io.set_target_to_stderr()

	svm.set_epsilon(epsilon)
	svm.parallel.set_num_threads(1)
	if weight != 1.0:
		svm.set_C(svmC, svmC*weight)
	svm.train()

	return svm


def main(argv = sys.argv):
	usage = "Usage: %prog [options] POSITIVE_SEQ NEGATIVE_SEQ"
	desc  = "1. take two files(FASTA format) as input, 2. train an SVM and store the trained SVM weights"
	parser = optparse.OptionParser(usage=usage, description=desc)                                                                              
	parser.add_option("-C", dest="svmC", type="float", default=1, \
			help="set the regularization parameter svmC (default=1)")

	parser.add_option("-e", dest="epsilon", type="float", default=0.00001, \
			help="set the precision parameter epsilon (default=0.00001)")

	parser.add_option("-w", dest="weight", type="float", default=1.0, \
			help="set the weight for positive set (default=1.0)")

	parser.add_option("-k", dest="kmerlen", type="int",default=6, \
			help="set the length of k-mer (default = 6)")

	parser.add_option("-o", dest="output", default="svm_enh_weights.out", \
  			help="set the name of output svm weight file (default=svm_enh_weights.out)")

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 2:
		parser.error("incorrect number of arguments")
		sys.exit(0)

	svmC = options.svmC
	epsilon = options.epsilon
	weight = options.weight
	kmerlen = options.kmerlen
	output = options.output

	"""
	set global variable
	!!!IMPORTANT!!! order of generation should be preserved
	"""
	global g_kmers_list
	g_kmers_list = generate_kmers_list(kmerlen)

	global g_revcomp_mapping
	g_revcomp_mapping = generate_revcomp_mapping_table(kmerlen)
	
	posf = args[0]
	negf = args[1]
	
	sys.stderr.write('SVM parameters:\n')
	sys.stderr.write('  svm-C: ' + str(svmC) + '\n')
	sys.stderr.write('  epsilon: ' + str(epsilon) + '\n')
	sys.stderr.write('  weight: ' + str(weight) + '\n')
	sys.stderr.write('\n')

	sys.stderr.write('Other options:\n')
	sys.stderr.write('  kmerlen: ' + str(kmerlen) + '\n')
	sys.stderr.write('  output: ' + output + '\n')
	sys.stderr.write('\n')

	sys.stderr.write('Input args:\n')
 	sys.stderr.write('  positive sequence file: ' + posf + '\n')
	sys.stderr.write('  negative sequence file: ' + negf + '\n')
	sys.stderr.write('\n')

	seqs_pos, sids_pos = read_fastafile(posf)
	seqs_neg, sids_neg = read_fastafile(negf)
	npos = len(seqs_pos)
	nneg = len(seqs_neg)
	seqs = seqs_pos + seqs_neg
	sids = sids_pos + sids_neg
	sys.stderr.write('numer of total positive seqs: ' + str(npos) + '\n')
	sys.stderr.write('numer of total negative seqs: ' + str(nneg) + '\n')

	#generate labels
	labels = [numpy.double(1)]*npos + [numpy.double(-1)]*nneg
	char_feats_tr = StringCharFeatures(seqs, DNA)

	sys.stderr.write('..kernel building..\n')
	kernel = []
	if kmerlen <= 8:
		feats_tr = StringWordFeatures(DNA)
		feats_tr.obtain_from_char(char_feats_tr, kmerlen-1, kmerlen, 0, False)
		feats_tr = non_redundant_word_features(feats_tr, kmerlen)
		kernel = CommWordStringKernel(feats_tr, feats_tr)
	else:
		feats_tr = StringUlongFeatures(DNA)
		feats_tr.obtain_from_char(char_feats_tr, kmerlen-1, kmerlen, 0, False)
		feats_tr = non_redundant_ulong_features(feats_tr, kmerlen)
		kernel = CommUlongStringKernel(feats_tr, feats_tr)

	sys.stderr.write('..svm learning..\n')
	svm = svm_learn(kernel, numpy.array(labels), svmC, epsilon, weight)

	sys.stderr.write('..calculating SVM weights..\n')
	alphas = svm.get_alphas()
	support_vector_ids = svm.get_support_vectors()

	w = numpy.array([0]*(2**(2*kmerlen)), numpy.double)

	for i in xrange(len(alphas)):
		x = numpy.array([0]*(2**(2*kmerlen)), numpy.double)
		for kmerid in feats_tr.get_feature_vector(int(support_vector_ids[i])):
			x[kmerid] += 1
		x = x/numpy.sqrt(numpy.sum(x**2))
		w += (alphas[i] * x)

	sys.stderr.write('..writing a result..\n')

	f = open(output, 'w')
	f.write("bias\t" + str(svm.get_bias()) + "\n")

	rc = {'A':'T', 'G':'C', 'C':'G', 'T':'A'}
	revcomp = lambda s: ''.join([rc[s[i]] for i in xrange(len(s)-1, -1, -1)])

	for i in xrange(len(w)): 
		if i == g_revcomp_mapping[i]:
			f.write(g_kmers_list[i] + "\t" + revcomp(g_kmers_list[i]) + "\t" \
					+ str(w[i]) + "\n")

	f.close()

if __name__=='__main__': main()
