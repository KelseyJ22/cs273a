#!/usr/local/bin/python
"""
	svm_enh_classify.py; classify sequences using shogun toolbox
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


def read_svmwfile(filename, kmerlen):
	try:
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()

	except IOError, (errno, strerror):
		print "I/O error(%d): %s" % (errno, strerror)
		sys.exit(0)

	kmers = g_kmers_list

	kmer_id_dict = {}
	for i in xrange(len(kmers)): 
		kmer_id_dict[kmers[i]] = i

	bias = float(lines[0].split()[1])

	svmw = [0]*(2**(2*kmerlen))
	for line in lines[1:]:
		s = line.split()
		if len(s[0]) != kmerlen:
			print "length of k-mer in the SVM weight file is different from the kmerlen: ", len(s[0]), "<>", kmerlen
			sys.exit()

		kid = kmer_id_dict[s[0]]
		svmw[kid] = s[2]

	return bias, numpy.array(svmw, numpy.double)


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


def main(argv = sys.argv):
	usage = "Usage: %prog [options] SVM_WEIGHTS TEST_SEQ"
	desc  = "1. take two files(one is in FASTA format to score, the other is SVM weight file generated from svm_enh_train.py) as input, 2. score each sequence in the given file"
	parser = optparse.OptionParser(usage=usage, description=desc)                                                                              

	parser.add_option("-k", dest="kmerlen", type="int",default=6, \
			help="set the length of k-mer (default = 6)")

	parser.add_option("-o", dest="output", default="svm_enh_scores.out", \
  			help="set the name of output score file (default=svm_enh_scores.out)")

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 2:
		parser.error("incorrect number of arguments")
		sys.exit(0)

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
	
	svmwf = args[0]
	seqf = args[1]
	
	sys.stderr.write('Options:\n')
	sys.stderr.write('  kmerlen: ' + str(kmerlen) + '\n')
	sys.stderr.write('  output: ' + output + '\n')
	sys.stderr.write('\n')

	sys.stderr.write('Input args:\n')
	sys.stderr.write('  SVM weights file: ' + svmwf + '\n')
 	sys.stderr.write('  sequence file: ' + seqf + '\n')
	sys.stderr.write('\n')

	seqs_te, sids_te = read_fastafile(seqf)
	bias, svmw = read_svmwfile(svmwf, kmerlen)

	char_feats_te = StringCharFeatures(seqs_te, DNA)
	feats_te = []

	if kmerlen <= 8:
		feats_te = StringWordFeatures(DNA)
		feats_te.obtain_from_char(char_feats_te, kmerlen-1, kmerlen, 0, False)
		feats_te = non_redundant_word_features(feats_te, kmerlen)
	else:
		feats_te = StringUlongFeatures(DNA)
		feats_te.obtain_from_char(char_feats_te, kmerlen-1, kmerlen, 0, False)
		feats_te = non_redundant_ulong_features(feats_te, kmerlen)

	sys.stderr.write('..scoring..\n')
	f = open(output, 'w')
	for i in xrange(len(seqs_te)):
		x = numpy.array([0]*(2**(2*kmerlen)), numpy.double)
		for kmerid in feats_te.get_feature_vector(i):
			x[kmerid] += 1
		x = x/numpy.sqrt(numpy.sum(x**2))
		s = numpy.dot(svmw, x) + bias
		f.write(sids_te[i] + "\t" + str(s) + "\n")

	f.close()

if __name__=='__main__': main()
