filenames = ['H3K27me3.tab', 'H3K36me3.tab', 'H3K4me1.tab', 'H3K4me3.tab', 'H3K9ac.tab']

def process_veb():
	res = dict()
	f = open('vistaEnhancerBrowser.txt', 'r')
	line = f.readline()
	sequence = ''
	key = ''
	while(line):
		if line[0] == '>':
			pieces = line[1:].split('|')
			if len(pieces) > 1:
				location = pieces[1]
				chrom = location.split(':')[0]
				start_end = location.split(':')[1].split('-')
				key = chrom + ' '+ start_end[0] + ' ' + start_end[1].rstrip()

		elif len(line.rstrip()) > 0:
			sequence += line.rstrip()

		else:
			res[key] = sequence
			sequence = ''
			key = ''
		line = f.readline()
	return res


def process_veb_bed():
	res = dict()
	f = open('vistaEnhancerBrowser.bed', 'r')
	line = f.readline()
	while(line):
		pieces = line.split('\t')
		location = pieces[0] + ' ' + pieces[1] + ' ' + pieces[2].rstrip()
		number = pieces[3]
		res[int(number)] = location
		line = f.readline()

	return res


def display(dictionary):
	for entry in dictionary:
		print str(entry) + ' ' + str(dictionary[entry])


def parse():
	global filenames
	veb_bed_map = process_veb_bed()
	loc_to_seq_map = process_veb()
	avg_to_seq = dict()
	for fname in filenames:
		f = open(fname, 'r')
		line = f.readline()
		while(line):
			if line[0] != '#':
				pieces = line.split('\t')
				avg = pieces[4]
				if int(pieces[0]) in veb_bed_map:
					loc = veb_bed_map[int(pieces[0])]
					if loc in loc_to_seq_map:
						seq = loc_to_seq_map[loc]
						avg_to_seq[avg] = seq # add to the final map
			line = f.readline()
		f.close()
	display(avg_to_seq)


parse()