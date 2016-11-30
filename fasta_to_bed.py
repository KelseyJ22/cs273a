def change_format(line, count):
	pieces = line.split('|')
	location = pieces[1].split(':')
	chrom = location[0]
	indexes = location[1].split('-')
	return str(chrom) + '\t' + str(indexes[0]) + '\t' + str(indexes[1]) + '\t' + str(count) + '\n'


def process():
	f = open('vistaEnhancerBrowser.txt', 'r')
	o = open('vistaEnhancerBrowser.bed', 'w')
	line = f.readline()
	count = 0
	while(line):
		if line[0] == '>':
			o.write(change_format(line, count))
			count += 1
		line = f.readline()
	f.close()
	o.close()

process()