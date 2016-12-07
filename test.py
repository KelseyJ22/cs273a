f = open('H3K36me3.tab', 'r')
line = f.readline()
count = 0
while(line):
	print line
	print count
	line = f.readline()
	count += 1