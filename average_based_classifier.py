# reads through the provided vistaEnhancerBrowser.txt file of train/test
# data and stores it in a map for use by train and test
def read():
	f = open('vistaEnhancerBrowser.txt', 'r')
	line = f.readline()
	info = ''
	sequence = ''
	parsed = dict()
	while(line):
		if line[0] == '>':
			info = line[1:]
			if len(sequence) > 0 and len(info) > 0:
				parsed[sequence] = info
				info = ''
				sequence = ''
		elif len(line) > 0:
			sequence += line
			
		line = f.readline()
	f.close()
	return parsed


# reads through the training data, parses out the labels provided and given scores,
# and maps k-mers of length 4, 5, 6, 7, and 8 to the labels for the sequence that
# contains them
def train(training_data):
	k_mer_to_tag_scores = dict() # will be a map from k-mer to list of dicts mapping label to score
	for sequence in training_data:

		tags = training_data[sequence]
		labels = dict() # create new for tagged each sequence
		pieces = tags.split('|')
		for entry in pieces:
			if entry.find('[') != -1:
				type = entry[0:entry.find('[')]
				score = entry[entry.find('['):]
				labels[type] = score

		for i in range(0, len(sequence)):
			# map 4-mers, 5-mers, 6-mers, 7-mers, and 8-mers in this sequence to all
			# provided with the full sequence 
			if i+4 < len(sequence):
				seq = sequence[i:i+4]
				if seq in k_mer_to_tag_scores:
					k_mer_to_tag_scores[seq].append(labels)
				else:
					k_mer_to_tag_scores[seq] = [labels]
			elif i+5 < len(sequence):
				seq = sequence[i:i+5]
				if seq in k_mer_to_tag_scores:
					k_mer_to_tag_scores[seq].append(labels)
				else:
					k_mer_to_tag_scores[seq] = [labels]			
			elif i+6 < len(sequence):
				seq = sequence[i:i+6]
				if seq in k_mer_to_tag_scores:
					k_mer_to_tag_scores[seq].append(labels)
				else:
					k_mer_to_tag_scores[seq] = [labels]
			elif i+7 < len(sequence):
				seq = sequence[i:i+7]
				if seq in k_mer_to_tag_scores:
					k_mer_to_tag_scores[seq].append(labels)
				else:
					k_mer_to_tag_scores[seq] = [labels]
			elif i+8 < len(sequence):
				seq = sequence[i:i+8]
				if seq in k_mer_to_tag_scores:
					k_mer_to_tag_scores[seq].append(labels)
				else:
					k_mer_to_tag_scores[seq] = [labels]

	return k_mer_to_tag_scores


# sets up the final label dict with all of the labels that will be predicted
def init_tags():
	tags = dict()
	# these will all be filled in with some value out of 12
	tags['neural tube'] = 0
	tags['hindbrain'] = 0
	tags['midbrain'] = 0
	tags['limb'] = 0
	tags['cranial nerve'] = 0
	tags['mesenchyme'] = 0
	tags['heart'] = 0
	tags['branchial arch'] = 0
	tags['nose'] = 0
	tags['dorsal root ganglion'] = 0
	tags['trigeminal V'] = 0
	tags['neural_tube'] = 0

	# will either be 1 or 0
	tags['positive'] = 0
	return tags


# divides the summed score for each label by the number of times the label
# appeared to come up with the predicted enhancement score
def compute_average(tags, counts, num_positive):
	results = init_tags()
	for label in tags:
		if tags[label] != 0 and counts[label] != 0:
			average = float(tags[label])/float(counts[label])
			results[label] = average

	if num_positive > 0:
		results['positive'] = 1
	print results
	return results


# updates the score tally in the tags and counts dicts
def add_score(name, tags, counts, label, numerical_score):
	if label.find(name) != -1:
		score = numerical_score[numerical_score.find('[')+1:numerical_score.find('/')]
		tags[name] += int(score)
		counts[name] += 1
	return [tags, counts]


# works through the list of labels suggested by the k-mers found in
# some sequence in the test data and then computes the most likely
# numeric value for each of the different enhancement types and
# positive/negative (enhanced/not enhanced) boolean value
def find_likeliness(labels):
	tags = init_tags()
	counts = init_tags()
	num_positive = 0
	for entry in labels:
		for label_list in entry:
			for label in label_list:
				numerical_score = label_list[label]
				dicts = add_score('neural tube', tags, counts, label, numerical_score)
				dicts = add_score('hindbrain', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('midbrain', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('limb', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('cranial nerve', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('mesenchyme', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('heart', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('branchial arch', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('heart', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('nose', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('dorsal root ganglion', dicts[0], dicts[1], label, numerical_score)
				dicts = add_score('trigeminal V', dicts[0], dicts[1], label, numerical_score)

				if label.find('positive') != -1:
					num_positive += 1
				if label.find('negative') != -1:
					num_positive -= 1

	return compute_average(dicts[0], dicts[1], num_positive)


# for each sequence in the test data, compiles a list of all of the labels from
# k-mers that were encountered in the training data and computes the most likely
# label for the given sequence using a weighted average
def test(encountered, test_data):
	classifications = dict()
	for sequence in test_data: # ignore the tags (treat the dict as a list)

		# will be a list of all labels encountered by sequences stored in the map constructed with training data
		labels = list()
		for i in range(0, len(sequence)):
			if i+3 < len(sequence):
				seq = sequence[i:i+3]
				if seq in encountered:
					labels.append(encountered[seq])
			if i+4 < len(sequence):
				seq = sequence[i:i+4]
				if seq in encountered:
					labels.append(encountered[seq])
			elif i+5 < len(sequence):
				seq = sequence[i:i+5]
				if seq in encountered:
					labels.append(encountered[seq])		
			elif i+6 < len(sequence):
				seq = sequence[i:i+6]
				if seq in encountered:
					labels.append(encountered[seq])
			elif i+7 < len(sequence):
				seq = sequence[i:i+7]
				if seq in encountered:
					labels.append(encountered[seq])
			elif i+8 < len(sequence):
				seq = sequence[i:i+8]
				if seq in encountered:
					labels.append(encountered[seq])
			elif i+9 < len(sequence):
				seq = sequence[i:i+9]
				if seq in encountered:
					labels.append(encountered[seq])
			elif i+10 < len(sequence):
				seq = sequence[i:i+10]
				if seq in encountered:
					labels.append(encountered[seq])

		label = find_likeliness(labels)
		classifications[sequence] = label
	return classifications


# compares the classification to the golden label provided with the data
# to assess its accuracy
def display(classifications, test_data):
	for sequence in classifications:
		print classifications[sequence]
		print test_data[sequence]
		print '\n\n\n'


# high-level function to parse the provided data, train, and then test
def predict_enhancements():
	data = read()
	training_data = dict(data.items()[len(data)/3:])
	testing_data = dict(data.items()[:len(data)/3])
	probabilities = train(training_data)
	classifications = test(probabilities, testing_data)
	display(classifications, testing_data)

predict_enhancements()
