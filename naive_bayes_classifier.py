import pandas as pd
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import classification_report

LABEL_SIZE = 8

def classify(data):
	vec_size = len(data.columns) - LABEL_SIZE
	features = data[data.columns[0:vec_size]]
	labels_matrix = data[data.columns[vec_size:vec_size + LABEL_SIZE]]
	enhancer_labels = labels_matrix[labels_matrix.columns[0]]
	classifier = GaussianNB()
	classifier.fit(features, enhancer_labels.tolist())
	predictions = classifier.predict(features)
	print classification_report(enhancer_labels.tolist(), predictions)


def read_pickle():
	k = 6
	filename = "pickled_vista_data_k_{0}.pkl".format(k)
	data = pd.read_pickle(filename)
	return data


def process():
	data = read_pickle()
	classify(data)

process()