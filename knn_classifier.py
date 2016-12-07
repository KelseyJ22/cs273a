from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report
#from sklearn.model_selection import cross_val_score
#if this import complains, comment it out and uncomment the import above
from sklearn.cross_validation import cross_val_score
import pandas as pd


LABEL_SIZE = 8

def classify_enhancer(df):
	vec_size = len(df.columns) - LABEL_SIZE
	features = df[df.columns[0:vec_size]]
	labels_matrix = df[df.columns[vec_size:vec_size + LABEL_SIZE]]
	enhancer_labels = labels_matrix[labels_matrix.columns[0]]
	m = KNeighborsClassifier(n_neighbors=3)
	m.fit(features, enhancer_labels.tolist())
	#scores = cross_val_score(m, features, enhancer_labels.tolist(), cv=5)
	#print scores
	#print "Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2)
	predictions = m.predict(features)
	print classification_report(enhancer_labels.tolist(), predictions)

def read_pickle():
	k = 6
	filename = "pickled_vista_data_k_{0}.pkl".format(k)
	data = pd.read_pickle(filename)
	return data

def main():
	df = read_pickle()
	classify_enhancer(df)

main()
