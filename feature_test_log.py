from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
#from sklearn.model_selection import cross_val_score
#if this import complains, comment it out and uncomment the import above
from sklearn.cross_validation import cross_val_score
from sklearn import metrics
import pandas as pd


LABEL_SIZE = 8

def show_most_informative_features(clf, n=20):
    feature_names = [i for i in range(45)]
    coefs_with_fns = sorted(zip(clf.coef_[0], feature_names))
    top = zip(coefs_with_fns[:n], coefs_with_fns[:-(n + 1):-1])
    for (coef_1, fn_1), (coef_2, fn_2) in top:
        print "\t%.4f\t%-15s\t\t%.4f\t%-15s" % (coef_1, fn_1, coef_2, fn_2)

def classify_enhancer(df):
	vec_size = len(df.columns) - LABEL_SIZE
	features = df[df.columns[0:vec_size]]
	labels_matrix = df[df.columns[vec_size:vec_size + LABEL_SIZE]]
	enhancer_labels = labels_matrix[labels_matrix.columns[0]]
	m = LogisticRegression(class_weight="balanced")
	m.fit(features, enhancer_labels.tolist())
	predictions = m.predict(features)
	print classification_report(enhancer_labels.tolist(), predictions)
	fpr, tpr, thresholds = metrics.roc_curve(enhancer_labels.tolist(), predictions)
	print metrics.auc(fpr, tpr)
	#m.fit(features, enhancer_labels)
	scores = cross_val_score(m, features, enhancer_labels.tolist(), cv=5)
	print scores
	print "Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2)
	show_most_informative_features(m)

def read_pickle():
	k = 4
	filename = "noepi_pickled_vista_data_k_{0}.pkl".format(k)
	data = pd.read_pickle(filename)
	return data

def main():
	df = read_pickle()
	classify_enhancer(df)

main()
