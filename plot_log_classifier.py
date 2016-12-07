from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, roc_curve
from sklearn import metrics
#from sklearn.model_selection import cross_val_score
#if this import complains, comment it out and uncomment the import above
from sklearn.cross_validation import cross_val_score, KFold
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import interp


# LABEL_SIZE = 8
LABEL_SIZE = 1

def classify_enhancer(df):
	vec_size = len(df.columns) - LABEL_SIZE
	features = df[df.columns[0:vec_size]]
	labels_matrix = df[df.columns[vec_size:vec_size + LABEL_SIZE]]
	enhancer_labels = labels_matrix[labels_matrix.columns[0]]

	X = features.values
	y = list(enhancer_labels.values)

	kf = KFold(n=len(y), n_folds=5)

	tprs = []
	base_fpr = np.linspace(0, 1, 101)

	plt.figure(figsize=(5, 5))

	for i, (train, test) in enumerate(kf):
		model = LogisticRegression(C=500.0, class_weight='balanced').fit(X[train], [y[i] for i in train])
		y_score = model.predict_proba(X[test])
		fpr, tpr, _ = roc_curve([y[i] for i in test], y_score[:, 1])

		plt.plot(fpr, tpr, 'b', alpha=0.15)
		tpr = interp(base_fpr, fpr, tpr)
		tpr[0] = 0.0
		tprs.append(tpr)

	tprs = np.array(tprs)
	mean_tprs = tprs.mean(axis=0)
	std = tprs.std(axis=0)

	tprs_upper = np.minimum(mean_tprs + std, 1)
	tprs_lower = mean_tprs - std


	plt.plot(base_fpr, mean_tprs, 'b')
	plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)

	plt.plot([0, 1], [0, 1],'r--')
	plt.xlim([-0.01, 1.01])
	plt.ylim([-0.01, 1.01])
	plt.ylabel('True Positive Rate')
	plt.xlabel('False Positive Rate')
	plt.axes().set_aspect('equal', 'datalim')
	plt.savefig('enh_fb_log_random4000_k_3.png')
	plt.close()

	print 'auROC Curve for 5-fold CV:'
	print metrics.auc(base_fpr, mean_tprs)

def read_pickle():
	# k = 3
	# filename = "pickled_vista_data_k_{0}.pkl".format(k)
	filename = 'enh_fb_random4000_k_3.pkl'
	data = pd.read_pickle(filename)
	return data

def main():
	df = read_pickle()
	classify_enhancer(df)

main()
