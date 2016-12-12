from sklearn.svm import SVC
from sklearn.metrics import roc_curve, accuracy_score, auc, classification_report, average_precision_score
from sklearn.cross_validation import KFold
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import interp


# Must change this to correspond to label length! Check the dataset!
LABEL_SIZE = 8
# LABEL_SIZE = 1

# Much thanks to Alexey Grigorev at 
# https://stats.stackexchange.com/questions/186337/average-roc-for-repeated-10-fold-cross-validation-with-probability-estimates/
# for the bulk of this code!
def classify_enhancer(df):
	vec_size = len(df.columns) - LABEL_SIZE
	features = df[df.columns[0:vec_size]]
	labels_matrix = df[df.columns[vec_size:vec_size + LABEL_SIZE]]
	# labels are [enhancer, brain, forebrain, midbrain, hindbrain, limb, neural tube, heart] for vista data,
	# just [enhancer] for Beer lab data
	enhancer_labels = labels_matrix[labels_matrix.columns[2]]

	X = features.values
	y = list(enhancer_labels.values)

	kf = KFold(n=len(y), n_folds=5)

	tprs = []
	base_fpr = np.linspace(0, 1, 101)

	plt.figure(figsize=(5, 5))

	accuracies = []
	preds = []
	scores = []

	for i, (train, test) in enumerate(kf):
		print 'Evaluating fold ', i
		model = SVC(C=10000.0, class_weight='balanced', probability=True).fit(X[train], [y[i] for i in train])
		y_score = model.predict_proba(X[test])
		y_pred = model.predict(X[test])
		preds.extend(y_pred)
		scores.extend(y_score[:, 1])
		accuracies.append(accuracy_score([y[i] for i in test], y_pred))
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
	plt.savefig('pickled_vista_data_fb_k_6.png')
	plt.close()

	print 'Classification report (5-fold CV):'
	print classification_report(y, preds)

	print 'Mean accuracy (5-fold CV):'
	print np.mean(accuracies)

	print 'auPRC curve for 5-fold CV:'
	print average_precision_score(y, scores)

	print 'auROC Curve for 5-fold CV:'
	print auc(base_fpr, mean_tprs)

def read_pickle():
	filename = 'pickled_vista_data_k_6.pkl'
	data = pd.read_pickle(filename)
	return data

def main():
	df = read_pickle()
	classify_enhancer(df)

main()
