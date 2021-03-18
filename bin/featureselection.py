import warnings
warnings.filterwarnings("ignore")
import time
import numpy as np
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
import joblib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from utils import *

def print_header():
	print("\n------------------------------------------\n")
	sys.stdout.write("selecting sample features... ")
	sys.stdout.flush()

def L1_regression(rank, st, X_, y_, a):
	acc_list=[]
	max_coef=[]
	max_acc=.0
	max_model=''
	ncross=10
	nzfeats_comm=set()
	topfeats_comm=set()
	top_n=max(int(rank*0.1), 1) # select top 10% features, select at least 1 feature
	for i in range(ncross):
		crossratio=float(1)/ncross
		X_train, X_test, y_train, y_test = train_test_split(X_, y_, test_size=crossratio)

		reg = linear_model.Lasso(alpha = a, max_iter=5000)	# L1 regularization
		reg.fit(X_train, y_train)
		
		pred=reg.predict(X_test)
		pred_label=np.copy(pred)
		pred_label[pred_label<0.5]=0
		pred_label[pred_label>=0.5]=1
		test_label=np.asarray(y_test)

		acc=np.equal(pred_label, test_label)
		ratio=np.sum(acc)/acc.size
		acc_list.append(ratio)

		if ratio > max_acc: 
			max_acc=ratio
			max_coef=reg.coef_
			max_model=reg
			selfeat_n=np.count_nonzero(reg.coef_)

		nzfeats,=np.nonzero(reg.coef_)
		for f in list(nzfeats): nzfeats_comm.add(f)

		top_f=np.argsort(np.abs(reg.coef_))[-top_n:]
		for f in list(top_f): topfeats_comm.add(f)

	return max_coef, max_acc, max_model, selfeat_n, nzfeats_comm, topfeats_comm

def select_features(info, dat_info):
	continue_check(info, 'feature selection')
	print_header()
	strout=[]

	rank=info.rank
	
	tddata=np.load("%s/components/%s"%(info.outdir, info.decompfile), allow_pickle=True)

	X=tddata[-1]	# sample slice

	label_counts=[0]
	count=0
	for labkey, labval in dat_info.samp_labels.items():
		count+=labval
		label_counts.append(count)

	rep_n=0
	while True:
		subtype_coef=np.zeros(shape=(len(dat_info.samp_labels),rank))	# empty subtype coefficient array
		
		subtype_models=[]
		subtype_features=[]
		subtype_topfeatures=[]

		subtype_labels=list(dat_info.samp_labels.keys())
		min_feat_n=rank
		alpha=info.alpha	# L1 alpha penalty
		for i in range(len(subtype_labels)):
			y=np.zeros(label_counts[-1])
			y[label_counts[i]:label_counts[i+1]]=1

			coef_, acc_, max_model_, feat_n_, nzfeats_, topn_=L1_regression(info.rank, subtype_labels[i], X, y, alpha)
			subtype_coef[i]=coef_		
			subtype_models.append(max_model_)
			subtype_features.append(nzfeats_)
			subtype_topfeatures.append(topn_)
			
			if feat_n_ < min_feat_n: 
				min_feat_n=feat_n_

		unionfeats=set()
		for i in range(0,len(subtype_labels)):
			for j in range(0,len(subtype_labels)):
				if j>i:
					cf=subtype_topfeatures[i].intersection(subtype_topfeatures[j])
					uf=subtype_topfeatures[i].union(subtype_topfeatures[j])
					ratio=len(cf)/len(uf)
					for f in uf: unionfeats.add(f)
		selfeats=unionfeats

		X_sel=tddata[-1][:,list(selfeats)]	# sample slice
		strout.append("Selected subtype features:\t%d"%(len(selfeats)))

		y=[]
		for samplab, sampval in data_info.samp_labels.items():
			y+=[samplab]*sampval

		# Build multilayer perceptron classifier
		ncross=10
		acc_list=[]
		for i in range(ncross):
			crossratio=float(1)/ncross
			X_train, X_test, y_train, y_test = train_test_split(X_sel, y, test_size=crossratio)
			clf_mlp = MLPClassifier(solver='adam', alpha=1e-5, max_iter=5000).fit(X_train, y_train)
			score_mlp=clf_mlp.score(X_test, y_test)
			acc_list.append(score_mlp)
			
			filename = "%s/patient_models/MLP_r%d_i%d.mod"%(info.outdir, info.rank, i)
			joblib.dump(clf_mlp, filename)
		acc_avg=np.mean(acc_list)

		# save results
		feat_f=open("%s/sample_features_r%d.txt"%(info.outdir, info.rank), "w")
		for idx, i in enumerate(subtype_topfeatures):
			for fidx in sorted(i):
				feat_f.write("%d\t%s\n"%(fidx, subtype_labels[idx]))
		feat_f.close()

		acc_outf=open("%s/accuracy_patients_r%d.txt"%(info.outdir, info.rank), "w")
		for idx, score in enumerate(acc_list):
			acc_outf.write("%d\t%f\n"%(idx, score))
		strout.append("Classification accuracy:\t%f"%(acc_avg))
		acc_outf.write("average\t%f\n"%(acc_avg))
		acc_outf.close()

		print("done\n")
		for s in strout: print("%s"%(s))

		break;


