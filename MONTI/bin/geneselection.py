import sys
import warnings
warnings.filterwarnings("ignore")
import joblib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.manifold import TSNE

import survival_mod as surv
from utils import *

def print_header():
	print("\n------------------------------------------\n")
	sys.stdout.write("selecting feature genes... ")
	sys.stdout.flush()

def get_featuregenes(info, dat_info):
	continue_check(info, 'gene selection')
	print_header()
	strout=[]

	# load selected features
	st_feats={}
	feats=set()
	for line in open("%s/sample_features_r%d.txt"%(info.outdir, info.rank)):
		tok=line.strip().split()
		fidx=int(tok[0])
		st=tok[1]
		if st not in st_feats: st_feats[st]=[]
		st_feats[st].append(fidx)
		feats.add(fidx)
	feats=sorted(list(feats))

	# load raw tensor data
	tddata=np.load("%s/components/%s"%(info.outdir, info.decompfile), allow_pickle=True)
	
	label_counts=[0]
	count=0
	for labkey, labval in dat_info.samp_labels.items():
		count+=labval
		label_counts.append(count)

	# select gene with max feat value
	G=tddata[-2]	# gene slice
	feat_genes={}
	st_genes={}
	sel_genes=set()
	for gid in range(G.shape[0]):
		maxf=np.argmax(G[gid,:])
		if maxf in feats:
			if maxf not in feat_genes:
				feat_genes[maxf]=[]
			feat_genes[maxf].append(gid)
			sel_genes.add(gid)

		for st in st_feats:
			if maxf in st_feats[st]:
				if st not in st_genes:
					st_genes[st]=[]
				st_genes[st].append(gid)

	sel_genes=list(sorted(sel_genes))
	strout.append("Total %d genes selected."%(len(sel_genes)))
	
	gout=open("%s/feature_genes_r%d.txt"%(info.outdir, info.rank),"w")
	for i in dat_info.samp_labels:
		for g in st_genes[i]:
			gout.write("%s\t%d\t%s\n"%(dat_info.gene_ensid[g], g, i))
		strout.append("  - %s: %s genes"%(i, len(st_genes[i])))
	gout.close()
	rawdata=np.load(info.fin)
	seldat=rawdata[:,sel_genes,:]

	X=[]
	for p in range(dat_info.data_n):
		d=seldat[:,:,p].flatten()
		X.append(d)
	
	y=[]
	for samplab, sampval in data_info.samp_labels.items():
		y+=[samplab]*sampval

	ncross=10
	acc_list=[]
	for i in range(ncross):
		crossratio=float(1)/ncross
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=crossratio)
		clf_mlp = MLPClassifier(solver='sgd', alpha=1e-5, max_iter=5000).fit(X_train, y_train)
		score_mlp=clf_mlp.score(X_test, y_test)
		acc_list.append(score_mlp)
		filename = "%s/gene_models/MLP_r%d_i%d.mod"%(info.outdir, info.rank, i)
		joblib.dump(clf_mlp, filename)
	acc_avg=np.mean(acc_list)

	acc_outf=open("%s/accuracy_genes_r%d.txt"%(info.outdir, info.rank), "w")
	for idx, score in enumerate(acc_list):
		acc_outf.write("%d\t%f\n"%(idx, score))
	strout.append("Classification accuracy:\t%f"%(acc_avg))
	acc_outf.write("average\t%f\n"%(acc_avg))
	acc_outf.close()

	print('done')
	for s in strout: print("%s"%(s))


