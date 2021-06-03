import sys
import numpy as np
import pandas as pd
import qnorm
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn import linear_model
from sklearn import metrics

__version__="1.0"

class data_object:
	samples=[]
	groups={}	
	data_n=0
	omics=[]
	genelist=[]


# get label information of samples
def get_sample_info(motypes, train_labels, gene_l):
	data_info = data_object()
	data_info.samples=pd.DataFrame(train_labels, columns=['Sampid', 'Class'])
	data_info.data_n=train_labels.shape[0]
	(group, counts) = np.unique(train_labels.iloc[:,1], return_counts=True)
	data_info.groups=dict(zip(group, counts))
	data_info.omics=list(motypes)
	data_info.genelist=gene_l
	return data_info

## log2 & quantile normalization
def proc_qnorm(d):	
	d_ln=np.log2(d+1)
	d_q=qnorm.quantile_normalize(d_ln, axis=1, ncpus=8)
	return d_q

## min-max scaling
def proc_scale(d):
	scaler = MinMaxScaler()
	d_s = MinMaxScaler(feature_range=(0, 1)).fit(d.T).transform(d.T).T
	df_s = pd.DataFrame(d_s, columns = list(d.columns))
	df_s=df_s.round(6)
	return df_s

def make_mir_gcentric(mirf, mrnaf, mirtarget):
	sep=','
	print('reading mirna file...')
	rf = open(mirf)
	mi_barcodes = rf.readline().strip().split(sep)[1:]
	miExprs = dict()
	while(True):
		line = rf.readline()
		if not line:
			break
		token = line.strip().split(sep)
		miName = token[0].strip()
		miExprs[miName] = list(map(lambda x: float(x), token[1:]))
	rf.close()

	rf = open(mrnaf)
	gene_barcodes = rf.readline().strip().split(sep)[1:]
	geneExprs = dict()
	while(True):
		line = rf.readline()
		if not line:
			break
		token = line.strip().split(sep)
		geneName = token[0].strip()
		geneExprs[geneName] = list(map(lambda x: float(x), token[1:]))
	rf.close()

	miNames = miExprs.keys()
	geneNames = geneExprs.keys()

	for miName in miNames:
		newLst = []
		oldLst = miExprs[miName]
		for barcode in mi_barcodes:
			newLst.append(oldLst[mi_barcodes.index(barcode)])
		miExprs[miName] = newLst

	for geneName in geneNames:
		newLst = []
		oldLst = geneExprs[geneName]
		for barcode in mi_barcodes:
			if barcode not in gene_barcodes: continue
			newLst.append(oldLst[gene_barcodes.index(barcode)])
		geneExprs[geneName] = newLst

	print('reading mirna gene target file...')
	rf = open(mirtarget)
	targetRevDict = dict()
	while(True):
		line = rf.readline()
		if not line:
			break
		token = line.strip().split('\t')
		miName = token[0].strip()
		geneName = token[1].strip()
		if geneName not in targetRevDict:
			targetRevDict[geneName] = [miName]
		else:
			targetRevDict[geneName].append(miName)
	rf.close()

	geneNames = geneExprs.keys()
	mir_data=[]

	def writeList(lst):
		row_data=[]
		for i in range(len(lst)):
			row_data.append(float(lst[i]))
		mir_data.append(row_data)

	print('converting miRNAs to gene level data...')
	gene_list=[]
	for geneNameIdx, geneName in enumerate(geneNames):
		gene_list.append(geneName)
		if geneName not in targetRevDict:
			writeList([0.0 for _ in range(len(mi_barcodes))])
			continue

		miNames = targetRevDict[geneName]
		miNames = list(filter(lambda x: x in miExprs, miNames))
		lst = []
		for i in range(len(mi_barcodes)):
			miAvg = 0.0
			num = 0.0
			for j in range(len(miNames)):
				num += 1.0
				miAvg = miAvg * (1.0 - 1.0/num) + miExprs[miNames[j]][i] / num  # running average

			lst.append(miAvg)
		arr = np.asarray(lst)
		lst = list(arr)
		writeList(lst)

	return pd.DataFrame(mir_data, columns=mi_barcodes, index=gene_list)
	


def make_methylation_gcentric(methf, mrnaf, gtid, methpromprob):
	sep=','
	rf = open(methf)
	meth_barcodes = rf.readline().strip().split(sep)[1:]
	methBetaValues = dict()
	print('reading methylation data...')
	while(True):
		line = rf.readline()
		if not line:
			break
		tok = line.strip().split(sep)
		pid = tok[0].strip()	# probe id

		if 'NA' in tok[1:]:
			methBetaValues[pid]=[-1.0 if x=='NA' else float(x) for x in tok[1:]]
		else:
			methBetaValues[pid]=[float(x) for x in tok[1:]]

	rf.close()

	rf = open(mrnaf)
	gene_barcodes = rf.readline().strip().split(sep)[1:]
	geneExprs = dict()
	while(True):
		line = rf.readline()
		if not line:
			break
		token = line.strip().split(sep)
		geneName = token[0].strip()
		geneExprs[geneName] = list(map(lambda x: float(x), token[1:]))
	rf.close()

	probe_ids = methBetaValues.keys()
	geneNames = geneExprs.keys()
	for pid in probe_ids:
		newLst = []
		oldLst = methBetaValues[pid]
		for barcode in meth_barcodes:
			newLst.append(oldLst[meth_barcodes.index(barcode)])
		methBetaValues[pid] = newLst

	for geneName in geneNames:
		newLst = []
		oldLst = geneExprs[geneName]
		for barcode in meth_barcodes:
			if barcode not in gene_barcodes: continue
			newLst.append(oldLst[gene_barcodes.index(barcode)])
		geneExprs[geneName] = newLst

	# tid to gid converstion dictionary
	gid_dict={}
	for line in open(gtid):
		gid, tid=line.strip().split('\t')
		gid_dict[tid]=gid

	print('reading TSS promoter probe annotation...')
	rf = open(methpromprob)
	targetRevDict = dict()
	while(True):
		line = rf.readline()
		if not line:
			break
		tok = line.strip().split('\t')
		pid = tok[0].strip()
		tid=tok[4].strip()
		if tid not in gid_dict:
			continue
		geneName=gid_dict[tid]
		if geneName not in targetRevDict:
			targetRevDict[geneName] = [pid]
		else:
			targetRevDict[geneName].append(pid)
	rf.close()

	geneNames = geneExprs.keys()
	meth_data=[]
	def writeList(lst):
		row_data=[]
		for i in range(len(lst)):
			row_data.append(float(lst[i]))
		meth_data.append(row_data)


	print('converting methylation to gene level data...')
	gene_n=len(geneExprs.keys())
	gene_list=[]
	for geneNameIdx, geneName in enumerate(geneNames):
		gene_list.append(geneName)
		if geneName not in targetRevDict:
			writeList([0.0 for _ in range(len(meth_barcodes))])
			continue

		pid = targetRevDict[geneName]
		pid = list(filter(lambda x: x in methBetaValues, pid))
		lst = []
		for i in range(len(meth_barcodes)):
			methAvg = 0.0
			num = 0.0
			for j in range(len(pid)):
				if methBetaValues[pid[j]][i] == -1: continue	# skip NA values
				num += 1.0
				methAvg = methAvg * (1.0 - 1.0/num) + methBetaValues[pid[j]][i] / num 	# running average

			lst.append(float(methAvg))

		writeList(lst)

	return pd.DataFrame(meth_data, columns=meth_barcodes, index=gene_list)


# decomposing omics tensor using PARAFAC
def decompose_tensor(data, r, maxi, tolv):
	import tensorly as tl
	from tensorly.decomposition import parafac
	tl.set_backend('numpy')

	##
	## To choose a subset of omics types, the data can be modified as below.
	## The order of omics are based on our examples, hence please select the indexes according to the customized omics order
	##

	## pair of omics as data
	# data=data[[0,1],:,:]	# GE, ME
	# data=data[[0,2],:,:]	# GE, MI
	# data=data[[1,2],:,:]	# ME, MI

	## single omics as data
	# data=data[0,:,:]	# GE
	# data=data[1,:,:]	# ME
	# data=data[2,:,:]	# MI
	# print(data.shape)

	# performing tensor decomposition (PARAFAC)
	components = tl.decomposition.non_negative_parafac(data, rank=r, n_iter_max=maxi, init='random', tol=10e-8)

	# save decomposed components	
	# np.save(outf, components.factors)

	return components.factors



def L1_regression(rank, st, X_, y_, a):
	acc_list=[]
	max_coef=[]
	max_acc=.0
	max_model=''
	ncross=10
	nzfeats_comm=set()
	topfeats_comm=set()
	top_n=max(int(rank*0.05), 1) # select top 5% features, select at least 1 feature
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

def select_features(P, alpha, dat_info, odir, cvidx):	
	rank=P.shape[-1]
	label_counts=[0]
	count=0
	for labkey, labval in dat_info.groups.items():
		count+=labval
		label_counts.append(count)

	rep_n=0
	while True:
		subtype_coef=np.zeros(shape=(len(dat_info.groups),rank))	# empty subtype coefficient array
		
		subtype_models=[]
		subtype_features=[]
		subtype_topfeatures=[]

		subtype_labels=list(dat_info.groups.keys())
		min_feat_n=rank
		# alpha=info.alpha	# L1 alpha penalty
		for i in range(len(subtype_labels)):
			y=np.zeros(label_counts[-1])
			y[label_counts[i]:label_counts[i+1]]=1

			coef_, acc_, max_model_, feat_n_, nzfeats_, topn_=L1_regression(rank, subtype_labels[i], P, y, alpha)
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

		print("Selected features:\t%d\n"%(len(selfeats)))

		# write selected features to file
		stfeatures=[]
		feat_f=open("%s/cvs/cv_%d/out/sample_features.txt"%(odir, cvidx), "w")
		for idx, i in enumerate(subtype_topfeatures):
			for fidx in sorted(i):
				feat_f.write("%d\t%s\n"%(fidx, subtype_labels[idx]))
				stfeatures.append([fidx, subtype_labels[idx]])
		feat_f.close()
		stdf=pd.DataFrame(stfeatures, columns=['FeatureIndex', 'Group'])

		return list(selfeats), stdf


def get_featuregenes(G, cv, dat_info, odir, cvidx, rank):
	# load selected features
	st_feats={}
	feats=set()
	for line in open("%s/cvs/cv_%d/out/sample_features.txt"%(odir, cvidx)):
		tok=line.strip().split()
		fidx=int(tok[0])
		st=tok[1]
		if st not in st_feats: st_feats[st]=[]
		st_feats[st].append(fidx)
		feats.add(fidx)
	feats=sorted(list(feats))

	label_counts=[0]
	count=0
	for labkey, labval in dat_info.groups.items():
		count+=labval
		label_counts.append(count)

	# select gene with max feat value
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
	print('Selected genes:	%d'%(len(sel_genes)))
	
	feature_genes=[]
	gout=open("%s/cvs/cv_%d/out/feature_genes.txt"%(odir, cvidx),"w")
	for i in dat_info.groups:
		for g in st_genes[i]:
			gout.write('%s,%s\t%d\t%s\n'%(dat_info.genelist.iloc[g,0], dat_info.genelist.iloc[g,1], g, i))
			feature_genes.append([dat_info.genelist.iloc[g,0], dat_info.genelist.iloc[g,1], g, i])
		print("  - %s: %s genes"%(i, len(st_genes[i])))
	gout.close()
	ftg=pd.DataFrame(feature_genes, columns=['ENSGID', 'GeneSymbol', 'GeneIndex', 'Group'])

	return sel_genes, ftg


def measure_accuracy(cv, sel_genes, dat_info, odir, cvidx):

	# load train data
	rawdata=cv.train_data
	seldat=rawdata[:,sel_genes,:]
	X_train=[]
	for p in range(dat_info.data_n):
		if seldat.ndim>=3:
			d=seldat[:,:,p].flatten()
		else:
			d=seldat[:,p].flatten()
		X_train.append(d)
	X_train=np.asarray(X_train)
	y_train=cv.train_labels['Class'].tolist()

	# load test data
	test_np=cv.test_data
	test_np=test_np[:,sel_genes,:]
	test_labs=cv.test_labels['Class'].tolist()
	
	X_test=[]
	for p in range(test_np.shape[-1]):
		if test_np.ndim>=3:
			d=test_np[:,:,p].flatten()
		else:
			d=test_np[:,p].flatten()
		X_test.append(d)
	X_test=np.asarray(X_test)
	
	y_test=test_labs

	clf_mlp = MLPClassifier(solver='adam', alpha=1e-5, max_iter=5000).fit(X_train, y_train)
	y_pred = clf_mlp.predict(X_test)
	crep=metrics.classification_report(y_test, y_pred)
	
	acc_f=open('%s/cvs/cv_%d/out/accuracy_result.txt'%(odir, cvidx), 'w')
	acc_f.write(crep)
	acc_f.close()

	return crep

