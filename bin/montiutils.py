import numpy as np
import pandas as pd
from pathlib import Path
import qnorm
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from sklearn.preprocessing import MinMaxScaler
from sklearn.manifold import TSNE


class cv_object:
	train_data=[]
	test_data=[]
	train_labels=[]
	test_labels=[]
	samp_info=None
	tddat=None
	feature_samples=None
	feature_sample_labels=None
	feature_genes=None
	gene_labels=None
	predres=None


def make_env(odir, cv_n):
	for i in range(cv_n):
		Path('%s/cvs/cv_%d/out/components'%(odir, i)).mkdir(parents=True, exist_ok=True)
		Path('%s/cvs/cv_%d/out/plots'%(odir, i)).mkdir(parents=True, exist_ok=True)
		Path('%s/cvs/cv_%d/out/gene_models'%(odir, i)).mkdir(parents=True, exist_ok=True)
		Path('%s/cvs/cv_%d/out/patient_models'%(odir, i)).mkdir(parents=True, exist_ok=True)


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

# splitting samples into train and test sets
def split_data(modata, class_labels, genelist, odir, cv_n):
	samp_info=class_labels
	
	cv_lst=[]
	for i in range(cv_n):
		print('generating cross validation data set %d...'%(i))
		cvdat=cv_object()

		tratio=0.10	# ratio of test sample
		test_idx=sorted(np.random.randint(samp_info.shape[0], size=int(samp_info.shape[0]*tratio)))	# index of test samples
		test=samp_info.iloc[test_idx,:].sort_values('Class')
		train=samp_info.iloc[~samp_info.index.isin(test_idx),:].sort_values('Class')

		test.to_csv('%s/cvs/cv_%d/test.list'%(odir, i), index=False)
		train.to_csv('%s/cvs/cv_%d/train.list'%(odir, i), index=False)

		omics=['gene', 'methylation', 'mirna']
		df_train_norm_set=[]
		df_test_norm_set=[]
		for o in omics:
			df=modata[o]

			# gene filter
			# df=df[df.index.isin(genelist['ENSGID'])]
			df.loc[genelist['ENSGID']]
			df_train=df.loc[:, train.Sampid]

			# log2 quantile normalization
			if o in ['gene', 'mirna']:
				df_train=np.log2(df_train+1)
				df=np.log2(df+1)

			df_train_norm=qnorm.quantile_normalize(df_train, axis=1, ncpus=8)
			df_norm=qnorm.quantile_normalize(df, axis=1, ncpus=8)

			# scaling to 0~1
			df_train_norm=proc_scale(df_train_norm)
			df_train_norm.index=df_train.index
			df_train_norm_set.append(df_train_norm.values)
			
			df_norm=proc_scale(df_norm)
			df_norm.index=df.index
			
			df_test_norm=df_norm.loc[:, test['Sampid'].tolist()]			
			df_test_norm_set.append(df_test_norm.values)

			df_train_norm.to_csv('%s/cvs/cv_%d/%s_train.txt'%(odir, i, o))
			df_test_norm.to_csv('%s/cvs/cv_%d/%s_test.txt'%(odir, i, o))


		# merge the omics file into a single tensor object
		cvdat.train_data=np.asarray(df_train_norm_set)
		cvdat.test_data=np.asarray(df_test_norm_set)

		cvdat.train_labels=train.reset_index(drop=True)
		cvdat.test_labels=test.reset_index(drop=True)

		cv_lst.append(cvdat)

	return cv_lst


##### plotting functions

color_palette=["darkcyan", "blue", "crimson", "orange", "yellow", "saddlebrown", "turquoise"]

def plot_sample_features(P, fts, cv_train, dat_info):
	selfeats=fts['FeatureIndex']
	X_sel=P[:,selfeats]
	tsne = TSNE(n_components=2, verbose=0, n_iter=1000, early_exaggeration=5)
	transformed = tsne.fit_transform(X_sel)
	xs = transformed[:,0]
	ys = transformed[:,1]

	plt.title("t-SNE")
	for gidx, g in enumerate(dat_info.groups):
		idx=cv_train[cv_train['Class']==g].index.tolist()
		plt.scatter(x=xs[idx],y=ys[idx], c=color_palette[gidx], s=20, label=g)

	leg=plt.legend(fancybox=False, bbox_to_anchor=(1, 1), loc="upper right")
	leg.get_frame().set_edgecolor('black')

	plt.xlabel("Component 1")
	plt.ylabel("Component 2")

	return plt

def plot_venn(gl, dat_info):
	from venn import venn
	gdict={}
	for gidx, g in enumerate(dat_info.groups):
		gdict[g]=set(gl.GeneSymbol[gl['Group']==g].tolist())

	v=venn(gdict, cmap=color_palette, fontsize=10, legend_loc="upper right")
	plt.title("Gene Venn diagram")

	return v

# heatmaps
def plot_sample_hmap(P, fs, dat_info):
	dat=P[:,fs]
	group_palette= dict(zip(dat_info.groups.keys(), color_palette[:len(dat_info.groups.keys())]))
	rcols=list(map(group_palette.get, dat_info.samples.Class.tolist()))

	g=sns.clustermap(dat, method="complete", cmap="vlag", standard_scale=0, row_colors=rcols)
	handles = [Patch(facecolor=group_palette[g]) for g in group_palette]
	plt.legend(handles, group_palette, title=None, bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right', frameon=False)

	ax = g.ax_heatmap
	ax.set_xlabel("Features")
	ax.set_ylabel("Samples")
	
	return plt


def plot_gene_hmap(G, O, fs, fg, dat_info):
	dat=G[fg, :]
	dat=dat[:, fs]

	omics_dict = dict(zip(np.arange(len(dat_info.omics)), dat_info.omics))
	omics_palette= dict(zip(omics_dict.keys(), color_palette[:len(dat_info.omics)]))
	omics_col=dict(zip(dat_info.omics, color_palette[:len(dat_info.omics)]))

	ccols=[]
	for i in range(O.shape[1]):
		maxf=np.argmax(O[:,i])
		ccols.append(omics_palette[maxf])

	g=sns.clustermap(dat, method="average", cmap="vlag", standard_scale=0, col_colors=ccols)
	handles = [Patch(facecolor=omics_col[g]) for g in omics_col]
	plt.legend(handles, omics_col, title=None, bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right', frameon=False)

	ax = g.ax_heatmap
	ax.set_xlabel("Features")
	ax.set_ylabel("Genes")
	
	return plt



def plot_gene(modata, cv, gidx, dat_info):

	ensgid=modata['gene'].index[gidx]
	gsym=dat_info.genelist[dat_info.genelist['ENSGID']==ensgid]['GeneSymbol'].values[0]

	x=np.arange(dat_info.data_n)
	# y=modata['gene'].iloc[gidx,:]	# pair 1
	# z=modata['methylation'].iloc[gidx,:]	# pair 2
	y=cv.train_data[0, gidx, :]	# pair 1
	z=cv.train_data[1, gidx, :]	# pair 2
	
	row=1
	col=3
	fig=plt.figure(figsize=(8,8))

	ax=fig.add_subplot(1, 1, 1, projection= '3d')

	# first pair
	ax.plot(x, y, '.', c='green', zdir='z', markersize=1)
	p2 = np.poly1d(np.polyfit(x, y, 2))
	ax.plot(x, p2(x), '-', c='black', zdir='z', linewidth=1)

	# second pair
	ax.plot(x, z, '.', c='blue', zdir='y', zs=1.0, markersize=1)
	p2 = np.poly1d(np.polyfit(x, z, 3))
	ax.plot(x, p2(x), '-', c='black', zdir='y', linewidth=1.2, zs=1.0)

	# annotate subtypes with line bars
	s_ind=0
	for gidx, g in enumerate(dat_info.groups):
		idx=cv.train_labels[cv.train_labels['Class']==g].index.tolist()
		ax.plot([idx[0], idx[-1]], [1,1], '-', c=color_palette[gidx], zdir='y', linewidth=3, zs=1.03, label=g) # subtype bars

	ax.set_ylim([0, 1])
	ax.set_zlim([0, 1])

	ax.set_xlabel('Samples')
	ax.set_ylabel('%s'%(dat_info.omics[0]))
	ax.set_zlabel('%s'%(dat_info.omics[1]))

	leg=plt.legend(fancybox=False, bbox_to_anchor=(1, 1), loc="upper right")
	leg.get_frame().set_edgecolor('black')

	plt.title('%s (%s)'%(gsym, ensgid))

	plt.show()

	return

	# pcorstr="%s (%s)\nr = %g\np = %g"%(gsym, g_st_dict[gidx], pcorgenes[g][ridx+5], pcorgenes[g][ridx+1])

	# ax.text2D(0.91, -0.02, s=pcorstr, fontsize=10, transform=ax.transAxes)

	# annotate subtypes with line bars
	# s_ind=0
	# for stidx, st in enumerate(dat_info.samp_labels):
	# 	e_ind=s_ind+dat_info.samp_labels[st]
	# 	ax.plot([s_ind,e_ind], [1,1], '-', c=color_palette[stidx], zdir='y', linewidth=3, zs=1.03) # subtype bars
	# 	s_ind+=dat_info.samp_labels[st]

	# ax.set_xlabel('Samples')
	# ax.set_ylabel('%s'%(omics_label[ridx][0]))
	# ax.set_zlabel('%s'%(omics_label[ridx][1]))

	# ax.set_ylim([0, 1])
	# ax.set_zlim([0, 1])



# survival plots



