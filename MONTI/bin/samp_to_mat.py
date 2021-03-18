#!/usr/bin/env python
import sys
import subprocess

import numpy as np
import pandas as pd
import qnorm

from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA

from arg_parser import *

delim=','

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

# check genelist
def get_genelist(i):
	dat=pd.read_csv(info.fin[0], header=0, index_col=0, sep=delim)	# genes in data	
	pcgenelist=pd.read_csv(info.geneinfo, header=0, sep='\t')		# restrict to protein coding genes
	gl=pcgenelist[pcgenelist['Gene stable ID'].isin(dat.index)]		# final genelist

	return gl

if __name__=='__main__':
	info = samp2mat()	# parsing input arguements

	# create output directory
	cmd="mkdir -p %s"%(info.outdir)
	proc=subprocess.Popen(cmd, shell=True, executable='/bin/bash')
	proc.wait()

	# read sampleinfo data
	sampdat=pd.read_csv(info.sampinfo, header=None, sep=delim)

	# check and cleanup gene list (restrict to protein coding genes)
	geneinfo=get_genelist(info)	

	# check data integrity - filter out samples not in sampinfo file
	print('checking data integrity...')
	for f in info.fin:
		fname=f.split('/')[-1]
		matdat=pd.read_csv(f, header=0, index_col=0, sep=delim)
		indatsamp=[]
		for idx, i in sampdat.iterrows():
			if i[0] in matdat.columns:
				indatsamp.append(i[0])
				
		## gene selection
		dat=matdat[indatsamp]
		dat=dat.loc[geneinfo['Gene stable ID'],]	# select gene subset by geneannotationdat file

		print('normalizing data... %s, %s'%(f, dat.shape))
		## log2 & quantile normalization
		dat_q=proc_qnorm(dat)

		## scaling data
		df_scaled=proc_scale(dat_q)
		df_scaled.index=dat.index
		df_scaled.to_csv('%s/%s.%s'%(info.outdir, fname, info.group), sep=',')

	# write final sample info
	out_sinfo=open('%s/sampinfo_%s.txt'%(info.outdir, info.group), 'w')
	for idx, i in sampdat.iterrows():
		if i[0] in indatsamp:
			out_sinfo.write('%s\n'%(delim.join([str(x) for x in i])))
	out_sinfo.close()

	# write gene info
	geneinfo.to_csv('%s/geneinfo_%s.txt'%(info.outdir, info.group), index=False, sep=',')
	
	# merge the omics file into a single tensor object
	print('merging final omics matrices into a tensor...')
	tensor_dat=[]
	for f in info.fin:
		fname=f.split('/')[-1]
		dat=pd.read_csv('%s/%s.%s'%(info.outdir, fname, info.group), header=0, index_col=0, sep=',')
		tensor_dat.append(dat.values)
	np.save('%s/tensor.%s.npy'%(info.outdir, info.group), tensor_dat)   # write tensor to file

