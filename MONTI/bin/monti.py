#!/usr/bin/env python
import sys

from arg_parser import *
from utils import *
from decompose_tensor import *
from featureselection import *
from geneselection import *
from drawplots import *
from preprocess_data import *

if __name__=='__main__':

	info = parse_args()	# parsing input arguements

	valid=sanity_check(info)	# check input data
	if not valid: sys.exit();	# ERROR, some input files are incorrect

	progress_check(info)

	preprocess(info)	# preprocess data (log2, quantile normalization, tensor formation)

	make_env(info)	# create output directory

	dat_info=get_sample_info(info)	# get information on sample labels
	print_init(info, dat_info)	# print init messages

	decompose_tensor(info)	# run PARAFAC

	select_features(info, dat_info)	# select cancer subtype associated features

	get_featuregenes(info, dat_info)	# select cancer subtype associated genes

	if info.plot:
		plot_genes(info, dat_info)		# draw plots
		plot_sample_tSNE(info, dat_info)	# draw sample t-SNE plot

	



