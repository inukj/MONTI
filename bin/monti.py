#!/usr/bin/env python
import sys

from arg_parser import *
from utils import *
from decompose_tensor import *
from featureselection import *
from geneselection import *
from drawplots import *

if __name__=='__main__':

	info = parse_args()	# parsing input arguements

	valid=sanity_check(info)	# check input data
	if not valid: sys.exit();

	progress_check(info)
	make_env(info)

	dat_info=get_sample_info(info)	# get information on sample labels
	print_init(info, dat_info)		# print init messages

	decompose_tensor(info)			# run PARAFAC
	select_features(info, dat_info)	# selecting features
	get_featuregenes(info, dat_info) # selecting feature associated genes

	if info.plot:
		plot_genes(info, dat_info)			# draw plots
		plot_sample_tSNE(info, dat_info)	# draw sample t-SNE plot

	



