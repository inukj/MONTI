import sys
import numpy as np
np.set_printoptions(suppress=True)
import tensorly as tl
from tensorly.decomposition import parafac
tl.set_backend('numpy')

from utils import *

# decomposing tensor file using PARAFAC
def decompose_tensor(info):
	continue_check(info, 'tensor decomposition')
	if not info.decomp: 
		sys.stdout.write('Skipping tensor decomposition.\n')
		return;
	sys.stdout.write("decomposing tensor (may take some time)... ")
	sys.stdout.flush()

	data=np.load(info.fin)
	
	##
	## To choose a subset of omics types, the data can be modified as below.
	## The order of omics are based on our examples, hence please select the indexes according to the customized omics order
	##

	## pair of omics as data
	# data=data[[0,1],:,:]	GE, ME
	# data=data[[0,2],:,:]	# GE, MI
	# data=data[[1,2],:,:]	# ME, MI

	## single omics as data
	# data=data[0,:,:]	# GE
	# data=data[1,:,:]	# ME
	# data=data[2,:,:]	# MI
	print(data.shape)

	# performing tensor decomposition (PARAFAC)
	components = tl.decomposition.non_negative_parafac(data, rank=info.rank, n_iter_max=info.dmaxiter, init='random', tol=10e-8)

	# save decomposed components	
	f="%s/components/%s"%(info.outdir, info.decompfile)
	np.save(f, components.factors)

	print('done')

