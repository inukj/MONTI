import sys
import argparse

class input_info:
	args = {}
	fin = ''
	indir = ''
	outdir = ''
	rank = ''
	sampinfo = ''
	geneinfo=''
	omics_n=0
	survival = ''
	plot = ''
	dmaxiter = ''
	alpha = ''
	pre = ''
	decomp = 1
	prog = 0
	decompfile=''
	group = ''
	intput_files=[]
	top_n=0

def parse_args():
	info = input_info()
	parser = argparse.ArgumentParser(prog="monti.py")

	# required arguements
	parser.add_argument('-f', '--input_file', type=str, default=100, help='input tensor file', required=True)
	parser.add_argument('-r', '--rank', type=int, default=50, help='number of rank features', required=True)
	parser.add_argument('-s', '--sample_info', type=str, help='label information of samples', required=True)
	parser.add_argument('-g', '--gene_info', type=str, help='ordered list of genes', required=True)

	# optional arguements
	parser.add_argument('--alpha', type=float, default=0.01, help='alpha value for L1 regularization',required=False)
	parser.add_argument('-o', '--outdir', type=str, default='output', help='output directory', required=False)
	parser.add_argument('--plot', action='store_true', help='plot genes', required=False)
	parser.add_argument('--dmax_iter', type=int, default=300, help='maximum iteration for tensor decomposition',required=False)
	
	args = vars(parser.parse_args())
	info.args = args
	
	info.fin = args['input_file']
	info.indir = info.fin.rsplit("/")[0]
	info.outdir = args['outdir']
	info.rank = args['rank']
	info.sampinfo = args['sample_info']
	info.geneinfo = args['gene_info']
	info.survival = args['survival_info']
	info.plot = args['plot']
	info.dmaxiter = args['dmax_iter']
	info.alpha = args['alpha']
	info.pre = args['preprocess_dir']
	info.input_files=[args['sample_info'], args['gene_info'], args['survival_info']]

	return info


def samp2mat():
	info = input_info()
	parser = argparse.ArgumentParser(prog="samp_to_mat.py")

	# required arguements	
	parser.add_argument('-i', '--exprmat', type=str, nargs='+', help='list of all available omics samples', required=True)
	parser.add_argument('-s', '--sample_info', type=str, help='label information of samples', required=True)
	parser.add_argument('-g', '--gene_info', type=str, help='list of genes to be used', required=True)
	parser.add_argument('-r', '--group_name', type=str, help='name of group', required=True)

	# optional arguements
	parser.add_argument('-o', '--outdir', type=str, default='input_data', help='output directory', required=False)

	args = vars(parser.parse_args())
	info.args = args
	
	info.fin = args['exprmat']
	info.outdir = args['outdir']
	info.sampinfo = args['sample_info']
	info.geneinfo = args['gene_info']
	info.group = args['group_name']
	
	return info
