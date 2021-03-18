import sys
import subprocess
import numpy as np
from pathlib import Path

__version__="1.0"

class data_info:
	samp_labels={}
	samp_caseid={}
	samp_caseidL=[]
	gene_ensid={}
	data_n=0

def print_init(info, dat_info):
	print("\n------------------------------------------\n")
	print("Starting MONTI (v%s)\n"%(__version__))
	
	print("Sample classes: [%s]"%(', '.join(dat_info.samp_labels)))
	print("Samples: %d"%(dat_info.data_n))
	print("Omics: %d"%(info.omics_n))
	print("Genes: %d"%(len(dat_info.gene_ensid)))
	
	print()
	print("Rank: %d"%(info.rank))
	print("alpha: %g"%(info.alpha))

	print("Output directory: \"%s\""%(info.outdir))
	print("\n------------------------------------------\n")

# create environment for storing results
def make_env(info):
	# # clear output directory
	# if info.decomp:
	# 	cmd="rm %s"%(info.outdir)
	# 	proc=subprocess.Popen(cmd, shell=True, executable='/bin/bash')
	# 	proc.wait()

	# make output directories
	for d in ['components', 'patient_models', 'gene_models', 'plots']:
		cmd="mkdir -p %s/%s"%(info.outdir, d)
		proc=subprocess.Popen(cmd, shell=True, executable='/bin/bash')
		proc.wait()

# get label information of samples
def get_sample_info(info):
	dat_info = data_info()
	sampdat=np.genfromtxt(info.sampinfo, delimiter=',', dtype=str)
	label_count=np.unique(sampdat[:, 1], return_counts=True)

	for i, j in zip(label_count[0], label_count[1]): 
		dat_info.samp_labels[i] = j

	# load sample case ids
	for ldx, line in enumerate(open("%s"%(info.sampinfo))):
		tok=line.strip().split(',')
		caseid=tok[0]
		st=tok[1]
		dat_info.samp_caseid[ldx]=(caseid, st)
		dat_info.samp_caseidL.append(caseid)
	dat_info.data_n=len(dat_info.samp_caseid)

	# load gene annotations
	for idx, ensid in enumerate(open('%s'%(info.geneinfo))):
		tok=ensid.strip().split('\t')
		dat_info.gene_ensid[idx]=tok[0]
	
	return dat_info

def sanity_check(info):

	print("\n------------------------------------------\n")
	sys.stdout.write("checking input files... ")
	sys.stdout.flush()

	strout=[]
	valid=1
	for fl in info.input_files:	# read metadata files
		if fl==None: 
			continue
		try:
			with open(fl, 'r') as fh:
				pass;
		except FileNotFoundError:
			valid=0
			strout.append("[ERROR] File does not exist: %s"%(fl))

	# try loading input tensor file
	fname=info.fin
	indata = Path(fname)
	if indata.is_file():
		dat=np.load(indata)
		info.omics_n=dat.shape[0]

	if valid:
		print('done')
		print('Everything looks ok.')
	else:
		print('ERROR')
		for e in strout:
			print(e)

	return valid

def continue_check(info, phase):
	if not info.prog: 
		print("Something went wrong before \'%s\'. Exiting..."%(phase))
		sys.exit()	# something is wrong

def progress_check(info):
	info.prog=0
	# check decomposition file (npy)
	flag=1	
	info.decompfile='r%d_td.npy'%(info.rank)
	fname="%s/components/%s"%(info.outdir, info.decompfile)
	tdfile = Path(fname)
	if tdfile.is_file():
		print("\n------------------------------------------\n")
		info.decomp=1	# default: override
		override=input("Previous tensor decomposition results exist. Override? [y/n]: ")
		if override!='y':
			info.decomp=0
	else:
		flag=0

	info.prog=1













