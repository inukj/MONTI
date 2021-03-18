import sys
import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import pearsonr
from sklearn.manifold import TSNE

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from mpl_toolkits.mplot3d import Axes3D
from utils import *
from itertools import chain, combinations

color_palette=["darkcyan", "blue", "crimson", "orange", "yellow", "saddlebrown", "turquoise"]

def print_header():
	print("\n------------------------------------------\n")

def powerset(iterable):
	# powerset of list s
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def print_progress (iteration, total, prefix='', suffix='', decimals=1, bar_length=50):
	str_format = "{0:." + str(decimals) + "f}"
	percents = str_format.format(100 * (iteration / float(total)))
	filled_length = int(round(bar_length * iteration / float(total)))
	bar = '#' * filled_length + '-' * (bar_length - filled_length)
	sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
	if iteration == total:
		sys.stdout.write('\n')
	sys.stdout.flush()

def plot_genes(info, dat_info):
	continue_check(info, 'plotting')
	print_header()

	rawdata=np.load(info.fin)
	genelist=np.loadtxt("%s"%(info.geneinfo), delimiter=',', skiprows=1, dtype=str)

	st_dict={}
	gset=set()
	for line in open("%s/feature_genes_r%d.txt"%(info.outdir, info.rank)):
		tok=line.strip().split('\t')
		gidx=int(tok[1])
		st=(tok[2])
		if st not in st_dict:
			st_dict[st]=[]
		st_dict[st].append(gidx)
		gset.add(gidx)

	st_assoc_genes={}
	g_st_dict={}
	for gidx in list(gset):
		st_label=""
		for st in dat_info.samp_labels:
			if gidx in st_dict[st]: 
				st_label+="%s,"%(st)
		st_label=st_label.strip(",")

		if st_label not in st_assoc_genes:
			st_assoc_genes[st_label]=[]
		st_assoc_genes[st_label].append(gidx)
		g_st_dict[gidx]=st_label

	st_assoc_labels=[]
	st_pset=list(powerset(dat_info.samp_labels.keys()))
	for i in st_pset:
		if len(i)>0: 
			st_assoc_labels.append(','.join(i))

	for stlab in st_assoc_labels:
		if stlab not in st_assoc_genes: continue
		pcorgenes=[]
		pval=0.05/len(st_assoc_genes[stlab])	# Bonferroni correction
		
		for idx, i in enumerate(st_assoc_genes[stlab]):
			ge=rawdata[0,i,:]
			me=rawdata[1,i,:]
			mi=rawdata[2,i,:]

			ge_me=pearsonr(ge, me)
			ge_mi=pearsonr(ge, mi)
			me_mi=pearsonr(me, mi)

			# select genes with at least one omics relation with significant p-value
			if ge_me[1]<pval or ge_mi[1]<pval:
				min_p=np.min([ge_me[1], ge_mi[1], me_mi[1]])
				max_r=np.max([ge_me[0], ge_mi[0], me_mi[0]])
				min_p_oid=np.argmin([ge_me[1], ge_mi[1], me_mi[1]])
				pcorgenes.append([i, ge_me[1], ge_mi[1], me_mi[1], min_p, ge_me[0], ge_mi[0], me_mi[0], max_r, min_p_oid])


		# NO PCORGENES -- pass
		if len(pcorgenes)==0 : continue 	# pass if no pcorgenes

		pcorgenes=np.asarray(pcorgenes)
		gsorted=np.argsort(pcorgenes[:,4])
		print("plotting %s genes... "%(stlab))

		omics_pair=[[0,1],[0,2]]
		omics_label=[["Gene","Methylation"],["Gene","miRNA"]]
		omics_color=["red", "blue", "green"]

		# plot figures
		row=4
		col=3
		pos=1
		fsx=20
		fsy=20
		fig=plt.figure(figsize=(fsx,fsy))
		pdf_out = matplotlib.backends.backend_pdf.PdfPages("%s/plots/gene_plots_%s.pdf"%(info.outdir, stlab))

		print_progress(0, len(gsorted), prefix = 'Progress:', suffix = 'Complete')
		for gx, g in enumerate(gsorted):
			if pos==0:
				pdf_out.savefig(fig)
				fig.clf()
				pos=1
				fig=plt.figure(figsize=(fsx,fsy))

			gidx=int(pcorgenes[g][0])
			ridx=np.argmin(pcorgenes[g][1:3])
			x=np.asarray(range(rawdata.shape[2]))
			y=rawdata[omics_pair[ridx][0],gidx,:]	# pair 1
			z=rawdata[omics_pair[ridx][1],gidx,:]	# pair 2

			ax=fig.add_subplot(row, col, pos, projection= '3d')

			# first pair
			ax.plot(x, y, '.', c=omics_color[omics_pair[ridx][0]], zdir='z', markersize=1)
			p2 = np.poly1d(np.polyfit(x, y, 2))
			ax.plot(x, p2(x), '-', c='black', zdir='z', linewidth=1.2)

			# second pair
			ax.plot(x, z, '.', c=omics_color[omics_pair[ridx][1]], zdir='y', zs=1.0, markersize=1)
			p2 = np.poly1d(np.polyfit(x, z, 3))
			ax.plot(x, p2(x), '-', c='black', zdir='y', linewidth=1.2, zs=1.0)

			if genelist[gidx].size>1:
				gsym=genelist[gidx][1]
			else:
				gsym=genelist[gidx][0]
			
			pcorstr="%s (%s)\nr = %g\np = %g"%(gsym, g_st_dict[gidx], pcorgenes[g][ridx+5], pcorgenes[g][ridx+1])

			ax.text2D(0.91, -0.02, s=pcorstr, fontsize=10, transform=ax.transAxes)

			# annotate subtypes with line bars
			s_ind=0
			for stidx, st in enumerate(dat_info.samp_labels):
				e_ind=s_ind+dat_info.samp_labels[st]
				ax.plot([s_ind,e_ind], [1,1], '-', c=color_palette[stidx], zdir='y', linewidth=3, zs=1.03) # subtype bars
				s_ind+=dat_info.samp_labels[st]

			ax.set_xlabel('Samples')
			ax.set_ylabel('%s'%(omics_label[ridx][0]))
			ax.set_zlabel('%s'%(omics_label[ridx][1]))

			ax.set_ylim([0, 1])
			ax.set_zlim([0, 1])

			pos=(pos+1)%((col*row)+1)

			print_progress(gx + 1, len(gsorted), prefix = 'Progress:', suffix = 'Complete')

		pdf_out.savefig(fig)
		pdf_out.close()

		plt.close('all')


def plot_sample_tSNE(info, dat_info):
	continue_check(info, 'plotting')
	sys.stdout.write("drawing sample tSNE plots... ")
	sys.stdout.flush()

	selfeats=[]
	for line in open('%s/sample_features_r%s.txt'%(info.outdir, info.rank)):
		tok=line.strip().split('\t')
		selfeats.append(int(tok[0]))

	tddata=np.load("%s/components/r%d_td.npy"%(info.outdir, info.rank), allow_pickle=True)
	X=tddata[2]

	X_sel=X[:,selfeats]
	tsne = TSNE(n_components=2, verbose=0, n_iter=1000, early_exaggeration=5)
	transformed = tsne.fit_transform(X_sel)
	xs = transformed[:,0]
	ys = transformed[:,1]

	plt.title("t-SNE")
	s_ind=0
	for stidx, st in enumerate(dat_info.samp_labels):
		e_ind=s_ind+dat_info.samp_labels[st]
		plt.scatter(x=xs[s_ind:e_ind],y=ys[s_ind:e_ind],c=color_palette[stidx], s=20, label=st)
		s_ind+=dat_info.samp_labels[st]

	lim_val=np.max(np.abs(reject_outliers(transformed)))*2	# removing outliers for better visual of plot
	
	plt.xlim(-lim_val,lim_val)
	plt.ylim(-lim_val,lim_val)
	leg=plt.legend(fancybox=False, bbox_to_anchor=(1, 1), loc="upper right")
	leg.get_frame().set_edgecolor('black')
	plt.savefig('%s/plots/sample_tSNE.pdf'%(info.outdir))
	plt.close('all')

	print('done')


