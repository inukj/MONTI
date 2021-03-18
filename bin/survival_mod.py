import warnings
warnings.filterwarnings("ignore")
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from sklearn.cluster import AgglomerativeClustering
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from sklearn.preprocessing import LabelEncoder


font = {'family': 'sans-serif',
	        'color':  'black',
	        'weight': 'normal',
	        'size': 11}

class Patient(object):
	def __init__(self, id):
		self.pid=id
		self.stat=""
		self.dead=0
		self.time=-1
		self.group=-1

def survival_analysis(X, component, info, dinfo):
	tmpdat={}	# patient survival data
	for ldx, line in enumerate(open(info.survival)):
		tok=line.strip().split(',')
		pid=tok[0]
		tmpdat[pid]=Patient(pid)
		tmpdat[pid].stat=int(tok[1])
		tmpdat[pid].time=int(tok[2])
		# patients[pid].group=cluster_labels[ldx]
	
	# match with sample list
	clinical_data={}	# patient survival data
	pat_sel=[]
	for pat in dinfo.samp_caseid:
		bid=dinfo.samp_caseid[pat][0]
		print(pat, bid)
		if bid in tmpdat:
			clinical_data[bid]=Patient(bid)
			clinical_data[bid].stat=tmpdat[pid].stat
			clinical_data[bid].time=tmpdat[pid].time
			# pat_sel.append(pat)
			pat_sel.append(dinfo.samp_caseid[pat][1])

	# X_=X[pat_sel]
	# print(X)
	print(pat_sel)
	pat_sel_encode = LabelEncoder().fit_transform(pat_sel)

	# perform hclust on patients with survival info
	hclust = AgglomerativeClustering(n_clusters=4, affinity='euclidean',linkage='ward', compute_full_tree=True)
	cluster_labels = hclust.fit_predict(X_)

	print(cluster_labels)

	sys.exit()

	T=[]
	E=[]
	groups=[]
	for p, v in patients.items():
		T.append(v.time)
		E.append(v.stat)
		groups.append(v.group)
	T=np.array(T)
	E=np.array(E)
	groups=np.array(groups)
	c1=(groups == 0)
	c2=(groups == 1)

	fig, ax = plt.subplots(figsize=(8,8), nrows=1, ncols=1, squeeze=True)
	kmf_c1 = KaplanMeierFitter()
	kmf_c1.fit(T[c1], E[c1], label='c1').plot(ax=ax, ci_show=False, c='lime')

	kmf_c2 = KaplanMeierFitter()
	kmf_c2.fit(T[c2], E[c2], label='c2').plot(ax=ax, ci_show=False, c='red')

	ax.legend(edgecolor='black', fancybox=False)
	add_at_risk_counts(kmf_c1, kmf_c2, ax=ax)
	results=logrank_test(T[c1], T[c2], E[c1], E[c2], alpha=0.99)

	plt.text(.1, .1, "p=%g"%(results.p_value), transform=ax.transAxes, fontdict=font)
	plt.ylim(0, 1)
	plt.savefig("%s/survival_%s_plot.pdf"%(info.outdir, component))

	return results.p_value
