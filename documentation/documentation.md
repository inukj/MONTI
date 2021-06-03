# MONTI v1.2.1

The MONTI is a python program for analyzing multi-omics data and extracting features associated with certain clinical attributes or phenotypes.

This document describes the functions and objects used by MONTI.

## Functions

###make_env(@outdir)
> creates the output directory where results will be stored
- *@outdir* (STRING) is the path name of the directory to store the results
> returns NULL


###make_methylation_gcentric(@methmat, @genemat, @gid_tid, @methprobes)
> converts the methylation raw data into gene-level data using methylation probes in the gene promoter regions
- *@methmat* (STRING) is the path to the raw methylation value file
- *@genemat* (STRING) is the path to the raw gene expression (mRNA) file
- *@gid_tid* (STRING) is the path to the mapping information between transcript ID and Gene ID
> returns a gene-level methylation dataframe


###make_mir_gcentric(@mirmat, @genemat, @mir_gene)
> converts the miRNA raw expression data into gene-level data using the miRNA-gene target annotation data
- *@mirmat* (STRING) is the path to the raw miRNA expression file
- *@genemat* (STRING) is the path to the raw gene expression (mRNA) file
- *@mir_gene* (STRING) is the path to the miRNA-gene target file
> returns a gene-level miRNA dataframe


###split_data(@modata, @label_data, @genelist, @outdir, @cv_n, @test_data_ratio)
> creates train and test data sets
- *@modata* (DICT) is a dictionary object including all the raw omics data
- *@label_data* (DF) is a dataframe containing the clincal data per sample
- *@genelist* (DF) is a dataframe containing the genes to used for the analysis
- *@outdir* (STRING) the path to the resul directory
- *@cv_n* (INT) the number of train/test data sets to be generated (e.g., if to perform 10-cross validation cv_n=10)
- *@test_data_ratio* (FLOAT) the ratio size of the test data to be used, which is usually 1/@cv_n
> returns a *cv_object* object that contains the entire cross validation dataset


###decompose_tensor(@train_data, @rank, @maxiter, @tol_val)
> decomposes the input multi-omics tensor data
- *@train_data* (np.array) is a numpy object containing the multi-omics train data
- *@rank* (INT) the number of ranks to decompose the tensor
- *@maxiter* (INT) the maximum number of iterations for tensor decomposition (higher the less will be the residual)
- *@tol_val* (float) the minimum tolerance value for convergence during decomposition (lower the less will be the residual)
> returns the decompose omics *O*, gene *G* and sample *P* components, which are np.array objects


###get_sample_info(@omics, @train_labels, @genelist)
> fetches the information of train data
- *@omics* (LIST) is a list describing the type of omics in order (e.g., ['gene', 'methylation', 'mirna'])
- *@train_labels* (DF) is a dataframe including the labels of each sample in the train data
- *@genelist* (DF) is a dataframe containing the genes to used for the analysis
> returns a *data_object* that contains various information of the train data


###select_features(@P, @alpha, @sample_info, @oudir, @i)
> selects features from the sample *P* component using the train data information
- *@P* (np.array) is an array of the decomposed sample component
- *@alpha* (FLOAT) is the L1 alpha penalty value used during the feature selection (a larger @alpha will penalize more and thus select less features)
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function
- *@outdir* (STRING) the path to the resul directory
- *@i* (INT) the @i'th train data set from which features will be extracted
> returns a list of selected features


###get_featuregenes(@G, @[i], @sample_info, @oudir, @i, @rank)
> selectes feature genes based on the selected features from *select_features*
- *@G* (np.array) is an array of the decomposed gene component
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function
- *@outdir* (STRING) the path to the resul directory
- *@i* (INT) the @i'th train data set from which features will be extracted
- *@rank* (INT) the number of ranks to decompose the tensor
> returns a list of the feature associated genes


###measure_accuracy(@test_data, @feature_genes, @sample_info, @outdir, @i)
> measure several accuracy metrics of classyfing the sample clinical features using the selected features genes 
- *@test_data* (np.array) is a numpy object containing the multi-omics test data
- *@feature_genes* (LIST) is a list of feature associated genes
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function
- *@outdir* (STRING) the path to the resul directory
- *@i* (INT) the @i'th train data set from which features will be extracted
> returns the measured accuracy metrics


###plot_sample_features(@P, @feature_labels, @train_labels, @sample_info)
> plots the tSNE plot on the sample *P* component only using the selected features
- *@P* (np.array) is an array of the decomposed sample component
- *@feature_labels* (DF) is a dataframe containing the labels of each selected feature
- *@train_labels* (DF) is a dataframe containing the labels of the train data
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function


###plot_venn(@gene_labels, @sample_info)
> plots the Venn diagram showing the status of feature genes
- *@gene_labels* (DF) is a dataframe containing the labels of the genes
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function


###plot_sample_hmap(@P, @feature_samples, @sample_info)
> plots the heatmap on the sample *P* component using only the selected features
- *@P* (np.array) is an array of the decomposed sample component
- *@feature_samples* (DF) is a dataframe containing the labels of each selected feature
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function


###plot_gene_hmap(@G, @O, @feature_samples, @feature_genes, @sample_info)
> plots the heatmap on the gene *G* component using only the feature genes
- *@G* (np.array) is an array of the decomposed gene component
- *@O* (np.array) is an array of the decomposed omics component
- *@feature_genes* (DF) is a dataframe containing the labels of each selected feature gene
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function


###plot_gene(@modata, @train_data, @gidx, @sample_info)
> plots the raw expression data of two omics for visualizing correlation
- *@modata* (DICT) is a dictionary object including all the raw omics data
- *@train_data* (np.array) is a numpy object containing the multi-omics train data
- *@gidx* (INT) the index of a gene in genelist
- *@sample_info* (data_object) contains information of the train samples, which is returned by the *get_sample_info* function



## Objects

class cv_object:
	train_data=[]
	test_data=[]
	train_labels=[]
	test_labels=[]
	samp_info=None
	tddat=None
	feature_samples=None
	feature_labels=None
	feature_genes=None
	gene_labels=None
	predres=None


class data_object:
	samples=[]
	groups={}	
	data_n=0
	omics=[]
	genelist=[]


