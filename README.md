# MONTI: A Multi-Omics Non-negative Tensor decomposition framework for the Integrated analysis of large cancer cohorts

Multi-omics data is frequently measured to characterize biological mechanisms underlying phenotypes. Complex relationships in multi-omics data, if mined, can lead to more accurate classification of patient samples according to the phenotypes.

MONTI (Multi-Omics Non-negative Tensor decomposition for Integrative analysis) is a tool that can be used to integrate and analyze large sets of multi-omics data. MONTI identifies gene regulatory multi-omics features specific to a group of samples that share a common biological trait.

Below is an illustration of the analysis workflow of MONTI.
![workflow](./images/monti_workflow.jpg)

The output of MONTI is a simple gene list with information of their associated subtypes, which can be used for further downstream analysis. For example, the Venn diagram below shows the genes that are found to be associated to colorectal cancer subtypes CMS1, CMS2, CMS3 and CMS4. These genes showed to be informative in separating the four subtypes as shown in the t-SNE plot.
<!--![example output](./images/monti_outputexample.png =250x)-->
<img src="./images/monti_outputexample.png" alt="example output" width="600"/>

---

## Download MONTI
```bash
git clone https://github.com/inukj/MONTI.git
```

## Installing MONTI
* Python version>=3.6 is required
* The python modules below are required which can be installed manually or using the `install_monti.py` script
  * `tensorly`, `argparse`, `joblib`, `matplotlib`, `lifelines`, `seaborn`, `qnorm`

## Running MONTI
The input to MONTI are gene-centric omics matrices. The omics matrices must be of the same dimensions, which are matrices of _m x n_. Here, _m_ is the number of genes and _n_ the number of patients. Hence, every omics data are bundled into units of genes to generate a gene-centric matrix.

Assuming that the omics are processed into gene-centric format, each omics data are further pre-processed. This can be done using the `samp_to_mat.py` code, which compiles omics into patient matched data and normalizes data accordingly. For example, if mRNA, methylation and miRNA omics data are given, patients with all the three omics data are selected, which will serve as the data used by MONTI. The output of `samp_to_mat.py` is the input data used by MONTI. Users may prepare their own omics data.


```bash
usage: samp_to_mat.py [-h] -i OMICS1_FILE OMICS2_FILE OMICS3_FILE ... -s SAMPLE_INFO -g GENE_INFO
	[-l GROUP_LABEL] [-o OUTDIR]

# mandatory arguements
-i: the omics matrices in CSV format
-s: a two column text file that contains sample IDs and their clinical features (e.g., subtype)
-g: a two column text file including the genes to be used (by default a list or 14K protein coding genes are provided)

# optional arguements
-l: the label of the group (just an indicator)
-o: the output path
```

Now MONTI can be run using the input data as follows:
```bash
usage: monti.py [-h] -f INPUT_FILE -r RANK -s SAMPLE_INFO -g GENE_INFO
	[-surv SURVIVAL_INFO] [-o OUTDIR] [--plot]
	[--dmax_iter DMAX_ITER] [--alpha ALPHA]
	[-pre PREPROCESS_DIR]

# mandatory arguements
-f: the input tensor data (a numpy ndarray)
-r: the number of ranks which the tensor is to be decomposed (type: int)
-s: a two column text file that contains sample IDs and its associated breast cancer subtype
-g: a two column text file including the genes to be used (by default a list or 14K protein coding genes are provided)

# optional arguements
--damx_iter: the number of maximum iterations during tensor decomposition (default: 300)
--alpha: the L1 penalty weight (type: float, default: 0.01)
--plot:	 if set tSNE and omics correlation plots are drawn
-o:	 the output directory name (default: 'output')
```

## Output of MONTI
In the output directory the following output files can be found. 
* sample_features_<rank>.txt: the index of selected features and their labels
* feature_genes_<rank>.txt: the feature associated genes
* accuracy_patients_r<rank>.txt: the classification accuracy using patient features
* accuracy_genes_r<rank>.txt: the classification accuracy using fea
* patient_models/: the classification models generated using the patient features
* gene_models/: the classification models generated using the gene features
* plots/: 
	gene_plots_<subtype>.pdf: multi-omics scatter plot of the feature associated genes
	sample_tSNE.pdf: the t-SNE plot of the patient features

---

## Brief MONTI tutorial using COAD omics data (mRNA, methylaion, miRNA)
Below is a brief tutorial using the COAD mRNA, methylation and miRNA omics data, which will reproduce the results in our study.

1. Get the data
```bash
cd <path to MONTI>/dataset
tar -xzvf COAD_data.tar.gz
```

2. Generate input data
``` bash
cd <path to MONTI>
monti_input_dir="dataset/COAD/inputdata/"
bin/samp_to_mat.py -i dataset/COAD/data/omics_COAD_gene_genecentric.csv dataset/COAD/data/omics_COAD_meth450_genecentric.csv dataset/COAD/data/omics_COAD_mirna_genecentric.csv -s dataset/COAD/data/subtype_info.txt -l subtype -g annotation/gene_info_withheader.txt -o $monti_input_dir
```

3. Run MONTI
>The tensor decomposition takes time, hence, we will use a pre-decomposed data. The pre-decomposed data is located in `dataset/COAD/output/components/r150_td.npy`.
``` bash
rank=150
monti_outputdir="dataset/COAD/output"
bin/monti.py -f $monti_input_dir/tensor.subtype.npy -s $monti_input_dir/sampinfo_subtype.txt -g $monti_input_dir/geneinfo_subtype.txt -r $rank -o $monti_outputdir --plot
```

