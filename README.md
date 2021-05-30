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

## Install MONTI
* Python version>=3.6 is required
* The python modules below are required which can be installed manually or using the `install_monti.py` script
  * `tensorly`, `argparse`, `joblib`, `matplotlib`, `lifelines`, `seaborn`, `qnorm`

## Running MONTI
The input to MONTI are gene-level omics matrices. The omics matrices must be of the same dimensions, which are matrices of _m x n_. Here, _m_ is the number of genes and _n_ the number of patients. Hence, every omics data are bundled into units of genes to generate a gene-level matrix.

A well documented jupyter notebook is available in the link below.
https://github.com/inukj/MONTI/blob/main/monti_tutorial_coad.ipynb


The tutorial includes the full procedure of analyzing the colon cancer subtype data (TCGA-COAD), including raw data processing, gene-level transofrmation and subtype associated gene selection.
Also, various plotting functions are available for visualizing the result.






