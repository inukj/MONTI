# MONTI: A Multi-Omics Non-negative Tensor decomposition framework for the Integrated analysis of large cancer cohorts

Multi-omics data is frequently measured to characterize biological mechanisms underlying phenotypes. Complex relationships in multi-omics data, if mined, can lead to more accurate classification of patient samples according to the phenotypes.

MONTI (Multi-Omics Non-negative Tensor decomposition for Integrative analysis) is a tool that can be used to integrate and analyze large sets of multi-omics data. MONTI identifies gene regulatory multi-omics features specific to a group of samples that share a common biological trait.

Below is an illustration of the analysis workflow of MONTI.
![workflow](./images/monti_workflow.jpg)

The output of MONTI is a simple gene list with information of their associated subtypes, which can be used for further downstream analysis. For example, the Venn diagram below shows the genes that are found to be associated to colorectal cancer subtypes CMS1, CMS2, CMS3 and CMS4. These genes showed to be informative in separating the four subtypes as shown in the t-SNE plot.
<!--![example output](./images/monti_outputexample.png =250x)-->
<img src="./images/monti_outputexample.png" alt="example output" width="600"/>

---

<!-- ## Download MONTI
```bash
git clone https://github.com/inukj/MONTI.git
``` -->

## Install MONTI
MONTI is developed in python3 and can be installed as below
```bash
pip install monti
```

## Tutorial using colon cancer data (TCGA-COAD)
A brief tutorial for using MONTI can be found under the ['tutorial'](https://github.com/inukj/MONTI/tree/main/tutorial) directory.

Before starting the tutorial, the dataset should be downloaded.
After download decompress data by
```bash
cd <download_path>
tar -xzvf tutorial_data_coad.tar.gz
```

The *<download_path>* should also be used as the tutorial directory, or you can simply move the data to another directory to be used for the tutorial.

The data includes three omics data, 1) gene expression (mRNA), 2) methylation level and 3) miRNA expression.
They are raw data directly collected from the TCGA portal.

In the [jupyter notebook](https://github.com/inukj/MONTI/blob/main/tutorial/tutorial_coad.ipynb) shows an example of how to integrate multi-omics data in a gene-level manner and extract features that can well distinguish the molecular subtypes of COAD.

The tutorial includes the below analysis procedures:
* gene-level transformation
* normalization
* feature selection
* classification accuracy measurement and
* plotting of the results










