# MONTI: A Multi-Omics Non-negative Tensor decomposition framework for the Integrated analysis of large cancer cohorts

Multi-omics data is frequently measured to characterize biological mechanisms underlying phenotypes. Complex relationships in multi-omics data, if mined, can lead to more accurate classification of patient samples according to the phenotypes.

MONTI (Multi-Omics Non-negative Tensor decomposition for Integrative analysis) is a tool that can be used to integrate and analyze large sets of multi-omics data. MONTI identifies gene regulatory multi-omics features specific to a group of samples that share a common biological trait.

Below is an illustration of the analysis workflow of MONTI.
![workflow](./images/monti_workflow.jpg)

The output of MONTI is a simple gene list with information of their associated subtypes, which can be used for further downstream analysis. For example, the Venn diagram below shows the genes that are found to be associated to colorectal cancer subtypes CMS1, CMS2, CMS3 and CMS4. These genes showed to be informative in separating the four subtypes as shown in the t-SNE plot.
<!--![example output](./images/monti_outputexample.png =250x)-->
<div style="text-align:center"><img src="./images/monti_outputexample.png" alt="example output" width="600"/></div>

## Using the STAD data

### Generate input data for MONTI
```bash
monti_input_dir="dataset/STAD/inputdata"
bin/samp_to_mat.py -i dataset/STAD/data/omics_STAD_gene_genecentric.csv dataset/STAD/data/omics_STAD_meth450_genecentric.csv dataset/STAD/data/omics_STAD_mirna_genecentric.csv -s dataset/STAD/subtype_info.txt -r subtype -g dataset/gene_info_withheader.txt -o $monti_input_dir
```

### running MONTI
```bash
rank=150	# number of ranks (features)
monti_outputdir="dataset/STAD/output"
bin/monti.py -f $monti_input_dir/tensor.subtype.npy -s $monti_input_dir/sampinfo_subtype.txt -g $monti_input_dir/geneinfo_subtype.txt -r $rank -o $monti_outputdir --plot
```

## Using the COAD data

### Generate input data for MONTI
``` bash
monti_input_dir="dataset/COAD/inputdata"
bin/samp_to_mat.py -i dataset/COAD/data/omics_COAD_gene_genecentric.csv dataset/COAD/data/omics_COAD_meth450_genecentric.csv dataset/COAD/data/omics_COAD_mirna_genecentric.csv -s dataset/COAD/data/subtype_info.txt -r subtype -g dataset/gene_info_withheader.txt  -o $monti_input_dir
```

### running MONTI
``` bash
rank=150
monti_outputdir="dataset/COAD/output"
bin/monti.py -f $monti_input_dir/tensor.subtype.npy -s $monti_input_dir/sampinfo_subtype.txt -g $monti_input_dir/geneinfo_subtype.txt -r $rank -o $monti_outputdir --plot
```
