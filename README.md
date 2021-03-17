# MONTI tutorial


## Using the STAD data

### Generate input data for MONTI
` bash
monti_input_dir="dataset/STAD/inputdata"
bin/samp_to_mat.py -i dataset/STAD/data/omics_STAD_gene_genecentric.csv dataset/STAD/data/omics_STAD_meth450_genecentric.csv dataset/STAD/data/omics_STAD_mirna_genecentric.csv -s dataset/STAD/subtype_info.txt -r subtype -g dataset/gene_info_withheader.txt  -o $monti_input_dir
`

### running MONTI
``` bash
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
time bin/monti.py -f $monti_input_dir/tensor.subtype.npy -s $monti_input_dir/sampinfo_subtype.txt -g $monti_input_dir/geneinfo_subtype.txt -r $rank -o $monti_outputdir --plot
```
