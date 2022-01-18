import os
import pandas as pd
from glob import glob

## change sample_set to pick source to download MC3 data from:
## TCGAmutations is faster, but not verified (pre-compiled by PoisonAlien on github)
sample_set = "MC3/GDC" #"MC3/PoisonAlien"
#type_file = 'ref_files/cancer_types.txt'
type_file = 'foo'
data_dir = os.path.join("data", sample_set)
mutsigcv_dir =  os.path.join('mutsigcv_results', sample_set)

def get_jobnames(type_file):
	return pd.read_csv(type_file, sep = '\t')['abbrev'].to_list()

## rules 

rule all:	
	input: expand(os.path.join(mutsigcv_dir, "{job}.sig_genes.txt"), job = get_jobnames(type_file))


rule save_mafs_TCGAmutations:
	output: "data/MC3/PoisonAlien/{job}.mutSig.maf"
	params: "{job}", data_dir
	shell: 
		"Rscript download_mc3.R {params}"

rule download_mc3_gdc:
	output: "data/MC3/GDC/mc3.v0.2.8.PUBLIC.maf.gz"
	shell:
		"wget -O {output} https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc"

rule save_mafs_gdc:
	input: "data/MC3/GDC/mc3.v0.2.8.PUBLIC.maf.gz"
	output: expand("data/MC3/GDC/{job}.mutSig.maf", job = get_jobnames(type_file))
	params: data_dir
	shell:
		"Rscript download_mc3.R {input} {params}"

rule run_mutsigcv:
	input: os.path.join(data_dir, "{job}.mutSig.maf")
	output: multiext(os.path.join(mutsigcv_dir, "{job}"), ".sig_genes.txt")
	params: os.path.join(mutsigcv_dir, "{job}")
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: 
		"""
		taskset -c 27,28,29,30 MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} \
		ref_files/exome_full192.coverage.txt ref_files/gene.covariates.txt {params} \
		ref_files/mutation_type_dictionary_file.txt ref_files/chr_files_hg19 > {log[0]} 2> {log[1]}
		"""


