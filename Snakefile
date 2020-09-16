import os
from glob import glob

sample_set = "harmonized" #"MC3"

data_dir = os.path.join("GDCdata", sample_set)
results_dir = os.path.join("results", sample_set)
summary_dir = os.path.join("summary", sample_set)

def get_jobnames(data_dir):
	maf_location = os.path.join(data_dir, "*.mutSig.maf")
	return [os.path.basename(x)[:-11] for x in glob(maf_location)]

rule all:	
	input: os.path.join(summary_dir, "index.html")

rule knit_dash:
	#input: expand("results/{set}/{job}.sig_genes.txt", set = sample_sets)
	input: expand(os.path.join(results_dir, "{job}.sig_genes.txt"), job = get_jobnames(data_dir))
	output: os.path.join(summary_dir, "index.html")
	params: summary_dir, os.path.join("..", data_dir), os.path.join("..", results_dir)
	shell: 
		"""
		Rscript -e "rmarkdown::render('scripts/dash.Rmd', output_file = '{output}', output_dir = '{params[0]}', params = list(maf_dir = '{params[1]}', res_dir = '{params[2]}'))"
		"""

rule run_mutsigcv:
	input: os.path.join(data_dir, "{job}.mutSig.maf")
	output: multiext(os.path.join(results_dir, "{job}"), ".sig_genes.txt", ".mutcateg_discovery.txt", ".categs.txt", ".coverage.txt", ".mutations.txt")
	params: os.path.join(results_dir, "{job}")
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: 
		"""
		taskset -c 27,28,29,30 MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} \
		ref_files/exome_full192.coverage.txt ref_files/gene.covariates.txt {params} \
		ref_files/mutation_type_dictionary_file.txt ref_files/chr_files_hg19 > {log[0]} 2> {log[1]}
		"""
