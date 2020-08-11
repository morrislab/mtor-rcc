import glob
import os

maf_pattern = "GDCdata/MC3/*.mutSig.maf"
jobs = [os.path.basename(x) for x in glob.glob(maf_pattern)]
jobs = [x[:-11] for x in jobs]


rule all:
	#input: expand("summary/MC3/{job}.html", job = jobs)	
	input: "summary/MC3/index.html"

rule knit_dash:
	input: expand("results/MC3/{job}.sig_genes.txt", job = jobs)
	output: "summary/MC3/index.html"
	params: "summary/MC3", "../GDCdata/MC3", "../results/MC3"
	shell: 
		"""
		Rscript -e "rmarkdown::render('scripts/dash.Rmd', output_file = '{output}', output_dir = '{params[0]}', params = list(maf_dir = '{params[1]}', res_dir = '{params[2]}'))"
		"""

rule run_mutsigcv:
	input: "GDCdata/MC3/{job}.mutSig.maf"
	output: multiext("results/MC3/{job}", ".sig_genes.txt", ".mutcateg_discovery.txt", ".categs.txt", ".coverage.txt", ".mutations.txt")
	params: "results/MC3/{job}"
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: 
		"""
		taskset -c 27,28,29,30 MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} \
		ref_files/exome_full192.coverage.txt ref_files/gene.covariates.txt {params} \
		ref_files/mutation_type_dictionary_file.txt ref_files/chr_files_hg19 > {log[0]} 2> {log[1]}
		"""
