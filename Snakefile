import glob
import os

maf_pattern = "GDCdata/harmonized/*.mutSig.maf"
jobs = [os.path.basename(x) for x in glob.glob(maf_pattern)]
jobs = [x[:-11] for x in jobs]


rule all:
	#input: expand("summary/harmonized/{job}.html", job = jobs)	
	input: "summary/harmonized/index.html"

rule knit_dash:
	input: expand("results/harmonized/{job}.sig_genes.txt", job = jobs)
	output: "summary/harmonized/index.html"
	params: "summary/harmonized", "../GDCdata/harmonized", "../results/harmonized"
	shell: 
		"""
		Rscript -e "rmarkdown::render('scripts/dash.Rmd', output_file = '{output}', output_dir = '{params[0]}', params = list(maf_dir = '{params[1]}', res_dir = '{params[2]}'))"
		"""

#rule knit_plots:
#	input: "results_harmonized/{job}.sig_genes.txt"
#	output: "summary/{job}.html"
#	params: "{job}"
#	log: multiext("logs/{job}", ".stdout", ".stderr")
#	shell: 
#		"""
#		Rscript -e "rmarkdown::render('scripts/plot_maftools.Rmd', params = list(id = '{params}'), output_file = 'summary/{params}.html', output_dir = 'summary')" > {log[0]} 2> {log[1]}
#		"""

rule run_mutsigcv:
	input: "GDCdata/harmonized/{job}.mutSig.maf"
	output: multiext("results/harmonized/{job}", ".sig_genes.txt", ".mutcateg_discovery.txt", ".categs.txt", ".coverage.txt", ".mutations.txt")
	params: "results/harmonized/{job}"
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: 
		"""
		taskset -c 27,28,29,30 MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} \
		ref_files/exome_full192.coverage.txt ref_files/gene.covariates.txt {params} \
		ref_files/mutation_type_dictionary_file.txt ref_files/chr_files_GRC38 > {log[0]} 2> {log[1]}
		"""
