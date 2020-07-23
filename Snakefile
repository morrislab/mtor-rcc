import glob

with open('joblist.txt') as f:
	jobs = f.read().splitlines()

maf_dict = {x: glob.glob("data/{}/*.maf".format(x))[0] for x in jobs}

rule all:
	input: "summary/index.html"

rule knit_summary:
	input: expand("summary/{job}.html", job = jobs)
	output: "summary/index.html"
	log: multiext("logs/all", ".stdout", ".stderr")
	shell: 
		"""
		Rscript -e "rmarkdown::render('summary.Rmd', output_file = '{output}', output_dir = 'summary')" > {log[0]} 2> {log[1]}
		"""

rule knit_plots:
	input: "results/{job}.sig_genes.txt"
	output: "summary/{job}.html"
	params: "{job}"
	log: multiext("logs/{job}", ".stdout", ".stderr")
	shell: 
		"""
		Rscript -e "rmarkdown::render('plot_maftools.Rmd', params = list(id = '{params}'), output_file = 'summary/{params}.html', output_dir = 'summary')" > {log[0]} 2> {log[1]}
		"""

rule run_mutsigcv:
	input: lambda wildcards: maf_dict[wildcards.job]
	output: multiext("results/{job}", ".sig_genes.txt", ".mutcateg_discovery.txt", ".categs.txt", ".coverage.txt", ".mutations.txt")
	params: "results/{job}"
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: "MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} exome_full192.coverage.txt gene.covariates.txt {params} mutation_type_dictionary_file.txt chr_files_GDC > {log[0]} 2> {log[1]}"

	

