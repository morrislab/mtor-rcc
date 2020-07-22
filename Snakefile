import glob

with open('joblist.txt') as f:
	jobs = f.read().splitlines()

maf_dict = {x: glob.glob("data/{}/*.maf".format(x))[0] for x in jobs}

rule all:
	input: expand("results/{job}.sig_genes.txt", job = jobs)

rule run_mutsig:
	input: lambda wildcards: maf_dict[wildcards.job]
	output: multiext("results/{job}", ".sig_genes.txt", ".mutcateg_discovery.txt", ".categs.txt", ".coverage.txt", ".mutations.txt")
	params: "results/{job}"
	log: multiext("logs/{job}", ".stdout", ".stderr")
	benchmark: "logs/{job}.time"
	shell: "MutSigCV_1.41/run_MutSigCV.sh /home/cait/Desktop/mtor/v901 {input} exome_full192.coverage.txt gene.covariates.txt {params} mutation_type_dictionary_file.txt chr_files_GDC > {log[0]} 2> {log[1]}"


