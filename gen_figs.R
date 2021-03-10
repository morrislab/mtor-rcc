# gen_figs.R
# prep data for all figures derived from PCAWG drivers working group analysis

## MC3 maf data
# retrieved with R package TCGAmutations

data_dir <- "~/data/mc3-mafs"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}
abbrevs <- gsub("TCGA-", "", (read.table('cancer_types.txt', header = T, sep = "\t")$abbrev))

mafs <- TCGAmutations::tcga_load(abbrevs)
for (m in names(mafs)){
    maftools::prepareMutSig(mafs[[m]], file.path(data_dir, paste0('TCGA-',m)))
}

## PCAWG driver data
# retrieved from icgc api

data_dir <- "~/data/pcawg-drivers"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}

data_source <- "https://dcc.icgc.org/api/v1/download?fn=/PCAWG/networks/final_integration_results_2017_03_16.tar.gz"
data_dest <- file.path(data_dir, basename(data_source))

# retrieve data 
if(!file.exists(data_dest)){
	options(timeout=1e6)
	download.file(data_source, dest = data_dest)
}

# select CDS files only
fns <- untar(data_dest, list=TRUE)
sel <- grepl(fns, pattern = "CDS.combined_p_values.automatic_method_removal.txt")
#sel <- sel & grepl(fns, pattern = "Kidney-RCC|Panc-Endocrine|Uterus-AdenoCa|ColoRect-AdenoCa|Stomach-AdenoCa")
untar(data_dest, files = fns[sel], exdir = path.expand(data_dir))

# plot for each study result
if(!file.exists('mutsigp.Rdata')){source('prep_mutsigp.R')}
if(!file.exists('pcawgp.Rdata')){source('prep_pcawgp.R')}

source('figs.R')

pdf('plots/cait_mutsig_p.pdf')
p_plots(mutsig_p, cancer_names, main = 'MutSigCV p-values', pq = p)
dev.off()

pdf('plots/cait_mutsig_q.pdf')
p_plots(mutsig_q, cancer_names, main = 'MutSigCV q-values', pq = q)
dev.off()

pdf('plots/pcawg_mutsig_p.pdf')
p_plots(pcawg_p, data_frame(abbrev = pcawg_names, brief = pcawg_names), 
	main = 'PCAWG-Driver Working Group MutSig p-values', pq = p)
dev.off()

pdf('plots/pcawg_trimmed-brown_p.pdf')
p_plots(pcawg_b, data_frame(abbrev = pcawg_names, brief = pcawg_names),
	main = 'PCAWG-Driver Working Group Brown Test q-values', pq = q)
dev.off()

pdf('plots/tmb_scatter.pdf')
tmb_plot(cancer_names)
dev.off()

