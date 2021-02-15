# gen_figs.R
# prep data for all figures derived from PCAWG drivers working group analysis

## MC3 maf data
# retrieved with R package TCGAmutations

data_dir <- "~/data/mc3-mafs"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}

mafs <- TCGAmutations::tcga_load(c('KIRC', 'UCEC', 'STAD', 'COAD'))
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

# select target files
fns <- untar(data_dest, list=TRUE)
sel <- grepl(fns, pattern = "CDS.combined_p_values.automatic_method_removal.txt")
sel <- sel & grepl(fns, pattern = "Kidney-RCC|Panc-Endocrine|Uterus-AdenoCa|ColoRect-AdenoCa|Stomach-AdenoCa")
fns <- fns[sel]

untar(data_dest, files = fns, exdir = path.expand(data_dir))

# plot for each study result
for(fn in file.path(data_dir, fns)){
	res <- read.table(fn, header = T)	
}

