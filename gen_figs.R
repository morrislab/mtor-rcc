# gen_figs.R
# prep data for all figures derived from PCAWG drivers working group analysis

## MC3 maf data
# retrieved with R package TCGAmutations

data_dir <- "~/data/mc3-mafs"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}
abbrevs <- gsub("TCGA-", "", (read.table('cancer_types.txt', header = T, sep = "\t")$abbrev))

mafs <- TCGAmutations::tcga_load(abbrevs)
names(mafs) <- paste0("TCGA-", names(mafs))

##for (m in names(mafs)){
##    maftools::prepareMutSig(mafs[[m]], file.path(data_dir, paste0('TCGA-',m)))
##}

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
for(f in fns[sel]){
	if (!file.exists(file.path(data_dir,f))){untar(data_dest, files = f, exdir = path.expand(data_dir)) }
}

# cancer types
cancer_names <- read.table('cancer_types.txt', header = T, sep = "\t")

# genes of interest
goi <- c("TTN", "TP53", "KRAS", "BRCA1", "BRCA2", "ATM", "WNT1", "AKT1", "AKT2", "AKT3", "TSC1", "TSC2", "EGFR", "C7orf60", "GATSL3", 
         "FGFR1", "ERBB2", "ERBB3", "ERBB4", "ROS1", "MET", "ALK", "FLT1", "PDGFRA", "FLT3", "FLT4", "RET", "FGFR2", 
         "FGFR3", "DEPDC5", "NPRL2", "NPRL3", "MIOS", "SEH1L", "SEC13", "WDR24", "WDR59", "SLC38A9", "MTOR", "PTEN", "PIK3CA", "PIK3CB", "PIK3CG", "PIK3CD")
goi <- factor(goi, levels = goi, ordered = T)


# plot for each study result
source('prep_mutsigp.R')
source('prep_pcawgp.R')

source('figs.R')

pdf('plots/cait_mutsig_p.pdf')
p_plots(mutsig_p, cancer_names, main = 'MutSigCV p-values', pq = 'p')
dev.off()

pdf('plots/cait_mutsig_q.pdf')
p_plots(mutsig_q, cancer_names, main = 'MutSigCV q-values', pq = 'q')
dev.off()

pdf('plots/pcawg_mutsig_p.pdf')
p_plots(pcawg_p, data_frame(abbrev = pcawg_names, brief = pcawg_names), 
	main = 'PCAWG-Driver Working Group MutSig p-values', pq = 'p')
dev.off()

pdf('plots/pcawg_trimmed-brown_p.pdf')
p_plots(pcawg_b, data_frame(abbrev = pcawg_names, brief = pcawg_names),
	main = 'PCAWG-Driver Working Group Brown Test q-values', pq = 'q')
dev.off()




pdf('plots/tmb_scatter.pdf')
tmb_plot(mutRates, cancer_names, main = "MTOR Alterations not Linear with Absolute TMB")
tmb_plot(mutRates, cancer_names, main = "MTOR Alterations not Linear with Absolute TMB", do_labels = c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC','TCGA-SKCM', 'TCGA-LUSC', 'TCGA-LUAD', 'TCGA-BLCA', 'TCGA-STAD'))
tmb_plot(mutRates[rownames(mutRates) != 'TCGA-SKCM',], cancer_names[cancer_names$abbrev != 'TCGA-SKCM',], main = "MTOR Alterations not Linear with Absolute TMB (Melanoma excluded)")
tmb_plot(mutRates[rownames(mutRates) != 'TCGA-SKCM',], cancer_names[cancer_names$abbrev != 'TCGA-SKCM',], main = "MTOR Alterations not Linear with Absolute TMB (Melanoma excluded)", , do_labels = c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC','TCGA-LUSC', 'TCGA-LUAD', 'TCGA-BLCA'))
dev.off()


pdf('plots/co-mut_venns.pdf')
venn_pik3ca_pten_mtor(mafs[['TCGA-KIRC']], main = 'Kidney RCC Mutation Co-occurance (n samples in cohort)')
grid.newpage()
venn_pik3ca_pten_mtor(mafs[['TCGA-UCEC']], main = 'Kidney RCC Mutation Co-occurance (n samples in cohort)')
#somaticInteractions(maf = mafs[['TCGA-KIRC']], top = 25, pvalue = c(0.05, 0.1))
#somaticInteractions(maf = mafs[['TCGA-UCEC']], top = 25, pvalue = c(0.05, 0.1))
#somaticInteractions(maf = mafs[['TCGA-KIRC']], genes = c('PTEN', 'PIK3CA', 'MTOR'), pvalue = c(0.05, 0.1))
#somaticInteractions(maf = mafs[['TCGA-UCEC']], genes = c('PTEN', 'PIK3CA', 'MTOR'), pvalue = c(0.05, 0.1))
dev.off()


#prep_mutex_mat(mafs[['TCGA-KIRC']]) %>% write.csv("kirc_comut.csv")
#prep_mutex_mat(mafs[['TCGA-UCEC']]) %>% write.csv("ucec_comut.csv")

