library(maftools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(VennDiagram)
library(viridis)

# taken from rmd

# genes of interest
goi <- c("TTN", "TP53", "KRAS", "BRCA1", "BRCA2", "ATM", "WNT1", "AKT1", "AKT2", "AKT3", "TSC1", "TSC2", "EGFR", 
         "FGFR1", "ERBB2", "ERBB3", "ERBB4", "ROS1", "MET", "ALK", "FLT1", "PDGFRA", "FLT3", "FLT4", "RET", "FGFR2", 
         "FGFR3", "DEPDC5", "NPRL2", "NPRL3", "MIOS", "SEH1L", "SEC13", "WDR24", "WDR59", "SLC38A9", "MTOR", "PTEN", "PIK3CA", "PIK3CB", "PIK3CG", "PIK3CD")
goi <- factor(goi, levels = goi, ordered = T)

# TMB calculations
med_tmb <- function(m){ median(tmb(m)$total_perMB) }

# get Tumor_Sample_Barcode
tsb <- function(m){ m@data$Tumor_Sample_Barcode }


# MC3
params = list(maf_dir = '../GDCdata/MC3',
			  res_dir = '../results/MC3',
			  names_tr = "../ref_files/cancer_types.txt"
			  )

# mafs
maf_dir <- params$maf_dir
maf_fns <- list.files(maf_dir, pattern = "*.maf$")
mafs <- sapply(paste0(maf_dir, "/", maf_fns), read.maf)
names(mafs) <- sapply(maf_fns, function(x){sub('\\.mutSig\\.maf$', '', x) })

# cancer types
cancer_names <- read.table(params$names_tr, header = T, sep = "\t")

# mutsigCV results
res_dir <- params$res_dir
res_fns <- list.files(res_dir, pattern = "*.sig_genes.txt")
mutsig_p <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(mafs)) {
  ord <- match(c, substr(res_fns, 1, nchar(res_fns)-14))
  mutsig_res <- read.table(paste0(res_dir, "/", res_fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
}

sampleSizes <- sapply(mafs, function(x){as.numeric(x@summary[ID %in% 'Samples', summary])})
mutRates <- sapply(mafs, function(x){getGeneSummary(x)[match(goi, getGeneSummary(x)$Hugo_Symbol), MutatedSamples]})
mutRates <- as.data.frame(apply(mutRates, denom = sampleSizes, FUN = function(x, denom){ x / denom * 100}, MAR = 1))
colnames(mutRates) <- goi


# renal subtyping

renal_params = list(maf_dir = '../GDCdata/MC3-subtyped',
			  		res_dir = '../results/MC3-subtyped',
			  		names_tr = "../ref_files/cancer_subtypes.txt"
			  		)

maf_dir <- renal_params$maf_dir
maf_fns <- list.files(maf_dir, pattern = "^TCGA-KIRC*")
renal_mafs <- sapply(paste0(maf_dir, "/", maf_fns), read.maf)
names(renal_mafs) <- sapply(maf_fns, function(x){sub('\\.mutSig\\.maf$', '', x) })

# cancer types
renal_names <- read.table(renal_params$names_tr, header = T, sep = "\t")

# mutsigCV results
res_dir <- renal_params$res_dir
res_fns <- list.files(res_dir, pattern = "*.sig_genes.txt")
mutsig_p <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(renal_mafs)) {
  ord <- match(c, substr(res_fns, 1, nchar(res_fns)-14))
  mutsig_res <- read.table(paste0(res_dir, "/", res_fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
}

renal_sampleSizes <- sapply(renal_mafs, function(x){as.numeric(x@summary[ID %in% 'Samples', summary])})
renal_mutRates <- sapply(renal_mafs, function(x){getGeneSummary(x)[match(goi, getGeneSummary(x)$Hugo_Symbol), MutatedSamples]})
renal_mutRates <- as.data.frame(apply(renal_mutRates, denom = renal_sampleSizes, FUN = function(x, denom){ x / denom * 100}, MAR = 1))
colnames(renal_mutRates) <- goi

#endometrial subtyping

endo_params = list(maf_dir = '../GDCdata/MC3-subtyped',
			  		res_dir = '../results/MC3-subtyped',
			  		names_tr = "../ref_files/cancer_subtypes.txt"
			  		)

maf_dir <- endo_params$maf_dir
maf_fns <- list.files(maf_dir, pattern = "^TCGA-UCEC*")
endo_mafs <- sapply(paste0(maf_dir, "/", maf_fns), read.maf)
names(endo_mafs) <- sapply(maf_fns, function(x){sub('\\.mutSig\\.maf$', '', x) })

# cancer types
endo_names <- read.table(endo_params$names_tr, header = T, sep = "\t")

# mutsigCV results
res_dir <- endo_params$res_dir
res_fns <- list.files(res_dir, pattern = "*.sig_genes.txt")
mutsig_p <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(endo_mafs)) {
  ord <- match(c, substr(res_fns, 1, nchar(res_fns)-14))
  mutsig_res <- read.table(paste0(res_dir, "/", res_fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
}

endo_sampleSizes <- sapply(endo_mafs, function(x){as.numeric(x@summary[ID %in% 'Samples', summary])})
endo_mutRates <- sapply(endo_mafs, function(x){getGeneSummary(x)[match(goi, getGeneSummary(x)$Hugo_Symbol), MutatedSamples]})
endo_mutRates <- as.data.frame(apply(mutRates, denom = sampleSizes, FUN = function(x, denom){ x / denom * 100}, MAR = 1))
colnames(endo_mutRates) <- goi

