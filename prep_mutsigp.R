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
med_tmb <- function(maf, captureSize = 50, logScale = TRUE){
  # modified from maftools

  if(!is(object = maf, class2 = "MAF")){
    stop("Input must be an MAF object")
  }

  maf.mutload = getSampleSummary(maf)[,.(Tumor_Sample_Barcode, total)]
  maf.mutload[,total_perMB := total/captureSize]
  maf.mutload[,total_perMB_log := log10(total_perMB)]
  maf.mutload = maf.mutload[order(total_perMB, decreasing = FALSE)]

  medload = median(maf.mutload[,total_perMB], na.rm = TRUE)
  return(medload)
}

# get Tumor_Sample_Barcode
tsb <- function(m){ m@data$Tumor_Sample_Barcode }


# MC3
params = list(maf_dir = '~/data/mc3-mafs',
              res_dir = '~/data/mutsig-mc3',
              names_tr = "cancer_types.txt"
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
mutsig_p <- mutsig_q <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(mafs)) {
  ord <- match(c, substr(res_fns, 1, nchar(res_fns)-14))
  mutsig_res <- read.table(paste0(res_dir, "/", res_fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14,15)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
  mutsig_q[[c]] <- poi$q[match(goi, poi$gene)]
}

sampleSizes <- sapply(mafs, function(x){as.numeric(x@summary[ID %in% 'Samples', summary])})
mutRates <- sapply(mafs, function(x){getGeneSummary(x)[match(goi, getGeneSummary(x)$Hugo_Symbol), MutatedSamples]})
mutRates <- as.data.frame(apply(mutRates, denom = sampleSizes, FUN = function(x, denom){ x / denom * 100}, MAR = 1))
colnames(mutRates) <- goi

save(file = 'mutsigp.Rdata', mutsig_p, mutsig_q)