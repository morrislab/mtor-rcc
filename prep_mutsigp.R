library(maftools)
library(dplyr)
library(reshape2)
if(!file.exists('mutsigp.Rdata')){

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


save(file = 'mutsigp.Rdata', mutsig_p, mutsig_q, mutRates)

}else{load('mutsigp.Rdata')}
