if(!file.exists('pcawgp.Rdata')){

data_dir <- '~/data/pcawg-drivers/xchip/cga_home/gtiao/PCAWG/Oct_2016/final_integration_results_2017_03_16'
fns <- list.files(data_dir)

pcawg_names <- gsub('.CDS.combined_p_values.automatic_method_removal.txt', '', fns)
pcawg_names <- pcawg_names[!grepl('meta', pcawg_names) & !grepl('Pancan', pcawg_names)]
pcawg_p <- pcawg_b <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in pcawg_names) {
  #ord <- match(c, gsub('.CDS.combined_p_values.automatic_method_removal.txt', '', fns))
  res <- read.table(file.path(data_dir, paste0(c, '.CDS.combined_p_values.automatic_method_removal.txt')), header = T)
  res <- res[!is.na(res$ID),]
  res$gene <- unlist(stringr::str_split(res$ID, "::"))[c(F, F, T, F)]
  poi <- res[res$gene %in% goi,]
  pcawg_p[[c]] <- poi$MutSig[match(goi, poi$gene)]
  pcawg_b[[c]] <- poi$Brown_observed_trimmed[match(goi, poi$gene)]
}

save(file = 'pcawgp.Rdata', pcawg_p, pcawg_b)

} else{load('pcawgp.Rdata')}

