library(maftools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(viridis)

# taken from rmd

# genes of interest
goi <- c("TTN", "TP53", "KRAS", "BRCA1", "BRCA2", "ATM", "WNT1", "AKT1", "AKT2", "AKT3", "TSC1", "TSC2", "EGFR", 
         "FGFR1", "ERBB2", "ERBB3", "ERBB4", "ROS1", "MET", "ALK", "FLT1", "PDGFRA", "FLT3", "FLT4", "RET", "FGFR2", 
         "FGFR3", "DEPDC5", "NPRL2", "NPRL3", "MIOS", "SEH1L", "SEC13", "WDR24", "WDR59", "SLC38A9", "MTOR", "PTEN", "PIK3CA", "PIK3CB", "PIK3CG", "PIK3CD")
goi <- factor(goi, levels = goi, ordered = T)


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

mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = cancer_names$brief[match(variable, cancer_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

png(file = '../plots_tmp/1c_mutsig.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.75, label = 'p = 0.05', size = 2) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()

png(file = '../plots_tmp/s1_mutsig_drivers.png')
df %>% 
	subset(Gene %in% c("TP53", "KRAS", "BRCA1", "BRCA2", "TTN")) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.75, label = 'p = 0.05', size = 2) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()

png(file = '../plots_tmp/s2_mutsig_rtks.png')
df %>% 
	subset(Gene %in% c("ERBB2", "EGFR")) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.75, label = 'p = 0.05', size = 2) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()

png(file = '../plots_tmp/s2_mutsig_gator1.png')
df %>% 
	subset(Gene %in% c("DEPDC5", "NPRL2", "NPRL3")) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.1, label = 'p = 0.05', size = 2) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()

png(file = '../plots_tmp/s2_mutsig_gator2.png')
df %>% 
	subset(Gene %in% c("MIOS", "SEH1L", "SEC13", "WDR24", "WDR59"))%>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.1, label = 'p = 0.05', size = 2) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()

# Waterfall plots

for (m in c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-MESO', 'TCGA-STAD')){
	png(file = paste0('../plots_tmp/2c_waterfall_', m, '.png'))
	PlotOncogenicPathways(mafs[[m]], pathways = "PI3K")
	title(main = cancer_names$brief[match(m, cancer_names$abbrev)])
	dev.off()
}
		

# Subtyped renal

params = list(maf_dir = '../GDCdata/MC3-subtyped',
			  res_dir = '../results/MC3-subtyped',
			  names_tr = "../ref_files/cancer_subtypes.txt"
			  )

maf_dir <- params$maf_dir
maf_fns <- list.files(maf_dir, pattern = "^TCGA-KIRC*")
renal_mafs <- sapply(paste0(maf_dir, "/", maf_fns), read.maf)
names(renal_mafs) <- sapply(maf_fns, function(x){sub('\\.mutSig\\.maf$', '', x) })

# cancer types
cancer_names <- read.table(params$names_tr, header = T, sep = "\t")

# mutsigCV results
res_dir <- params$res_dir
res_fns <- list.files(res_dir, pattern = "*.sig_genes.txt")
mutsig_p <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(renal_mafs)) {
  ord <- match(c, substr(res_fns, 1, nchar(res_fns)-14))
  mutsig_res <- read.table(paste0(res_dir, "/", res_fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
}

sampleSizes <- sapply(renal_mafs, function(x){as.numeric(x@summary[ID %in% 'Samples', summary])})
mutRates <- sapply(renal_mafs, function(x){getGeneSummary(x)[match(goi, getGeneSummary(x)$Hugo_Symbol), MutatedSamples]})
mutRates <- as.data.frame(apply(mutRates, denom = sampleSizes, FUN = function(x, denom){ x / denom * 100}, MAR = 1))
colnames(mutRates) <- goi

mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = cancer_names$name[match(variable, cancer_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

png(file = '../plots_tmp/s3_mutsig_renal_subtypes.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 1.75, label = 'p = 0.05', size = 3) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()



