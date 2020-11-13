source('figs-preprocess.R')

mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = cancer_names$brief[match(variable, cancer_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

png(file = '../plots_tmp/mutsig_mtor.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, x=value, y=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 2) +
		geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 3.4, label = 'p = 0.001', size = 2) +
		labs(title= 'Mutation frequency significance', y = '', x = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		theme_classic()
dev.off()

png(file = '../plots_tmp/mutsig_drivers.png')
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

png(file = '../plots_tmp/mutsig_rtks.png')
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

png(file = '../plots_tmp/mutsig_gator1.png')
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

png(file = '../plots_tmp/mutsig_gator2.png')
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
	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
	lil_m <- subsetMaf(mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
	title(main = cancer_names$brief[match(m, cancer_names$abbrev)])
	dev.off()
}
		

# Venn (option)
for (m in c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-MESO', 'TCGA-STAD')){
	venn.diagram(
	  x = list(tsb(subsetMaf(mafs[[m]], genes = c("PIK3CA"))), tsb(subsetMaf(mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(mafs[[m]], genes = c("MTOR"))) ),
	  category.names = c("PIK3CA", "PTEN", "MTOR"),
	  filename = paste0('../plots_tmp/venn_', m, '.png'),
	  output = TRUE ,
	          imagetype="png" ,
	          height = 480 , 
	          width = 480 , 
	          resolution = 300,
	          compression = "lzw",
	          lwd = 1,
	          col=viridis(3),
	          fill = alpha( viridis(3), 0.3 ),
	          cex = 0.5,
	          fontfamily = "sans",
	          cat.cex = 0.3,
	          cat.default.pos = "outer",
	          cat.pos = c(-27, 27, 135),
	          cat.dist = c(0.04, 0.04, 0.04),
	          cat.fontfamily = "sans",
	          #cat.col = viridis(3),
	          rotation = 1
	        )
}


# Subtyped renal

renal_mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = renal_names$name[match(variable, renal_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

png(file = '../plots_tmp/mutsig_renal_subtypes.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 1.75, label = 'p = 0.05', size = 3) +
		geom_hline(yintercept = -log10(0.01), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 2, label = 'p = 0.01', size = 3) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()


# Waterfall plots

for (m in names(renal_mafs)){
	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
	lil_m <- subsetMaf(renal_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
	title(main = cancer_names$brief[match(m, cancer_names$abbrev)])
	dev.off()
}


# Venn (option)
for (m in names(renal_mafs)){
	venn.diagram(
	  x = list(tsb(subsetMaf(renal_mafs[[m]], genes = c("PIK3CA"))), tsb(subsetMaf(renal_mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(renal_mafs[[m]], genes = c("MTOR"))) ),
	  category.names = c("PIK3CA", "PTEN", "MTOR"),
	  filename = paste0('../plots_tmp/venn_', m, '.png'),
	  output = TRUE ,
	          imagetype="png" ,
	          height = 480 , 
	          width = 480 , 
	          resolution = 300,
	          compression = "lzw",
	          lwd = 1,
	          col=viridis(3),
	          fill = alpha( viridis(3), 0.3 ),
	          cex = 0.5,
	          fontfamily = "sans",
	          cat.cex = 0.3,
	          cat.default.pos = "outer",
	          cat.pos = c(-27, 27, 135),
	          cat.dist = c(0.04, 0.04, 0.04),
	          cat.fontfamily = "sans",
	          #cat.col = viridis(3),
	          rotation = 1
	        )
}

# do-over for m4 (no PI3KCA)

m = names(renal_mafs)[[4]]
venn.diagram(
  x = list(tsb(subsetMaf(renal_mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(renal_mafs[[m]], genes = c("MTOR"))) ),
  category.names = c("PTEN", "MTOR"),
  filename = paste0('../plots_tmp/venn_', m, '.png'),
  output = TRUE ,
          imagetype="png" ,
          height = 480 , 
          width = 480 , 
          resolution = 300,
          compression = "lzw",
          lwd = 1,
          col=viridis(3)[2:3],
          fill = alpha( viridis(3), 0.3 )[2:3],
          cex = 0.5,
          fontfamily = "sans",
          cat.cex = 0.3,
          cat.default.pos = "outer",
          #cat.pos = c(-27, 27, 135),
          cat.dist = c(0.04, 0.04),
          cat.fontfamily = "sans",
          #cat.col = viridis(3),
          #rotation = 1
        )

# Subtyped endometrial

endo_mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = endo_names$name[match(variable, endo_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

png(file = '../plots_tmp/mutsig_endo_subtypes.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 1.75, label = 'p = 0.05', size = 3) +
		geom_hline(yintercept = -log10(0.01), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 2, label = 'p = 0.01', size = 3) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()


# Waterfall plots

for (m in names(endo_mafs)){
	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
	lil_m <- subsetMaf(endo_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
	title(main = cancer_names$brief[match(m, cancer_names$abbrev)])
	dev.off()
}
		
# Venn (option)
for (m in names(endo_mafs)){
	venn.diagram(
	  x = list(tsb(subsetMaf(endo_mafs[[m]], genes = c("PIK3CA"))), tsb(subsetMaf(endo_mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(endo_mafs[[m]], genes = c("MTOR"))) ),
	  category.names = c("PIK3CA", "PTEN", "MTOR"),
	  filename = paste0('../plots_tmp/venn_', m, '.png'),
	  output = TRUE ,
	          imagetype="png" ,
	          height = 480 , 
	          width = 480 , 
	          resolution = 300,
	          compression = "lzw",
	          lwd = 1,
	          col=viridis(3),
	          fill = alpha( viridis(3), 0.3 ),
	          cex = 0.5,
	          fontfamily = "sans",
	          cat.cex = 0.3,
	          cat.default.pos = "outer",
	          cat.pos = c(-27, 27, 135),
	          cat.dist = c(0.04, 0.04, 0.04),
	          cat.fontfamily = "sans",
	          #cat.col = viridis(3),
	          rotation = 1
	        )
}


