source('helpers.R')

mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
		   cancer = cancer_names$brief[match(variable, cancer_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	arrange(desc(value), desc(Gene)) %>%
	select(cancer) %>% unique() %>% unlist() -> lvls

# vanilla
png(file = '../plots_tmp/mutsig_mtor.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
	ggplot(aes(x=value, y=cancer, fill = Gene)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 2.5) +
		geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 3.5, label = 'p = 0.001', size = 2.5) +
		labs(title= 'Mutation frequency significance', y = '', x = expression(-log[10](p-value)) ) +
		scale_fill_viridis(discrete = T, direction = -1) +
		#facet_wrap( ~ cancer, ncol = 1) + 
		theme_classic()  -> p1
p1
dev.off()

# with pies
png(file = '../plots_tmp/mutsig_mtor_pies.png', width = 800)
mutRates[c('MTOR', "PTEN", "PIK3CA")] %>% 
	tibble::rownames_to_column(var = "abbrev") %>% 
	melt() %>% rename(Gene = variable, frac = value) %>%
	mutate(cancer = ordered(cancer_names$brief[match(abbrev, cancer_names$abbrev)], levels = lvls), antifrac = 100-frac) %>%
	melt() %>%tidyr::replace_na(list(value = 0)) %>%
	mutate(variable = factor(ifelse(variable == "antifrac", "antifrac", Gene))) %>%
	ggplot(aes(x = 0.5, y=value, fill=variable)) + 
		geom_bar(stat = "identity", width = 1, position = position_fill()) + 
		coord_polar("y", start=0) + 
		scale_fill_manual(values=c(viridis(3)[3:1], "#d3d3d3")) + 
		facet_wrap(~cancer+Gene, ncol = 3) +
		labs(title= 'Fraction of samples mutated') +
		theme_minimal() + labs(x = "", y = "") + 
		theme(strip.background = element_blank(), strip.text = element_blank(), 
			  legend.position = "none", panel.grid  = element_blank(),
			  panel.spacing = unit(-2, "lines"), 
        	  axis.ticks = element_blank(), axis.text = element_blank()) +
        scale_x_continuous(limits = c(0, 3)) -> p2
egg::ggarrange(p1, p2, ncol =2, widths = c(6,1))
dev.off()

# other mutsig plots

png(file = '../plots_tmp/mutsig_drivers.png')
df %>% 
	subset(Gene %in% c("TP53", "KRAS", "BRCA1", "BRCA2", "TTN")) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, x=value, y=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 2.5) +
		geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 3.5, label = 'p = 0.001', size = 2.5) +
		labs(title= 'Mutation frequency significance', y = '', x = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		theme_classic()
dev.off()

png(file = '../plots_tmp/mutsig_rtks.png')
df %>% 
	subset(Gene %in% c("ERBB2", "EGFR")) %>%
	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', x = 1, y = 1.75, label = 'p = 0.05', size = 3) +
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
		annotate('text', x = 1, y = 1.2, label = 'p = 0.05', size = 3) +
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
		annotate('text', x = 1, y = 1.2, label = 'p = 0.05', size = 3) +
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


# TMB vs mtor alt frequency

png(file = '../plots_tmp/tmb_mtor_scatter.png', width = 600)
cancer_names %>% 
	mutate(med_TMB = unlist(lapply(mafs, med_tmb))[cancer_names$abbrev], 
		   MTOR_alt_frac = mutRates['MTOR'][match(cancer_names$abbrev, rownames(mutRates)),]) %>%
	tidyr::replace_na(list(med_TMB = 0, MTOR_alt_frac = 0))%>%
	ggplot(aes(x = MTOR_alt_frac, y = med_TMB, col = brief)) + 
		geom_point() +  
		geom_text(aes(hjust=1.02,vjust=0.01,label=ifelse(cancer_names$abbrev %in% c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-SKCM', 'TCGA-STAD'),cancer_names$brief,''))) +
		labs(title= "MTOR alterations not linear with absolute TMB", x = 'MTOR fraction of samples altered (%)', y = 'Median TMB/MB') +
		theme_classic() + theme(legend.position = "none")
dev.off()


# Subtyped renal

renal_mutsig_p %>%
	melt() %>%
	mutate(Gene = rep(factor(goi, ordered = T), ncol(renal_mutsig_p)),
		   cancer = renal_names$name[match(variable, renal_names$abbrev)],
		   value = pmin(-log10(value), 5)) -> df

df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	arrange(desc(value), desc(Gene)) %>%
	select(cancer) %>% unique() %>% unlist() -> lvls

# vanilla
png(file = '../plots_tmp/mutsig_renal_subtypes.png')
df %>% 
	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
	mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
	ggplot(aes(x=value, y=cancer, fill = Gene)) + 
		geom_bar(position="dodge", stat="identity") +
		geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 2.5) +
		geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
		annotate('text', y = 1, x = 3.5, label = 'p = 0.001', size = 2.5) +
		labs(title= 'Mutation frequency significance', y = '', x = expression(-log[10](p-value)) ) +
		scale_fill_viridis(discrete = T, direction = -1) +
		#facet_wrap( ~ cancer, ncol = 1) + 
		theme_classic()  -> p1
p1
dev.off()

# with pies
png(file = '../plots_tmp/mutsig_renal_subtypes_pies.png', width = 800)
renal_mutRates[c('MTOR', "PTEN", "PIK3CA")] %>% 
	tibble::rownames_to_column(var = "abbrev") %>% 
	melt() %>% rename(Gene = variable, frac = value) %>%
	mutate(cancer = ordered(renal_names$brief[match(abbrev, renal_names$abbrev)], levels = lvls), antifrac = 100-frac) %>%
	melt() %>%tidyr::replace_na(list(value = 0)) %>%
	mutate(variable = factor(ifelse(variable == "antifrac", "antifrac", Gene))) %>%
	ggplot(aes(x = 0.5, y=value, fill=variable)) + 
		geom_bar(stat = "identity", width = 1, position = position_fill()) + 
		coord_polar("y", start=0) + 
		scale_fill_manual(values=c(viridis(3)[3:1], "#d3d3d3")) + 
		facet_wrap(~cancer+Gene, ncol = 3) +
		labs(title= 'Fraction of samples mutated') +
		theme_minimal() + labs(x = "", y = "") + 
		theme(strip.background = element_blank(), strip.text = element_blank(), 
			  legend.position = "none", panel.grid  = element_blank(),
			  panel.spacing = unit(-2, "lines"), 
        	  axis.ticks = element_blank(), axis.text = element_blank()) +
        scale_x_continuous(limits = c(0, 3)) -> p2
egg::ggarrange(p1, p2, ncol =2)
dev.off()

# Waterfall plots

for (m in names(renal_mafs)){
	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
	lil_m <- subsetMaf(renal_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
	title(main = renal_names$brief[match(m, renal_names$abbrev)])
	dev.off()
}


# Venn (option)
for (m in names(renal_mafs)){
	venn.diagram(
	  x = list(tsb(subsetMaf(renal_mafs[[m]], genes = c("PIK3CA"))), 
	  		   tsb(subsetMaf(renal_mafs[[m]], genes = c("PTEN"))), 
	  		   tsb(subsetMaf(renal_mafs[[m]], genes = c("MTOR"))) ),
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

# do-over for m4 (no PIK3CA)

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
	mutate(Gene = rep(factor(goi, ordered = T), ncol(endo_mutsig_p)),
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
		geom_hline(yintercept = -log10(0.001), size=0.2, linetype="dashed") +
		annotate('text', x = 0.5, y = 3.5, label = 'p = 0.001', size = 3) +
		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
		scale_fill_viridis(discrete = T, direction = -1) +
		coord_flip() + theme_classic()
dev.off()


# Waterfall plots

for (m in names(endo_mafs)){
	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
	lil_m <- subsetMaf(endo_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
	title(main = endo_names$brief[match(m, endo_names$abbrev)])
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


