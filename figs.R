

# https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2
ch_palettes <- list(
  `flat` = cividis(30)[1:20],
  `hl` = c("#d11141", cividis(50)[5:20])
)

ch_pal <- function(palette = "hl", reverse = FALSE, ...) {
  pal <- ch_palettes[[palette]]
  if (reverse) pal <- rev(pal)
  colorRampPalette(pal, ...)
}


p_plots <- function(mutsig_p, cancer_names, main = 'Mutation frequency significance', pq = p){ 

	if(pq == 'p'){xlab = expression(-log[10](p-value))
	}else{xlab = expression(-log[10](q-value))}

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
	#png(file = '../plots_tmp/mutsig_mtor.png')
	df %>% 
		subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(x=value, y=cancer, fill = Gene)) + 
			geom_bar(position="dodge", stat="identity") +
			geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
			annotate('text', y = 1, x = 1.75, label = paste0(pq, ' = 0.05'), size = 3.5) +
			geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
			annotate('text', y = 1, x = 3.5, label = paste0(pq, ' = 0.001'), size = 3.5) +
			labs(title= main, y = '', x = xlab ) +
			#facet_wrap( ~ cancer, ncol = 1) + 
			#scale_fill_viridis(discrete = T, direction = -1, option = "D", end = 0.7) +
			discrete_scale("fill", 'hl', palette = ch_pal()) +
			theme_classic() + 
			theme(axis.text.y = element_text(colour = "black", size = 10)) -> p1

	#dev.off()

	# with pies
	#png(file = '../plots_tmp/mutsig_mtor_pies.png', width = 800)
	#mutRates[c('MTOR', "PTEN", "PIK3CA")] %>% 
	#	tibble::rownames_to_column(var = "abbrev") %>% 
	#	melt() %>% rename(Gene = variable, frac = value) %>%
	#	mutate(cancer = ordered(cancer_names$brief[match(abbrev, cancer_names$abbrev)], levels = lvls), antifrac = 100-frac) %>%
	#	melt() %>%tidyr::replace_na(list(value = 0)) %>%
	#	mutate(variable = factor(ifelse(variable == "antifrac", "antifrac", Gene))) %>%
	#	ggplot(aes(x = 0.5, y=value, fill=variable)) + 
	#		geom_bar(stat = "identity", width = 1, position = position_fill()) + 
	#		coord_polar("y", start=0) + 
	#		scale_fill_manual(values=c(viridis(3)[3:1], "#d3d3d3")) + 
	#		facet_wrap(~cancer+Gene, ncol = 3) +
	#		labs(title= 'Fraction of samples mutated') +
	#		theme_minimal() + labs(x = "", y = "") + 
	#		theme(strip.background = element_blank(), strip.text = element_blank(), 
	#			  legend.position = "none", panel.grid  = element_blank(),
	#			  panel.spacing = unit(-2, "lines"), 
	#        	  axis.ticks = element_blank(), axis.text = element_blank()) +
	#        scale_x_continuous(limits = c(0, 3)) -> p2
	#egg::ggarrange(p1, p2, ncol =2, widths = c(6,1))
	#dev.off()
	
	# other mutsig plots
	
	#png(file = '../plots_tmp/mutsig_drivers.png')
	df %>% 
		subset(Gene %in% c("TP53", "KRAS", "BRCA1", "BRCA2", "TTN")) %>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(fill=Gene, x=value, y=cancer)) + 
			geom_bar(position="dodge", stat="identity") +
			geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
			annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 3.5) +
			geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
			annotate('text', y = 1, x = 3.5, label = 'p = 0.001', size = 3.5) +
			labs(title= main, y = '', x = xlab ) + 
			#scale_fill_viridis(discrete = T, direction = -1, option = "E", end = 0.8) +
			discrete_scale("fill", 'hl', palette = ch_pal('flat')) +
			theme_classic()+ 
			theme(axis.text.y = element_text(colour = "black", size = 10)) -> p2

	#dev.off()

	viz <- list(geom_bar(position="dodge", stat="identity"),
				geom_hline(yintercept = -log10(0.05), size=0.3, linetype="dashed"),
				annotate('text', x = 1, y = 1.75, label = 'p = 0.05', size = 3.5),
				geom_hline(yintercept = -log10(0.001), size=0.3, linetype="dashed"),
				annotate('text', x = 1, y = 3.5, label = 'p = 0.001', size = 3.5),
				labs(title= main, x = '', y = xlab ),
				discrete_scale("fill", 'hl', palette = ch_pal('flat')),
				coord_flip(), theme_classic(), ylim(0, 5),
				theme(axis.text.y = element_text(colour = "black", size = 10))
			   )

	
	#png(file = '../plots_tmp/mutsig_rtks.png')
	df %>% 
		subset(Gene %in% c("ERBB2", "EGFR")) %>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(fill=Gene, y=value, x=cancer)) + viz -> p3
			
	#dev.off()
	
	
	#png(file = '../plots_tmp/mutsig_gator1.png')
	df %>% 
		subset(Gene %in% c("DEPDC5", "NPRL2", "NPRL3")) %>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(fill=Gene, y=value, x=cancer)) + viz  -> p4
	#dev.off()
	
	#png(file = '../plots_tmp/mutsig_gator2.png')
	df %>% 
		subset(Gene %in% c("MIOS", "SEH1L", "SEC13", "WDR24", "WDR59"))%>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(fill=Gene, y=value, x=cancer)) + viz  -> p5
	#dev.off()


	print(p1)
	print(p2)
	print(p3)
	print(p4)
	print(p5)

}


# # Waterfall plots
# 
# for (m in c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-MESO', 'TCGA-STAD')){
# 	#png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
# 	lil_m <- subsetMaf(mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
# 	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
# 	title(main = cancer_names$brief[match(m, cancer_names$abbrev)])
# 	#dev.off()
# }
# 		
# 
# # Venn (option)
# for (m in c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-STAD')){
# 	venn.diagram(
# 	  x = list(tsb(subsetMaf(mafs[[m]], genes = c("PIK3CA"))), tsb(subsetMaf(mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(mafs[[m]], genes = c("MTOR"))) ),
# 	  category.names = c("PIK3CA", "PTEN", "MTOR"),
# 	  filename = NULL,
# 	  output = TRUE ,
# 	  col=viridis(3),
# 	  fill = alpha( viridis(3), 0.3 ),
# 	  cex = 1,
# 	  fontfamily = "sans",
# 	  cat.cex = 0.3,
# 	  cat.default.pos = "outer",
# 	  cat.pos = c(-27, 27, 135),
# 	  cat.dist = c(0.04, 0.04, 0.04),
# 	  cat.fontfamily = "sans",
# 	  #cat.col = viridis(3),
# 	  rotation = 1
# 	  ) %>% grid.draw()
# }
# 

# TMB vs mtor alt frequency
tmb_plot <- function(cancer_names){
	cancer_names %>% 
		mutate(med_TMB = unlist(lapply(mafs, med_tmb))[cancer_names$abbrev], 
			   MTOR_alt_frac = mutRates['MTOR'][match(cancer_names$abbrev, rownames(mutRates)),]) %>%
		tidyr::replace_na(list(med_TMB = 0, MTOR_alt_frac = 0)) %>%
		ggplot(aes(x = MTOR_alt_frac, y = med_TMB)) + 
			geom_point() +  
			geom_text(aes(hjust=1.03,vjust=0.02,label=ifelse(cancer_names$abbrev %in% c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-SKCM'),cancer_names$brief,''))) +
			labs(title= "MTOR alterations not linear with absolute TMB", x = 'MTOR fraction of samples altered (%)', y = 'Median TMB/MB') +
			geom_smooth(method = "lm", se=FALSE, color="black", formula = y~x) +
			ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
			theme_classic() + theme(legend.position = "none")
}



# # Subtyped renal
# 
# renal_mutsig_p %>%
# 	melt() %>%
# 	mutate(Gene = rep(factor(goi, ordered = T), ncol(renal_mutsig_p)),
# 		   cancer = renal_names$name[match(variable, renal_names$abbrev)],
# 		   value = pmin(-log10(value), 5)) -> df
# 
# df %>% 
# 	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
# 	arrange(desc(value), desc(Gene)) %>%
# 	select(cancer) %>% unique() %>% unlist() -> lvls
# 
# # vanilla
# png(file = '../plots_tmp/mutsig_renal_subtypes.png')
# df %>% 
# 	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
# 	mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
# 	ggplot(aes(x=value, y=cancer, fill = Gene)) + 
# 		geom_bar(position="dodge", stat="identity") +
# 		geom_vline(xintercept = -log10(0.05), size=0.2, linetype="dashed") +
# 		annotate('text', y = 1, x = 1.7, label = 'p = 0.05', size = 2.5) +
# 		geom_vline(xintercept = -log10(0.001), size=0.2, linetype="dashed") +
# 		annotate('text', y = 1, x = 3.5, label = 'p = 0.001', size = 2.5) +
# 		labs(title= 'Mutation frequency significance', y = '', x = expression(-log[10](p-value)) ) +
# 		scale_fill_viridis(discrete = T, direction = -1) +
# 		#facet_wrap( ~ cancer, ncol = 1) + 
# 		theme_classic()  -> p1
# p1
# dev.off()
# 
# # with pies
# png(file = '../plots_tmp/mutsig_renal_subtypes_pies.png', width = 800)
# renal_mutRates[c('MTOR', "PTEN", "PIK3CA")] %>% 
# 	tibble::rownames_to_column(var = "abbrev") %>% 
# 	melt() %>% rename(Gene = variable, frac = value) %>%
# 	mutate(cancer = ordered(renal_names$brief[match(abbrev, renal_names$abbrev)], levels = lvls), antifrac = 100-frac) %>%
# 	melt() %>%tidyr::replace_na(list(value = 0)) %>%
# 	mutate(variable = factor(ifelse(variable == "antifrac", "antifrac", Gene))) %>%
# 	ggplot(aes(x = 0.5, y=value, fill=variable)) + 
# 		geom_bar(stat = "identity", width = 1, position = position_fill()) + 
# 		coord_polar("y", start=0) + 
# 		scale_fill_manual(values=c(viridis(3)[3:1], "#d3d3d3")) + 
# 		facet_wrap(~cancer+Gene, ncol = 3) +
# 		labs(title= 'Fraction of samples mutated') +
# 		theme_minimal() + labs(x = "", y = "") + 
# 		theme(strip.background = element_blank(), strip.text = element_blank(), 
# 			  legend.position = "none", panel.grid  = element_blank(),
# 			  panel.spacing = unit(-2, "lines"), 
#         	  axis.ticks = element_blank(), axis.text = element_blank()) +
#         scale_x_continuous(limits = c(0, 3)) -> p2
# egg::ggarrange(p1, p2, ncol =2)
# dev.off()
# 
# # Waterfall plots
# 
# for (m in names(renal_mafs)){
# 	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
# 	lil_m <- subsetMaf(renal_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
# 	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
# 	title(main = renal_names$brief[match(m, renal_names$abbrev)])
# 	dev.off()
# }
# 
# 
# # Venn (option)
# for (m in names(renal_mafs)){
# 	venn.diagram(
# 	  x = list(tsb(subsetMaf(renal_mafs[[m]], genes = c("PIK3CA"))), 
# 	  		   tsb(subsetMaf(renal_mafs[[m]], genes = c("PTEN"))), 
# 	  		   tsb(subsetMaf(renal_mafs[[m]], genes = c("MTOR"))) ),
# 	  category.names = c("PIK3CA", "PTEN", "MTOR"),
# 	  filename = paste0('../plots_tmp/venn_', m, '.png'),
# 	  output = TRUE ,
# 	          imagetype="png" ,
# 	          height = 480 , 
# 	          width = 480 , 
# 	          resolution = 300,
# 	          compression = "lzw",
# 	          lwd = 1,
# 	          col=viridis(3),
# 	          fill = alpha( viridis(3), 0.3 ),
# 	          cex = 0.5,
# 	          fontfamily = "sans",
# 	          cat.cex = 0.3,
# 	          cat.default.pos = "outer",
# 	          cat.pos = c(-27, 27, 135),
# 	          cat.dist = c(0.04, 0.04, 0.04),
# 	          cat.fontfamily = "sans",
# 	          #cat.col = viridis(3),
# 	          rotation = 1
# 	        )
# }
# 
# # do-over for m4 (no PIK3CA)
# 
# m = names(renal_mafs)[[4]]
# venn.diagram(
#   x = list(tsb(subsetMaf(renal_mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(renal_mafs[[m]], genes = c("MTOR"))) ),
#   category.names = c("PTEN", "MTOR"),
#   filename = paste0('../plots_tmp/venn_', m, '.png'),
#   output = TRUE ,
#           imagetype="png" ,
#           height = 480 , 
#           width = 480 , 
#           resolution = 300,
#           compression = "lzw",
#           lwd = 1,
#           col=viridis(3)[2:3],
#           fill = alpha( viridis(3), 0.3 )[2:3],
#           cex = 0.5,
#           fontfamily = "sans",
#           cat.cex = 0.3,
#           cat.default.pos = "outer",
#           #cat.pos = c(-27, 27, 135),
#           cat.dist = c(0.04, 0.04),
#           cat.fontfamily = "sans",
#           #cat.col = viridis(3),
#           #rotation = 1
#         )
# 
# # Subtyped endometrial
# 
# endo_mutsig_p %>%
# 	melt() %>%
# 	mutate(Gene = rep(factor(goi, ordered = T), ncol(endo_mutsig_p)),
# 		   cancer = endo_names$name[match(variable, endo_names$abbrev)],
# 		   value = pmin(-log10(value), 5)) -> df
# 
# png(file = '../plots_tmp/mutsig_endo_subtypes.png')
# df %>% 
# 	subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
# 	mutate(cancer = ordered(cancer, levels = rev(unique(cancer[order(desc(value), desc(Gene))])))) %>%
# 	ggplot(aes(fill=Gene, y=value, x=cancer)) + 
# 		geom_bar(position="dodge", stat="identity") +
# 		geom_hline(yintercept = -log10(0.05), size=0.2, linetype="dashed") +
# 		annotate('text', x = 0.5, y = 1.75, label = 'p = 0.05', size = 3) +
# 		geom_hline(yintercept = -log10(0.001), size=0.2, linetype="dashed") +
# 		annotate('text', x = 0.5, y = 3.5, label = 'p = 0.001', size = 3) +
# 		labs(title= 'Mutation frequency significance', x = '', y = expression(-log[10](p-value)) ) + 
# 		scale_fill_viridis(discrete = T, direction = -1) +
# 		coord_flip() + theme_classic()
# dev.off()
# 
# 
# # Waterfall plots
# 
# for (m in names(endo_mafs)){
# 	png(file = paste0('../plots_tmp/waterfall_', m, '.png'))
# 	lil_m <- subsetMaf(endo_mafs[[m]], genes = c("MTOR", "PTEN", "PIK3CA"))
# 	PlotOncogenicPathways(lil_m, pathways = "PI3K", fontSize = 1)
# 	title(main = endo_names$brief[match(m, endo_names$abbrev)])
# 	dev.off()
# }
# 		
# # Venn (option)
# for (m in names(endo_mafs)){
# 	venn.diagram(
# 	  x = list(tsb(subsetMaf(endo_mafs[[m]], genes = c("PIK3CA"))), tsb(subsetMaf(endo_mafs[[m]], genes = c("PTEN"))), tsb(subsetMaf(endo_mafs[[m]], genes = c("MTOR"))) ),
# 	  category.names = c("PIK3CA", "PTEN", "MTOR"),
# 	  filename = paste0('../plots_tmp/venn_', m, '.png'),
# 	  output = TRUE ,
# 	          imagetype="png" ,
# 	          height = 480 , 
# 	          width = 480 , 
# 	          resolution = 300,
# 	          compression = "lzw",
# 	          lwd = 1,
# 	          col=viridis(3),
# 	          fill = alpha( viridis(3), 0.3 ),
# 	          cex = 0.5,
# 	          fontfamily = "sans",
# 	          cat.cex = 0.3,
# 	          cat.default.pos = "outer",
# 	          cat.pos = c(-27, 27, 135),
# 	          cat.dist = c(0.04, 0.04, 0.04),
# 	          cat.fontfamily = "sans",
# 	          #cat.col = viridis(3),
# 	          rotation = 1
# 	        )
# }


