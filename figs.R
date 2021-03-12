library(dplyr)
library(maftools)
library(VennDiagram)
library(reshape2)
library(ggplot2)

# https://drsimonj.svbtle.com/creating-corporate-colour-palettes-for-ggplot2
ch_palettes <- list(
  `flat` = c("#5380a2", "#15203a", "#327138", "#d05c39", "#d79e35"),
  `hl` = c("#d11141", "#858585", "#15304A")
)

ch_pal <- function(palette = "hl", reverse = FALSE, ...) {
  pal <- ch_palettes[[palette]]
  if (reverse) pal <- rev(pal)
  colorRampPalette(pal, ...)
}

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

	#png(file = '../plots_tmp/other_nutrient_sensing.png')
        df %>% 
		subset(Gene %in% c("GATSL3", "GATSL2","C7orf60", "SESN2"))%>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(fill=Gene, y=value, x=cancer)) + viz +
		discrete_scale("fill", 'hl', palette = ch_pal('flat'), labels = c('CASTOR1', 'CASTOR2', 'SAMTOR', 'SESN2')) -> p6
	#dev.off()               

	print(p1)
	print(p2)
	print(p3)
	print(p4)
	print(p6)

}


# Venn 
venn_pik3ca_pten_mtor <- function(m, main = 'Mutation Co-occurance (n samples in cohort)'){
	venn.diagram(
	  x = list(tsb(subsetMaf(m, genes = c("PIK3CA"))), tsb(subsetMaf(m, genes = c("PTEN"))), tsb(subsetMaf(m, genes = c("MTOR"))) ),
	  category.names = c("MTOR", "PTEN", "PIK3CA"),
	  filename = NULL,
	  output = TRUE ,
	  col=ch_pal('hl')(3),
	  fill = alpha( ch_pal('hl')(3), 0.3 ),
	  cex = 1,
	  fontfamily = "sans",
	  cat.cex = 1,
	  cat.default.pos = "outer",
	  cat.pos = c(-27, 27, 135),
	  cat.dist = c(0.04, 0.04, 0.04),
	  cat.fontfamily = "sans",
	  rotation = 1
	) %>% grid.draw()
}


# TMB vs mtor alt frequency
tmb_plot <- function(mutRates, cancer_names, show_na = T, main = "MTOR alterations not linear with absolute TMB", do_labels = c('TCGA-KIRC', 'TCGA-COAD', 'TCGA-UCEC', 'TCGA-SKCM')){
	cancer_names %>% 
		mutate(med_TMB = unlist(lapply(mafs, med_tmb))[cancer_names$abbrev], 
			   MTOR_alt_frac = mutRates['MTOR'][match(cancer_names$abbrev, rownames(mutRates)),]) -> df
		if(show_na){df %>% tidyr::replace_na(list(med_TMB = 0, MTOR_alt_frac = 0)) -> df}
		
		df %>%
		ggplot(aes(x = MTOR_alt_frac, y = med_TMB)) + 
			geom_point() +  
			geom_text(aes(hjust=1.03,vjust=0.02,label=ifelse(cancer_names$abbrev %in% do_labels,cancer_names$brief,''))) +
			labs(title= main, x = 'MTOR fraction of samples altered (%)', y = 'Median TMB/MB') +
			geom_smooth(method = "lm", formula = y~x, se=FALSE, colour="black", weight = 0.7) +
			ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
			theme_classic() + theme(legend.position = "none")
}



prep_mutex_mat <- function(maf, outfn){

  subset(maf@data, select = c('Hugo_Symbol', 'Tumor_Sample_Barcode')) %>% 
  	mutate(Tumor_Sample_Barcode = as.character(Tumor_Sample_Barcode)) %>%
  	dcast(Hugo_Symbol ~ Tumor_Sample_Barcode)  %>%
  	tibble() %>% 
  	tibble::column_to_rownames('Hugo_Symbol') %>%
  	replace(((.)>1),1)

}



