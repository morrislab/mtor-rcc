library(dplyr)
library(maftools)
library(VennDiagram)
library(reshape2)
library(ggplot2)
library(readr)
library(tibble)
library(stringr)

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

# get Tumor_Sample_Barcode
tsb <- function(m){ m@data$Tumor_Sample_Barcode }


p_plots <- function(mutsig_p, cancer_names, gene_names, main = 'MutSigCV sigificance', pal_name = 'hl', pq = 'p', lvls = NULL){ 

	if(pq == 'p'){xlab = expression(-log[10](p-value))
	}else{xlab = expression(-log[10](q-value))}

	mutsig_p %>%
		melt() %>%
		mutate(Gene = rep(factor(goi, ordered = T), ncol(mutsig_p)),
			   cancer = cancer_names$brief[match(variable, cancer_names$abbrev)],
			   value = pmin(-log10(value), 5)) -> df
    
    if (is.null(lvls)){
        df %>% 
		subset(Gene %in% c('PIK3CA', 'PTEN', 'MTOR')) %>%
		arrange(desc(value), desc(Gene)) %>%
		select(cancer) %>% unique() %>% unlist() -> lvls
    }
	
	df %>% 
		subset(Gene %in% gene_names) %>%
		mutate(cancer = ordered(cancer, levels = rev(lvls))) %>%
		ggplot(aes(x=value, y=cancer, fill = Gene)) + 
			geom_bar(position="dodge", stat="identity") +
			geom_vline(xintercept = -log10(0.05/length(lvls)), size=0.2, linetype="dashed") +
			annotate('text', y = 1, x = -log10(0.05/length(lvls)) + 0.65, 
                     label = paste0(pq, ' = 0.05 / ', length(lvls) ), size = 3.5) +
            #geom_vline(xintercept = -log10(0.1/length(lvls)), size=0.2, linetype="dashed") +
			#annotate('text', y = 1, x = -log10(0.1/length(lvls)) + 0.5, 
            #         label = paste0(pq, ' = 0.1 / ', length(lvls) ), size = 3.5) +
			labs(title= main, y = '', x = xlab ) +
			discrete_scale("fill", '', palette = ch_pal(pal_name)) +
			theme_classic() + 
            xlim(0, 5) +
			theme(axis.text.y = element_text(colour = "black", size = 10)) %>%
    return()

}

# bonferroni correct pvalues
bonferroni <- function(p_mat, axis=1){
    p_mat <- p_mat * dim(p_mat)[axis]
    p_mat[p_mat >1] <- 1
    return(p_mat)
}



is.driver <- function(str){
    return ( !is.na(str) & (grepl('driver', str) | grepl('_rec', str)) )
}

mut_venn <- function(df, cex = 2, cat.dist = c(0.06, 0.06, 0.06), cat.pos = c(-27, 27, 135),
                     main = ''){
  venn.diagram(
      x = list(colnames(df)[as.logical(colSums(df[grepl('PIK3CA', rownames(df)),])) ],
               colnames(df)[as.logical(colSums(df[grepl('PTEN', rownames(df)),])) ],
               colnames(df)[as.logical(colSums(df[grepl('MTOR', rownames(df)),])) ]),
      category.names = c("PIK3CA", "PTEN", "MTOR"),
      filename = NULL,
      output = TRUE ,
      col= rev(ch_pal('hl')(3)),
      fill = alpha( rev(ch_pal('hl')(3)), 0.3 ),
      cex = cex,
      ext.text = T,
      fontfamily = "sans",
      cat.cex = 2,
      cat.default.pos = "outer",
      cat.pos = cat.pos,
      cat.dist = cat.dist,
      cat.fontfamily = "sans",
      rotation = 1,
      main = paste0(main, " (n=", sum(colSums(df)>0), ")"),
      main.cex = 1,
      main.fontfamily = 'sans'
  ) %>% grid.draw()
}


