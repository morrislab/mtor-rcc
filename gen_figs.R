###### gen_figs.R ######
# create figures from TCGA, PCAWG data

cancer_names <- read.table('cancer_names.txt', header = T, sep = "\t")

## MC3 maf data
# MC3 data can be retrieved with R package TCGAmutations
abbrevs <- gsub("TCGA-", "", (cancer_names$abbrev))
mafs <- TCGAmutations::tcga_load(abbrevs)
names(mafs) <- paste0("TCGA-", names(mafs))

## uncomment following if desired for mutsig analysis 
#data_dir <- "~/data/mc3-mafs"
#if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}
#for (m in names(mafs)){
#    maftools::prepareMutSig(mafs[[m]], file.path(data_dir, paste0('TCGA-',m)))
#}

## PCAWG driver data
# retrieved from icgc 

data_dir <- "~/data/pcawg-drivers"
if(!dir.exists(data_dir)){dir.create(data_dir, recursive = TRUE)}

data_source <- "https://dcc.icgc.org/api/v1/download?fn=/PCAWG/networks/final_integration_results_2017_03_16.tar.gz"
data_dest <- file.path(data_dir, basename(data_source))

# retrieve data 
if(!file.exists(data_dest)){
	options(timeout=1e6)
	download.file(data_source, dest = data_dest)
}

# select CDS files only
fns <- untar(data_dest, list=TRUE)
sel <- grepl(fns, pattern = "CDS.combined_p_values.automatic_method_removal.txt")
#sel <- sel & grepl(fns, pattern = "Kidney-RCC|Panc-Endocrine|Uterus-AdenoCa|ColoRect-AdenoCa|Stomach-AdenoCa")
for(f in fns[sel]){
	if (!file.exists(file.path(data_dir,f))){untar(data_dest, files = f, exdir = path.expand(data_dir)) }
}

## genes of interest
goi <- c("MTOR", "PTEN", "PIK3CA", "TTN", "TP53", "KRAS", "BRCA1", "BRCA2", 
         "DEPDC5", "NPRL2", "NPRL3", "MIOS", "SEH1L", "SEC13", "WDR24", "WDR59")
goi <- factor(goi, levels = goi, ordered = T)

## bar order
lvls <- c("Bladder Carcinoma", "Breast Carcinoma", "Cervical SCC", "Colon Adenocarcinoma", "Head & Neck SCC", "Brain Glioma (LG)", "Rectal Adenocarcinoma", "Stomach Adenocarcinoma", "Endometrial Carcinoma",  "Uterine Carcinosarcoma",  "Glioblastoma Multiforme",  "Kidney Chromophobe",  "Kidney RCCC",  "Lung SCC",  "Prostate Adenocarcinoma",  "Cutaneous Melanoma",  "Renal Papillary Cell",  "Sarcoma",  "Esophageal Carcinoma",  "Ovarian Carcinoma",  "Hepatocellular Carcinoma",  "Lung Adenocarcinoma",  "Thyroid Carcinoma",  "Mesothelioma",  "Pancreatic Adenocarcinoma",  "Germ Cell Tumors",  "DLBCL",  "Cholangiocarcinoma",  "Adrenocortical Carcinoma",  "Uveal Melanoma",  "Thymoma",  "Acute Myeloid Leukemia",  "Pheochromocytoma & Para")



###### prep data #####

## TCGA MC3
data_dir <- '~/data/mutsig-mc3'
fns <- list.files(data_dir, pattern = "*.sig_genes.txt")
mutsig_p <- data.frame(row.names = goi)

# set up mutsig dataframe for genes of interest
for (c in names(mafs)) {
  ord <- match(c, substr(fns, 1, nchar(fns)-14))
  mutsig_res <- read.table(paste0(data_dir, "/", fns[ord]), header=T)
  poi <- mutsig_res[mutsig_res$gene %in% goi, c(1,14,15)]
  mutsig_p[[c]] <- poi$p[match(goi, poi$gene)]
}

## PCAWG 
data_dir <- '~/data/pcawg-drivers/xchip/cga_home/gtiao/PCAWG/Oct_2016/final_integration_results_2017_03_16'
fns <- list.files(data_dir)

pcawg_names <- gsub('.CDS.combined_p_values.automatic_method_removal.txt', '', fns)
pcawg_names <- pcawg_names[!grepl('meta', pcawg_names) & !grepl('Pancan', pcawg_names)]

pcawg_p <- data.frame(row.names = goi)
# set up mutsig dataframe for genes of interest
for (c in pcawg_names) {
  res <- read.table(file.path(data_dir, paste0(c, '.CDS.combined_p_values.automatic_method_removal.txt')), header = T)
  res <- res[!is.na(res$ID),]
  res$gene <- unlist(stringr::str_split(res$ID, "::"))[c(F, F, T, F)]
  poi <- res[res$gene %in% goi,]
  pcawg_p[[c]] <- poi$MutSig[match(goi, poi$gene)]
}

## Venn diagram data is downloaded from cBioPortal oncoprint tab in tabular format
## for Kidney Renal Clear Cell Carcinoma (TCGA, PanCancer Atlas) 
## and Uterine Corpus Endometrial Carcinoma (TCGA, PanCancer Atlas)

source('figs.R')

###### plot ######

pdf('plots.pdf')

# fig 1D
p_plots(mutsig_p, cancer_names, gene_names = c("MTOR", "PTEN", "PIK3CA"), main = 'TCGA MC3 MutSigCV p-values: MTOR and Upstream', pal_name = 'hl', lvls = lvls)

# fig 4C
p_plots(mutsig_p, cancer_names, gene_names = c("DEPDC5", "NPRL2", "NPRL3"), main = 'TCGA MC3 MutSigCV p-values: GATOR1 complex', pal_name = 'flat', lvls = lvls)

# fig 4C
p_plots(mutsig_p, cancer_names, gene_names = c("MIOS", "SEH1L", "SEC13", "WDR24", "WDR59"), main = 'TCGA MC3 MutSigCV p-values: GATOR2 complex', pal_name = 'flat', lvls = lvls)

# fig S1D
p_plots(pcawg_p, data_frame(abbrev = pcawg_names, brief = pcawg_names), main = 'PCAWG-Driver Working Group MutSig p-values', gene_names = c("MTOR", "PTEN", "PIK3CA"), pal_name = 'hl')

# fig S1E
p_plots(mutsig_p, cancer_names, gene_names = c("TTN", "TP53", "KRAS", "BRCA1", "BRCA2"), main = 'TCGA MC3 MutSigCV p-values', pal_name = 'flat', lvls = lvls)


# Venn for all mutations
overrideTriple = T


grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

read_tsv('kirc_PATIENT_DATA_oncoprint.tsv', skip_empty_rows = F) %>%
    subset(track_type != 'STRUCTURAL_VARIANT') %>%
    mutate(track = str_c(track_name, '_', track_type)) %>%
    select(-c(track_name, track_type)) %>%
    column_to_rownames('track') %>%
    mutate(across(everything(), ~!is.na(.x))) %>%
    mut_venn(main = 'RCCC patients with any alteration', cat.dist = c(0.06, 0.06, 0.03), cat.pos = c(0, 0, 0))

grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

read_tsv('ucec_PATIENT_DATA_oncoprint.tsv', skip_empty_rows = F) %>%
    subset(track_type != 'STRUCTURAL_VARIANT') %>%
    mutate(track = str_c(track_name, '_', track_type)) %>%
    select(-c(track_name, track_type)) %>%
    column_to_rownames('track') %>%
    mutate(across(everything(), ~!is.na(.x))) %>%
    mut_venn(main = 'Endometrial carcinoma patients with any alteration', cex = c(2,2,2,0.9,1.7,2,2), cat.pos = c(-27, 27, 180))

grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

rm(overrideTriple)
read_tsv('ucec_PATIENT_DATA_oncoprint.tsv', skip_empty_rows = F) %>%
    subset(track_type != 'STRUCTURAL_VARIANT') %>%
    mutate(track = str_c(track_name, '_', track_type)) %>%
    select(-c(track_name, track_type)) %>%
    column_to_rownames('track') %>%
    mutate(across(everything(), ~!is.na(.x))) %>%
    mut_venn(main = 'Endometrial carcinoma patients with any alteration', cat.pos = c(-27, 27, 180))


# Venn for putative drivers only
grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

read_tsv('kirc_PATIENT_DATA_oncoprint.tsv', skip_empty_rows = F) %>%
    subset(track_type != 'STRUCTURAL_VARIANT') %>%
    mutate(track = str_c(track_name, '_', track_type)) %>%
    select(-c(track_name, track_type)) %>%
    column_to_rownames('track') %>%
    mutate(across(everything(), ~is.driver(.x))) %>%
    mut_venn(main = 'RCCC patients with driver alterations', cat.dist = c(0.06, 0.06, 0.06), cat.pos = c(-27, 70, 180))

grid.newpage()
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc")))

read_tsv('ucec_PATIENT_DATA_oncoprint.tsv', skip_empty_rows = F) %>%
    subset(track_type != 'STRUCTURAL_VARIANT') %>%
    mutate(track = str_c(track_name, '_', track_type)) %>%
    select(-c(track_name, track_type)) %>%
    column_to_rownames('track') %>%
    mutate(across(everything(), ~is.driver(.x))) %>%
    mut_venn(main = 'Endometrial carcinoma patients with driver alterations', 
             cat.dist = c(0.08, 0.06, 0.06), cat.pos = c(-27, 27, 0))

dev.off()

