## BiocManager::install("PoisonAlien/TCGAmutations")
library(TCGAmutations)
library(dplyr)

## helpers to download data

save_poisonalien_mc3 <- function(cancer_type, data_dir){
    # download from compiled data made available by PoisonAlien on github

    if(!file.exists(data_dir)){
        dir.create(data_dir, recursive = T)
    }

    cancer_type %>% 
        gsub(pattern='TCGA-', replacement='') %>%
        tcga_load() %>%
        prepareMutSig(fn = file.path(data_dir, cancer_type))

}

save_gdc_mc3 <- function(maf_path, data_dir){
    # process the mc3.v0.2.8.PUBLIC.maf.gz MC3 file from GDC
    if(!file.exists(data_dir)){
        dir.create(data_dir, recursive = T)
    }

    message("o Reading MAF")
    mc3 <- readr::read_tsv(maf_path,progress = TRUE, col_types = readr::cols())
    
    message("o Adding project_id information")
    clinical <- TCGAbiolinks::colDataPrepare(unique(mc3$Tumor_Sample_Barcode))[c('project_id')] %>%
        tibble::rownames_to_column(var = "Tumor_Sample_Barcode") %>%
        filter(complete.cases(.))

    # read into maftools maf class
    maf <- read.maf(mc3, clinicalData = clinical, removeDuplicatedVariants = F)

    # generate per-cancer mafs
    for (project in unique(clinical$project_id)) {
            message(paste0('o Preparing ', project, '.mutSig.maf'))
            sel <- clinical$Tumor_Sample_Barcode[clinical$project_id == project]
            subsetMaf(maf, tsb=sel) %>%
                prepareMutSig(fn = file.path(data_dir, project))
    }

}


# download MAF and prepare for mutsigCV analysis for cmd arguments
args <- commandArgs(trailingOnly=TRUE)
data_dir <- args[2]

if (grepl("PoisonAlien", data_dir)){
    print('retrieving MAF from TCGAmutations package')
    save_poisonalien_mc3(args[1], args[2])

} else if (grepl("GDC", data_dir)){
    print('retrieving MAFs from GDC MC3 file')
    save_gdc_mc3(args[1], args[2])

} else {stop('download selection not understood')}
