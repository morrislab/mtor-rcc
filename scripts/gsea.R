library(TCGAbiolinks)
library(DESeq2)
library(maftools)

library(BiocParallel)
register(MulticoreParam(3))

library(clusterProfiler)
libarary(enrichplot)
organism <- "org.Hs.eg.db"
library(organism, character.only = TRUE)

proj_name = "TCGA-KIRC"
maf_fn = sprintf("GDCdata/harmonized/%s.mutSig.maf", proj_name) 
counts_fn = sprintf("%s_counts.rda", proj_name)

if(!file.exists(counts_fn)){
	counts_query <- GDCquery(project = proj_name, 
	                  data.category = "Transcriptome Profiling", 
	                  data.type = "Gene Expression Quantification", 
	                  workflow.type = "HTSeq - Counts")
	GDCdownload(counts_query)
	counts <- GDCprepare(query = counts_query, save = TRUE, 
         	         	 save.filename = counts_fn)
} else {
	load(counts_fn)
	counts <- data
	rm(data)
}

maf <- read.maf(maf_fn)

# leave silent mutaitons to the side for now. 
counts@colData$mtor_mut <- (counts@colData$bcr_patient_barcode %in% 
							substr(maf@data$Tumor_Sample_Barcode,1,16)[maf@data$Hugo_Symbol == "MTOR"])

dds <- DESeqDataSet(counts[,600:610], design = ~ mtor_mut)
dds <- DESeq(dds, parallel = T)
res <- results(dds)

gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = T)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")



