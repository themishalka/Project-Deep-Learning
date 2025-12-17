############################################################
### Downloading TCGA expression matrix and clinical data ###
############################################################

### Loading packages ---------------------------------------

#BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SingleCellExperiment)
library(jsonlite)


### Querying TCGA expression data from GDC -----------------

query <- GDCquery(
  project = "TCGA-GBM", # cancer type: glioblastoma
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
getResults(query) #visualising data retrieved

### Downloading raw files ----------------------------------

GDCdownload(query, directory = "../data/raw_gbm_files")

### Parsing the data into a SummarizedExperiment object ----

gbm_exp <- GDCprepare(query,directory = "../data/raw_gbm_files")

gbm_exp

counts <- assay(gbm_exp) # extracting raw counts


### Querying TCGA clinical data from GDC -------------------

clinical_gbm <- GDCquery_clinic(
  project = "TCGA-GBM",
  type = "clinical"
)


### Downloading processed files as CSV ---------------------

# Expression matrix (genes X samples)

expr_df <- as.data.frame(counts)
expr_df$gene_id <- rownames(expr_df)
expr_df <- expr_df[, c("gene_id", setdiff(colnames(expr_df), "gene_id"))]

# Sample metadata (samples X clinical/covariates)

sample_df <- as.data.frame(colData(gbm_exp))
sample_df$sample_id <- rownames(sample_df)

# Gene annotations (genes X annotations)

gene_df <- as.data.frame(rowData(gbm_exp))
gene_df$gene_id <- rownames(gene_df)

# Saving to data folder 

write.csv(expr_df, "../data/gbm_expression.csv", row.names = FALSE)
write.csv(gene_df, "../data/gbm_genes.csv", row.names = FALSE)

#samples_df has listed columns so we apply transformation to save and load in python
df <- as.data.frame(sample_df, stringsAsFactors = FALSE)
is_list <- vapply(df, is.list, logical(1))

df[is_list] <- lapply(df[is_list], function(col) {
  vapply(col, function(el) toJSON(el, auto_unbox = TRUE, null = "null"), character(1))
})

# Now df has only scalar columns (strings/numbers) 
saveRDS(df, "C:/Users/joann/Desktop/M2/Deep_Learning/data/gbm_samples.rds", version=2)



















