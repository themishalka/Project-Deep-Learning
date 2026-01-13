The data for this projected includes bulk RNA-seq and corresponding metadata from low-grade gliomas and glioblastomas as they appear in the TCGA (Cancer Genome Atlas) dataset.
They were downloaded and saved using the TCGA_file_download.R script and initial preprocessing was conducted within preprocessing.Rmd which resulted in files that were saved locally and reloaded in the next steps of the analysis.
Regarding the supervised approach, a Multi-Layer Perceptron model was built and trained to predict tumor grade from gene expression data. All steps of this analysis, from hyperparameter tuning to feature importance
can be found in the mlp_tumor_grade.ipynb folder. The final model run, including the hyperparameter tuning history and the final configuration, as well as primary results
is available in our results folder.
Functional enrichment analysis on the identified important genes for the MLP model was performed within the enrichment_important_genes.R script
