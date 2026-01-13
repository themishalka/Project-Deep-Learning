<<<<<<< HEAD
The project uses public TCGA RNA-seq gene expression counts and clinical metadata(patient and sample data) retrieved from the Genomic Data Commons (GDC) using the R package TCGAbiolinks, within  TCGA_file_download.R. This includes expression matrices for glioblastoma and low-grade gliomas (gbm_expression.csv and lgg_expression.csv). Gene annotations were also saved in two separate files (gbm_genes.csv and lgg_genes.csv). The sample metadata were saved in two rds files: gbm_samples.rds and lgg_samples.rds. The aforementioned can be found in our data folder. 

Preprocessing was conducted within the preprocessing.Rmd, which expects gbm_expression.csv, gbm_genes.csv, gbm_samples,rds, lgg_expression.csv, lgg_genes.csv, lgg_samples.rds as input. The preprocessing script produces merged_expression.csv, a normalized and filtered expression matrix (samples x genes), merged_genes.csv (gene annotations for the retained genes) and merged_samples.rds, the merged sample metadata (GBM+LGG, aligned to expression). The aforementioned can be found in our data folder.

For our supervised approach, a Multi-Layer Perceptron model was built and trained to predict tumor grade from the preprocessed TCGA gene expression matrix. The model expects the outputs of the preprocessing step. Each run of the model creates a timestamped folder, e.g.: 20260105_181536 in the indicated output folder which needs to be hardcoded by the user. This includes a csv file with the sample IDs in the held-out test set (test_ids.csv), the validation accuracy and macro-F1 per cross-validation fold (cv_results.json) and the best hyoer-parameter configuration chosen by hyperband (best_hyperparameters.json). The final selected model is also saved (best_model.keras), as well as the final scaler used for test evaluation (scaler.joblib) and the final test metrics (test_metrics.json). Configuration and package versions can be found in config.json and versions.json for reproducibility. The aforementioned can be found in our results/MLP_final_results folder.

The produced plots presenting the model's evaluation and hyper-parameter selection are found in results/MLP_final_results/Tuning_and_Evaluation_Plots. Feature importance analysis csv and figures are stored in results/MLP_final_results/Feature_Importance. 

Some downstream analysis (e.g loading a trained model/reading tuner or integrated gradient outputs to produce the graphs) requires hardcoding run_dir to point to the specific timestamped run folder you want to inspect. 
=======
The data for this projected includes bulk RNA-seq and corresponding metadata from low-grade gliomas and glioblastomas as they appear in the TCGA (Cancer Genome Atlas) dataset.
They were downloaded and saved using the TCGA_file_download.R script and initial preprocessing was conducted within preprocessing.Rmd which resulted in files that were saved locally and reloaded in the next steps of the analysis.
Regarding the supervised approach, a Multi-Layer Perceptron model was built and trained to predict tumor grade from gene expression data. All steps of this analysis, from hyperparameter tuning to feature importance
can be found in the mlp_tumor_grade.ipynb folder. The final model run, including the hyperparameter tuning history and the final configuration, as well as primary results
is available in our results script.
Functional enrichment analysis on the identified important genes for the MLP model was performed within the enrichment_important_genes.R script
>>>>>>> 27ff31e1c94bfd8d4d5e52589f345a650a9b214e
