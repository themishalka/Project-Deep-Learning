# Investigating the transcriptomic diversity of gliomas 

## Project overview 

Gliomas are heterogeneous tumors originating from glial cells and represent the most common form of tumors in the central nervous system. They are divided into two main types, Low-grade gliomas (LGG, grade 1, 2 and 3) which contain several subtypes (e.g. astrocytomas or oligodendrogliomas), and Glioblastoma Multiforme, the most aggressive type (GBM, grade 4). Those different glioma grades and types were characterised because of specific associated phenotypes and degree of malignancy. In this project, we investigated whether the diversity of gliomas could be observed at a transcriptomic level, using two deep learning approaches, an autoencoder and a multi-layer perceptron (MLP) classifier. 

## Pipeline overview 

1. Download GBM and LGG Transcriptome profiling data from TCGA
2. Preprocessing (cleaning, normalisation & feature selection)
3. Unsupervised approach (autoencoder)
5. Supervised approach (MLP)

## Requirements 

The analysis pipeline was implemented on: 
- R version 4.5.1 (Data download, preprocessing & functional enrichment)
- Python version 3.13.5. (Model training and evaluation, hyperparameter tuning, analysis)

Deep learning models were built as Functional APIs in the TensorFlow (v2.20.0)-integrated Keras (v3.11.3) framework in Python. Scikit-learn v1.7.1 was used for complementary machine learning approaches (PCA, UMAP, t-SNE, kNN...). 

## Detailed workflow 

#### Data preparation

The project uses public TCGA bulkRNA-seq gene expression counts and clinical metadata (patient and sample data) retrieved from the Genomic Data Commons (GDC) using the R package `TCGAbiolinks`, within `TCGA_file_download.R`. This includes expression matrices for glioblastoma and low-grade gliomas (`gbm_expression.csv` and `lgg_expression.csv`). Gene annotations were also saved in two separate files (`gbm_genes.csv` and `lgg_genes.csv`). The sample metadata were saved in two rds files: `gbm_samples.rds` and `lgg_samples.rds.` The aforementioned can be found in the `data` folder. 

Preprocessing was conducted within the `preprocessing.Rmd` notebook, which expects `gbm_expression.csv`, `gbm_genes.csv`, `gbm_samples.rds`, `lgg_expression.csv`, `lgg_genes.csv`, `lgg_samples.rds` as input. The preprocessing script produces `merged_expression.csv`, a normalized and filtered expression matrix (samples x genes), `merged_genes.csv` (gene annotations for the retained genes) and `merged_samples.rds`, the merged sample metadata (GBM+LGG, aligned to expression). The aforementioned can be found in the `data` folder.

#### Unsupervised approach (autoencoder)

An autoencoder was built and trained to reconstruct the input preprocessed TCGA gene expression matrix following compression in the encoder. Since autoencoders are unsupervised methods, the `data/merged_expression.csv` was used as an input without label. The script `script/model_ae.py` enables the training of the autoencoder and outputs the validation loss (MSE & MAE) for a given set of parameters. The hyperparameter tuning pipeline was implemented in the notebook `script/tuning_ae.ipynb`. All parameter combinations tested and their associated MSE can be visualised for each stage in the `scripts/tuning_model` folder, as well as the top model(s) at each tuning stage. The final model was trained on the optimal set of parameters and saved alongside the latent representations `Z_train` and `Z_test` and the reconstruction of the input `X_train_rec` and `X_test_rec`; they are available in the `scripts/final_model` folder. The final model, test latent representation and train and test reconstructions were analysed in the `script/analysis_ae.ipynb` notebook. 

The produced plots presenting the hyperparameter tuning, quality of the reconstruction, informativeness of the latent representation and feature importance/functional enrichment are found in `results/tuning_ae`, `results/analysis_ae` and `results/enrichment` respectively. 

Useful custom functions were stored in the `scripts/utils.py`, necessary import for the `tuning_ae.ipynb` and `analysis_ae.ipynb` scripts. 

#### Supervised approach (MLP)

For our supervised approach, a Multi-Layer Perceptron model was built and trained to predict tumor grade from the preprocessed TCGA gene expression matrix. The model expects the outputs of the preprocessing step. Each run of the model creates a timestamped folder, e.g.: 20260105_181536 in the indicated output folder which needs to be hardcoded by the user. This includes a csv file with the sample IDs in the held-out test set (test_ids.csv), the validation accuracy and macro-F1 per cross-validation fold (cv_results.json) and the best hyoer-parameter configuration chosen by hyperband (best_hyperparameters.json). The final selected model is also saved (best_model.keras), as well as the final scaler used for test evaluation (scaler.joblib) and the final test metrics (test_metrics.json). Configuration and package versions can be found in config.json and versions.json for reproducibility. The aforementioned can be found in our results/MLP_final_results folder.

The produced plots presenting the model's evaluation and hyper-parameter selection are found in results/MLP_final_results/Tuning_and_Evaluation_Plots. Feature importance analysis csv and figures are stored in results/MLP_final_results/Feature_Importance. 

Some downstream analysis (e.g loading a trained model/reading tuner or integrated gradient outputs to produce the graphs) requires hardcoding run_dir to point to the specific timestamped run folder you want to inspect. 
