
# GO ORA + GSEA for Glioma "Important Genes" (based on Integrated Gradients)


# ---- packages ----
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(forcats)
library(readr)

# ---- USER: set paths ----
ig_csv <- "C:/Users/joann/Desktop/M2/Deep_Learning/MLP_run/20260105_181536/final/grad_importance/grad_importance_per_class_gene_names.csv"  # <-- your file
out_dir <- "C:/Users/joann/Desktop/M2/Deep_Learning/MLP_run/20260105_181536/final"      # <-- output folder

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load IG table ----
ig <- readr::read_csv(ig_csv, show_col_types = FALSE)

# Expect columns:
# gene_name, mean_abs_ig_G2, mean_abs_ig_G3, mean_abs_ig_G4
required <- c("gene_name", "mean_abs_ig_G2", "mean_abs_ig_G3", "mean_abs_ig_G4")
missing <- setdiff(required, colnames(ig))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

# clean gene symbols
ig <- ig %>%
  mutate(gene_name = as.character(gene_name)) %>%
  mutate(gene_name = str_trim(gene_name)) %>%
  filter(!is.na(gene_name), gene_name != "")

# Universe:
# all genes in the expression matrix

universe_genes <- unique(ig$gene_name)

# ---- MSigDB gene sets for GSEA ----
# 1: GO BP 
m_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP")
pathways_bp <- split(m_bp$gene_symbol, m_bp$gs_name)

# 2: Hallmark
m_h <- msigdbr(species = "Homo sapiens", collection = "H")
pathways_h <- split(m_h$gene_symbol, m_h$gs_name)

# helper: wrap long pathway names for plots
wrap_name <- function(x, width=50) gsub("(.{1,50})(\\s|$)", "\\1\n", x)

# helper: save clusterProfiler barplot as png
save_cp_barplot <- function(enrich_obj, file, title, showCategory=15) {
  png(file, width=1200, height=800, res=150)
  print(barplot(enrich_obj, showCategory=showCategory, title=title))
  dev.off()
}

# helper: fgsea NES barplot (top +/-)
save_fgsea_barplot <- function(fg, file, title, top_n=20) {
  # fg: fgsea result data.frame
  fg <- fg %>% filter(!is.na(padj))
  if (nrow(fg) == 0) return(NULL)
  
  top_pos <- fg %>% filter(padj < 0.05, NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = top_n)
  top_neg <- fg %>% filter(padj < 0.05, NES < 0) %>% arrange(NES) %>% slice_head(n = top_n)
  top_both <- bind_rows(top_pos, top_neg)
  
  if (nrow(top_both) == 0) {
    # If nothing passes padj<0.05, we take top by padj for visualization
    top_both <- fg %>% arrange(padj) %>% slice_head(n = top_n)
  }
  
  top_both <- top_both %>%
    mutate(pathway_label = wrap_name(pathway, 50),
           pathway_label = fct_reorder(pathway_label, NES))
  
  p <- ggplot(top_both, aes(x=pathway_label, y=NES, fill=NES > 0)) +
    geom_col() +
    coord_flip() +
    labs(x=NULL, y="NES", title=title, fill="NES > 0") +
    theme_minimal() +
    theme(axis.text.y = element_text(size=8))
  
  ggsave(file, p, width=10, height=7, dpi=200)
}

# Convert list-columns (e.g., fgsea leadingEdge) into strings so write.csv works
flatten_list_cols <- function(df) {
  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    for (cn in names(df)[is_list_col]) {
      df[[cn]] <- vapply(df[[cn]], function(x) paste(x, collapse=";"), character(1))
    }
  }
  df
}


# Main loop over classes

classes <- c("G2", "G3", "G4")

for (cl in classes) {
  message("Running GO ORA + GSEA for class: ", cl)
  
  score_col <- paste0("mean_abs_ig_", cl)
  
  # ---- TOP 50 for GO ORA ----
  top50 <- ig %>%
    dplyr::select(gene_name, all_of(score_col)) %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    arrange(desc(.data[[score_col]])) %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::pull(gene_name) %>%
    unique()
  
  # ---- GO ORA (Biological Process) ----
  
  # testing whether important genes are over-represented in GO terms.
  ego_bp <- enrichGO(
    gene          = top50,
    universe      = universe_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  ego_df <- as.data.frame(ego_bp)
  write.csv(ego_df, file=file.path(out_dir, paste0("GO_ORA_BP_top50_", cl, ".csv")), row.names = FALSE)
  
  if (!is.null(ego_bp) && nrow(ego_df) > 0) {
    save_cp_barplot(
      ego_bp,
      file=file.path(out_dir, paste0("GO_ORA_BP_barplot_top50_", cl, ".png")),
      title=paste0("GO BP enrichment (ORA) — Top 50 important genes (", cl, ")"),
      showCategory=15
    )
  }
  
  # ---- GSEA (fgsea) using full ranking (GO BP) ----
  # importance scores are non-negative (abs IG),
  # so GSEA here means: which pathways are enriched among the most "important" genes.
  ranks <- ig %>%
    dplyr::select(gene_name, all_of(score_col)) %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::group_by(gene_name) %>%                 # handle duplicates
    dplyr::summarize(score = max(.data[[score_col]]), .groups="drop") %>%
    dplyr::arrange(desc(score))
  
  stats <- ranks$score
  names(stats) <- ranks$gene_name
  stats <- sort(stats, decreasing = TRUE)
  
  set.seed(1)
  fg_bp <- fgsea(
    pathways = pathways_bp,
    stats    = stats,
    nperm    = 10000,
    minSize  = 10,
    maxSize  = 500
  ) %>% as.data.frame() %>% arrange(padj)
  
  fg_bp <- flatten_list_cols(fg_bp)
  write.csv(fg_bp, file=file.path(out_dir, paste0("GSEA_fGSEA_GO_BP_", cl, ".csv")), row.names = FALSE)
  
  save_fgsea_barplot <- function(fg, file, title, top_n=20) {
    fg <- fg %>% dplyr::filter(!is.na(padj), !is.na(NES))
    if (nrow(fg) == 0) return(NULL)
    
    # choose top terms for plotting
    top_pos <- fg %>% dplyr::filter(padj < 0.05, NES > 0) %>% dplyr::arrange(dplyr::desc(NES)) %>% dplyr::slice_head(n = top_n)
    top_neg <- fg %>% dplyr::filter(padj < 0.05, NES < 0) %>% dplyr::arrange(NES) %>% dplyr::slice_head(n = top_n)
    top_both <- dplyr::bind_rows(top_pos, top_neg)
    
    # fallback if nothing passes padj<0.05
    if (nrow(top_both) == 0) {
      top_both <- fg %>% dplyr::arrange(padj) %>% dplyr::slice_head(n = top_n)
    }
    
    # add -log10(padj) for display
    top_both <- top_both %>%
      dplyr::mutate(
        neglog10_padj = -log10(padj + 1e-300),  # avoid log(0)
        pathway_label = wrap_name(pathway, 50),
        pathway_label = forcats::fct_reorder(pathway_label, NES),
        nes_sign = NES > 0
      )
    
    # position label slightly beyond bar end
    top_both <- top_both %>%
      dplyr::mutate(
        label_x = ifelse(NES >= 0, NES + 0.03, NES - 0.03),
        label_hjust = ifelse(NES >= 0, 0, 1),
        label_txt = sprintf("-log10(FDR)=%.2f", neglog10_padj)
      )
    
    p <- ggplot(top_both, aes(x = pathway_label, y = NES, fill = nes_sign)) +
      geom_col() +
      coord_flip() +
      geom_text(aes(y = label_x, label = label_txt, hjust = label_hjust),
                size = 3) +
      labs(x = NULL, y = "NES", title = title) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8)) +
      guides(fill = "none")  # <-- removes TRUE/FALSE legend
    
    ggsave(file, p, width = 11, height = 7, dpi = 200)
  }
  
  
  set.seed(1)
  fg_h <- fgsea(
    pathways = pathways_h,
    stats    = stats,
    nperm    = 10000,
    minSize  = 5,
    maxSize  = 500
  ) %>% as.data.frame() %>% arrange(padj)
  
  fg_h <- flatten_list_cols(fg_h)
  write.csv(fg_h, file=file.path(out_dir, paste0("GSEA_fGSEA_HALLMARK_", cl, ".csv")), row.names = FALSE)
  
  save_fgsea_barplot(
    fg_h,
    file=file.path(out_dir, paste0("GSEA_HALLMARK_barplot_", cl, ".png")),
    title=paste0("GSEA (fgsea) — Hallmark ranked by IG importance (", cl, ")"),
    top_n=20
  )
} 

message("Done. Results saved in: ", out_dir)



