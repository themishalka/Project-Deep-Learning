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
library(ggridges)

# ---- USER: set paths ----
ig_csv <- "../results/FI_ae.csv"      # <-- your file
out_dir <- "../results/enrichment/"   # <-- output folder
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load IG table ----
ig <- readr::read_csv(ig_csv, show_col_types = FALSE, col_select = -1)

# Expect columns: gene_name, FI
required <- c("gene_name", "FI")
missing <- setdiff(required, colnames(ig))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

# ---- clean gene symbols ----
ig$ensembl_id <- sub("\\..*$", "", ig$gene_name) # remove version
map <- bitr(
  ig$ensembl_id,
  fromType = "ENSEMBL",
  toType   = "SYMBOL",
  OrgDb    = org.Hs.eg.db
)
ig <- ig %>%
  dplyr::left_join(map, by = c("ensembl_id" = "ENSEMBL")) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::mutate(gene_name = str_trim(as.character(SYMBOL))) %>%
  dplyr::filter(!is.na(gene_name), gene_name != "")

# Universe of genes
universe_genes <- unique(ig$gene_name)

# ---- MSigDB gene sets ----
m_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP")
pathways_bp <- split(m_bp$gene_symbol, m_bp$gs_name)

m_h <- msigdbr(species = "Homo sapiens", collection = "H")
pathways_h <- split(m_h$gene_symbol, m_h$gs_name)

# ---- helper functions ----
wrap_name <- function(x, width=50) gsub("(.{1,50})(\\s|$)", "\\1\n", x)

flatten_list_cols <- function(df) {
  is_list_col <- vapply(df, is.list, logical(1))
  if (any(is_list_col)) {
    for (cn in names(df)[is_list_col]) {
      df[[cn]] <- vapply(df[[cn]], function(x) paste(x, collapse=";"), character(1))
    }
  }
  df
}

# GO ORA barplot (clusterProfiler object)
save_cp_barplot <- function(enrich_obj, file, title, showCategory=15) {
  png(file, width=1200, height=800, res=150)
  print(barplot(enrich_obj, showCategory=showCategory, title=title))
  dev.off()
}

# Ridge plot for fgsea results
save_fgsea_ridgeplot <- function(fg, file, title, top_n = 20) {
  
  fg <- fg %>% 
    dplyr::filter(!is.na(padj), !is.na(NES))
  
  if (nrow(fg) == 0) return(NULL)
  
  # Select top pathways by significance
  top_pos <- fg %>% 
    dplyr::filter(NES > 0) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::slice_head(n = top_n)
  
  top_neg <- fg %>% 
    dplyr::filter(NES < 0) %>% 
    dplyr::arrange(padj) %>% 
    dplyr::slice_head(n = top_n)
  
  top_both <- dplyr::bind_rows(top_pos, top_neg)
  
  # Clean pathway names + labels
  top_both <- top_both %>%
    dplyr::mutate(
      pathway_clean = pathway %>%
        stringr::str_remove("^GOBP_") %>%
        stringr::str_remove("^HALLMARK_") %>%
        stringr::str_replace_all("_", " "),
      pathway_clean = forcats::fct_reorder(pathway_clean, NES),
      fdr_label = paste0("FDR = ", formatC(padj, format = "e", digits = 2))
    )
  
  p <- ggplot(top_both, aes(y = pathway_clean)) +
    geom_segment(
      aes(x = 0, xend = NES, yend = pathway_clean, color = NES > 0),
      linewidth = 1
    ) +
    geom_point(
      aes(
        x = NES,
        color = NES > 0,
        size = -log10(padj)
      )
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = title,
      color = "NES > 0",
      size  = expression(-log[10]~FDR)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  
  ggsave(file, p, width = 11, height = 8, dpi = 200)
}



# ---- Main loop ----
classes <- c("ae")

for (cl in classes) {
  message("Running GO ORA + GSEA for class: ", cl)
  
  score_col <- "FI"
  
  # ---- TOP 50 for ORA ----
  top50 <- ig %>%
    dplyr::select(gene_name, all_of(score_col)) %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::arrange(dplyr::desc(.data[[score_col]])) %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::pull(gene_name) %>%
    unique()
  
  # ---- GO ORA ----
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
  
  # Save GO ORA barplot
  if (!is.null(ego_bp) && nrow(ego_df) > 0) {
    save_cp_barplot(
      ego_bp,
      file=file.path(out_dir, paste0("GO_ORA_BP_barplot_top50_", cl, ".png")),
      title=paste0("GO BP enrichment (ORA) — Top 50 important genes (", cl, ")"),
      showCategory=15
    )
  }
  
  # ---- GSEA ranks ----
  ranks <- ig %>%
    dplyr::select(gene_name, all_of(score_col)) %>%
    dplyr::filter(!is.na(.data[[score_col]])) %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(score = max(.data[[score_col]]), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(score))
  
  stats <- ranks$score
  names(stats) <- ranks$gene_name
  stats <- sort(stats, decreasing = TRUE)
  
  # ---- GSEA GO BP ----
  set.seed(1)
  fg_bp <- fgsea(
    pathways = pathways_bp,
    stats    = stats,
    nperm    = 10000,
    minSize  = 10,
    maxSize  = 500
  ) %>% as.data.frame() %>% dplyr::arrange(padj)
  fg_bp <- flatten_list_cols(fg_bp)
  write.csv(fg_bp, file=file.path(out_dir, paste0("GSEA_fGSEA_GO_BP_", cl, ".csv")), row.names = FALSE)
  
  save_fgsea_ridgeplot(
    fg_bp,
    file = file.path(out_dir, paste0("GSEA_GO_BP_ridgeplot_", cl, ".png")),
    title = paste0("GSEA — GO Biological Process (", cl, ")"),
    top_n = 20
  )
  
  # ---- GSEA Hallmark ----
  set.seed(1)
  fg_h <- fgsea(
    pathways = pathways_h,
    stats    = stats,
    nperm    = 10000,
    minSize  = 5,
    maxSize  = 500
  ) %>% as.data.frame() %>% dplyr::arrange(padj)
  fg_h <- flatten_list_cols(fg_h)
  write.csv(fg_h, file=file.path(out_dir, paste0("GSEA_fGSEA_HALLMARK_", cl, ".csv")), row.names = FALSE)
  
  save_fgsea_ridgeplot(
    fg_h,
    file = file.path(out_dir, paste0("GSEA_HALLMARK_ridgeplot_", cl, ".png")),
    title = paste0("GSEA — Hallmark (", cl, ")"),
    top_n = 20
  )
}

message("Done. Results saved in: ", out_dir)
