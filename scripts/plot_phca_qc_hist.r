 
  # =============================
 # 5) QC plot: max-weight + entropy distributions
 # =============================
 
 plot_pcha_qc_hist <- function(qc_df, value_col, title = NULL, bins = 50) {
   # -----------------------------
   # Inputs
   # -----------------------------
   if (!is.data.frame(qc_df)) stop("qc_df must be a data.frame.", call. = FALSE)
   if (!(value_col %in% colnames(qc_df))) stop("value_col not found in qc_df: ", value_col, call. = FALSE)
   if (is.null(title)) title <- paste0("QC: ", value_col)
   
   ggplot(qc_df, aes(x = .data[[value_col]])) +
     geom_histogram(bins = bins) +
     theme_classic(base_size = 13) +
     labs(title = title, x = value_col, y = "Cells") +
     theme(plot.title = element_text(face = "bold"),
           axis.title = element_text(face = "bold"))
 }


  qcA5 <- compute_pcha_weight_qc_from_W(W5)
 
 p_qc_max <- plot_pcha_qc_hist(qcA5, "max_weight", title = "Archetype Max Weight Per Cell")
 print(p_qc_max)
 ggsave("A5_QC_max_weight_hist.png", p_qc_max, width = 6, height = 4, dpi = 300)
 
 p_qc_ent <- plot_pcha_qc_hist(qcA5, "entropy_norm", title = "Archetype Entropy")
 print(p_qc_ent)
 ggsave("A5_QC_entropy_hist.png", p_qc_ent, width = 6, height = 4, dpi = 300)