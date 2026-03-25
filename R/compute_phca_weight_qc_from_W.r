 # =============================
 # 5) QC plot: max-weight + entropy distributions
 # =============================
 compute_pcha_weight_qc_from_W <- function(W, near_zero = 1e-3) {
   # -----------------------------
   # Inputs
   # -----------------------------
   if (!is.matrix(W)) stop("W must be a matrix (cells x k).", call. = FALSE)
   if (is.null(rownames(W))) stop("W must have rownames (cell IDs).", call. = FALSE)
   
   rs <- rowSums(W)
   mx <- apply(W, 1L, max)
   
   eps <- 1e-12
   Wp <- pmax(W, 0)
   Wn <- Wp / pmax(rowSums(Wp), eps)
   ent <- -rowSums(Wn * log(pmax(Wn, eps)))
   ent_norm <- ent / log(ncol(W))
   spars <- rowMeans(W < near_zero)
   
   data.frame(
     cell_id = rownames(W),
     row_sum = rs,
     max_weight = mx,
     entropy_norm = ent_norm,
     sparsity_near_zero = spars,
     stringsAsFactors = FALSE
   )
 }
 

