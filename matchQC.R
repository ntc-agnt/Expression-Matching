# matchQC.R
#
# Quality control for a matched control-gene set. Reports covariate balance and
# draws distribution plots comparing the genes of interest ("input") to the
# matched control genes (and, optionally, to all genes as a genome-wide
# reference), so a biologist can judge at a glance whether the match is good.
#
# Balance metric: standardized mean difference (SMD) per covariate,
#   SMD = (mean_matched - mean_input) / sd_input.
# By convention |SMD| < 0.1 means the matched set is well balanced on that
# covariate; larger values mean it differs systematically from the input.
#
# Arguments:
#   geneset     genes of interest (the input/focal set).
#   matchedset  the matched control genes to QC.
#   dat         data.frame: first column GeneID, remaining columns the covariates
#               (same frame used for matching).
#   covars      covariate column names (default: all but the first).
#   all_genes   optional vector of gene IDs for a genome-wide reference line.
#   out_prefix  output file prefix; writes <prefix>_balance.txt and
#               <prefix>_QC.png.
#
# Returns (invisibly) the balance table (data.frame).

matchQC <- function(geneset, matchedset, dat, covars = NULL, all_genes = NULL,
                    out_prefix = "match_qc") {
  ids <- as.character(dat[[1]])
  rownames(dat) <- ids
  covar_cols <- if (is.null(covars)) colnames(dat)[-1] else covars
  getv <- function(g, col) as.numeric(dat[intersect(as.character(g), ids), col])

  input_v   <- lapply(covar_cols, function(c) getv(geneset, c))
  matched_v <- lapply(covar_cols, function(c) getv(matchedset, c))
  names(input_v) <- names(matched_v) <- covar_cols

  # --- covariate-balance table ---
  smd <- function(m, o) {
    if (length(m) < 1 || length(o) < 2 || sd(o) == 0) return(NA_real_)
    (mean(m) - mean(o)) / sd(o)
  }
  bal <- vapply(covar_cols, function(c) smd(matched_v[[c]], input_v[[c]]), numeric(1))
  balance <- data.frame(covariate = covar_cols,
                        SMD = round(bal, 4),
                        balanced = ifelse(abs(bal) <= 0.1, "yes", "NO"),
                        stringsAsFactors = FALSE)
  write.table(balance, paste0(out_prefix, "_balance.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Covariate balance of matched set vs input (|SMD| <= 0.1 = balanced):")
  for (i in seq_len(nrow(balance))) {
    message(sprintf("  %-16s SMD = %+.3f   %s", balance$covariate[i], balance$SMD[i],
                    ifelse(balance$balanced[i] == "yes", "OK", "IMBALANCED")))
  }
  if (any(balance$balanced == "NO")) {
    message("  NOTE: at least one covariate is imbalanced (|SMD| > 0.1); the ",
            "matched set may not be a fair control on that covariate.")
  }

  # --- distribution plots: one quantile panel per covariate ---
  # Each panel shows the value at each quantile for the input genes (black) and
  # the matched controls (red dashed); when the two curves overlie, the match is
  # good. The grey "all genes" line shows the genome-wide background the controls
  # are drawn from.
  if (!capabilities("png")) {
    warning("no PNG device available; skipping distribution plots.", call. = FALSE)
    return(invisible(balance))
  }
  qg <- seq(0, 1, 0.02)
  qq <- function(v) stats::quantile(v, qg, na.rm = TRUE)
  png(paste0(out_prefix, "_QC.png"),
      width = 540 * length(covar_cols), height = 540, res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(1, length(covar_cols)), mar = c(4.2, 4.6, 3.2, 1))
  for (col in covar_cols) {
    iv <- input_v[[col]]; mv <- matched_v[[col]]
    av <- if (!is.null(all_genes)) getv(all_genes, col) else NULL
    yl <- range(c(iv, mv, av), finite = TRUE)
    plot(qg, qq(iv), type = "l", lwd = 2.6, col = "black", ylim = yl,
         xlab = "quantile", ylab = col,
         main = sprintf("%s  (SMD = %+.3f)", col, bal[[col]]))
    if (!is.null(av)) lines(qg, qq(av), lwd = 3, col = "grey70")
    lines(qg, qq(iv), lwd = 2.6, col = "black")        # redraw input on top of grey
    lines(qg, qq(mv), lwd = 2,   col = "red2", lty = 2)
    legend("bottomright", bty = "n", lwd = 2, lty = c(1, 1, 2),
           col = c("grey70", "black", "red2"),
           legend = c("all genes", "genes of interest (input)", "matched control"),
           cex = 0.85)
  }
  invisible(balance)
}
