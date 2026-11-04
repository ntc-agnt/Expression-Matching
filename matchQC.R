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

  # --- distribution plots: quantile, density, and ECDF per covariate ---
  # Each row is one covariate; the three views (value-by-quantile, kernel
  # density, and cumulative distribution / ECDF) all compare the genes of
  # interest (black) with the matched controls (red dashed), against all genes
  # (grey) as the background. When the red and black curves overlie, that
  # covariate is well matched; a gap from the grey curve shows the matched set
  # is a targeted control rather than a random genome-wide sample.
  if (!capabilities("png")) {
    warning("no PNG device available; skipping distribution plots.", call. = FALSE)
    return(invisible(balance))
  }
  qg <- seq(0, 1, 0.02)
  qq <- function(v) stats::quantile(v, qg, na.rm = TRUE)
  pal <- c(all = "grey70", input = "black", matched = "red2")
  png(paste0(out_prefix, "_QC.png"),
      width = 1500, height = 430 * length(covar_cols), res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(length(covar_cols), 3), mar = c(4.2, 4.6, 3.2, 1))
  first <- TRUE
  for (col in covar_cols) {
    iv <- input_v[[col]]; mv <- matched_v[[col]]
    av <- if (!is.null(all_genes)) getv(all_genes, col) else NULL

    # (1) value by quantile
    yl <- range(c(iv, mv, av), finite = TRUE)
    plot(qg, qq(iv), type = "l", lwd = 2.6, col = pal["input"], ylim = yl,
         xlab = "quantile", ylab = col,
         main = sprintf("%s  (SMD = %+.3f) - quantile", col, bal[[col]]))
    if (!is.null(av)) lines(qg, qq(av), lwd = 3, col = pal["all"])
    lines(qg, qq(iv), lwd = 2.6, col = pal["input"])     # redraw input on top of grey
    lines(qg, qq(mv), lwd = 2,   col = pal["matched"], lty = 2)

    # (2) kernel density (needs >= 2 finite values per set)
    plot_density <- length(iv[is.finite(iv)]) >= 2 && length(mv[is.finite(mv)]) >= 2
    if (plot_density) {
      dI <- density(iv, na.rm = TRUE); dM <- density(mv, na.rm = TRUE)
      dA <- if (!is.null(av) && length(av[is.finite(av)]) >= 2) density(av, na.rm = TRUE) else NULL
      ym <- max(c(dI$y, dM$y, if (!is.null(dA)) dA$y))
      plot(dI, col = pal["input"], lwd = 2.6, ylim = c(0, ym),
           xlab = col, ylab = "density", main = paste0(col, " - density"))
      if (!is.null(dA)) { polygon(dA, col = "grey92", border = NA); lines(dA, col = pal["all"], lwd = 3) }
      lines(dI, col = pal["input"], lwd = 2.6)
      lines(dM, col = pal["matched"], lwd = 2, lty = 2)
    } else {
      plot.new(); title(main = paste0(col, " - density (too few points)"))
    }
    if (first) legend("topright", bty = "n", lwd = 2, lty = c(1, 1, 2), col = pal,
                      legend = c("all genes", "genes of interest (input)", "matched control"),
                      cex = 0.8)

    # (3) empirical cumulative distribution (ECDF)
    plot(stats::ecdf(iv), col = pal["input"], lwd = 2.6, do.points = FALSE, verticals = TRUE,
         xlab = col, ylab = "cumulative fraction", main = paste0(col, " - ECDF"))
    if (!is.null(av)) plot(stats::ecdf(av), add = TRUE, col = pal["all"], lwd = 3,
                           do.points = FALSE, verticals = TRUE)
    plot(stats::ecdf(iv), add = TRUE, col = pal["input"], lwd = 2.6, do.points = FALSE, verticals = TRUE)
    plot(stats::ecdf(mv), add = TRUE, col = pal["matched"], lwd = 2, lty = 2,
         do.points = FALSE, verticals = TRUE)
    first <- FALSE
  }
  invisible(balance)
}
