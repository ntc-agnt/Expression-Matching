# functions.R
#
# Shared function library for the Expression-Matching tools. Everything the
# drivers source lives here, in one place:
#   expressionMatched()  original decile (quantile) expression matcher (legacy;
#                         used by main.R and by matchit_match.R --method expression)
#   matchitMatched()     nearest-neighbour Mahalanobis covariate matcher (MatchIt)
#   matchQC()            match-quality QC (balance table + distribution plots)
#
# Sourced by main.R, matchit_match.R, and tests/benchmark_matchers.R.

## ===========================================================================
## expressionMatched() - legacy decile expression matcher
## ===========================================================================
expressionMatched <- function(geneset,ctrlgeneset,dat,numQuantGroups,CtrlSetRatio=1,setSeedNum=1){  
  #pull out the relevant columns from the dataframe given
  geneIDs <- as.character(dat[,1])
  expressionValues <- as.numeric(dat[,2])
  
  
  #pull out geneset quantile values
  de_set_quantiles <- quantile(expressionValues[geneIDs %in% as.character(geneset)],seq(0,1,1/numQuantGroups))
  
  #calculate the number of ctrl genes to pull per quantile group
  num_per_quantile <- CtrlSetRatio*round(length(geneset)/numQuantGroups,0)
  
  #create an empty character object to add expression matched genes too
  ctrlGenes <- character()
  
  #loop to pull out genes from each quantile group
  for (i in 1:(length(de_set_quantiles)-1)){
    #use set.seed to make the random sampling reproducible
    set.seed(setSeedNum)
    #Randomly select "num_per_quantile" genes from the ctrlgeneset for each quantile group of geneset expression
    x <- sample(geneIDs[geneIDs %in% as.character(ctrlgeneset) & 
                          expressionValues >= de_set_quantiles[i] & 
                          expressionValues < de_set_quantiles[(i+1)]],num_per_quantile)
    #add the genes to the growing vector of ctrlGenes                              
    ctrlGenes <- c(ctrlGenes,x)
    
  }
  #return the finished list of ctrlGenes
  return(ctrlGenes)
}

## ===========================================================================
## matchitMatched() - nearest-neighbour Mahalanobis covariate matcher
## ===========================================================================
#
# Covariate-matched control-gene selection via 1:1 nearest-neighbour matching
# (the CRAN package MatchIt). This is an alternative to expressionMatched() that
# can match on MORE THAN ONE covariate at once -- e.g. length-normalized
# expression AND gene length.
#
# Why length too? Decile expression matching balances expression only and leaves
# gene length unbalanced, which is a real confounder for PRO-seq (length affects
# read counts and pausing metrics). Nearest-neighbour matching on the Mahalanobis
# distance pairs each gene of interest with its closest *non-interest* background
# gene in the joint covariate space, so every supplied covariate is balanced
# together. The matching is deterministic (no random seeds) and stays accurate
# for small gene sets, where quantile binning and propensity-rejection sampling
# both struggle.
#
# Requires the CRAN package 'MatchIt' (install.packages("MatchIt")).
#
# Arguments:
#   geneset     character vector of genes of interest (the "focal" set).
#   ctrlgeneset character vector of candidate control genes (the background pool).
#   dat         data.frame whose first column is the GeneID and whose remaining
#               columns are the numeric covariates to match on. Covariate column
#               NAMES are used in the matching formula.
#   covars      optional names of the covariate columns to use (default: every
#               column of 'dat' except the first).
#   ratio       number of control genes to match per gene of interest (default 1).
#
# Returns: character vector of matched control gene IDs (a subset of ctrlgeneset;
#          never contains a gene from geneset).

matchitMatched <- function(geneset, ctrlgeneset, dat, covars = NULL, ratio = 1) {
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("matchitMatched() needs the 'MatchIt' package; install it with ",
         "install.packages(\"MatchIt\").", call. = FALSE)
  }
  if (ncol(dat) < 2) {
    stop("'dat' must have a GeneID column plus at least one covariate column; ",
         "got ", ncol(dat), " column(s).", call. = FALSE)
  }
  ids <- as.character(dat[[1]])
  covar_cols <- if (is.null(covars)) colnames(dat)[-1] else covars
  missing_cols <- setdiff(covar_cols, colnames(dat))
  if (length(missing_cols)) {
    stop("covariate column(s) not found in 'dat': ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  rownames(dat) <- ids

  geneset <- intersect(as.character(geneset), ids)
  # Controls must be disjoint from the genes of interest, so a gene of interest
  # can never be returned as its own control.
  ctrl <- setdiff(intersect(as.character(ctrlgeneset), ids), geneset)
  if (length(geneset) < 1) {
    stop("none of the genes of interest are present in 'dat'.", call. = FALSE)
  }

  # Build the matching frame and drop any gene with a non-finite covariate
  # (e.g. zero expression after a log, or a negative length after the window
  # subtraction); matchit cannot use those rows.
  use <- c(geneset, ctrl)
  cov <- as.data.frame(lapply(dat[use, covar_cols, drop = FALSE], as.numeric))
  rownames(cov) <- use
  keep <- apply(cov, 1, function(r) all(is.finite(r)))
  cov <- cov[keep, , drop = FALSE]
  treat <- as.integer(rownames(cov) %in% geneset)
  n_focal <- sum(treat == 1L)
  n_ctrl  <- sum(treat == 0L)
  if (n_focal < 1) {
    stop("no genes of interest remain after dropping non-finite covariates.",
         call. = FALSE)
  }
  if (n_ctrl < n_focal * ratio) {
    stop("control pool too small: ", n_ctrl, " usable candidate(s) for ",
         n_focal, " gene(s) of interest at ratio ", ratio,
         ". Loosen the pool definition or lower 'ratio'.", call. = FALSE)
  }

  mf <- data.frame(treat = treat, cov, check.names = FALSE)
  matching_formula <- stats::as.formula(
    paste("treat ~", paste(covar_cols, collapse = " + ")))

  # 1:1 nearest-neighbour matching on the Mahalanobis distance. Mahalanobis
  # matches directly in covariate space (scaled by the covariate covariance),
  # rather than collapsing covariates to a single propensity score, so each
  # covariate is balanced jointly. Deterministic apart from rare exact-distance
  # ties; set the seed so those ties resolve reproducibly.
  set.seed(1)
  fit <- MatchIt::matchit(matching_formula, data = mf, method = "nearest",
                          distance = "mahalanobis", replace = FALSE, ratio = ratio)
  matched <- MatchIt::match.data(fit)
  control_ids <- rownames(matched)[matched$treat == 0L]
  if (length(control_ids) < 1) {
    stop("matching returned no control genes.", call. = FALSE)
  }
  control_ids
}


## ===========================================================================
## matchQC() - match-quality QC
## ===========================================================================
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
