# matchitMatched.R
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
