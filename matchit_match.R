# matchit_match.R
#
# Reusable control-gene matching workflow built on MatchIt. Given a set of genes
# of interest and the standard NTC DESeq2 + GGA outputs, it selects an
# expression- (and optionally length-, GC-, ...) matched set of unchanged control
# genes by 1:1 nearest-neighbour Mahalanobis matching, and writes the matched
# list plus the full match-quality QC (per-covariate quantile/density/ECDF plots
# and a standardized-mean-difference balance table).
#
# METHOD (--method, default "matchit"):
#   matchit     nearest-neighbour Mahalanobis matching on the --covars (recommended)
#   expression  legacy decile matcher (expressionMatched, expression-only, as the
#               deployed main.R); alias "legacy". Kept for reproducing older control
#               sets; ignores any covariates other than expression.
#
# COVARIATES (--covars, comma-separated; default "expression"):
#   expression   length-normalized counts/kb (baseline condition), log2(x+1)   [from makeHeatmap.rds]
#   length       gene length (DominantTES-TSS+1), log10                        [from GGA WORKING_FILE]
#   <gga column> any numeric column of the GGA WORKING_FILE, by name
#   <file column>any numeric column of --covar-file (GeneID + numeric columns)
#   log2:<name> / log10:<name>  request a transform on a GGA/file covariate
#   e.g.  --covars expression              (default: expression only)
#         --covars expression,length
#         --covars expression,length,gc    (gc supplied via --covar-file)
#
# COVARIATE SOURCES: makeHeatmap.rds (expression) -> GGA WORKING_FILE (length +
# its numeric columns) -> --covar-file (anything else, e.g. precomputed GC).
#
# Usage:
#   Rscript matchit_match.R --genes genes.txt --deseq-dir DIR --gga GGA.WORKING_FILE.txt \
#     --manifest DESeq_manifest.txt --baseline WT [--window 0] [--covars expression] \
#     [--covar-file extra.txt] [--ratio 1] [--outdir matchit_output]
#
# Requires: R (base) + the CRAN package 'MatchIt'.

p_thresh  <- 0.1    # padj cutoff: genes with padj > p_thresh are "not significant"
fc_thresh <- 1.3    # linear fold-change band [1/fc_thresh, fc_thresh] = "unchanged"

## --- locate this script and load the repo's matching + QC libraries ---
cArgs <- commandArgs()
runningScript <- normalizePath(sapply(strsplit(cArgs[which(startsWith(cArgs, "--file="))], "--file="), "[", 2))
here <- function(f) sub("matchit_match.R", f, runningScript)
source(here("matchitMatched.R"))
source(here("matchQC.R"))
source(here("expressionMatched.R"))   # legacy decile matcher (--method expression)
`%||%` <- function(a, b) if (is.null(a)) b else a   # default for an unset flag

## --- flag parsing + validation ---
usage <- paste(
  "Usage: Rscript matchit_match.R --genes FILE --deseq-dir DIR --gga FILE",
  "--manifest FILE --baseline NAME [--method matchit] [--window 0] [--covars expression]",
  "[--covar-file FILE] [--ratio 1] [--outdir matchit_output]")
argsArg <- which(startsWith(cArgs, "--args"))
args <- if (length(argsArg)) cArgs[(argsArg + 1):length(cArgs)] else character(0)
parse_flags <- function(a) {
  out <- list(); i <- 1
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("unexpected argument '", a[i], "'.\n", usage, call. = FALSE)
    key <- sub("^--", "", a[i])
    if (key == "help") { out$help <- TRUE; i <- i + 1; next }
    if (i == length(a)) stop("missing value for --", key, ".\n", usage, call. = FALSE)
    out[[key]] <- a[i + 1]; i <- i + 2
  }
  out
}
opt <- parse_flags(args)
if (isTRUE(opt$help) || !length(args)) { message(usage); quit(save = "no", status = 0) }

required <- c("genes", "deseq-dir", "gga", "manifest", "baseline")
missing_req <- setdiff(required, names(opt))
if (length(missing_req)) stop("missing required argument(s): --",
                              paste(missing_req, collapse = ", --"), "\n", usage, call. = FALSE)
window   <- suppressWarnings(as.numeric(opt[["window"]] %||% "0"))
ratio    <- suppressWarnings(as.integer(opt[["ratio"]] %||% "1"))
covspec  <- opt[["covars"]] %||% "expression"
out_dir  <- opt[["outdir"]] %||% "matchit_output"
# Matching algorithm: "matchit" (default) = nearest-neighbour Mahalanobis on the
# --covars; "expression" (alias "legacy") = the original decile expressionMatched,
# expression-only, matching the deployed main.R behaviour.
method <- tolower(opt[["method"]] %||% "matchit")
if (method == "legacy") method <- "expression"
if (!method %in% c("matchit", "expression")) {
  stop("--method must be 'matchit' (default) or 'expression'/'legacy'; got '", method, "'.", call. = FALSE)
}
if (is.na(window)) stop("--window must be numeric.", call. = FALSE)
if (is.na(ratio) || ratio < 1) stop("--ratio must be an integer >= 1.", call. = FALSE)
for (f in c(opt$genes, opt$gga, opt$manifest)) if (!file.exists(f)) stop("file not found: ", f, call. = FALSE)
if (!dir.exists(opt[["deseq-dir"]])) stop("DESeq2 dir not found: ", opt[["deseq-dir"]], call. = FALSE)
if (!is.null(opt[["covar-file"]]) && !file.exists(opt[["covar-file"]])) stop("covar-file not found: ", opt[["covar-file"]], call. = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes_of_interest <- scan(opt$genes, character(), quiet = TRUE)
if (!length(genes_of_interest)) stop("no gene IDs read from ", opt$genes, call. = FALSE)

## --- unchanged control pool: not significant AND within fold-change band, in EVERY
##     contrast (same rule as main.R; DESeq column 6 = padj, 7 = linear foldchange) ---
results_tables <- Sys.glob(paste(opt[["deseq-dir"]], "*", "*_results_table.txt", sep = "/"))
if (!length(results_tables)) stop("no */*_results_table.txt under ", opt[["deseq-dir"]], call. = FALSE)
results_list <- lapply(results_tables, read.table, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
unchanged <- Reduce(intersect, lapply(results_list, function(x)
  rownames(x)[!is.na(x[, 6]) & !is.na(x[, 7]) & x[, 6] > p_thresh & x[, 7] <= fc_thresh & x[, 7] >= fc_thresh^-1]))

## --- GGA: length, type filter, and source of named covariates ---
gga <- read.table(opt$gga, header = TRUE, sep = "\t", row.names = 1, quote = "", stringsAsFactors = FALSE)
for (col in c("GeneType", "DominantTSS", "DominantTES")) if (!col %in% colnames(gga)) stop("GGA missing column: ", col, call. = FALSE)
gene_len_all <- abs(gga$DominantTES - gga$DominantTSS) + 1; names(gene_len_all) <- rownames(gga)
# control pool: protein-coding, >= 1 kb (full TSS-TES length), minus genes of interest
pool_cand <- intersect(unchanged, rownames(gga))
pool <- setdiff(pool_cand[gga[pool_cand, "GeneType"] %in% "protein_coding" & gene_len_all[pool_cand] >= 1000],
                genes_of_interest)
if (!length(pool)) stop("no unchanged control genes passed the filters.", call. = FALSE)

## --- expression (counts/kb over baseline), for all annotated genes ---
mh_rds <- Sys.glob(paste(opt[["deseq-dir"]], "matrix", "*.makeHeatmap.rds", sep = "/"))
if (length(mh_rds) != 1) stop("expected exactly one matrix/*.makeHeatmap.rds; found ", length(mh_rds), ".", call. = FALSE)
mh <- readRDS(mh_rds)
meta <- read.table(opt$manifest, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
if (!"condition" %in% colnames(meta)) stop("manifest needs a 'condition' column.", call. = FALSE)
ref_samples <- rownames(meta)[meta$condition %in% opt$baseline]
if (!length(ref_samples)) stop("baseline '", opt$baseline, "' not in manifest 'condition' column.", call. = FALSE)
ref_cols <- make.names(sapply(strsplit(ref_samples, "/"), tail, 1), allow_ = FALSE)
if (length(setdiff(ref_cols, colnames(mh)))) stop("baseline sample column(s) missing from makeHeatmap matrix: ",
                                                  paste(setdiff(ref_cols, colnames(mh)), collapse = ", "), call. = FALSE)
common <- intersect(rownames(mh), rownames(gga))
counts_kb <- rowMeans(mh[common, ref_cols, drop = FALSE]) / ((gene_len_all[common] - window) / 1000)

## --- optional user covariate file (GeneID + numeric columns) ---
covfile <- NULL
if (!is.null(opt[["covar-file"]])) {
  covfile <- read.table(opt[["covar-file"]], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
}

## --- resolve the requested covariates into one matching frame ---
gga_num <- gga[, vapply(gga, is.numeric, logical(1)), drop = FALSE]   # numeric GGA columns
apply_tr <- function(v, tr) switch(tr, none = v, log2 = log2(v), log10 = log10(v))
resolve_one <- function(token) {
  tr <- "none"; nm <- token
  if (startsWith(token, "log2:"))  { tr <- "log2";  nm <- sub("^log2:",  "", token) }
  if (startsWith(token, "log10:")) { tr <- "log10"; nm <- sub("^log10:", "", token) }
  if (nm == "expression") return(list(name = "expression", value = log2(counts_kb + 1)))
  if (nm == "length")     return(list(name = "length",     value = log10(gene_len_all[common])))
  if (nm %in% colnames(gga_num))                 return(list(name = nm, value = apply_tr(gga_num[common, nm], tr)))
  if (!is.null(covfile) && nm %in% colnames(covfile)) return(list(name = nm, value = apply_tr(covfile[common, nm], tr)))
  stop("covariate '", nm, "' not found as a keyword, a numeric GGA column, or a --covar-file column.", call. = FALSE)
}
tokens <- trimws(strsplit(covspec, ",")[[1]])
resolved <- lapply(tokens, resolve_one)
dat <- data.frame(GeneID = common, stringsAsFactors = FALSE, row.names = common)
for (r in resolved) dat[[r$name]] <- as.numeric(r$value)
covar_names <- vapply(resolved, function(x) x$name, character(1))

## --- select the matched control set by the chosen method ---
if (method == "matchit") {
  message(sprintf("Method: matchit (nearest-neighbour Mahalanobis) on: %s",
                  paste(covar_names, collapse = ", ")))
  matched <- matchitMatched(genes_of_interest, pool, dat, covars = covar_names, ratio = ratio)
} else {
  # Legacy decile matcher: expression-only, 10 quantile groups, as in main.R.
  if (!identical(covar_names, "expression")) {
    message("NOTE: --method expression matches on expression ONLY; ignoring other ",
            "covariate(s): ", paste(setdiff(covar_names, "expression"), collapse = ", "))
  }
  message("Method: expression (legacy decile expressionMatched, expression only)")
  exp_dat <- data.frame(GeneID = common, Expression = as.numeric(counts_kb))  # raw counts/kb, as main.R
  matched <- tryCatch(
    unique(stats::na.omit(expressionMatched(genes_of_interest, pool, exp_dat, 10, setSeedNum = 1))),
    error = function(e) stop("legacy expression matching failed (", conditionMessage(e),
      "); a quantile bin likely has too few control genes - use --method matchit.", call. = FALSE))
}

## --- write the control list (filename stamped with the method for provenance) ---
out_list <- file.path(out_dir, sprintf("matched_controls.%s.txt", method))
writeLines(matched, out_list)
writeLines(c(sprintf("method\t%s", method),
             sprintf("covariates\t%s", if (method == "expression") "expression" else paste(covar_names, collapse = ",")),
             sprintf("baseline\t%s", opt$baseline),
             sprintf("n_genes_of_interest\t%d", length(genes_of_interest)),
             sprintf("n_matched\t%d", length(matched)),
             sprintf("generated\t%s", format(Sys.time()))),
           file.path(out_dir, sprintf("run_info.%s.txt", method)))
message(sprintf("Wrote %d matched control genes -> %s", length(matched), out_list))

## --- full match-quality QC (quantile/density/ECDF per covariate + balance table) ---
# QC always runs on the requested covariates, so a legacy run visibly shows any
# covariate (e.g. length) it left unbalanced.
qc_prefix <- file.path(out_dir, sprintf("match_qc.%s", method))
matchQC(genes_of_interest, matched, dat, covars = covar_names, all_genes = common,
        out_prefix = qc_prefix)
message(sprintf("QC written to %s_QC.png and %s_balance.txt", qc_prefix, qc_prefix))
