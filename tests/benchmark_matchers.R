# benchmark_matchers.R
#
# Stress-test and compare control-gene MATCHING methods on contrived-but-plausible
# synthetic data. Builds a realistic gene population (lognormal lengths, expression
# with a weak length dependence, and GC content as a third covariate), then defines
# several focal "genes of interest" sets each confounded in a way that exploits one
# method's blind spot, and reports how well each method balances every covariate.
#
# Methods compared (a method is skipped with a note if its package is absent):
#   - Expression : the repo decile expression matcher (functions.R: expressionMatched()), expr only
#   - MatchIt    : 1:1 nearest-neighbour Mahalanobis matching (functions.R: matchitMatched()),
#                  on expr+length and on expr+length+GC
#   - matchRanges: nullranges propensity-rejection matching (only if 'nullranges'
#                  is installed), on expr+length and expr+length+GC
#
# Scenarios (focal-set rules):
#   S0_null        random focal                         -> nobody should fail
#   S1_lengthconf  short-biased focal                   -> trips expr-only (length)
#   S2_interaction short AND highly expressed (corner)  -> trips expr-only + propensity
#   S3_sparse      top ~3% expression (thin pool)       -> no common support; all strain
#   S4_hiddenGC    high-GC focal                        -> expr+len leaves GC; +GC fixes
#   S5_tight_conf  large short-biased focal (tight pool)-> propensity fails; NN strains
#
# Outputs (written to <output_dir>, default "benchmark_out"):
#   scorecard.txt / scorecard.png   : |SMD| of each method vs input, per covariate
#   <scenario>_QC.png / _balance.txt: matchQC distribution plots for the MatchIt
#                                     (expr+len) match of each scenario
#
# Usage:
#   Rscript benchmark_matchers.R [output_dir] [n_genes] [seed]
#
# Requires: R (base) + MatchIt; optionally the Bioconductor package 'nullranges'.

## --- locate this script and load the repo's matching + QC libraries ---
cArgs <- commandArgs()
fileArg <- which(startsWith(cArgs, "--file="))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))
# This script lives in tests/; the matcher/QC libraries it sources are in the repo
# root one level up.
repo_root <- dirname(dirname(runningScript))
here <- function(f) file.path(repo_root, f)
source(here("functions.R"))   # expressionMatched(), matchitMatched(), matchQC()

if (!requireNamespace("MatchIt", quietly = TRUE)) {
  stop("benchmark_matchers.R needs the 'MatchIt' package; install.packages(\"MatchIt\").",
       call. = FALSE)
}
have_mr <- requireNamespace("nullranges", quietly = TRUE) &&
           requireNamespace("GenomicRanges", quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a   # null-coalescing for error messages

## --- arguments ---
argsArg <- which(startsWith(cArgs, "--args"))
args <- if (length(argsArg)) cArgs[(argsArg + 1):length(cArgs)] else character(0)
out_dir <- if (length(args) >= 1) args[1] else "benchmark_out"
n_genes <- if (length(args) >= 2) as.integer(args[2]) else 12000L
seed    <- if (length(args) >= 3) as.integer(args[3]) else 42L
if (is.na(n_genes) || n_genes < 500) stop("n_genes must be an integer >= 500.", call. = FALSE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) stop("could not create output directory: ", out_dir, call. = FALSE)

## --- simulate a realistic gene population ---
# Lengths ~ lognormal with the median (~25 kb) and spread seen in real PRO-seq gene
# sets; expression (log2 counts/kb) decreases weakly with length; GC is centred ~0.45
# and is mildly higher for highly expressed / shorter genes (housekeeping-like).
simulate_population <- function(N, seed) {
  set.seed(seed)
  loglen <- pmin(pmax(rnorm(N, log10(25000), 0.42), 3), log10(2.2e6))
  expr   <- 4.8 - 0.30 * (loglen - log10(25000)) + rnorm(N, 0, 1.5)
  gc     <- 0.44 + 0.02 * scale(expr)[, 1] - 0.015 * scale(loglen)[, 1] + rnorm(N, 0, 0.04)
  gc     <- pmin(pmax(gc, 0.25), 0.70)
  data.frame(gene = sprintf("G%06d", seq_len(N)),
             expr = expr, loglen = loglen, gc = gc, len = round(10^loglen),
             stringsAsFactors = FALSE)
}
pop <- simulate_population(n_genes, seed)
rownames(pop) <- pop$gene
pct <- function(x) rank(x) / length(x)

## --- scenario focal-set definitions (pool = all genes not in the focal set) ---
scenarios <- list(
  S0_null        = function(p) sample(p$gene, 400),
  S1_lengthconf  = function(p) sample(p$gene, 400, prob = (1 - pct(p$loglen))^3),
  S2_interaction = function(p) { c <- p$gene[pct(p$expr) > 0.70 & pct(p$loglen) < 0.30]
                                 sample(c, min(300, length(c))) },
  S3_sparse      = function(p) p$gene[pct(p$expr) > 0.97],
  S4_hiddenGC    = function(p) sample(p$gene, 400, prob = pct(p$gc)^3),
  S5_tight_conf  = function(p) sample(p$gene, round(0.40 * nrow(p)), prob = (1 - pct(p$loglen))^2)
)

## --- matchers: each returns matched control IDs, or an empty vector on failure ---
match_expression <- function(focal, pool) {
  ed <- data.frame(GeneID = c(focal, pool), Expression = pop[c(focal, pool), "expr"])
  tryCatch(unique(stats::na.omit(expressionMatched(focal, pool, ed, 10, setSeedNum = 1))),
           error = function(e) structure(character(0), err = conditionMessage(e)))
}
match_matchit <- function(focal, pool, covars) {
  tryCatch(matchitMatched(focal, pool, pop[, c("gene", covars)], covars = covars),
           error = function(e) structure(character(0), err = conditionMessage(e)))
}
match_matchranges <- function(focal, pool, covars) {
  if (!have_mr) return(structure(character(0), err = "nullranges not installed"))
  mk <- function(ids) { m <- pop[ids, ]
    g <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, width = pmax(1, m$len)), gene = ids)
    for (cv in covars) GenomicRanges::mcols(g)[[cv]] <- m[[cv]]; g }
  f <- stats::as.formula(paste("~", paste(covars, collapse = "+")))
  tryCatch(unique(nullranges::matchRanges(mk(focal), mk(pool), covar = f,
                                          method = "rejection", replace = FALSE)$gene),
           error = function(e) structure(character(0), err = conditionMessage(e)))
}

methods <- list(
  `Expression (expr)`  = function(f, p) match_expression(f, p),
  `MatchIt (expr+len)` = function(f, p) match_matchit(f, p, c("expr", "loglen")),
  `MatchIt (+GC)`      = function(f, p) match_matchit(f, p, c("expr", "loglen", "gc")))
if (have_mr) {
  methods[["matchRanges (expr+len)"]] <- function(f, p) match_matchranges(f, p, c("expr", "loglen"))
  methods[["matchRanges (+GC)"]]      <- function(f, p) match_matchranges(f, p, c("expr", "loglen", "gc"))
} else {
  message("NOTE: 'nullranges' not installed; skipping matchRanges methods ",
          "(install via BiocManager::install(\"nullranges\") to include them).")
}

## --- run every scenario x method; record covariate balance and run status ---
smd <- function(m, o) if (length(m) < 1 || length(o) < 2 || sd(o) == 0) NA_real_ else (mean(m) - mean(o)) / sd(o)
rows <- list()
for (sn in names(scenarios)) {
  focal <- intersect(scenarios[[sn]](pop), pop$gene)
  pool  <- setdiff(pop$gene, focal)
  oE <- pop[focal, "expr"]; oL <- pop[focal, "loglen"]; oG <- pop[focal, "gc"]
  message(sprintf("\n=== %s === focal=%d pool=%d (ratio %.1fx)", sn, length(focal),
                  length(pool), length(pool) / length(focal)))
  for (mn in names(methods)) {
    res <- methods[[mn]](focal, pool)
    if (length(res) < 1) {
      message(sprintf("  %-24s FAILED: %s", mn, attr(res, "err") %||% "empty"))
      rows[[length(rows) + 1]] <- data.frame(scenario = sn, method = mn,
        SMD_expr = NA, SMD_len = NA, SMD_gc = NA, status = "FAIL", stringsAsFactors = FALSE)
    } else {
      se <- smd(pop[res, "expr"], oE); sl <- smd(pop[res, "loglen"], oL); sg <- smd(pop[res, "gc"], oG)
      message(sprintf("  %-24s expr=%+.3f len=%+.3f gc=%+.3f (n=%d)", mn, se, sl, sg, length(res)))
      rows[[length(rows) + 1]] <- data.frame(scenario = sn, method = mn,
        SMD_expr = se, SMD_len = sl, SMD_gc = sg, status = "ok", stringsAsFactors = FALSE)
    }
  }
  # Per-scenario distribution QC of the MatchIt (expr+len) match, checked on ALL
  # three covariates so a left-over confounder (e.g. GC) is visible.
  mi <- match_matchit(focal, pool, c("expr", "loglen"))
  if (length(mi) >= 2) {
    matchQC(focal, mi, pop, covars = c("expr", "loglen", "gc"), all_genes = pop$gene,
            out_prefix = file.path(out_dir, sn))
  }
}
res <- do.call(rbind, rows)
write.table(res, file.path(out_dir, "scorecard.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## --- scorecard figure: one heatmap per covariate (scenario x method) ---
colf <- function(x) if (is.na(x)) "grey75" else if (abs(x) <= 0.1) "#1a9850" else if (abs(x) <= 0.25) "#fdae61" else "#d73027"
scen <- names(scenarios); meth <- names(methods)
if (capabilities("png")) {
  png(file.path(out_dir, "scorecard.png"),
      width = 480 * 3, height = 90 + 70 * length(scen), res = 120)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(1, 3), mar = c(9, 9, 3, 1))
  covlabs <- c(SMD_expr = "expression SMD", SMD_len = "length SMD", SMD_gc = "GC SMD")
  for (cn in names(covlabs)) {
    plot(NA, xlim = c(0.5, length(meth) + 0.5), ylim = c(0.5, length(scen) + 0.5),
         axes = FALSE, xlab = "", ylab = "", main = covlabs[[cn]])
    for (i in seq_along(scen)) for (j in seq_along(meth)) {
      r <- res[res$scenario == scen[i] & res$method == meth[j], ]
      y <- length(scen) - i + 1
      fail <- nrow(r) && r$status == "FAIL"
      val <- if (nrow(r)) r[[cn]] else NA
      rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = if (fail) "grey30" else colf(val), border = "white")
      text(j, y, if (fail) "FAIL" else if (is.na(val)) "-" else sprintf("%+.2f", val),
           cex = 0.8, col = "white", font = 2)
    }
    axis(2, at = length(scen):1, labels = scen, las = 1, tick = FALSE, cex.axis = 0.8)
    axis(1, at = seq_along(meth), labels = meth, las = 2, tick = FALSE, cex.axis = 0.72)
  }
}
message(sprintf("\nBenchmark complete. Scorecard and per-scenario QC written to '%s/'.", out_dir))
