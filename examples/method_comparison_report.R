#!/usr/bin/env Rscript
# method_comparison_report.R
# -----------------------------------------------------------------------------
# Worked example: compare three control-gene matchers on a REAL DESeq2 dataset
# and emit a multi-page PDF report plus a summary table.
#
# The three matchers (all from the repo's functions.R, all drawing controls from
# the SAME unchanged-gene pool) are:
#   * expressionMatched (expr)   - legacy decile matcher, expression only
#   * MatchIt (expr)             - nearest-neighbour Mahalanobis, expression only
#   * MatchIt (expr+len)         - nearest-neighbour Mahalanobis, expression + gene length
#
# For every contrast and both directions (Up / Dn) it reports how well each
# matcher balances expression and gene length (standardized mean difference,
# |SMD| <= 0.1 = balanced), and it draws the joint expression x length "mechanism"
# view that shows what MatchIt actually does (which pool gene each input gene is
# paired with). See examples/README.md for full usage.
#
# Inputs (flags):
#   --dds  FILE   a DESeq2 DESeqDataSet saved with saveRDS() (must have a
#                 'condition' factor in colData and run through DESeq()).
#   --gga  FILE   the GGA Dominant.TSS.TES WORKING_FILE (tab-delimited, gene IDs
#                 in column 1) providing GeneType, GeneName, DominantTSS/TES.
#   --contrasts   comma-separated num/den pairs over the 'condition' levels,
#                 e.g. "tam/ctrl_003,siTTF2/ctrl_002". Baseline (den) expression
#                 is used as the matching covariate.
#   --out  DIR    output directory (default: matching_report_output).
#   --mech        which contrast drives the mechanism pages (default: last one).
# Thresholds (flags, with defaults matching the repo's conventions):
#   --up-fc 1.5  --padj 0.01            (differential call: |linear FC| and padj)
#   --unch-fc 1.3  --unch-padj 0.1      (unchanged pool: linear FC band and padj)
#   --min-len 1000  --genetype protein_coding  --deciles 10
# -----------------------------------------------------------------------------

suppressMessages({ library(DESeq2) })

cArgs <- commandArgs()
runningScript <- normalizePath(sapply(strsplit(cArgs[which(startsWith(cArgs, "--file="))], "--file="), "[", 2))
# This script lives in examples/; functions.R is at the repo root (one level up).
repo_root <- dirname(dirname(runningScript))
source(file.path(repo_root, "functions.R"))   # expressionMatched(), matchitMatched(), matchQC()
`%||%` <- function(a, b) if (is.null(a)) b else a

## --- flag parsing ---
usage <- paste(
  "Usage: Rscript examples/method_comparison_report.R --dds FILE --gga FILE",
  "--contrasts num/den[,num/den...] [--out DIR] [--mech num/den]",
  "[--up-fc 1.5] [--padj 0.01] [--unch-fc 1.3] [--unch-padj 0.1]",
  "[--min-len 1000] [--genetype protein_coding] [--deciles 10]")
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

required <- c("dds", "gga", "contrasts")
missing_req <- setdiff(required, names(opt))
if (length(missing_req)) stop("missing required argument(s): --",
                              paste(missing_req, collapse = ", --"), "\n", usage, call. = FALSE)
for (f in c(opt$dds, opt$gga)) if (!file.exists(f)) stop("file not found: ", f, call. = FALSE)

out_dir     <- opt[["out"]] %||% "matching_report_output"
up_fc       <- as.numeric(opt[["up-fc"]]     %||% "1.5")
padj_de     <- as.numeric(opt[["padj"]]      %||% "0.01")
unch_fc     <- as.numeric(opt[["unch-fc"]]   %||% "1.3")
unch_padj   <- as.numeric(opt[["unch-padj"]] %||% "0.1")
min_len     <- as.numeric(opt[["min-len"]]   %||% "1000")
gene_type   <- opt[["genetype"]] %||% "protein_coding"
deciles     <- as.integer(opt[["deciles"]]   %||% "10")
lfc_up      <- log2(up_fc)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

parse_contrasts <- function(s) {
  parts <- strsplit(strsplit(s, ",")[[1]], "/")
  lapply(parts, function(p) {
    if (length(p) != 2) stop("each contrast must be num/den; got '", paste(p, collapse="/"), "'.", call. = FALSE)
    trimws(p)
  })
}
cands  <- parse_contrasts(opt$contrasts)
mech_c <- if (!is.null(opt[["mech"]])) trimws(strsplit(opt[["mech"]], "/")[[1]]) else cands[[length(cands)]]

## --- load data ---
dds <- readRDS(opt$dds)
if (!"condition" %in% colnames(as.data.frame(colData(dds))))
  stop("the dds colData has no 'condition' column.", call. = FALSE)
levs <- levels(dds$condition)
for (ct in c(cands, list(mech_c))) if (!all(ct %in% levs))
  stop("condition level(s) not found: ", paste(setdiff(ct, levs), collapse=", "),
       "\navailable: ", paste(levs, collapse=", "), call. = FALSE)
gga <- read.table(opt$gga, header = TRUE, sep = "\t", row.names = 1, quote = "",
                  stringsAsFactors = FALSE)
mh  <- counts(dds, normalized = TRUE)
len <- abs(gga$DominantTES - gga$DominantTSS) + 1; names(len) <- rownames(gga)
common <- intersect(rownames(mh), rownames(gga))
smd <- function(m, o) if (length(m) < 1 || length(o) < 2 || sd(o) == 0) NA_real_ else (mean(m) - mean(o)) / sd(o)

## For a contrast, return the Up/Dn input lists, the shared unchanged pool, and a
## per-gene covariate frame (expression = log2 counts/kb at baseline; length).
prep_contrast <- function(num, den) {
  res <- as.data.frame(results(dds, contrast = c("condition", num, den)))
  bs  <- colnames(dds)[dds$condition == den]
  ck  <- rowMeans(mh[common, bs, drop = FALSE]) / (len[common] / 1000)          # counts/kb (baseline)
  gm  <- data.frame(gene = common, counts_kb = ck, expr = log2(ck + 1),
                    loglen = log10(len[common]), row.names = common, stringsAsFactors = FALSE)
  gm  <- gm[is.finite(gm$expr) & is.finite(gm$loglen), ]
  elig <- rownames(gm)[gga[rownames(gm), "GeneType"] %in% gene_type & len[rownames(gm)] >= min_len]
  up_  <- intersect(rownames(res)[!is.na(res$padj) & res$padj <= padj_de & res$log2FoldChange >=  lfc_up], elig)
  dn_  <- intersect(rownames(res)[!is.na(res$padj) & res$padj <= padj_de & res$log2FoldChange <= -lfc_up], elig)
  lfc_lin <- 2^res$log2FoldChange
  unch <- rownames(res)[!is.na(res$padj) & !is.na(lfc_lin) & res$padj > unch_padj &
                        lfc_lin <= unch_fc & lfc_lin >= 1/unch_fc]
  pool <- setdiff(intersect(unch, elig), c(up_, dn_))                            # unchanged, minus either tail
  list(gm = gm, up = up_, dn = dn_, pool = pool, num = num, den = den)
}

## Three matched sets for an arbitrary input set (up or dn), all from the SAME pool.
run_methods <- function(P, input) {
  gm <- P$gm; pool <- P$pool
  ed  <- data.frame(GeneID = c(input, pool), Expression = gm[c(input, pool), "counts_kb"])  # raw counts/kb
  em  <- tryCatch(unique(na.omit(expressionMatched(input, pool, ed, deciles, setSeedNum = 1))),
                  error = function(e) character(0))
  mie <- tryCatch(matchitMatched(input, pool, gm[, c("gene", "expr")], covars = "expr"),
                  error = function(e) character(0))
  mil <- tryCatch(matchitMatched(input, pool, gm[, c("gene", "expr", "loglen")], covars = c("expr", "loglen")),
                  error = function(e) character(0))
  list(`expressionMatched (expr)` = em, `MatchIt (expr)` = mie, `MatchIt (expr+len)` = mil)
}

## Return the actual 1:1 pairs MatchIt forms (which pool gene each input gene is matched to),
## replicating the repo's matchitMatched() setup so the pairs are faithful to the report.
matchit_pairs <- function(gm, input, pool, covars) {
  use  <- c(input, pool)
  mf   <- data.frame(treat = as.integer(use %in% input),
                     gm[use, covars, drop = FALSE], row.names = use, check.names = FALSE)
  mf   <- mf[apply(mf[, covars, drop = FALSE], 1, function(r) all(is.finite(r))), , drop = FALSE]
  set.seed(1)
  f    <- stats::as.formula(paste("treat ~", paste(covars, collapse = " + ")))
  fit  <- MatchIt::matchit(f, data = mf, method = "nearest", distance = "mahalanobis",
                           replace = FALSE, ratio = 1)
  mm   <- fit$match.matrix
  data.frame(input_gene = rownames(mm), ctrl = mm[, 1], row.names = NULL, stringsAsFactors = FALSE)
}

## ---- run all contrasts x {Up, Dn} ----
preps <- list()
for (ct in cands) preps[[paste(ct[1], "vs", ct[2])]] <- prep_contrast(ct[1], ct[2])
summary_tab <- data.frame(); figdata <- list()
dir_defs <- c(Up = "up", Dn = "dn")   # Up: linear FC >= up_fc; Dn: linear FC <= 1/up_fc; both padj <= padj_de
for (nm in names(preps)) {
  P <- preps[[nm]]
  for (dir in names(dir_defs)) {
    input <- P[[dir_defs[[dir]]]]
    sets  <- run_methods(P, input)
    for (mth in names(sets)) {
      ids <- sets[[mth]]
      summary_tab <- rbind(summary_tab, data.frame(
        Contrast = nm, Dir = dir, Method = mth, n_input = length(input), n_matched = length(ids),
        SMD_expr = round(smd(P$gm[ids,"expr"], P$gm[input,"expr"]), 3),
        SMD_length = round(smd(P$gm[ids,"loglen"], P$gm[input,"loglen"]), 3)))
    }
    figdata[[sprintf("%s [%s]", nm, dir)]] <- list(P = P, input = input, sets = sets, nm = nm, dir = dir)
  }
}
summary_tab$length_balanced <- with(summary_tab,
  ifelse(n_input < 10,           "n<10",     # too few input genes for a meaningful SMD
  ifelse(n_matched == 0,         "FAILED",   # matcher could not return a set (see note)
  ifelse(is.na(SMD_length),      "-",
  ifelse(abs(SMD_length) <= 0.1, "yes", "NO")))))
csv_path <- file.path(out_dir, "summary_table.csv")
write.csv(summary_tab, csv_path, row.names = FALSE)
cat("\n=== Summary table ===\n"); print(summary_tab, row.names = FALSE)

## ---- mechanism worked example (the --mech contrast) ----
mech_nm <- paste(mech_c[1], "vs", mech_c[2])
mechP <- preps[[mech_nm]] %||% prep_contrast(mech_c[1], mech_c[2])
g <- mechP$gm
mech_in <- mechP$dn                                          # Dn: often the length-skewed, decile-starved case
pe <- matchit_pairs(g, mech_in, mechP$pool, "expr")
pb <- matchit_pairs(g, mech_in, mechP$pool, c("expr", "loglen"))
mech <- merge(pe, pb, by = "input_gene", suffixes = c("_E", "_EL"))
mech$in_expr  <- g[mech$input_gene,"expr"]; mech$in_loglen  <- g[mech$input_gene,"loglen"]
mech$cE_expr  <- g[mech$ctrl_E,    "expr"]; mech$cE_loglen  <- g[mech$ctrl_E,    "loglen"]
mech$cEL_expr <- g[mech$ctrl_EL,   "expr"]; mech$cEL_loglen <- g[mech$ctrl_EL,   "loglen"]
mech$lenerrE  <- abs(mech$cE_loglen  - mech$in_loglen)
mech$lenerrEL <- abs(mech$cEL_loglen - mech$in_loglen)
sym <- function(ids) { s <- gga[ids, "GeneName"]; s[is.na(s) | s == ""] <- ids[is.na(s) | s == ""]; s }
md   <- mech[mech$ctrl_E != mech$ctrl_EL, ]
rep4 <- head(md[order(md$lenerrE - md$lenerrEL, decreasing = TRUE), ], 4)
med_fold <- if (nrow(mech)) c(E = 10^median(mech$lenerrE), EL = 10^median(mech$lenerrEL)) else c(E = NA, EL = NA)
# all-three view on the Up direction of the mech contrast
up_in <- mechP$up
edU   <- data.frame(GeneID = c(up_in, mechP$pool), Expression = g[c(up_in, mechP$pool), "counts_kb"])
emU   <- intersect(unique(na.omit(tryCatch(expressionMatched(up_in, mechP$pool, edU, deciles, setSeedNum = 1),
                                           error = function(e) character(0)))), rownames(g))
peU   <- matchit_pairs(g, up_in, mechP$pool, "expr")
pbU   <- matchit_pairs(g, up_in, mechP$pool, c("expr", "loglen"))

## ---- PDF report ----
cols <- c(all="grey70", input="black", `expressionMatched (expr)`="dodgerblue3",
          `MatchIt (expr)`="darkorange2", `MatchIt (expr+len)`="red2")
ecdf_panel <- function(gm, input, sets, col, lab) {
  plot(ecdf(gm[, col]), col = cols["all"], lwd = 3, do.points = FALSE, verticals = TRUE,
       main = lab, xlab = lab, ylab = "cumulative fraction")
  if (length(input) >= 1)
    plot(ecdf(gm[input, col]), add = TRUE, col = cols["input"], lwd = 2.6, do.points = FALSE, verticals = TRUE)
  for (mth in names(sets)) if (length(sets[[mth]]) > 1)
    plot(ecdf(gm[sets[[mth]], col]), add = TRUE, col = cols[mth], lwd = 1.8,
         lty = if (mth == "MatchIt (expr+len)") 1 else 2, do.points = FALSE, verticals = TRUE)
}
pdf_path <- file.path(out_dir, "matching_report.pdf")
pdf(pdf_path, width = 11, height = 8.5)

# page 1: title + summary table
par(mar = c(1,1,1,1)); plot.new()
text(0.5, 0.96, "Expression-matched control sets: method comparison on real data", font = 2, cex = 1.3)
text(0.5, 0.91, sprintf("Up (linear FC >= %.2f) and Dn (linear FC <= 1/%.2f) sets, both padj <= %g", up_fc, up_fc, padj_de), cex = 0.9)
text(0.02, 0.85, "Summary", adj = 0, font = 2, cex = 1.1)
tt <- summary_tab
hdr <- c("Contrast","Dir","Method","n_input","n_matched","SMD_expr","SMD_length","len_bal")
x   <- c(0.02, 0.26, 0.34, 0.56, 0.66, 0.75, 0.84, 0.93)
y <- 0.825
text(x, y, hdr, adj = 0, font = 2, cex = 0.60, family = "mono")
row_dy <- min(0.0275, (0.70) / max(nrow(tt), 1))
for (i in seq_len(nrow(tt))) {
  y <- y - row_dy
  if (i > 1 && tt$Contrast[i] != tt$Contrast[i-1]) y <- y - row_dy * 0.3
  row <- c(tt$Contrast[i], tt$Dir[i], tt$Method[i], tt$n_input[i], tt$n_matched[i],
           sprintf("%+.3f", tt$SMD_expr[i]), sprintf("%+.3f", tt$SMD_length[i]), tt$length_balanced[i])
  colr <- if (tt$length_balanced[i] %in% c("NO","FAILED")) "red3" else
          if (tt$length_balanced[i] == "n<10") "darkorange3" else "black"
  text(x, rep(y, length(x)), row, adj = 0, cex = 0.58, family = "mono", col = colr)
}
text(0.02, 0.095, "len_bal: black/green = balanced (|SMD|<=0.1); red = imbalanced (NO) or matcher returned nothing (FAILED); orange = n<10 (unreliable).", adj = 0, cex = 0.68)
text(0.02, 0.065, "Expression is matched to ~0 by all three methods; only MatchIt (expr+len) also balances gene length.", adj = 0, cex = 0.68)
text(0.02, 0.035, "FAILED: expressionMatched's 1:1 decile draw (no replacement) runs out of unchanged pool genes in the input's dense deciles; MatchIt still succeeds.", adj = 0, cex = 0.68, font = 3)

# page 2: how MatchIt works (text) + representative pairs
par(mar = c(1,1,1,1)); plot.new()
text(0.5, 0.955, "How MatchIt builds a control set — and why it differs from decile expression-matching", font = 2, cex = 1.25)
yy <- 0.925
line <- function(txt, dy = 0.040, cex = 0.86, font = 1, col = "black", xx = 0.035, mono = FALSE) {
  yy <<- yy - dy; text(xx, yy, txt, adj = 0, cex = cex, font = font, col = col, family = if (mono) "mono" else "")
}
gap <- function(d = 0.016) yy <<- yy - d
gap(0.02)
line("1.  Same pool, different rule.  Both matchers draw controls ONLY from the unchanged-gene pool (padj > unch-padj and linear FC within [1/unch-fc, unch-fc]).", font = 2)
line("     Every control MatchIt returns is a member of that pool — it never samples outside the unchanged set, and never relaxes what 'expression-matched' means.")
gap()
line("2.  expressionMatched (legacy) — rigid bins.  It sorts genes into expression deciles, then draws an equal number of pool genes from EACH decile,", font = 2)
line("     one-to-one, without reuse.  If one decile holds more input genes than the pool can supply there, the draw is impossible and the whole run aborts (FAILED).")
gap()
line("3.  MatchIt — nearest available neighbour.  No bins.  For each input gene it takes the closest still-unused pool gene by distance.  'Closest available'", font = 2)
line("     always exists (as long as the pool is >= the input set), so it degrades gracefully — a slightly worse match — instead of failing outright.")
line("     With covars = expr this is the SAME single yardstick (expression) as the legacy matcher, which is why it succeeds where the decile draw fails.")
gap()
line("4.  A second yardstick, not a new definition.  With covars = expr + length each gene becomes a point in 2-D (expression, length), and 'closest' means", font = 2)
line("     closest in BOTH at once (Mahalanobis distance, each axis scaled by its spread).  Expression is still matched to SMD ~ 0; length just stops being")
line(sprintf("     ignored.  This works because the unchanged pool (here %d genes) is large enough to contain genes that match on expression AND length together.", length(mechP$pool)))
gap(0.022)
line(sprintf("Worked example — %s, %s (%d genes).  Representative pairs where the two rules chose different controls (length in bp):",
             mech_nm, "Dn", length(mech_in)), font = 2, cex = 0.82)
gap(0.004)
line(sprintf("%-22s %-28s %-28s", "input gene (len)", "expr-only control (len)", "expr+len control (len)"), cex = 0.72, mono = TRUE, font = 2)
if (nrow(rep4)) for (i in seq_len(nrow(rep4))) {
  r <- rep4[i, ]
  line(sprintf("%-10s %9s   ->  %-10s %9s      %-10s %9s",
               sym(r$input_gene), format(len[r$input_gene], big.mark=","),
               sym(r$ctrl_E),  format(len[r$ctrl_E],  big.mark=","),
               sym(r$ctrl_EL), format(len[r$ctrl_EL], big.mark=",")),
       cex = 0.72, mono = TRUE, dy = 0.032)
}
gap(0.006)
if (is.finite(med_fold["E"]))
  line(sprintf("Median per-pair length error:  expr-only = %.1fx off,   expr+len = %.2fx off.   Expression is matched equally well by both; only expr+len also pins length.",
               med_fold["E"], med_fold["EL"]), cex = 0.8, font = 3)
line("The next pages show these pairs as a joint expression x length map — the view that actually captures what MatchIt does.", cex = 0.8, font = 3, dy = 0.032)

# page 3: mechanism figure — joint expression x length map with matched-pair links (mech contrast, Dn)
if (nrow(mech) > 1) {
  layout(matrix(c(1,2,3,3), nrow = 2, byrow = TRUE), heights = c(2.1, 1))
  par(mar = c(4.2,4.2,3,1), oma = c(2.2,0,2.4,0))
  scatter_links <- function(prefix, color, title) {
    plot(g$expr, g$loglen, pch = 16, cex = 0.2, col = adjustcolor("grey70", 0.35),
         xlab = "log2 expression (counts/kb)", ylab = "log10 gene length (bp)", main = title, cex.main = 0.95)
    xc <- mech[[paste0("c", prefix, "_expr")]]; yc <- mech[[paste0("c", prefix, "_loglen")]]
    segments(mech$in_expr, mech$in_loglen, xc, yc, col = adjustcolor(color, 0.15), lwd = 0.5)
    points(xc, yc, pch = 16, cex = 0.5, col = color)
    points(mech$in_expr, mech$in_loglen, pch = 16, cex = 0.5, col = "black")
    legend("topright", bty = "n", cex = 0.8, pch = 16, col = c("black", color, "grey70"),
           legend = c("input gene", "its matched control", "all genes"))
  }
  scatter_links("E",  cols["MatchIt (expr)"],     "MatchIt (expr only): length left uncontrolled")
  scatter_links("EL", cols["MatchIt (expr+len)"], "MatchIt (expr+len): length matched")
  xm <- quantile(c(mech$lenerrE, mech$lenerrEL), 0.995)
  plot(ecdf(mech$lenerrE), col = cols["MatchIt (expr)"], lwd = 2.4, do.points = FALSE, verticals = TRUE,
       xlim = c(0, xm), main = sprintf("Per-pair length mismatch (all %d pairs)", nrow(mech)),
       xlab = "|log10(control length / gene length)|   (0 = identical length; 0.30 = 2x; 1.0 = 10x)",
       ylab = "cumulative fraction of pairs")
  plot(ecdf(mech$lenerrEL), add = TRUE, col = cols["MatchIt (expr+len)"], lwd = 2.4, do.points = FALSE, verticals = TRUE)
  abline(v = log10(2), lty = 3, col = "grey50")
  legend("bottomright", bty = "n", cex = 0.9, lwd = 2.4, col = cols[c("MatchIt (expr)","MatchIt (expr+len)")],
         legend = c(sprintf("expr only  (median %.1fx off)", med_fold["E"]),
                    sprintf("expr+len  (median %.2fx off)", med_fold["EL"])))
  mtext(sprintf("Mechanism (%s, Dn, %d pairs): each segment links an input gene (black) to the unchanged control MatchIt paired it with", mech_nm, nrow(mech)),
        outer = TRUE, font = 2, cex = 1.0)
  if (length(mechP$dn) > 0 && length(figdata[[sprintf("%s [Dn]", mech_nm)]]$sets[["expressionMatched (expr)"]]) == 0)
    mtext("expressionMatched produced no set for this input (decile draw FAILED), so it has nothing to project here — see the Up page for all three methods.",
          outer = TRUE, side = 1, line = 0.7, font = 3, cex = 0.8)
}

# page 4: same idiom, all three methods on the mech contrast Up direction
if (length(up_in) >= 10) {
  layout(matrix(c(1,2,3,4,4,4), nrow = 2, byrow = TRUE), heights = c(2, 1))
  par(mar = c(4.2,4.2,3,1), oma = c(0,0,2.4,0))
  bg_scatter <- function(title) plot(g$expr, g$loglen, pch = 16, cex = 0.2, col = adjustcolor("grey70", 0.35),
         xlab = "log2 expression (counts/kb)", ylab = "log10 gene length (bp)", main = title, cex.main = 0.9)
  mech_leg <- function(color, ctrl_lab) legend("topright", bty = "n", cex = 0.72, pch = 16,
         col = c("black", color, "grey70"), legend = c("up gene (input)", ctrl_lab, "all genes"))
  if (length(emU) > 0) {
    bg_scatter("expressionMatched (expr): set-level draw (no 1:1 pairs)")
    points(g[emU, "expr"], g[emU, "loglen"], pch = 16, cex = 0.5, col = cols["expressionMatched (expr)"])
    points(g[up_in, "expr"], g[up_in, "loglen"], pch = 16, cex = 0.5, col = "black")
    mech_leg(cols["expressionMatched (expr)"], "its control set")
  } else {
    plot.new(); text(0.5, 0.5, "expressionMatched FAILED\n(no set for this input)", cex = 1.1, font = 3)
  }
  up_links <- function(pairs, color, title) {
    bg_scatter(title)
    segments(g[pairs$input_gene,"expr"], g[pairs$input_gene,"loglen"],
             g[pairs$ctrl,"expr"], g[pairs$ctrl,"loglen"], col = adjustcolor(color, 0.15), lwd = 0.5)
    points(g[pairs$ctrl,"expr"], g[pairs$ctrl,"loglen"], pch = 16, cex = 0.5, col = color)
    points(g[up_in,"expr"], g[up_in,"loglen"], pch = 16, cex = 0.5, col = "black")
    mech_leg(color, "its matched control")
  }
  up_links(peU, cols["MatchIt (expr)"],     "MatchIt (expr only): length left uncontrolled")
  up_links(pbU, cols["MatchIt (expr+len)"], "MatchIt (expr+len): length matched")
  plot(ecdf(g[up_in,"loglen"]), col = cols["input"], lwd = 2.8, do.points = FALSE, verticals = TRUE,
       main = "Gene-length distribution: up-regulated input vs each control set (marginal)",
       xlab = "log10 gene length (bp)", ylab = "cumulative fraction")
  if (length(emU) > 0)
    plot(ecdf(g[emU,"loglen"]), add = TRUE, col = cols["expressionMatched (expr)"], lwd = 2, lty = 2, do.points = FALSE, verticals = TRUE)
  plot(ecdf(g[peU$ctrl,"loglen"]), add = TRUE, col = cols["MatchIt (expr)"],    lwd = 2, lty = 2, do.points = FALSE, verticals = TRUE)
  plot(ecdf(g[pbU$ctrl,"loglen"]), add = TRUE, col = cols["MatchIt (expr+len)"], lwd = 2, lty = 1, do.points = FALSE, verticals = TRUE)
  legend("bottomright", bty = "n", cex = 0.85, lwd = 2, lty = c(1,2,2,1),
         col = cols[c("input","expressionMatched (expr)","MatchIt (expr)","MatchIt (expr+len)")],
         legend = c("up input (target)","expressionMatched","MatchIt (expr)","MatchIt (expr+len)"))
  mtext(sprintf("All three matchers, same idiom (%s, Up, n = %d): only MatchIt(expr+len) reproduces the input gene-length curve; the expression-only methods miss it",
        mech_nm, length(up_in)), outer = TRUE, font = 2, cex = 0.82)
}

# pages 5+: per-contrast, per-direction marginal ECDF panels (expression and length)
for (key in names(figdata)) {
  fd <- figdata[[key]]
  if (length(fd$input) < 10) next
  par(mfrow = c(1,2), mar = c(4.4,4.4,4,1), oma = c(0,0,2,0))
  ecdf_panel(fd$P$gm, fd$input, fd$sets, "expr", "log2 expression (counts/kb)")
  legend("bottomright", bty="n", cex=0.75, lwd=2, lty=c(1,1,2,2,1),
         col=cols[c("all","input","expressionMatched (expr)","MatchIt (expr)","MatchIt (expr+len)")],
         legend=c("all genes","input set","expressionMatched","MatchIt (expr)","MatchIt (expr+len)"))
  ecdf_panel(fd$P$gm, fd$input, fd$sets, "loglen", "log10 gene length")
  mtext(sprintf("%s  —  %s-regulated  (n = %d genes)", fd$nm, fd$dir, length(fd$input)),
        outer = TRUE, font = 2, cex = 1.1)
}
invisible(dev.off())
cat(sprintf("\nWrote:\n  %s\n  %s\n", csv_path, pdf_path))
