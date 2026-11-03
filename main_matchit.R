# main_matchit.R
#
# PROTOTYPE alternative driver to main.R. It selects the unchanged-gene control
# set by COVARIATE MATCHING on expression AND gene length (MatchIt nearest-
# neighbour), instead of decile expression-only matching, and writes the control
# gene list plus a QC report (covariate-balance table + distribution plots).
#
# Rationale: benchmarking on real NTC PRO-seq data showed decile expression
# matching leaves gene length badly imbalanced for most contrasts (a real
# confounder for PRO-seq), while joint nearest-neighbour matching balances
# expression and length together and is robust even for small gene sets. See the
# method-comparison report / repo skill for the full comparison.
#
# Usage (same 6 positional arguments as main.R):
#   Rscript main_matchit.R genes_ensg.txt DESeqOutputDir GGA.WORKING_FILE.txt \
#     <window_bp> DESeq_manifest.txt <baseline_condition>
#
# Requires: R (base) plus the CRAN package 'MatchIt'.

p_thresh  <- 0.1    # padj cutoff: genes with padj > p_thresh are "not significant"
fc_thresh <- 1.3    # linear fold-change band [1/fc_thresh, fc_thresh] = "unchanged"

# --- locate this script and load the matching + QC libraries (name-swap, as main.R) ---
cArgs <- commandArgs()
fileArg <- which(startsWith(cArgs, "--file="))
argsArg <- which(startsWith(cArgs, "--args"))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))
source(sub("main_matchit.R", "matchitMatched.R", runningScript))
source(sub("main_matchit.R", "matchQC.R", runningScript))

# --- parse and validate arguments (fail early with a clear message) ---
args <- if (length(argsArg)) cArgs[(argsArg + 1):length(cArgs)] else character(0)
usage <- paste0(
  "Usage: Rscript main_matchit.R genes.txt DESeqOutputDir GGA.WORKING_FILE.txt ",
  "<window_bp> DESeq_manifest.txt <baseline_condition>")
if (length(args) != 6) stop(usage, call. = FALSE)
genes_file       <- args[1]
deseqDir         <- args[2]
wf_path          <- args[3]
subtract_TSS_TES <- suppressWarnings(as.numeric(args[4]))  # bp subtracted from TSS-TES length
deseq_manifest   <- args[5]
ref_condition    <- args[6]

for (f in c(genes_file, wf_path, deseq_manifest)) {
  if (!file.exists(f)) stop("input file not found: ", f, call. = FALSE)
}
if (!dir.exists(deseqDir)) stop("DESeq2 output directory not found: ", deseqDir, call. = FALSE)
if (is.na(subtract_TSS_TES)) {
  stop("window_bp (argument 4) must be numeric; got '", args[4], "'.", call. = FALSE)
}

genes_of_interest <- scan(genes_file, character(), quiet = TRUE)
if (!length(genes_of_interest)) stop("no gene IDs read from ", genes_file, call. = FALSE)

# --- read DESeq2 results tables and define "unchanged" genes (same rule as main.R) ---
# DESeq column positions are hardcoded to the NTC results-table layout:
# column 6 = padj, column 7 = LINEAR foldchange. A gene is "unchanged" when it is
# not significant (padj > p_thresh) and within the fold-change band in EVERY
# contrast (intersection across all results tables).
results_tables <- Sys.glob(paste(deseqDir, "*", "*_results_table.txt", sep = "/"))
if (!length(results_tables)) {
  stop("no */*_results_table.txt files found under ", deseqDir, call. = FALSE)
}
results_list <- lapply(results_tables, read.table, header = TRUE, sep = "\t",
                       row.names = 1, stringsAsFactors = FALSE)
unchanged_all <- lapply(results_list, function(x)
  rownames(x)[!is.na(x[, 6]) & !is.na(x[, 7]) & x[, 6] > p_thresh &
              x[, 7] <= fc_thresh & x[, 7] >= fc_thresh^-1])
unchanged_intersect <- Reduce(intersect, unchanged_all)

# --- GGA filter: protein-coding, >= 1 kb (full TSS-TES length) ---
gga_wf <- read.table(wf_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
for (col in c("GeneType", "DominantTSS", "DominantTES")) {
  if (!col %in% colnames(gga_wf)) stop("GGA WORKING_FILE missing column: ", col, call. = FALSE)
}
gga_unchanged <- gga_wf[unchanged_intersect, ]
gga_unchanged <- gga_unchanged[gga_unchanged$GeneType %in% "protein_coding", ]
gga_unchanged <- gga_unchanged[(abs(gga_unchanged$DominantTES - gga_unchanged$DominantTSS) + 1) >= 1000, ]
unchanged_genes <- rownames(gga_unchanged)
if (!length(unchanged_genes)) stop("no unchanged control genes passed the filters.", call. = FALSE)

# --- length-normalized expression (counts/kb) over the baseline condition, for ALL genes ---
mh_rds <- Sys.glob(paste(deseqDir, "matrix", "*.makeHeatmap.rds", sep = "/"))
if (length(mh_rds) != 1) {
  stop("expected exactly one matrix/*.makeHeatmap.rds under ", deseqDir,
       "; found ", length(mh_rds), ".", call. = FALSE)
}
mh_counts <- readRDS(mh_rds)
deseq_meta <- read.table(deseq_manifest, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
if (!"condition" %in% colnames(deseq_meta)) stop("manifest needs a 'condition' column.", call. = FALSE)
ref_samples <- rownames(deseq_meta)[deseq_meta$condition %in% ref_condition]
if (!length(ref_samples)) {
  stop("baseline condition '", ref_condition,
       "' not found in the manifest 'condition' column.", call. = FALSE)
}
ref_cols <- make.names(sapply(strsplit(ref_samples, "/"), tail, 1), allow_ = FALSE)
missing_cols <- setdiff(ref_cols, colnames(mh_counts))
if (length(missing_cols)) {
  stop("baseline sample column(s) not found in makeHeatmap matrix: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

# Build covariates for every gene present in both the matrix and the annotation.
common   <- intersect(rownames(mh_counts), rownames(gga_wf))
gga_all  <- gga_wf[common, ]
gene_len <- abs(gga_all$DominantTES - gga_all$DominantTSS) + 1               # full length (bp)
# Expression model (same as main.R): mean baseline counts / effective length in
# kb. The window is subtracted only here, not in the >= 1 kb filter above.
counts_kb <- rowMeans(mh_counts[common, ref_cols, drop = FALSE]) /
             ((gene_len - subtract_TSS_TES) / 1000)
# log transforms put expression and length on comparable, roughly symmetric
# scales so the Mahalanobis distance is not dominated by either's raw spread.
dat <- data.frame(GeneID = common,
                  log2Expression = log2(counts_kb + 1),
                  log10Length    = log10(gene_len),
                  row.names = common, stringsAsFactors = FALSE)
covars <- c("log2Expression", "log10Length")

# --- match on expression + length and write the control list ---
matched <- matchitMatched(genes_of_interest, unchanged_genes, dat, covars = covars)
out_prefix <- paste0("matchit_matched-", ref_condition)
writeLines(matched, paste0(out_prefix, ".txt"))
message(sprintf("Wrote %d matched control genes -> %s.txt", length(matched), out_prefix))

# --- QC: balance table + distribution plots (input vs matched vs all genes) ---
matchQC(genes_of_interest, matched, dat, covars = covars,
        all_genes = common, out_prefix = out_prefix)
