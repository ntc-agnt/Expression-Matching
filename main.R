# probably only works for PRO-seq DE at the moment

cArgs <- commandArgs()
fileArg <- which(startsWith(cArgs, "--file="))
argsArg <- which(startsWith(cArgs, "--args"))
runningScript <- normalizePath(sapply(strsplit(cArgs[fileArg], "--file="), "[", 2))

source(sub("main.R", "expressionMatched.R", runningScript))
source(sub("main.R", "qc.R", runningScript))

p_thresh <- 0.1
fc_thresh <- 1.3

genes_of_interest <- scan(cArgs[argsArg+1], character()) # read in 1-column file of genes of interest
deseqDir <- cArgs[argsArg+2] # path to DESeq2 output dir
wf_path <- cArgs[argsArg+3] # path to GGA WORKING_FILE
subtract_TSS_TES <- as.numeric(cArgs[argsArg+4]) # bp to substract from TSS-TES length to match deseq window (e.g. 250 for +250-TES window)
deseq_manifest <- cArgs[argsArg+5] # path to deseq manifest
ref_condition <- cArgs[argsArg+6] # condition in deseq manifest to use for expression matching

# read in DESeq2 results tables
results_tables <- Sys.glob(paste(deseqDir, "*", "*_results_table.txt", sep = "/"))
results_list <- vector("list", length(results_tables))
for (i in 1:length(results_tables)){
  results_list[[i]] <- read.table(results_tables[i], header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
}
names(results_list) <- sapply(strsplit(sapply(strsplit(results_tables, "/"), tail, 1), "_results_table.txt"), "[", 1)

# We usually pick genes that are less than 1.3-fold changed and with p-values above 0.1
unchanged_all <- lapply(results_list, function(x) rownames(x)[!is.na(x[, 6]) & !is.na(x[, 7]) & x[, 6] > p_thresh & (x[, 7] <= fc_thresh & x[, 7] >= fc_thresh^-1)]) # x[, 6]=padj; x[, 7]=foldchange
unchanged_intersect <- Reduce(intersect, unchanged_all)

# use GGA WORKING_FILE to filter unchanged genes
gga_wf <- read.table(wf_path, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
gga_unchanged <- gga_wf[unchanged_intersect, ]
gga_unchanged <- gga_unchanged[gga_unchanged$GeneType %in% "protein_coding", ] # only keep protein-coding genes
gga_unchanged <- gga_unchanged[(abs(gga_unchanged$DominantTES - gga_unchanged$DominantTSS) + 1) >= 1000, ] # only keep genes at least 1 kb
unchanged_genes <- rownames(gga_unchanged)
gene_ids <- c(genes_of_interest, unchanged_genes)
gga_wf <- gga_wf[gene_ids, ]

# use length-normalized ref_condition counts as expression values
mh_rds <- Sys.glob(paste(deseqDir, "matrix", "*.makeHeatmap.rds", sep = "/"))
mh_counts <- readRDS(mh_rds)
deseq_meta <- read.table(deseq_manifest, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE)
mh_counts <- mh_counts[gene_ids, make.names(sapply(strsplit(rownames(deseq_meta)[deseq_meta$condition %in% ref_condition], "/"), tail, 1), allow_ = FALSE)]
mh_counts <- rowMeans(mh_counts) / ((abs(gga_wf$DominantTES - gga_wf$DominantTSS) + 1 - subtract_TSS_TES) / 1000) # counts/kb
exp_dat <- data.frame(GeneID=gene_ids, Expression=mh_counts)

# run expressionMatched() with 3 different seeds, writing the matched list and
# a QC report + violin plot for each
for (i in 1:3){
  exp_matched <- expressionMatched(genes_of_interest, unchanged_genes, exp_dat, 10, setSeedNum=i)
  out_prefix <- paste0("expression_matched-", ref_condition, "_", i)
  write(exp_matched, paste0(out_prefix, ".txt"))
  # QC: compare target / control / expression-matched expression distributions
  expressionMatchQC(genes_of_interest, unchanged_genes, exp_matched, exp_dat,
                    out_prefix = out_prefix,
                    plot_title = paste0(ref_condition, " (seed ", i, ")"))
}
