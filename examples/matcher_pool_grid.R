#!/usr/bin/env Rscript
# matcher_pool_grid.R
# -----------------------------------------------------------------------------
# Full factorial (3 matchers x 2 unchanged-pool definitions), visualized natively.
# For each contrast it draws:
#   * a 3x2 grid of the joint expression x length SCATTER (MatchIt shown with its native pair-links):
#       rows = unchanged-pool definition  (TOP: legacy hard cutoff   BOTTOM: lessAbs equivalence)
#       cols = matcher  (expMat, MatchIt:expr, MatchIt:expr+len)
#   * a 3x2 grid of gene-length ECDFs in the same layout
# plus one overall summary table (length SMD per cell).
#
# expMat = expressionMatched (decile, expression-only) -> point cloud (no 1:1 pairs)
# MatchIt:expr / MatchIt:expr+len -> nearest-neighbour Mahalanobis, drawn with pair-links.
# Expression is matched to ~0 by all three matchers in every cell; the summary table is the LENGTH axis.
#
# Flags:  --dds FILE  --gga FILE  --contrasts num/den[,num/den...]  [--out DIR]  [--dir up|dn]
#         [--up-fc 1.5] [--padj 0.01] [--unch-fc 1.3] [--unch-padj 0.1] [--lessabs-alpha 0.1]
#         [--min-len 1000] [--genetype protein_coding] [--deciles 10] [--bg 6000]
# See examples/README.md for details.
# -----------------------------------------------------------------------------
suppressMessages({ library(DESeq2) })
cArgs <- commandArgs()
runningScript <- normalizePath(sapply(strsplit(cArgs[which(startsWith(cArgs, "--file="))], "--file="), "[", 2))
source(file.path(dirname(dirname(runningScript)), "functions.R"))   # expressionMatched(), matchitMatched()
`%||%` <- function(a, b) if (is.null(a)) b else a

usage <- paste("Usage: Rscript examples/matcher_pool_grid.R --dds FILE --gga FILE",
  "--contrasts num/den[,num/den...] [--out DIR] [--dir up|dn] [--up-fc 1.5] [--padj 0.01]",
  "[--unch-fc 1.3] [--unch-padj 0.1] [--lessabs-alpha 0.1] [--min-len 1000] [--genetype protein_coding] [--deciles 10] [--bg 6000]")
argsArg <- which(startsWith(cArgs, "--args")); args <- if (length(argsArg)) cArgs[(argsArg+1):length(cArgs)] else character(0)
parse_flags <- function(a){ out<-list(); i<-1
  while (i<=length(a)){ if(!startsWith(a[i],"--")) stop("unexpected argument '",a[i],"'.\n",usage,call.=FALSE)
    key<-sub("^--","",a[i]); if(key=="help"){out$help<-TRUE;i<-i+1;next}
    if(i==length(a)) stop("missing value for --",key,".\n",usage,call.=FALSE); out[[key]]<-a[i+1]; i<-i+2 }; out }
opt <- parse_flags(args)
if (isTRUE(opt$help) || !length(args)) { message(usage); quit(save="no", status=0) }
for (r in c("dds","gga","contrasts")) if (is.null(opt[[r]])) stop("missing required --",r,"\n",usage,call.=FALSE)
for (f in c(opt$dds,opt$gga)) if (!file.exists(f)) stop("file not found: ",f,call.=FALSE)
outdir    <- opt[["out"]] %||% "matcher_pool_grid_output"
GRID_DIR  <- tolower(opt[["dir"]] %||% "up"); if (!GRID_DIR %in% c("up","dn")) stop("--dir must be up or dn", call.=FALSE)
up_fc     <- as.numeric(opt[["up-fc"]] %||% "1.5");  padj_de <- as.numeric(opt[["padj"]] %||% "0.01")
unch_fc   <- as.numeric(opt[["unch-fc"]] %||% "1.3"); unch_padj <- as.numeric(opt[["unch-padj"]] %||% "0.1")
min_len   <- as.numeric(opt[["min-len"]] %||% "1000"); gene_type <- opt[["genetype"]] %||% "protein_coding"
deciles   <- as.integer(opt[["deciles"]] %||% "10"); lessabs_alpha <- as.numeric(opt[["lessabs-alpha"]] %||% "0.1")
nbg       <- as.integer(opt[["bg"]] %||% "6000"); lfc_up <- log2(up_fc)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
cands <- lapply(strsplit(strsplit(opt$contrasts, ",")[[1]], "/"), function(p){ if(length(p)!=2) stop("each contrast must be num/den",call.=FALSE); trimws(p) })

dds <- readRDS(opt$dds)
if (!"condition" %in% colnames(as.data.frame(colData(dds)))) stop("dds colData has no 'condition' column.", call.=FALSE)
gga <- read.table(opt$gga, header=TRUE, sep="\t", row.names=1, quote="", stringsAsFactors=FALSE)
mh  <- counts(dds, normalized=TRUE)
len <- abs(gga$DominantTES - gga$DominantTSS)+1; names(len) <- rownames(gga)
common <- intersect(rownames(mh), rownames(gga))
smd <- function(m,o) if (length(m)<1||length(o)<2||sd(o)==0) NA_real_ else (mean(m)-mean(o))/sd(o)
bal <- function(n_in,n_ct,s) ifelse(n_in<10,"n<10",ifelse(n_ct==0,"FAIL",ifelse(is.na(s),"-",ifelse(abs(s)<=0.1,"yes","NO"))))

prep <- function(num,den){
  res <- as.data.frame(results(dds, contrast=c("condition",num,den)))
  bs  <- colnames(dds)[dds$condition==den]
  ck  <- rowMeans(mh[common,bs,drop=FALSE])/(len[common]/1000)
  gm  <- data.frame(gene=common, counts_kb=ck, expr=log2(ck+1), loglen=log10(len[common]), row.names=common, stringsAsFactors=FALSE)
  gm  <- gm[is.finite(gm$expr)&is.finite(gm$loglen),]
  elig <- rownames(gm)[gga[rownames(gm),"GeneType"] %in% gene_type & len[rownames(gm)]>=min_len]
  up_ <- intersect(rownames(res)[!is.na(res$padj)&res$padj<=padj_de&res$log2FoldChange>= lfc_up], elig)
  dn_ <- intersect(rownames(res)[!is.na(res$padj)&res$padj<=padj_de&res$log2FoldChange<=-lfc_up], elig)
  lfc_lin <- 2^res$log2FoldChange
  unchL <- rownames(res)[!is.na(res$padj)&!is.na(lfc_lin)&res$padj>unch_padj&lfc_lin<=unch_fc&lfc_lin>=1/unch_fc]
  legacy <- setdiff(intersect(unchL,elig), c(up_,dn_))
  resE <- as.data.frame(results(dds, contrast=c("condition",num,den), lfcThreshold=log2(unch_fc), altHypothesis="lessAbs", alpha=lessabs_alpha))
  lessabs <- setdiff(intersect(rownames(resE)[!is.na(resE$padj)&resE$padj<=lessabs_alpha],elig), c(up_,dn_))
  list(gm=gm, up=up_, dn=dn_, legacy=legacy, lessabs=lessabs)
}
em_set <- function(gm,input,pool) tryCatch(unique(na.omit(expressionMatched(input,pool,
            data.frame(GeneID=c(input,pool), Expression=gm[c(input,pool),"counts_kb"]), deciles, setSeedNum=1))), error=function(e) character(0))
mi_pairs <- function(gm,input,pool,covars){ use<-c(input,pool)
  mf<-data.frame(treat=as.integer(use %in% input), gm[use,covars,drop=FALSE], row.names=use, check.names=FALSE)
  mf<-mf[apply(mf[,covars,drop=FALSE],1,function(r) all(is.finite(r))),,drop=FALSE]
  tryCatch({set.seed(1); fit<-MatchIt::matchit(stats::as.formula(paste("treat ~",paste(covars,collapse=" + "))),
       data=mf, method="nearest", distance="mahalanobis", replace=FALSE, ratio=1)
       data.frame(input_gene=rownames(fit$match.matrix), ctrl=fit$match.matrix[,1], row.names=NULL, stringsAsFactors=FALSE)},
     error=function(e) data.frame(input_gene=character(0), ctrl=character(0), stringsAsFactors=FALSE)) }
methodf <- function(gm,input,pool,method){
  if (method=="expMat") list(set=em_set(gm,input,pool), pairs=NULL)
  else { cov<-if (method=="MatchIt:expr") "expr" else c("expr","loglen"); pr<-mi_pairs(gm,input,pool,cov); list(set=pr$ctrl, pairs=pr) } }
mlabs  <- c("expMat","MatchIt:expr","MatchIt:expr+len")
mcolor <- c("expMat"="dodgerblue3","MatchIt:expr"="darkorange2","MatchIt:expr+len"="red2")
pools  <- c(legacy="legacy", lessabs="lessAbs")

preps <- list(); summ <- data.frame()
for (ct in cands){ nm<-paste(ct[1],"vs",ct[2]); P<-prep(ct[1],ct[2]); preps[[nm]]<-P; g<-P$gm
  for (dd in c("Up","Dn")){ input<-P[[ifelse(dd=="Up","up","dn")]]
    for (pk in names(pools)) for (m in mlabs){ r<-methodf(g,input,P[[pk]],m)
      summ<-rbind(summ, data.frame(Contrast=nm, Dir=dd, pool=pools[pk], Method=m, n_input=length(input), n_ctrl=length(r$set),
        SMD_expr=round(smd(g[r$set,"expr"],g[input,"expr"]),3), SMD_len=round(smd(g[r$set,"loglen"],g[input,"loglen"]),3),
        flag=bal(length(input),length(r$set),smd(g[r$set,"loglen"],g[input,"loglen"])), stringsAsFactors=FALSE)) } } }
write.csv(summ, file.path(outdir,"grid_summary.csv"), row.names=FALSE)
flagcol <- function(f) if (f %in% c("NO","FAIL")) "red3" else if (f=="n<10") "darkorange3" else if (f=="yes") "forestgreen" else "grey40"

pdf(file.path(outdir,"matcher_pool_grid.pdf"), width=11, height=8.5)
## PAGE 1: overall summary table (length SMD per cell, both directions)
par(mar=c(1,1,1,1)); plot.new()
text(0.5,0.965,"Matcher x unchanged-pool: gene-length balance (SMD vs input)", font=2, cex=1.15)
text(0.5,0.932,"3 matchers x 2 pool definitions, each contrast and direction. Expression is matched to ~0 by all three matchers in every cell; this table is the LENGTH axis.", cex=0.62, font=3)
text(0.44,0.885,"— legacy pool (hard cutoff) —", font=2, cex=0.70, family="mono"); text(0.77,0.885,"— lessAbs pool (equivalence) —", font=2, cex=0.70, family="mono")
xc <- c(0.03,0.225, 0.30,0.395,0.49, 0.635,0.73,0.825)
text(xc, 0.86, c("Contrast","Dir","expMat","MI:expr","MI:e+l","expMat","MI:expr","MI:e+l"), adj=0, cex=0.60, family="mono", font=2)
cd <- unique(summ[,c("Contrast","Dir")]); cd<-cd[order(match(cd$Contrast,names(preps)),match(cd$Dir,c("Up","Dn"))),]
cellval <- function(nm,d,pl,m){ r<-summ[summ$Contrast==nm&summ$Dir==d&summ$pool==pl&summ$Method==m,]
  if(!nrow(r)) return(list(t="-",c="grey40")); f<-r$flag[1]
  list(t=if(f=="FAIL")"FAIL" else if(f=="n<10")"n<10" else sprintf("%+.2f", r$SMD_len[1]), c=flagcol(f)) }
y<-0.86
for (i in seq_len(nrow(cd))){ nm<-cd$Contrast[i]; d<-cd$Dir[i]; y<-y-0.03
  if (i>1 && nm!=cd$Contrast[i-1]) y<-y-0.008
  text(xc[1],y,nm,adj=0,cex=0.58,family="mono"); text(xc[2],y,d,adj=0,cex=0.58,family="mono")
  combos <- list(c("legacy","expMat"),c("legacy","MatchIt:expr"),c("legacy","MatchIt:expr+len"),c("lessAbs","expMat"),c("lessAbs","MatchIt:expr"),c("lessAbs","MatchIt:expr+len"))
  for (j in seq_along(combos)){ cv<-cellval(nm,d,combos[[j]][1],combos[[j]][2]); text(xc[j+2], y, cv$t, adj=0, cex=0.58, family="mono", col=cv$c) } }
y<-y-0.03
text(0.03,y,"green = length balanced (|SMD|<=0.1); red = imbalanced (NO) or matcher FAILED; orange = n<10 (too few).", adj=0, cex=0.66)
text(0.03,y-0.028,"Pattern: MatchIt:expr+len is green almost everywhere under BOTH pools. expMat and MatchIt:expr miss length and get WORSE under lessAbs", adj=0, cex=0.66)
text(0.03,y-0.052,"(length-biased pool); expMat also FAILs more often there (its decile draw starves on the smaller, denser pool).", adj=0, cex=0.66)
text(0.03,y-0.084,sprintf("The per-contrast scatter and ECDF grids that follow show the %s-regulated set (change with --dir).", GRID_DIR), adj=0, cex=0.66, font=3)

scatter_cell <- function(g, bg, input, set, pairs, color, title){
  plot(g[bg,"expr"], g[bg,"loglen"], pch=16, cex=0.2, col=adjustcolor("grey75",0.35), xlab="log2 expression", ylab="log10 length", main=title, cex.main=0.92, cex.lab=0.85)
  if (!is.null(pairs) && nrow(pairs)>0) segments(g[pairs$input_gene,"expr"], g[pairs$input_gene,"loglen"], g[pairs$ctrl,"expr"], g[pairs$ctrl,"loglen"], col=adjustcolor(color,0.18), lwd=0.5)
  if (length(set)>=1) points(g[set,"expr"], g[set,"loglen"], pch=16, cex=0.45, col=color)
  points(g[input,"expr"], g[input,"loglen"], pch=16, cex=0.5, col="black")
  if (length(set)<1) text(mean(range(g[bg,"expr"])), mean(range(g[bg,"loglen"])), "FAILED", col="red3", font=2, cex=1.3) }
ecdf_cell <- function(g, input, set, color, title){
  plot(ecdf(g[,"loglen"]), col="grey78", lwd=1.3, do.points=FALSE, verticals=TRUE, main=title, cex.main=0.92, xlab="log10 gene length", ylab="cum. fraction", cex.lab=0.85)
  if (length(set)>1) plot(ecdf(g[set,"loglen"]), add=TRUE, col=color, lwd=1.6, do.points=FALSE, verticals=TRUE)
  if (length(input)>=1) plot(ecdf(g[input,"loglen"]), add=TRUE, col="black", lwd=1.7, do.points=FALSE, verticals=TRUE) }
for (nm in names(preps)){ P<-preps[[nm]]; g<-P$gm; input<-P[[GRID_DIR]]
  if (length(input) < 3) next
  set.seed(7); bg <- sample(rownames(g), min(nbg, nrow(g)))
  cells <- list(); for (pk in names(pools)) for (m in mlabs) cells[[paste(pk,m)]] <- methodf(g,input,P[[pk]],m)
  par(mfrow=c(2,3), mar=c(4,4,2.6,1), oma=c(0,2.2,2.6,0))
  for (pk in names(pools)) for (m in mlabs){ r<-cells[[paste(pk,m)]]; scatter_cell(g, bg, input, r$set, r$pairs, mcolor[m], sprintf("%s  ·  %s", pools[pk], m)) }
  mtext(sprintf("%s  —  %s-regulated  (n = %d):  native joint expression × length view  (MatchIt shown with pair-links)", nm, GRID_DIR, length(input)), outer=TRUE, font=2, cex=0.95, line=0.6)
  mtext("lessAbs pool        legacy pool", outer=TRUE, side=2, cex=0.75, font=2, line=0.4)
  par(mfrow=c(2,3), mar=c(4,4,2.6,1), oma=c(0,2.2,2.6,0))
  for (pk in names(pools)) for (m in mlabs){ r<-cells[[paste(pk,m)]]
    ecdf_cell(g, input, r$set, mcolor[m], sprintf("%s  ·  %s", pools[pk], m))
    legend("bottomright", bty="n", cex=0.6, lwd=1.6, col=c("grey78","black",mcolor[m]), legend=c("all genes","input",m)) }
  mtext(sprintf("%s  —  %s-regulated  (n = %d):  gene-length ECDF  (control vs input vs all genes)", nm, GRID_DIR, length(input)), outer=TRUE, font=2, cex=0.95, line=0.6)
  mtext("lessAbs pool        legacy pool", outer=TRUE, side=2, cex=0.75, font=2, line=0.4)
}
invisible(dev.off())
cat(sprintf("\nWrote:\n  %s\n  %s\n", file.path(outdir,"matcher_pool_grid.pdf"), file.path(outdir,"grid_summary.csv")))
