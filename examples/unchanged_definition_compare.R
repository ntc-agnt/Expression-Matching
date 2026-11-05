#!/usr/bin/env Rscript
# unchanged_definition_compare.R
# -----------------------------------------------------------------------------
# Head-to-head of the two ways to DEFINE the "unchanged" control pool:
#   legacy  : hard cutoff      padj > unch-padj  AND  linear FC in [1/unch-fc, unch-fc]   ("not significant AND inside a band")
#   lessAbs : DESeq2 equivalence   results(lfcThreshold=log2(unch-fc), altHypothesis="lessAbs"), padj <= lessabs-alpha   ("CONFIDENTLY inside the band")
#
# Part A: how the two unchanged SETS differ (size, overlap, and their expression / length composition).
# Part B: downstream effect on matching (length balance achievable by each matcher drawing from each pool).
#
# Key finding this report surfaces: lessAbs is statistically cleaner (only genes the data can prove are
# unchanged) but its power-dependence biases the pool toward LONGER, higher-expressed genes, so it is a
# less representative background than the hard cutoff (which tracks all eligible genes).
#
# Flags:  --dds FILE  --gga FILE  --contrasts num/den[,num/den...]  [--out DIR]
#         [--up-fc 1.5] [--padj 0.01] [--unch-fc 1.3] [--unch-padj 0.1] [--lessabs-alpha 0.1]
#         [--min-len 1000] [--genetype protein_coding] [--deciles 10]
# See examples/README.md for details.
# -----------------------------------------------------------------------------
suppressMessages({ library(DESeq2) })
cArgs <- commandArgs()
runningScript <- normalizePath(sapply(strsplit(cArgs[which(startsWith(cArgs, "--file="))], "--file="), "[", 2))
source(file.path(dirname(dirname(runningScript)), "functions.R"))   # expressionMatched(), matchitMatched()
`%||%` <- function(a, b) if (is.null(a)) b else a

usage <- paste("Usage: Rscript examples/unchanged_definition_compare.R --dds FILE --gga FILE",
  "--contrasts num/den[,num/den...] [--out DIR] [--up-fc 1.5] [--padj 0.01]",
  "[--unch-fc 1.3] [--unch-padj 0.1] [--lessabs-alpha 0.1] [--min-len 1000] [--genetype protein_coding] [--deciles 10]")
argsArg <- which(startsWith(cArgs, "--args")); args <- if (length(argsArg)) cArgs[(argsArg+1):length(cArgs)] else character(0)
parse_flags <- function(a){ out<-list(); i<-1
  while (i<=length(a)){ if(!startsWith(a[i],"--")) stop("unexpected argument '",a[i],"'.\n",usage,call.=FALSE)
    key<-sub("^--","",a[i]); if(key=="help"){out$help<-TRUE;i<-i+1;next}
    if(i==length(a)) stop("missing value for --",key,".\n",usage,call.=FALSE); out[[key]]<-a[i+1]; i<-i+2 }; out }
opt <- parse_flags(args)
if (isTRUE(opt$help) || !length(args)) { message(usage); quit(save="no", status=0) }
for (r in c("dds","gga","contrasts")) if (is.null(opt[[r]])) stop("missing required --",r,"\n",usage,call.=FALSE)
for (f in c(opt$dds,opt$gga)) if (!file.exists(f)) stop("file not found: ",f,call.=FALSE)
outdir    <- opt[["out"]] %||% "unchanged_compare_output"
up_fc     <- as.numeric(opt[["up-fc"]] %||% "1.5");  padj_de <- as.numeric(opt[["padj"]] %||% "0.01")
unch_fc   <- as.numeric(opt[["unch-fc"]] %||% "1.3"); unch_padj <- as.numeric(opt[["unch-padj"]] %||% "0.1")
min_len   <- as.numeric(opt[["min-len"]] %||% "1000"); gene_type <- opt[["genetype"]] %||% "protein_coding"
deciles   <- as.integer(opt[["deciles"]] %||% "10"); lessabs_alpha <- as.numeric(opt[["lessabs-alpha"]] %||% "0.1")
lfc_up    <- log2(up_fc); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
cands <- lapply(strsplit(strsplit(opt$contrasts, ",")[[1]], "/"), function(p){ if(length(p)!=2) stop("each contrast must be num/den",call.=FALSE); trimws(p) })

dds <- readRDS(opt$dds)
if (!"condition" %in% colnames(as.data.frame(colData(dds)))) stop("dds colData has no 'condition' column.", call.=FALSE)
gga <- read.table(opt$gga, header=TRUE, sep="\t", row.names=1, quote="", stringsAsFactors=FALSE)
mh  <- counts(dds, normalized=TRUE)
len <- abs(gga$DominantTES - gga$DominantTSS)+1; names(len) <- rownames(gga)
common <- intersect(rownames(mh), rownames(gga))
# smd() (standardized mean difference) is provided by functions.R
bal <- function(n_in,n_ct,s) ifelse(n_in<10,"n<10",ifelse(n_ct==0,"FAILED",ifelse(is.na(s),"-",ifelse(abs(s)<=0.1,"yes","NO"))))

prep <- function(num,den){
  res <- as.data.frame(results(dds, contrast=c("condition",num,den)))
  bs  <- colnames(dds)[dds$condition==den]
  ck  <- rowMeans(mh[common,bs,drop=FALSE])/(len[common]/1000)
  gm  <- data.frame(gene=common, expr=log2(ck+1), loglen=log10(len[common]), len=len[common], counts_kb=ck, row.names=common, stringsAsFactors=FALSE)
  gm  <- gm[is.finite(gm$expr)&is.finite(gm$loglen),]
  elig <- rownames(gm)[gga[rownames(gm),"GeneType"] %in% gene_type & len[rownames(gm)]>=min_len]
  up_ <- intersect(rownames(res)[!is.na(res$padj)&res$padj<=padj_de&res$log2FoldChange>= lfc_up], elig)
  dn_ <- intersect(rownames(res)[!is.na(res$padj)&res$padj<=padj_de&res$log2FoldChange<=-lfc_up], elig)
  lfc_lin <- 2^res$log2FoldChange
  unchL <- rownames(res)[!is.na(res$padj)&!is.na(lfc_lin)&res$padj>unch_padj&lfc_lin<=unch_fc&lfc_lin>=1/unch_fc]
  legacy <- setdiff(intersect(unchL,elig), c(up_,dn_))
  resE <- as.data.frame(results(dds, contrast=c("condition",num,den), lfcThreshold=log2(unch_fc), altHypothesis="lessAbs", alpha=lessabs_alpha))
  lessabs <- setdiff(intersect(rownames(resE)[!is.na(resE$padj)&resE$padj<=lessabs_alpha],elig), c(up_,dn_))
  list(gm=gm, elig=elig, legacy=legacy, lessabs=lessabs, up=up_, dn=dn_)
}
preps <- list(); poolinfo <- data.frame()
for (ct in cands){ nm<-paste(ct[1],"vs",ct[2]); P<-prep(ct[1],ct[2]); preps[[nm]]<-P; g<-P$gm
  legonly<-setdiff(P$legacy,P$lessabs)
  poolinfo <- rbind(poolinfo, data.frame(Contrast=nm, n_legacy=length(P$legacy), n_lessAbs=length(P$lessabs),
    both=length(intersect(P$legacy,P$lessabs)), leg_only=length(legonly), lAbs_only=length(setdiff(P$lessabs,P$legacy)),
    medlen_leg=round(10^median(g[P$legacy,"loglen"])), medlen_lAbs=round(10^median(g[P$lessabs,"loglen"])),
    medlen_dropped=round(10^median(g[legonly,"loglen"])),
    medexpr_leg=round(median(g[P$legacy,"expr"]),2), medexpr_lAbs=round(median(g[P$lessabs,"expr"]),2), stringsAsFactors=FALSE)) }
write.csv(poolinfo, file.path(outdir,"pool_comparison.csv"), row.names=FALSE)

## Part B: matcher length balance from each pool
mcols <- c("expressionMatched (expr)","MatchIt (expr)","MatchIt (expr+len)")
run3 <- function(P, input, pool){ gm<-P$gm
  em <- tryCatch(unique(na.omit(expressionMatched(input, pool, data.frame(GeneID=c(input,pool), Expression=gm[c(input,pool),"counts_kb"]), deciles, setSeedNum=1))), error=function(e) character(0))
  mie<- tryCatch(matchitMatched(input, pool, gm[,c("gene","expr")], covars="expr"), error=function(e) character(0))
  mil<- tryCatch(matchitMatched(input, pool, gm[,c("gene","expr","loglen")], covars=c("expr","loglen")), error=function(e) character(0))
  list(`expressionMatched (expr)`=em, `MatchIt (expr)`=mie, `MatchIt (expr+len)`=mil) }
downstream <- data.frame()
for (nm in names(preps)){ P<-preps[[nm]]; g<-P$gm
  for (d in c(Up="up", Dn="dn")){ input<-P[[d]]
    for (poolnm in c("legacy","lessabs")){ sets<-run3(P,input,P[[poolnm]])
      for (m in mcols){ ids<-sets[[m]]
        downstream <- rbind(downstream, data.frame(Contrast=nm, Dir=names(which(c(Up="up",Dn="dn")==d)), pool=poolnm, Method=m,
          n_input=length(input), n_ctrl=length(ids), SMD_len=round(smd(g[ids,"loglen"],g[input,"loglen"]),3),
          len_bal=bal(length(input),length(ids),smd(g[ids,"loglen"],g[input,"loglen"])), stringsAsFactors=FALSE)) } } } }
write.csv(downstream, file.path(outdir,"downstream_length_balance.csv"), row.names=FALSE)

pdf(file.path(outdir,"unchanged_compare.pdf"), width=11, height=8.5)
## PAGE 1: definitions + sizes / overlap / covariate bias
par(mar=c(1,1,1,1)); plot.new()
text(0.5,0.965,"Defining the unchanged control pool: hard cutoff vs DESeq2 lessAbs equivalence", font=2, cex=1.12)
yy<-0.93
L<-function(t,dy=0.032,cex=0.74,font=1,col="black",xx=0.03,mono=FALSE){yy<<-yy-dy;text(xx,yy,t,adj=0,cex=cex,font=font,col=col,family=if(mono)"mono" else "")}
L(sprintf("legacy  (hard clip) : padj > %g  AND  linear FC in [1/%.2f, %.2f]       — 'not significant AND inside the band'", unch_padj, unch_fc, unch_fc), font=2)
L(sprintf("lessAbs (equivalence): results(lfcThreshold=log2(%.2f), altHypothesis='lessAbs'),  padj <= %g   — 'CONFIDENTLY inside the band'", unch_fc, lessabs_alpha), font=2)
L("Both are intersected with the same eligibility filter and exclude the up/down genes. Only the UNCHANGED rule differs.", cex=0.70, font=3)
yy<-yy-0.012
L("Set sizes and overlap", dy=0.034, font=2, cex=0.82)
x1<-c(0.03,0.30,0.42,0.54,0.65,0.78); text(x1, yy-0.028, c("Contrast","n_legacy","n_lessAbs","in_both","legacy_only","lessAbs_only"), adj=0, cex=0.70, family="mono", font=2); yy<-yy-0.028
for (i in seq_len(nrow(poolinfo))){ r<-poolinfo[i,]; yy<-yy-0.028; text(x1, rep(yy,6), c(r$Contrast,r$n_legacy,r$n_lessAbs,r$both,r$leg_only,r$lAbs_only), adj=0, cex=0.68, family="mono") }
yy<-yy-0.02
L("Covariate composition of each pool (the bias the definition introduces)", dy=0.034, font=2, cex=0.82)
x2<-c(0.03,0.30,0.44,0.58,0.74,0.87); text(x2, yy-0.028, c("Contrast","medLen_leg","medLen_lAbs","medLen_dropped","medExpr_leg","medExpr_lAbs"), adj=0, cex=0.68, family="mono", font=2); yy<-yy-0.028
for (i in seq_len(nrow(poolinfo))){ r<-poolinfo[i,]; yy<-yy-0.028
  text(x2, rep(yy,6), c(r$Contrast, paste0(format(r$medlen_leg,big.mark=","),"bp"), paste0(format(r$medlen_lAbs,big.mark=","),"bp"),
       paste0(format(r$medlen_dropped,big.mark=","),"bp"), r$medexpr_leg, r$medexpr_lAbs), adj=0, cex=0.66, family="mono") }
yy<-yy-0.030
text(0.03,yy,"Read it: the two sets overlap heavily but are NOT nested. lessAbs DROPS the low-power genes the band keeps (legacy_only) — short and lower-", adj=0, cex=0.68, font=3)
text(0.03,yy-0.026,"expressed (medLen_dropped) — and ADDS well-powered genes the hard cutoff called 'significant'. Both moves push the pool toward LONGER,", adj=0, cex=0.68, font=3)
text(0.03,yy-0.052,"higher-expressed genes (medLen_lAbs >> medLen_leg, which ≈ all eligible). lessAbs is a cleaner 'unchanged' definition but a LESS", adj=0, cex=0.68, font=3)
text(0.03,yy-0.078,"representative background — a length/expression bias the downstream matcher must correct.", adj=0, cex=0.68, font=3)

## PAGES 2-3: pool covariate distributions (ECDF), legacy vs lessAbs vs all-eligible
cpool <- c(all="grey70", legacy="dodgerblue3", lessAbs="red2")
pool_ecdf <- function(P, col, lab, title){ g<-P$gm
  plot(ecdf(g[P$elig,col]), col=cpool["all"], lwd=3, do.points=FALSE, verticals=TRUE, main=title, cex.main=0.95, xlab=lab, ylab="cumulative fraction")
  plot(ecdf(g[P$legacy,col]), add=TRUE, col=cpool["legacy"], lwd=2.2, lty=2, do.points=FALSE, verticals=TRUE)
  plot(ecdf(g[P$lessabs,col]), add=TRUE, col=cpool["lessAbs"], lwd=2.2, lty=1, do.points=FALSE, verticals=TRUE) }
nmv <- names(preps)
for (pg in split(nmv, ceiling(seq_along(nmv)/2))){
  par(mfrow=c(2,2), mar=c(4.3,4.3,3,1), oma=c(0,0,2,0))
  for (nm in pg){ P<-preps[[nm]]
    pool_ecdf(P,"expr","log2 expression (counts/kb)", sprintf("%s  —  expression", nm))
    legend("bottomright", bty="n", cex=0.78, lwd=2, lty=c(1,2,1), col=cpool[c("all","legacy","lessAbs")], legend=c("all eligible","legacy pool","lessAbs pool"))
    pool_ecdf(P,"loglen","log10 gene length (bp)", sprintf("%s  —  gene length", nm)) }
  if (length(pg)==1) { plot.new(); plot.new() }
  mtext("Pool composition: lessAbs vs legacy vs all eligible  (legacy tracks all-eligible; lessAbs is shifted toward longer, higher-expressed genes)", outer=TRUE, font=2, cex=0.92)
}

## PAGE 4: downstream length balance, legacy pool vs lessAbs pool
par(mfrow=c(1,1)); layout(1); par(mar=c(1,1,1,1)); plot.new()
text(0.5,0.965,"Downstream effect: length balance each matcher achieves, legacy pool vs lessAbs pool", font=2, cex=1.0)
text(0.5,0.932,"Same inputs and matchers; only the pool definition changes. yes = length balanced (|SMD|<=0.1); NO = imbalanced; FAILED = matcher returned nothing.", cex=0.66, font=3)
getf <- function(ct,dr,mth,pl){ r<-downstream[downstream$Contrast==ct&downstream$Dir==dr&downstream$Method==mth&downstream$pool==pl,]; if(!nrow(r))"-" else r$len_bal[1] }
cd <- unique(downstream[,c("Contrast","Dir")]); cd<-cd[order(match(cd$Contrast,names(preps)),match(cd$Dir,c("Up","Dn"))),]
xc <- c(0.03,0.30,0.37,0.56,0.76); yy<-0.885
text(xc, yy, c("Contrast","Dir","expressionMatched","MatchIt(expr)","MatchIt(expr+len)"), adj=0, cex=0.64, family="mono", font=2)
for (i in seq_len(nrow(cd))){ ct<-cd$Contrast[i]; dr<-cd$Dir[i]; yy<-yy-0.028
  if (i>1 && ct!=cd$Contrast[i-1]) yy<-yy-0.008
  text(xc[1],yy,ct,adj=0,cex=0.62,family="mono"); text(xc[2],yy,dr,adj=0,cex=0.62,family="mono")
  for (j in seq_along(mcols)){ fl<-getf(ct,dr,mcols[j],"legacy"); fe<-getf(ct,dr,mcols[j],"lessabs")
    text(xc[j+2], yy, sprintf("%s / %s", fl, fe), adj=0, cex=0.62, family="mono", col=if(fl!=fe)"red3" else "black") } }
yy<-yy-0.03
text(0.03,yy,"Each cell = length-balance verdict with the legacy pool / with the lessAbs pool. Red = the verdict changes.", adj=0, cex=0.68, font=2)
text(0.03,yy-0.030,"The expression-only matchers get WORSE on length under lessAbs (the pool's length bias shows through). MatchIt(expr+len) stays", adj=0, cex=0.66, font=3)
text(0.03,yy-0.054,"balanced except where the smaller lessAbs pool is too depleted at the input's length range. Net: lessAbs is a more principled", adj=0, cex=0.66, font=3)
text(0.03,yy-0.078,"'unchanged' definition, but its power-dependent length bias does not reduce — and can increase — the need to match on length explicitly.", adj=0, cex=0.66, font=3)
invisible(dev.off())
cat("\n=== Pool comparison ===\n"); print(poolinfo, row.names=FALSE)
cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n", file.path(outdir,"unchanged_compare.pdf"), file.path(outdir,"pool_comparison.csv"), file.path(outdir,"downstream_length_balance.csv")))
