# Expression-Matching — algorithm, arguments, and I/O schema

`main.R` is the whole program (plus the `expressionMatched.R` library it sources). It takes
**6 positional args** via `commandArgs()` (parsed relative to the `--args` marker), with no
usage message and no argument validation. Argument order is load-bearing.

## Positional arguments

| # | var in `main.R` | Meaning | Format |
|--|--|--|--|
| 1 | `genes_of_interest` | gene set of interest (e.g. up- or downregulated) | 1-column `.txt` of ENSG IDs (`scan(..., character())`) |
| 2 | `deseqDir` | DESeq2 output directory | dir globbed for `*/*_results_table.txt` and `matrix/*.makeHeatmap.rds` |
| 3 | `wf_path` | GGA WORKING_FILE | tab TSV, `row.names=1` (ENSG); needs `GeneType`, `DominantTSS`, `DominantTES` |
| 4 | `subtract_TSS_TES` | bp subtracted from TSS–TES length for the counts/kb normalization | numeric; `0` = full TSS–TES, `250` = +250→TES, etc. |
| 5 | `deseq_manifest` | DESeq2 sample manifest | tab TSV, `row.names=1` (sample path); needs a `condition` column |
| 6 | `ref_condition` | baseline condition for expression matching | string; must exactly equal a value in manifest `condition` (normally WT/control) |

`main.R` locates `expressionMatched.R` by `sub("main.R", "expressionMatched.R", <resolved
--file= path>)` — so the two files must sit together and keep those exact names.

## Algorithm

**1. Load DESeq2 results.** `Sys.glob(deseqDir/*/*_results_table.txt)`; each read with
`header=TRUE, sep="\t", row.names=1`. List names = basename stripped of `_results_table.txt`
(one entry per comparison/contrast directory).

**2. Select unchanged genes per comparison, then intersect.** Hardcoded thresholds
`p_thresh=0.1`, `fc_thresh=1.3`. A gene is "unchanged" in a comparison when, **by column
position**:
- `col 6` (padj) is non-NA and `> 0.1`, **and**
- `col 7` (foldChange, treated as **linear**) is non-NA and within `[1/1.3, 1.3]`
  (`x[,7] <= 1.3 & x[,7] >= 1.3^-1`).

`unchanged_intersect = Reduce(intersect, ...)` across all comparisons — a gene must be
unchanged in **every** results table.

**3. GGA filter.** Subset the WORKING_FILE to the intersected genes, then keep only:
- `GeneType == "protein_coding"`, **and**
- `|DominantTES − DominantTSS| + 1 ≥ 1000` (full TSS–TES length; the window arg is **not**
  applied here).

Result = `unchanged_genes`. Combined analysis set `gene_ids = c(genes_of_interest,
unchanged_genes)`.

**4. Build length-normalized expression values.** Read the single
`Sys.glob(deseqDir/matrix/*.makeHeatmap.rds)` via `readRDS`. Select rows = `gene_ids`,
columns = baseline-condition samples: manifest rows where `condition == ref_condition`, with
sample names taken as the last path component and normalized by
`make.names(..., allow_=FALSE)` to match the matrix column names. Then:
```
expr = rowMeans(baseline_counts) / ((|DominantTES − DominantTSS| + 1 − subtract_TSS_TES) / 1000)
```
i.e. **mean baseline counts per kb** of the (window-adjusted) gene length →
`exp_dat = data.frame(GeneID, Expression)`.

**5. Quantile-match, 3 seeds.** For `i in 1:3`, call
`expressionMatched(genes_of_interest, unchanged_genes, exp_dat, 10, setSeedNum=i)` and
`write()` the result to `expression_matched-<ref_condition>_<i>.txt` in the CWD.

### `expressionMatched(geneset, ctrlgeneset, dat, numQuantGroups, CtrlSetRatio=1, setSeedNum=1)`
From `expressionMatched.R`:
1. `geneIDs = dat[,1]`, `expressionValues = dat[,2]`.
2. Compute `numQuantGroups` (=10) quantile breakpoints of the **geneset's** expression:
   `quantile(expr[geneIDs %in% geneset], seq(0, 1, 1/numQuantGroups))`.
3. `num_per_quantile = CtrlSetRatio * round(length(geneset)/numQuantGroups, 0)`.
4. For each of the `numQuantGroups` bins: `set.seed(setSeedNum)` (reset **every** bin), then
   `sample(ctrl genes whose expr is in [break_i, break_{i+1}), num_per_quantile)`; append.
5. Return the concatenated control-gene vector.

## I/O summary

- **Inputs:** ENSG list `.txt`; DESeq2 dir (`*/*_results_table.txt` + `matrix/*.makeHeatmap.rds`);
  GGA WORKING_FILE (tab); DESeq2 manifest (tab, `condition` col); window int; baseline string.
- **Outputs:** `expression_matched-<baseline>_{1,2,3}.txt`, one ENSG per line, three
  alternative matched control lists (one per seed), written to the launch CWD. No out-dir
  option, no log file.

## Evolution, rationale & gotchas (distilled from history)
- Small and stable: ~6 commits. Core authored Feb 2022 (`gmnelson` / `benjaminmartin02`,
  Adelman lab); latest change (Sept 2025, `sgoldman101`) was README-only. Treat the logic as
  frozen; current maintainer likely NTC / `sgoldman101`.
- **PRO-seq only** — explicit code comment ("probably only works for PRO-seq DE at the
  moment"). The counts/kb length normalization is a nascent-transcription model.
- **DESeq column positions hardcoded** (`padj = col 6`, `FC = col 7`, FC linear not log2).
  This is the most fragile coupling — any change in the upstream DESeq results-table layout
  silently corrupts unchanged-gene selection. A design shortcut, not a feature.
- **`.makeHeatmap.rds` glob assumes exactly one match** — no length guard; multiple files
  break `readRDS`.
- **Matched set can be shorter than the gene set**: it pulls `round(n/10)` per bin, so
  `round(n/10)*10 ≠ n` in general, and `sample()` throws "cannot take a sample larger than
  the population" on a sparse bin. README remedy: lower `numQuantGroups` (hardcoded to 10) or
  the `CtrlSetRatio`. Running three seeds is the intended QC — compare all three lists.
- **`set.seed` is reset per bin**, so within one run each bin draws from the same seed; the
  three output files differ only by `setSeedNum` 1/2/3.
- **Everything is hardcoded** (thresholds, bin count, seed count, column mapping) — the
  README's "adjustable" knobs mean "clone the script and edit `main.R`."
