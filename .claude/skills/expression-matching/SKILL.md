---
name: expression-matching
description: >
  Use when building expression-matched unchanged-gene control/background lists for
  PRO-seq differential-expression analysis — matching a gene set of interest to
  non-DE "unchanged" genes by quantile-matched, length-normalized PRO-seq expression —
  or when a task in a related NTC repo touches this tool, its inputs, or its outputs.
  Consumes NTC DESeq2 output (results tables + makeHeatmap.rds matrix) and the GGA
  WORKING_FILE; its gene lists subset matrices downstream for heatmaps/metagenes
  (metagene_heatmap_tool) and as controls (pauseIndex, unverified). Covers running
  main.R, the unchanged-gene selection + expressionMatched() algorithm, and its I/O
  schema. Aliases: Expression-Matching, expressionMatched, PRO-seq expression matching.
---

# Expression-Matching

## Purpose
Builds **expression-matched control (background) gene lists** for PRO-seq DE analysis:
given a gene set of interest (e.g. upregulated genes), it selects an equal-size set of
**unchanged (non-DE) genes with a matched length-normalized expression distribution**, so
downstream heatmaps/metagenes compare against a fair background rather than all genes. It
wraps Ben Martin's (Adelman lab) `expressionMatched()` quantile-matching function for NTC
DESeq2 + GGA outputs. Repo: `Expression-Matching` (branch `master`); deployed at
`/n/data1/cores/ntc/scripts/NascentTranscriptionCore/Expression-Matching/`.

## Architecture
**Not a smartSlurm pipeline** — no `runAsPipeline`, no `#@` annotations, no `ssbatch`, no
flag folder, no Slurm. It is a **standalone R tool**: `main.R` is an `Rscript` CLI driver
that `source()`s the `expressionMatched.R` library (found by name-swapping `main.R` →
`expressionMatched.R` in its own resolved `--file=` path). Runs in seconds on a login node.
Keep both files together and named `main.R` / `expressionMatched.R`.

## I/O contract
**Inputs (6 positional args, order-sensitive, no usage message or validation):**
1. 1-column `.txt` of ENSG IDs of interest.
2. DESeq2 output dir (globs `*/*_results_table.txt` and `matrix/*.makeHeatmap.rds`).
3. GGA `WORKING_FILE` (tab, needs `GeneType` / `DominantTSS` / `DominantTES`).
4. Window bp to subtract from TSS–TES length for the expression normalization (0 = full
   TSS–TES; 250 = +250→TES, etc.).
5. DESeq2 manifest (tab, needs a `condition` column; row names are sample paths).
6. Baseline condition name (must exactly match a value in the manifest `condition` column;
   normally WT/control).

**Outputs:** `expression_matched-<baseline>_{1,2,3}.txt` — three ENSG-per-line control
lists (one per seed) written to the **current working directory** (no out-dir option).

Invocation:
```
Rscript main.R genes_ensg.txt DESeqOutputDir GGA.WORKING_FILE.txt <window_bp> \
  DESeq_manifest.txt <baseline_condition>
```
Full algorithm, column positions, and schemas: `reference/algorithm-and-io.md`.

## Dependencies & environment
- **R ≥4.x, base only** — no `library()` calls (uses `scan`, `read.table`, `Sys.glob`,
  `readRDS`, `quantile`, `sample`). TODO: confirm the exact R module / conda env on O2;
  likely the shared `ntc_basic` env (see `ntc-env` skill).
- Reference data: none of its own. Depends entirely on **GGA annotation** (WORKING_FILE;
  README: "currently only works with GGA annotation") and **NTC DESeq2 output**. No genome
  FASTA.

## Modifying it safely (gotchas & assumptions)
- **PRO-seq only** — code comment flags this; the length-normalized counts/kb model assumes
  nascent-transcription signal over gene bodies, not RNA-seq.
- **DESeq column order is hardcoded by position**: `padj = column 6`, `foldChange =
  column 7`, and FC is treated as **linear** (thresholded `[1/1.3, 1.3]`), not log2. A
  different DESeq results-table layout silently mis-selects.
- **`matrix/*.makeHeatmap.rds` glob has no length check** — multiple matching `.rds` files
  break `readRDS` (expects one path).
- **Matched list can be shorter than the gene set**: per bin it samples
  `round(length(geneset)/10)` genes, so `round(n/10)*10 ≠ n`, and `sample()` errors if a
  quantile bin is too sparse (README fix: lower `numQuantGroups`, currently hardcoded 10,
  or the control-set ratio).
- **`set.seed(setSeedNum)` is reset at the top of every bin loop** → the same seed drives
  each bin within a run; the 1/2/3 seeds give the three alternative lists.
- **All thresholds hardcoded** (`p_thresh=0.1`, `fc_thresh=1.3`, `numQuantGroups=10`,
  3 seeds): clone/edit `main.R` to change them.
- The **1 kb length filter uses the full TSS–TES** (`|TES−TSS|+1 ≥ 1000`); the window arg
  is subtracted only in the counts/kb normalization, not in the filter.

## Integration points
- **Upstream:** NTC PRO-seq **DESeq2** step (`*_results_table.txt` + `matrix/*.makeHeatmap.rds`),
  likely in `pipeline2` / `RNAseqMapping` (TODO confirm which emits `makeHeatmap.rds`); and
  **GGA** (`GeneAnnotationScripts` / Genome_Onboarding) for the WORKING_FILE.
- **Downstream:** its gene lists subset matrices for heatmaps/metagenes →
  `metagene_heatmap_tool`; plausibly `prepareBrowser`; `pauseIndex`/TSScall may use them as
  control sets (unverified).
- **Shared conventions:** ENSG IDs as row names; sample names normalized via
  `make.names(..., allow_=FALSE)` to match matrix columns.

## Reference files
- `reference/algorithm-and-io.md` — the 6 args, full unchanged-gene selection +
  quantile-matching algorithm, I/O schema, and evolution/gotchas.
- `.claude/reference/runAsPipeline-grammar.md` — smartSlurm grammar (this tool is *not*
  rAP, but most of the ecosystem is; useful context).
