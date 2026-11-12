---
name: expression-matching
description: >
  Use when building expression-matched (or covariate-matched) unchanged-gene control/
  background lists for PRO-seq differential-expression analysis — matching a gene set of
  interest to non-DE "unchanged" genes by expression (and optionally length / GC / other
  covariates) — or when a task in a related NTC repo touches this tool, its inputs, or its
  outputs. Consumes NTC DESeq2 output (results tables + makeHeatmap.rds matrix) and the GGA
  WORKING_FILE; its gene lists subset matrices downstream for heatmaps/metagenes
  (metagene_heatmap_tool) and as controls (pauseIndex, unverified). The primary entry point is
  now matchit_match.R (flag CLI, MatchIt matcher by default); the legacy decile matcher
  (expressionMatched()) is still reachable via --method and in archive/main.R. Aliases:
  Expression-Matching, expressionMatched, matchit_match, MatchIt matching, PRO-seq expression
  matching.
---

# Expression-Matching

## Purpose
Builds **matched control (background) gene lists** for PRO-seq DE analysis: given a gene set of
interest (e.g. upregulated genes), it selects a set of **unchanged (non-DE) genes whose
covariate distribution matches** — at minimum expression, optionally length/GC/other covariates
— so downstream heatmaps/metagenes compare against a fair background rather than all genes. It
offers two matchers: a **MatchIt nearest-neighbour covariate matcher** (default, recommended) and
the **legacy decile expression matcher** (Ben Martin, Adelman lab) for backward compatibility.
Repo: `Expression-Matching` (branch `main`; tag `v1.0`); deployed under
`/n/data1/cores/ntc/scripts/…/Expression-Matching/`. **Restructured 2026-11 — see "What changed".**

## Architecture
**Not a smartSlurm pipeline** — no `runAsPipeline`, no `#@` annotations, no `ssbatch`, no Slurm.
A **standalone R tool**, runs in seconds on a login node. Files:
- **`matchit_match.R`** — the **primary entry point**: a flag-driven `Rscript` CLI (~202 lines).
  Finds its library by string-substituting its own `--file=` path, so run it via `Rscript` and
  keep it next to `functions.R`.
- **`functions.R`** — the library (sourced, not run): `expressionMatched()` (legacy decile),
  `matchitMatched()` (MatchIt), and `matchQC()` (balance table + plots).
- **`archive/main.R`** — the **frozen legacy** positional-argument decile driver (sources
  `../functions.R`); kept for reproducibility of pre-2026-11 runs.
- **`tests/benchmark_matchers.R`** — synthetic stress-test suite (scenarios S0–S5) comparing the
  matchers; **`tests/qc.R`** is an orphaned older violin/Wilcoxon QC (nothing sources it).

## I/O contract
**`matchit_match.R` — flag CLI** (`--help` / no args prints usage):

Required: `--genes` (1-col ENSG text, no header), `--deseq-dir` (NTC DESeq2 output:
`*/*_results_table.txt` + `matrix/*.makeHeatmap.rds`), `--gga` (GGA `*WORKING_FILE.txt`, tab,
GeneID rownames, needs `GeneType`/`DominantTSS`/`DominantTES`), `--manifest` (tab, sample
rownames, `condition` column), `--baseline` (a value in `condition`).

Optional: `--method matchit|expression|legacy` (default `matchit`), `--covars` (comma list,
default `expression`), `--covar-file` (GeneID + numeric cols, tab), `--window` bp (subtracted
from length in counts/kb only, default 0), `--ratio` k (controls per focal, default 1),
`--outdir` (default `matchit_output`).

**Outputs** (in `--outdir`, stamped with method):
- `matched_controls.<method>.txt` — one gene ID per line (**single list; no seed replicates**).
- `run_info.<method>.txt` — method, covariates, baseline, n_genes_of_interest, n_matched, timestamp.
- `match_qc.<method>_balance.txt` — covariate, SMD, balanced? (threshold `|SMD| ≤ 0.1`,
  SMD = (mean_matched − mean_input)/sd_input). QC always runs on the requested covariates, even
  for a legacy run (so a covariate the legacy matcher leaves unbalanced is visible).
- `match_qc.<method>_QC.png` — per covariate: quantile, density, ECDF panels (input black,
  matched red-dashed, all-genes grey); skipped with a warning if PNG capability is absent.

```
Rscript matchit_match.R --genes genes.txt --deseq-dir DESeqOutputDir --gga GGA.WORKING_FILE.txt \
  --manifest DESeq_manifest.txt --baseline WT [--method matchit] [--covars expression,length] \
  [--ratio 1] [--window 0] [--outdir matchit_output]
```
**Legacy `archive/main.R`** (positional, decile, **no QC**): `Rscript archive/main.R genes.txt
DESeqDir GGA.WORKING_FILE.txt <window> manifest.txt <baseline>` → writes the old
`expression_matched-<baseline>_{1,2,3}.txt` (seeds 1–3) to the CWD.
Full algorithm, covariate resolution, column positions, and schemas: `reference/algorithm-and-io.md`.

## Dependencies & environment
- **R base + `MatchIt` (CRAN)** for the default matcher and the benchmark; the legacy decile
  method needs no extra package. Optional: `nullranges`/`GenomicRanges` (benchmark only),
  `vioplot` (only the orphaned `tests/qc.R`). No `library()` calls — uses `requireNamespace` + `::`.
- **Confirmed ntc_env gap (checked 2026-11-12, NTC-027a):** **`MatchIt` is in none of the 8
  `ntc_env` YAMLs** (`ntc`/`ntc_basic`, `proseq_env`, `rnaseq_env`, `gga` carry R-base + DESeq2/
  bioconductor but not MatchIt). EM is a downstream DESeq2-consuming R tool, so it most naturally
  runs under **`ntc`/`ntc_basic`** — add `r-matchit` (and, for the benchmark, optional
  `r-nullranges` + `r-genomicranges`) there. The legacy decile path still runs under bare R.
- Reference data: none of its own — depends entirely on **GGA annotation** (WORKING_FILE) and
  **NTC DESeq2 output**. No genome FASTA.

## Modifying it safely (gotchas & assumptions)
- **PRO-seq only** — the length-normalized counts/kb model assumes nascent-transcription signal
  over gene bodies, not RNA-seq.
- **DESeq column order is hardcoded by position:** `padj = col 6`, linear `foldChange = col 7`
  (treated linear, thresholded `[1/1.3, 1.3]`, not log2). "Unchanged" must hold in **every**
  contrast (intersection across all `*_results_table.txt`). A different results-table layout
  silently mis-selects. Thresholds `p_thresh=0.1`, `fc_thresh=1.3` are hardcoded.
- **Input guard (new):** stops if **<50 %** of `--genes` are found in the annotation (hints to
  check ENSG IDs / a stray header line); 50–99 % are **silently dropped** with only a printed count.
- **Pool now excludes the genes of interest** (old design did not) — pool = unchanged ∩ GGA,
  `GeneType == protein_coding`, length `|TES−TSS|+1 ≥ 1000`, minus the focal set.
- **Exactly one `matrix/*.makeHeatmap.rds`** is read (`readRDS`) — multiple matches break it.
- **MatchIt method:** `matchit(treat ~ covars, method="nearest", distance="mahalanobis",
  replace=FALSE, ratio=--ratio)`; deterministic (`set.seed(1)` for ties only); drops genes with
  non-finite covariates; **fails if the usable control pool < n_focal × ratio**.
- **Legacy method (`expression`/`legacy`):** raw counts/kb, 10 quantile groups, **seed 1 only**,
  expression-only (**extra `--covars` ignored with a NOTE**, `--ratio` not applied); an empty
  quantile bin makes `sample()` error — the new code catches it and suggests `--method matchit`.
  Its result may differ slightly from `archive/main.R` seed 1 because of the "minus focal" change.
- **Covariate names:** `expression` = log2(counts_kb+1), `length` = log10(len), any GGA numeric
  column by name, any `--covar-file` column; `log2:`/`log10:` prefixes add a transform; an unknown
  covariate name is a **fatal error**.
- **README is stale in places** — it still refers to a root `main.R` and auto-QC via `qc.R`
  (`_QC.txt`/`_violin.png`) that no current entry point produces; trust the code over the README.

## Integration points
- **Upstream (confirmed in code):** NTC **DESeq2** output (`*_results_table.txt` + `matrix/*.makeHeatmap.rds`)
  + the DESeq2 sample manifest, and **GGA** (`GeneAnnotationScripts` / Genome_Onboarding) WORKING_FILE.
  Which pipeline emits the DESeq2 output isn't named in code — it's the `pipeline2` / `RNAseqMapping`
  DESeq2 step (TODO confirm which emits `makeHeatmap.rds`). Optional new input: user `--covar-file`.
- **Downstream:** gene-ID lists subset matrices for heatmaps/metagenes → `metagene_heatmap_tool`;
  `pauseIndex` may use them as control sets (unverified; no explicit consumer named in code).
- **⚠ Naming caveat:** output filenames **changed** (`matched_controls.<method>.txt`, single list)
  vs the old `expression_matched-<baseline>_{1,2,3}.txt` (now only from `archive/main.R`). Any
  downstream tool/doc expecting the old names or the three replicates needs updating.
- **Shared conventions:** ENSG IDs as row names; sample names normalized via
  `make.names(..., allow_=FALSE)` to match matrix columns.

## What changed (2026-11 restructure)
- **New primary tool `matchit_match.R`**: flag CLI, a MatchIt covariate matcher (default), covariate
  spec (`--covars`), `--ratio`, `--covar-file`, SMD-balance + quantile/density/ECDF QC, an input
  guard, and the pool now excludes the focal genes.
- **Moved/merged:** `main.R` → `archive/main.R`; `expressionMatched.R` + QC merged into `functions.R`
  (`tests/qc.R` left orphaned).
- **Default behaviour:** MatchIt, expression-only unless `--covars expression,length` is given
  (expression-only matching leaves length unbalanced — the README recommends adding length).
- **Legacy decile still reachable** via `--method expression|legacy` (seed 1, no replicates) or the
  frozen `archive/main.R` (3 seeds, old filenames).
- **New dep MatchIt** and **new output filenames** (see caveats above).

## Reference files
- `reference/algorithm-and-io.md` — the flag CLI + legacy positional args, covariate resolution,
  the unchanged-gene selection, both matching algorithms (MatchIt + decile), QC, I/O schema, the
  benchmark suite, and evolution/gotchas.
- `.claude/reference/runAsPipeline-grammar.md` — smartSlurm grammar (this tool is *not* rAP, but
  most of the ecosystem is; useful context).
