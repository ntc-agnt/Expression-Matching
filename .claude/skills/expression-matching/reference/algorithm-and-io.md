# Expression-Matching — algorithm, arguments, and I/O schema

**Restructured 2026-11.** The program is now `matchit_match.R` (a flag-driven CLI) plus the
`functions.R` library it sources. The old positional driver `main.R` is frozen at
`archive/main.R`. `matchit_match.R` locates its library by `sub()`-ing its own resolved
`--file=` path to `functions.R`, so keep the two together and run via `Rscript`.

## `matchit_match.R` — flags

`--help` or no args prints usage. Required flags error out if missing or the file is absent.

| Flag | Required | Default | Meaning |
|--|--|--|--|
| `--genes` | yes | — | 1-column `.txt` of ENSG IDs of interest (`scan`, no header) |
| `--deseq-dir` | yes | — | NTC DESeq2 output dir (`*/*_results_table.txt` + `matrix/*.makeHeatmap.rds`) |
| `--gga` | yes | — | GGA `*WORKING_FILE.txt` (tab, GeneID rownames; needs `GeneType`/`DominantTSS`/`DominantTES`) |
| `--manifest` | yes | — | DESeq2 sample manifest (tab, sample rownames, `condition` column) |
| `--baseline` | yes | — | a value in `condition` (the baseline samples; normally WT/control) |
| `--method` | no | `matchit` | `matchit`, or `expression`/`legacy` (the decile matcher) |
| `--covars` | no | `expression` | comma list; see covariate resolution |
| `--covar-file` | no | none | tab TSV: `GeneID` + numeric columns (e.g. precomputed GC) |
| `--window` | no | `0` | bp subtracted from length **in counts/kb only** (not the length filter) |
| `--ratio` | no | `1` | controls per focal gene (MatchIt `ratio`); ignored by the legacy method |
| `--outdir` | no | `matchit_output` | output directory (created if absent) |

## Algorithm (matchit_match.R)

**1. Genes of interest.** `scan(--genes)`. Report how many are found in the annotation;
**`stop()` if < 50 %** (hint: check ENSG IDs / a stray header). 50–99 % are dropped silently
with a printed count.

**2. Unchanged pool.** `Sys.glob(--deseq-dir/*/*_results_table.txt)`, each read
`header=TRUE, sep="\t", row.names=1`. A gene is "unchanged" when **by column position**:
`col 6` (padj) non-NA and `> 0.1` (`p_thresh`), **and** `col 7` (linear foldChange) within
`[1/1.3, 1.3]` (`fc_thresh`). Must hold in **every** results table (`Reduce(intersect, …)`).

**3. GGA filter.** Keep intersected genes with `GeneType == "protein_coding"` and length
`|DominantTES − DominantTSS| + 1 ≥ 1000` (full TSS–TES; `--window` not applied here), **minus
the genes of interest** (new in 2026-11; the old design did not subtract them).

**4. Expression.** Read exactly one `Sys.glob(--deseq-dir/matrix/*.makeHeatmap.rds)` (`readRDS`).
Baseline samples = manifest rows where `condition == --baseline`; sample names taken as
`basename` and normalized `make.names(..., allow_=FALSE)` to match matrix columns.
`counts_kb = rowMeans(baseline cols) / ((len − --window)/1000)` over genes in the rds ∩ GGA.

**5. Covariate resolution** (`--covars`, comma list):
- `expression` → `log2(counts_kb + 1)`
- `length` → `log10(gene length)`
- any **GGA numeric column** by name; any **`--covar-file` column** by name
- `log2:`/`log10:` prefixes add a transform (GGA + file covariates)
- an **unknown covariate name is a fatal error**

**6. Match (`--method`):**

### `matchit` (default) — `matchitMatched()`
`MatchIt::matchit(treat ~ covars, method="nearest", distance="mahalanobis", replace=FALSE,
ratio=--ratio)`. 1:1 (or 1:k) nearest-neighbour matching in joint covariate (Mahalanobis) space,
without replacement. Deterministic (`set.seed(1)` for ties only). Drops genes with non-finite
covariates. **Fails if the usable control pool < n_focal × ratio.**

### `expression` / `legacy` — `expressionMatched()`
The original decile matcher, on **raw counts/kb**, 10 quantile groups, **seed 1 only**,
**expression-only** (extra `--covars` ignored with a NOTE; `--ratio` not applied). Result wrapped
`unique(na.omit(...))`. A `tryCatch` turns an empty-bin `sample()` error into a hint to use
`--method matchit`.

`expressionMatched(geneset, ctrlgeneset, dat, numQuantGroups, CtrlSetRatio=1, setSeedNum=1)`:
quantile breakpoints of the **geneset's** expression; `num_per_quantile =
CtrlSetRatio * round(length(geneset)/numQuantGroups)`; per bin `set.seed(setSeedNum)` then
`sample(ctrl genes in [break_i, break_{i+1}), num_per_quantile)`; return the concatenation.

**7. QC + write.** `matchQC()` always runs on the requested covariates (even for a legacy run, so
an unbalanced covariate is visible). Outputs below.

## Outputs (in `--outdir`, `<method>` = `matchit` | `expression`)
- `matched_controls.<method>.txt` — one gene ID per line (**single list, no seed replicates**).
- `run_info.<method>.txt` — key/value: method, covariates, baseline, n_genes_of_interest,
  n_matched, generated timestamp.
- `match_qc.<method>_balance.txt` — columns `covariate, SMD, balanced` (`|SMD| ≤ 0.1`;
  `SMD = (mean_matched − mean_input)/sd_input`).
- `match_qc.<method>_QC.png` — per covariate: quantile, density, ECDF panels (input black,
  matched red-dashed, all-genes grey). Skipped with a warning if PNG capability is missing.

## Legacy driver `archive/main.R` (frozen)
Positional args, in order: `genes deseqDir GGA.WORKING_FILE window manifest baseline`. Runs the
decile matcher for seeds 1–3 and writes `expression_matched-<baseline>_{1,2,3}.txt` to the CWD.
**No QC.** Kept for reproducing pre-2026-11 runs. Note: its control pool differs slightly from
`matchit_match.R --method expression` because the latter subtracts the focal genes (TODO: not
verified by running).

## Benchmark — `tests/benchmark_matchers.R`
`Rscript tests/benchmark_matchers.R [outdir=benchmark_out] [n_genes=12000] [seed=42]`. Synthetic
stress-test across scenarios **S0_null, S1_lengthconf, S2_interaction, S3_sparse, S4_hiddenGC,
S5_tight_conf**, comparing expression-only decile vs MatchIt on expr+len and expr+len+GC (and
optionally `nullranges::matchRanges`). Emits `scorecard.txt`/`.png` + per-scenario `_QC.png`/
`_balance.txt`. (`tests/qc.R` with `expressionMatchQC()` — the old violin/Wilcoxon QC — is
orphaned; nothing sources it.)

## Dependencies
- R base + **`MatchIt`** (CRAN) for `matchit` + the benchmark; legacy decile needs no extra pkg.
  Optional `nullranges`/`GenomicRanges` (benchmark), `vioplot` (orphaned `qc.R`). No `library()` —
  uses `requireNamespace` + `::`. **TODO: confirm which `ntc_env` env supplies R + MatchIt (likely
  a gap — MatchIt may not be in `ntc`/`ntc_basic`).**

## Evolution, rationale & gotchas
- Core decile matcher authored Feb 2022 (`benjaminmartin02`/`gmnelson`, Adelman lab); the NTC
  wrapper and the 2026-11 MatchIt rework are recent (`ntc-agnt`, Seth Goldman, `sgoldman101`).
  **Tool is in flux (almost all new work 2026-11-03/04).** Tag `v1.0`.
- **Why MatchIt:** expression-only matching leaves **length** (and other covariates) unbalanced;
  MatchIt matches jointly in covariate space. README recommends `--covars expression,length` (or
  more) for new work; the decile path is kept for backward compatibility.
- **DESeq column positions hardcoded** (`padj=col 6`, linear `FC=col 7`) — the most fragile
  coupling; any upstream results-table layout change silently corrupts selection.
- **Exactly one `.makeHeatmap.rds`** is read (no length guard).
- **Output filenames changed** — downstream tools/docs expecting `expression_matched-<baseline>_{1,2,3}.txt`
  or three replicates must update to `matched_controls.<method>.txt` (single list), or use the
  frozen `archive/main.R`.
- **README is stale** in places: it still refers to a root `main.R` and auto-QC via `qc.R`
  (`_QC.txt`/`_violin.png`) that no current entry point produces — trust the code.
- **PRO-seq only** — counts/kb is a nascent-transcription length-normalization model.
