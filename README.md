# PRO-seq Control-Gene Matching

Build a **control / background gene set** that is matched to a set of genes of
interest (e.g. up-regulated genes), for PRO-seq differential-expression analyses.
Controls are drawn from genes that are **unchanged** across the experiment and are
**matched to the query set's covariate distribution** — expression, and optionally
gene length and other covariates — so that downstream comparisons are not confounded
by those covariates.

Expression-only matching leaves **gene length** (a real PRO-seq confounder)
unbalanced; the recommended workflow matches on expression **and** length jointly.

---

## Tools in this repo

| Path | What it is |
|------|------------|
| **`matchit_match.R`** | **Recommended workflow.** Flagged CLI that selects controls by 1:1 nearest-neighbour Mahalanobis matching (`MatchIt`) on one or more covariates (expression, length, GC, …) and writes a full match-quality QC. Can also run the legacy expression-only algorithm. |
| `functions.R` | The library: `expressionMatched()` (legacy decile matcher), `matchitMatched()` (covariate matcher), `matchQC()` (balance table + plots). Sourced by the drivers; usable interactively. |
| `archive/main.R` | **Legacy** decile expression-matching driver (frozen, positional args). Kept for reproducing historical control sets. |
| `tests/` | `benchmark_matchers.R` — synthetic stress-test comparing the matchers across confound scenarios. |
| `examples/` | Reproducible reports on a **real** DESeq2 dataset comparing matchers and unchanged-pool definitions. See [`examples/README.md`](examples/README.md). |

---

## How it works

For a supplied list of genes of interest, the matcher:

1. **Defines the unchanged pool** — genes that are *not* differentially expressed in
   **every** contrast (`padj > 0.1` **and** linear fold change within `[1/1.3, 1.3]`),
   restricted to **protein-coding** genes **≥ 1 kb**. Those two thresholds are set per
   run with `--padj` and `--fc-thresh` and encode a judgement about what counts as
   "unchanged" (see *Defining "unchanged"* below). Your genes of interest are supplied
   by you and are not subjected to these filters.
2. **Computes per-gene covariates** from the **baseline** condition — length-normalized
   expression (counts/kb, `log2`), gene length (`log10`), plus any extra numeric
   covariates you provide.
3. **Selects matched controls** from the unchanged pool, 1:1 to the query set, by
   nearest-neighbour Mahalanobis matching on the chosen covariates (or, in legacy mode,
   by decile expression matching).
4. **Writes the matched list(s) and a match-quality QC** — a standardized-mean-difference
   (SMD) balance table (`|SMD| < 0.1` = balanced) and per-covariate distribution plots.

---

## Requirements

- **R ≥ 4.x**
- **`MatchIt`** (CRAN) — for the covariate matcher.
- **`DESeq2`** (Bioconductor) — only for the `examples/` reports, not for the matcher.

Inputs (NTC PRO-seq pipeline outputs):

- a plain-text file of **gene IDs of interest**, one per line;
- **DESeq2 outputs**: per-contrast `*_results_table.txt` files and a
  `matrix/*.makeHeatmap.rds` count matrix;
- a **GGA `WORKING_FILE`** annotation (provides `GeneType`, `GeneName`,
  `DominantTSS`/`DominantTES`; gene length = `|DominantTES − DominantTSS| + 1`);
- a tab-delimited **DESeq2 sample manifest** with a `condition` column.

---

## Quick start (recommended workflow)

```sh
Rscript matchit_match.R \
  --genes      genes_of_interest.txt \
  --deseq-dir  DESeqOutput/ \
  --gga        GGA.WORKING_FILE.txt \
  --manifest   DESeq_manifest.txt \
  --baseline   WT \
  --covars     expression,length \
  --outdir     matchit_output
```

| Flag | Default | Meaning |
|------|---------|---------|
| `--genes` | — | 1-column file of gene IDs of interest (e.g. up-regulated genes). |
| `--deseq-dir` | — | DESeq2 output folder (globbed for `*/*_results_table.txt` and `matrix/*.makeHeatmap.rds`). |
| `--gga` | — | GGA `WORKING_FILE` annotation. |
| `--manifest` | — | Tab-delimited DESeq2 manifest with a `condition` column. |
| `--baseline` | — | `condition` value whose samples define the baseline expression (normally the control/WT). |
| `--covars` | `expression` | Comma-separated covariates to match on (see below). |
| `--covar-file` | — | Optional `GeneID` + numeric columns table (e.g. precomputed GC). |
| `--method` | `matchit` | `matchit` = Mahalanobis covariate matching; `expression`/`legacy` = decile expression-only (reproduces old control sets). |
| `--window` | `0` | bp subtracted from each gene's length in the counts/kb normalization. |
| `--ratio` | `1` | Controls per gene of interest. |
| `--padj` | `0.1` | Unchanged pool: a gene is "not significant" if `padj > this` (in every contrast). |
| `--fc-thresh` | `1.3` | Unchanged pool: linear fold change must lie within `[1/this, this]`. |
| `--outdir` | `matchit_output` | Output directory. |

`--padj` and `--fc-thresh` set what counts as "unchanged" (see below); they are
recorded in `run_info.<method>.txt` for reproducibility.

**`--covars`** accepts the keywords `expression` (counts/kb, `log2`) and `length`
(`log10`), any numeric **GGA WORKING_FILE** column by name, and any column of
`--covar-file`; prefix with `log2:`/`log10:` to transform. To add GC, supply it via
`--covar-file` and pass `--covars expression,length,gc` — the Mahalanobis distance
generalizes to any number of numeric covariates.

**Outputs** (to `--outdir`, stamped with the method): `matched_controls.<method>.txt`,
`run_info.<method>.txt`, and `match_qc.<method>_QC.png` / `_balance.txt` (the SMD
balance table and per-covariate quantile/density/ECDF plots).

---

## Defining "unchanged"

The control pool is the set of genes called **unchanged** across all contrasts. Two
ways to define that were evaluated:

- **hard cutoff (used):** `padj > --padj` **and** linear fold change within a small band
  `[1/--fc-thresh, --fc-thresh]`, defaulting to `--padj 0.1` and `--fc-thresh 1.3`.
  Choosing these is a **user judgement** — the biologically meaningful "no change" effect
  size — which is the methodologically correct way to set it (there is no non-circular
  way to derive it from the same data). For this assay the `1.3×` band also sits at the
  measured replicate **noise floor**, so the default is not arbitrary. (The legacy
  `archive/main.R` keeps these as the `p_thresh` / `fc_thresh` script constants.)
- **`lessAbs` equivalence (evaluated, not used):** DESeq2's equivalence test. Rigorous for
  *stating* that genes are unchanged, but its power-dependence biases the pool toward
  long / high-expression genes and breaks the common support that nearest-neighbour
  matching needs. See the decision note and comparison reports in
  [`examples/`](examples/README.md).

Either way, **matching on the covariates (expression + length) is what controls the
confounders** — so an inclusive, user-thresholded unchanged pool is the right input.

---

## Match-quality QC

Every run writes the QC automatically. Read it as:

- **Balance table (`_balance.txt`):** standardized mean difference (SMD) between the
  matched controls and the query set for each covariate. **`|SMD| ≤ 0.1` = balanced.**
  Expression-only matching typically shows a balanced expression SMD but an imbalanced
  length SMD; joint `expression + length` matching balances both.
- **Distribution plots (`_QC.png`):** per-covariate quantile / density / ECDF curves for
  all genes, the query set, and the matched controls — the matched curve should sit on
  the query curve.
- **Fit overlap scatter (`_overlap.png`, ≥ 2 covariates):** the matched set (red) and the
  query set (black) in the *joint* covariate space (one panel per covariate pair), against
  all genes (grey). The marginal plots above can each look balanced while the joint
  distribution is displaced; here the red points should sit on top of the black. A matched
  cloud that is shifted or forms a band is a bad fit even when the marginals look fine.

---

## Legacy decile matcher (`archive/main.R`)

Frozen, positional-argument driver for the original expression-only decile matcher; kept
for reproducing historical control sets.

```sh
Rscript archive/main.R \
  genes_ensg.txt  DESeqOutput/  GGA.WORKING_FILE.txt  <window>  DESeq_manifest.txt  <baseline>
```

It writes three replicate control lists (seeds 1–3) — independent random draws with the
same intended properties, so inspecting all three checks that the matching is robust to
the draw. It emits the matched lists only (no QC); run the QC via the recommended
workflow, where the same decile algorithm is available as `matchit_match.R --method
expression`.

**Decile starvation:** if a run errors with *"cannot take a sample larger than the
population when 'replace = FALSE'"*, the unchanged pool is too thin in some expression
decile to fill the draw. Lower the number of quantile groups or the control ratio — or
use the Mahalanobis matcher (`--method matchit`), which degrades gracefully instead of
failing.

---

## Examples and benchmarks

- [`examples/`](examples/README.md) — reproducible reports on a real DESeq2 dataset:
  `method_comparison_report.R` (which matcher balances best), `matcher_pool_grid.R` (the
  full matcher × unchanged-pool grid with native MatchIt pair-link scatter views), and
  `unchanged_definition_compare.R` (hard cutoff vs `lessAbs`). They are the rationale
  behind the design choices above.
- `tests/benchmark_matchers.R` — stress-tests the matchers against synthetic confound
  scenarios and writes a covariate-balance scorecard.

---

## Interactive use

`functions.R` can be sourced directly. `matchitMatched(geneset, ctrlgeneset, dat,
covars, ratio)` returns matched control IDs from a covariate frame; `expressionMatched()`
is the legacy decile matcher; `matchQC()` writes the balance table and plots. See the
[Adelman lab Expression-Matching README](https://github.com/AdelmanLab/Expression-Matching)
for the original function documentation and a simulated example.

---

## Credit

- **Expression-matching algorithm:** Benjamin J. E. Martin (Karen Adelman Lab) —
  <https://github.com/AdelmanLab/Expression-Matching>.
- **NTC PRO-seq wrapper (legacy `main.R`):** Geoff Nelson.
- **Covariate matching workflow (`matchit_match.R`):** nearest-neighbour matching via
  `MatchIt` (Ho, Imai, King & Stuart).
- **QC workflow:** adapted from the Adelman lab README.
