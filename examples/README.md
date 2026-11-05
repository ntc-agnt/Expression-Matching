# Worked examples — matcher & pool comparisons on a real DESeq2 dataset

Reproducible reports that compare the control-gene matchers in this repo's
`functions.R` on a **real** DESeq2 dataset. Three scripts, all sharing the same
inputs (`--dds`, `--gga`, `--contrasts`) and conventions:

| Script | Question it answers |
|--------|---------------------|
| `method_comparison_report.R`        | Which **matcher** balances covariates best? (expressionMatched vs MatchIt expr vs MatchIt expr+len) |
| `unchanged_definition_compare.R`    | Which **unchanged-pool definition** should we use? (hard FC+padj cutoff vs DESeq2 `lessAbs` equivalence) |
| `matcher_pool_grid.R`               | The **full factorial**: 3 matchers × 2 pool definitions, visualized natively per contrast |

All three resolve `../functions.R` automatically (run them as files), take the
same thresholds, and are dataset-agnostic — `--contrasts` is the only
dataset-specific argument.

## `method_comparison_report.R`

Compares the three control-gene matchers in `functions.R`, writes a multi-page
PDF report plus a summary table. The reproducible version of the matcher report.

The three matchers, all drawing controls from the **same** unchanged-gene pool:

| Method                     | Algorithm                                             | Covariates        |
|----------------------------|-------------------------------------------------------|-------------------|
| `expressionMatched (expr)` | legacy decile matcher (`expressionMatched()`)         | expression only   |
| `MatchIt (expr)`           | nearest-neighbour Mahalanobis (`matchitMatched()`)    | expression only   |
| `MatchIt (expr+len)`       | nearest-neighbour Mahalanobis (`matchitMatched()`)    | expression + gene length |

The point of the comparison: all three match **expression** to a standardized
mean difference (SMD) of ~0, but expression-only matching leaves **gene length**
unbalanced whenever the input set is length-skewed, and the legacy decile matcher
can fail outright on such a set. Only `MatchIt (expr+len)` balances both.

## Requirements

- R with `DESeq2` and `MatchIt` installed (same as the rest of the repo).
- The script sources `../functions.R` automatically (it resolves its own path),
  so run it as a file — `Rscript examples/method_comparison_report.R ...` — not
  by pasting into a REPL.

## Inputs

| Flag          | Required | Meaning |
|---------------|----------|---------|
| `--dds FILE`  | yes      | A `DESeqDataSet` saved with `saveRDS()`, already run through `DESeq()`, with a `condition` factor in `colData`. |
| `--gga FILE`  | yes      | The GGA `Dominant.TSS.TES` **WORKING_FILE** (tab-delimited, gene IDs in column 1), providing `GeneType`, `GeneName`, `DominantTSS`, `DominantTES`. Gene length = `|DominantTES - DominantTSS| + 1`. |
| `--contrasts` | yes      | Comma-separated `num/den` pairs over the `condition` levels, e.g. `tam/ctrl_003,siTTF2/ctrl_002`. **Baseline (`den`) expression** is the matching covariate. These are **dataset-specific** — list the treatment-vs-control comparisons you care about. |
| `--out DIR`   | no       | Output directory (default `matching_report_output`). |
| `--mech`      | no       | Which contrast drives the two "mechanism" pages (default: the last one in `--contrasts`). |

Unchanged-pool definition (orthogonal to the matcher — see note below):

| Flag              | Default   | Meaning |
|-------------------|-----------|---------|
| `--unch-method`   | `legacy`  | How the **unchanged pool** is defined for the main report: `legacy` (not significant AND inside the FC band) or `lessabs` (DESeq2 equivalence test — *confidently* inside the band). |
| `--sensitivity`   | *(off)*   | Also run the whole comparison under **both** pool definitions and append a sensitivity section + `sensitivity_table.csv`. |
| `--lessabs-alpha` | `0.1`     | `padj` cutoff for the `lessAbs` equivalence call. |

Thresholds (defaults match the repo's conventions):

| Flag           | Default           | Meaning |
|----------------|-------------------|---------|
| `--up-fc`      | `1.5`             | Differential call: linear fold change. Up = `>= up-fc`, Dn = `<= 1/up-fc`. |
| `--padj`       | `0.01`            | Differential call: adjusted p-value cutoff. |
| `--unch-fc`    | `1.3`             | **Unchanged** pool: linear FC band `[1/unch-fc, unch-fc]`; also the `lessAbs` equivalence bound. |
| `--unch-padj`  | `0.1`             | **Unchanged** pool (legacy only): `padj > unch-padj`. |
| `--min-len`    | `1000`            | Minimum gene length (bp) for eligibility. |
| `--genetype`   | `protein_coding`  | `GeneType` value kept. |
| `--deciles`    | `10`              | Number of expression bins for the legacy decile matcher. |

The **unchanged pool is defined once per contrast and shared by all three
matchers**, intersected with the eligibility filters and minus the Up/Dn sets.
The matchers differ only in *how* they pick controls from that identical pool —
never in which genes are eligible.

### `lessAbs` is a pool definition, not a matcher

`lessAbs` (`results(dds, lfcThreshold=log2(unch-fc), altHypothesis="lessAbs")`)
is DESeq2's built-in **equivalence test**: it identifies genes *confidently*
within the fold-change band (reject |LFC| ≥ bound ⇒ unchanged). That answers
"*which genes are unchanged?*" — the **pool-definition** stage — not "*which
control matches this input gene?*" (the matcher stage). So it is **not** a fourth
matcher; it is an alternative way to build the pool all three matchers draw from.
It is compatible with `contrast=` (a thresholded Wald test on the MLE — no
shrinkage, unlike `apeglm`).

Use `--unch-method lessabs` to build the whole report on that pool, or
`--sensitivity` to compare both. In practice the `lessAbs` pool is **itself
length-biased** — proving equivalence needs statistical power, which favours
longer / well-expressed genes — so expression-only matching tends to get *worse*
under it, and even `MatchIt (expr+len)` can fall short when the stricter pool is
depleted at the input's length range. The matcher ranking is unchanged; the
takeaway is that the unchanged *definition* matters alongside the matcher.

## Example

```bash
Rscript examples/method_comparison_report.R \
  --dds   path/to/Lewis1_2_3_deseq.dds.rds \
  --gga   path/to/Lewis123_GGA_Dominant.TSS.TES.calls_WORKING_FILE.txt \
  --contrasts "sto1/dmso,sto1_noSark/dmso_noSark,siTTF2/ctrl_002,tam/ctrl_003" \
  --mech "tam/ctrl_003" \
  --sensitivity \
  --out   matching_report_output
```

## Outputs

- `summary_table.csv` — one row per (contrast × direction × method): `n_input`,
  `n_matched`, `SMD_expr`, `SMD_length`, and a `length_balanced` flag
  (`yes` / `NO` / `FAILED` / `n<10`).
- `sensitivity_table.csv` — only with `--sensitivity`: the same metrics under
  both the `legacy` and `lessAbs` pools side by side, with an `agree` column.
- `matching_report.pdf`:
  1. **Summary table** — color-coded (`NO`/`FAILED` red, `n<10` orange); the
     header names the unchanged-pool definition in force.
  2. **"How MatchIt builds a control set"** — plain-language explanation plus
     representative paired loci from the `--mech` contrast.
  3. **Mechanism (Dn)** — joint expression × length map with each input gene
     linked to the control MatchIt paired it with (expr vs expr+len), and the
     per-pair length-mismatch ECDF.
  4. **All three methods (Up), same idiom** — `expressionMatched` as a point
     cloud (it returns a *set*, not pairs), both MatchIt modes with pair links,
     and a marginal gene-length ECDF of every control set vs the input.
  5. **Per-contrast panels** — marginal expression and length ECDFs for each
     contrast and direction (`all genes` / `input` / three matched sets).
  6. **Unchanged-pool sensitivity** (only with `--sensitivity`) — pool sizes
     and a per-matcher balance verdict under `legacy` vs `lessAbs`, highlighting
     where the conclusion changes.

## Reading the flags in the summary

- `yes` — gene length balanced (`|SMD_length| <= 0.1`).
- `NO`  — length imbalanced; the matched set does not represent the input's
  length distribution.
- `FAILED` — the matcher returned nothing. For `expressionMatched` this is the
  decile-starvation case: its 1:1 without-replacement draw runs out of unchanged
  pool genes in the input's dense expression deciles. MatchIt (nearest-neighbour
  over the whole pool) degrades gracefully instead of failing.
- `n<10` — too few input genes for a meaningful SMD (ignore the number).

## Notes

- `--contrasts` are intrinsically dataset-specific; there is no sensible default.
  List the treatment-vs-control comparisons that are biologically coherent for
  your design (same control, same batch).
- Whichever `--unch-method` is chosen, all three matchers share the *same* pool,
  so the matcher comparison is always apples-to-apples. `lessAbs` changes the
  pool definition, not the matcher (see the section above); `--sensitivity`
  shows how much the pool choice matters.
- To extend MatchIt with a third covariate (e.g. GC content), add the column to
  the covariate frame and pass it through `matchitMatched(..., covars = c("expr",
  "loglen", "gc"))`; the Mahalanobis distance generalizes to any number of
  numeric covariates.

## `unchanged_definition_compare.R`

Head-to-head of the two ways to **define** the unchanged control pool, rather
than which matcher to use:

- `legacy` — hard cutoff: `padj > unch-padj` AND linear FC in `[1/unch-fc, unch-fc]`
  ("not significant AND inside a band").
- `lessAbs` — DESeq2's equivalence test
  (`results(lfcThreshold=log2(unch-fc), altHypothesis="lessAbs")`, `padj <= lessabs-alpha`):
  genes the data can **prove** are inside the band.

```bash
Rscript examples/unchanged_definition_compare.R \
  --dds path/to/deseq.dds.rds --gga path/to/GGA_WORKING_FILE.txt \
  --contrasts "sto1/dmso,siTTF2/ctrl_002,tam/ctrl_003" \
  --out unchanged_compare_output
```

Outputs `pool_comparison.csv`, `downstream_length_balance.csv`, and a 4-page PDF:
(1) set sizes, overlap and covariate composition of each pool; (2–3) pool
expression/length ECDFs vs all-eligible genes; (4) the length balance each
matcher achieves from each pool. **Finding:** the two sets overlap heavily but
aren't nested — `lessAbs` is statistically cleaner (only provably-unchanged
genes) but, because proving equivalence needs power (which scales with gene
length/expression), its pool is biased toward **longer, higher-expressed** genes.
So it's a less representative background, and expression-only matching gets worse
drawing from it.

## `matcher_pool_grid.R`

The full factorial, visualized natively. For each contrast it draws a **3×2 grid**
of the joint expression × length scatter — rows = pool definition (legacy top,
lessAbs bottom), columns = matcher (expMat, MatchIt:expr, MatchIt:expr+len), with
MatchIt shown using its native **pair-links** — and a matching **3×2 grid of
gene-length ECDFs**. Page 1 is an overall length-balance summary table.

```bash
Rscript examples/matcher_pool_grid.R \
  --dds path/to/deseq.dds.rds --gga path/to/GGA_WORKING_FILE.txt \
  --contrasts "sto1/dmso,siTTF2/ctrl_002,tam/ctrl_003" \
  --dir up --out matcher_pool_grid_output
```

Use `--dir up|dn` to pick the direction shown in the grids (summary table always
covers both); `--bg N` controls how many background genes are drawn in each
scatter (lower it to shrink the PDF). Outputs `grid_summary.csv` and the PDF.
