# PRO-seq Expression Matching

A command-line wrapper around Ben Martin's `expressionMatched()` function
(`expressionMatched.R`), tailored for outputs from the NTC PRO-seq pipeline.
It uses **quantile matching** to build a control gene set whose expression
distribution closely mirrors a target gene set of interest.

## How it works

`main.R` takes a list of genes of interest plus standard NTC pipeline outputs
and:

1. Identifies **unchanged** (non-differentially-expressed) genes from the DESeq2
   results.
2. Restricts those control candidates to **protein-coding genes ≥ 1 kb** using
   the GGA annotation.
3. Computes a **length-normalized expression value** (counts/kb) for every gene
   from the baseline condition.
4. Runs `expressionMatched()` to draw, for each expression quantile of the
   target set, a matching number of control genes — so the returned control set
   has a similar expression distribution to the target set.

The matching itself (quantile binning + reproducible random sampling) is
unchanged from the Adelman lab original; see `expressionMatched.R`.

## Requirements

- R (version ≥ 4.x recommended)
- A plain-text file with a **single column of gene IDs of interest** (e.g.
  upregulated genes), one ENSG ID per line
- **DESeq2 outputs** from the NTC pipeline:
  - per-contrast `*_results_table.txt` files (under `<DESeqOutput>/*/`)
  - a `matrix/*.makeHeatmap.rds` count matrix
- A **GGA `WORKING_FILE`** annotation table
- A tab-delimited **DESeq2 sample manifest** containing a `condition` column
- Optional (for QC violin plots): the `vioplot` R package. If it is not
  installed, QC falls back to base-R box plots automatically.

## Usage

```sh
# Takes only a few seconds with minimal resources.
# Argument order must be exactly as below.
Rscript main.R \
  gene_ensg.txt \
  path/to/DESeqOutput \
  path/to/GGA.WORKING_FILE.txt \
  <window> \
  path/to/DESeq_metadata.txt \
  <baseline condition name>
```

### Required arguments (in order)

| Arg # | Example | Description |
|-------|---------|-------------|
| 1 | `gene_ensg.txt` | 1-column `.txt` of ENSG IDs for the genes of interest (e.g. up- or down-regulated genes) |
| 2 | `DESeqOutput/` | Path to the DESeq2 output folder (globbed for `*/*_results_table.txt` and `matrix/*.makeHeatmap.rds`) |
| 3 | `gga_WORKING_FILE.tsv` | Path to the GGA file with `WORKING_FILE` in its name (currently requires GGA-style annotation) |
| 4 | `<window>` | Integer bp subtracted from each gene's TSS–TES length when computing length-normalized expression (e.g. `0` matches on TSS→TES; `250` matches on +250→TES) |
| 5 | `DESeq_metadata.txt` | Tab-delimited DESeq2 sample manifest; must include a `condition` column |
| 6 | `<baseline condition>` | A value from the manifest's `condition` column whose samples are used as the baseline for expression matching (normally the WT/control, to avoid matching on a perturbed condition) |

## Output

For each of three random seeds (1, 2, 3), `main.R` writes:

- `expression_matched-<baseline condition>_<seed>.txt` — the matched control
  gene list (one gene ID per line).
- `expression_matched-<baseline condition>_<seed>_QC.txt` — a QC report (see
  below).
- `expression_matched-<baseline condition>_<seed>_violin.png` — a QC violin
  plot (see below).

The three lists are **replicates**: independent random draws produced with
different seeds. Each list has the *same* intended properties — it is (1) not
differentially expressed, (2) approximately the same size as the input set, and
(3) matched to the target set's distribution of length-normalized expression.
Generating three lets you check that the matching is robust to the random draw;
it is suggested to inspect all three. Downstream analysis typically subsets
matrix files using the returned gene lists.

## QC workflow

The QC step is adapted from the Adelman lab README and runs automatically for
each seed (implemented in `qc.R`). For each matched list it compares three gene
sets — **control** (all candidate unchanged genes), **target** (your genes of
interest), and **expression-matched** (the drawn control subset) — using their
length-normalized expression values:

- **Summary statistics** (`summary()`) for each of the three sets.
- **Mann-Whitney / Wilcoxon rank-sum tests**:
  - control vs. target — typically *significant* (the full control set is not
    expression-matched), and
  - expression-matched vs. target — ideally *not significant* (a well-matched
    set should be statistically indistinguishable from the target).
- A **violin plot** of `log2` expression for control, expression-matched, and
  target sets, saved as a PNG.

Interpreting the QC: good matching looks like a non-significant
expression-matched-vs-target test and overlapping violins/summaries for the
expression-matched and target sets. If the expression-matched-vs-target test is
significant, the matching is poor — see Troubleshooting.

## Additional details

- **Unchanged-gene thresholds** (hardcoded in `main.R` as `p_thresh` and
  `fc_thresh`): `padj > 0.1` **and** fold change within `[1/1.3, 1.3]` (i.e. the
  linear fold-change column, not log2). A gene must satisfy this in **every**
  contrast (intersection across all `*_results_table.txt` files). Note: the
  thresholds are applied to fixed result-table **columns** (column 6 = `padj`,
  column 7 = fold change), so the DESeq2 tables must follow that column order.
- **Control-gene filtering:** control candidates are further restricted to
  `GeneType == "protein_coding"` and to genes whose annotated length
  (`DominantTES − DominantTSS + 1`) is **≥ 1 kb**. This length filter uses the
  full gene length and is **independent of the `<window>` argument** (the window
  only affects the counts/kb normalization). The **genes of interest are used as
  supplied** and are *not* subjected to the protein-coding or length filters.
- **Expression values:** mean length-normalized counts/kb across the baseline
  condition's replicates, taken from the `makeHeatmap.rds` matrix and divided by
  `(gene length + 1 − <window>) / 1000`.
- **Quantile groups:** 10, set in the `expressionMatched()` call in `main.R`.
  Adjustable only by editing the script.
- **Replicates:** 3 seeds (1, 2, 3).

## Interactive use of `expressionMatched()`

The `expressionMatched()` function in `expressionMatched.R` can also be sourced
and called directly in an interactive R session if you want full control over
the matching parameters (`numQuantGroups`, `CtrlSetRatio`, `setSeedNum`). See
the [Adelman lab's Expression-Matching
README](https://github.com/AdelmanLab/Expression-Matching) for the detailed
function documentation and a self-contained simulated example.

## Troubleshooting

If matching fails with an error such as:

```
Error in sample.int(length(x), size, replace, prob) :
  cannot take a sample larger than the population when 'replace = FALSE'
```

your control set does not contain enough genes to cover every expression
quantile of the target set. Reduce the number of quantile groups
(`numQuantGroups`, in the `expressionMatched()` call in `main.R`), or lower the
control-set ratio (`CtrlSetRatio`) when calling `expressionMatched()` directly.

## Customization (make your own clone/copy of the script)

- Change thresholds, quantile-group count, or seed count by editing `main.R`.
- Adapt the column mappings if using non-standard gene annotations or DESeq2
  table layouts.

## Credit

- **Expression-matching algorithm:** Benjamin J. E. Martin (Karen Adelman Lab) —
  see <https://github.com/AdelmanLab/Expression-Matching>.
- **NTC PRO-seq wrapper (`main.R`):** Geoff Nelson.
- **QC workflow:** adapted from the Adelman lab README.
