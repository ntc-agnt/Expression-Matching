# Expression-Matching — revision changelog (for the NTC ecosystem manager)

<!-- THIS FILE IS OWNED BY THIS REPO. The ecosystem manager (ntc-ecosystem-context) READS it to
track pipeline revisions and keep its overseer context current; it does NOT overwrite it.
Append a new entry at the TOP whenever you make a notable change (a revision / version bump /
results-affecting fix), so the manager has visibility. Keep entries to one or two lines. -->

Entry format (newest first):
```
## <YYYY-MM-DD> · <version/tag or claude/branch> · <NTC-### if any>
- <one-line summary>; affects results: <yes|no>
```

<!-- entries below, newest first -->

## 2026-11-04 · main
- matchit_match.R: add a gene-list sanity guard - report how many --genes IDs are found in the annotation vs dropped, and stop with a clear message if <50% match (wrong ID type or stray header). affects results: no

## 2026-11-04 · main
- Demote the legacy decile driver to archive/main.R (fixed its functions.R source path to the repo root). Updated README usage + overview. Algorithm/outputs unchanged. affects results: no

## 2026-11-04 · main
- Merge the function libraries into a single functions.R (expressionMatched, matchitMatched, matchQC); main.R, matchit_match.R and tests/benchmark_matchers.R now source it; removed expressionMatched.R / matchitMatched.R / matchQC.R. Resolves the matchit_match.R vs matchitMatched.R naming confusion. main.R change is the source path only; outputs identical. affects results: no

## 2026-11-04 · main
- Consolidate the MatchIt matcher work onto main: matchit_match.R (covariate matching workflow, --method default matchit, legacy expression mode), matchitMatched.R, matchQC.R. Moved benchmark_matchers.R and the legacy qc.R to tests/. Applied the corrected README and documented the new workflow. main.R unchanged (frozen legacy driver). affects results: no

## 2026-11-04 · claude/matchit-matcher
- matchit_match.R: add --method flag (default 'matchit'; 'expression'/'legacy' = original decile matcher, expression-only) so legacy control sets remain reproducible; outputs stamped with method (matched_controls.<method>.txt, run_info, QC) and shared QC exposes what legacy leaves unbalanced. Experimental branch, not merged; affects results: no

## 2026-11-04 · claude/matchit-matcher
- Add matchit_match.R: reusable flagged MatchIt matching workflow with a --covars spec (default expression; add length / any numeric GGA column / --covar-file column, with log2:/log10: transforms), sourcing covariates from makeHeatmap.rds + GGA WORKING_FILE + optional --covar-file, and emitting the full match-quality QC. Supersedes and removes the positional prototype main_matchit.R. Experimental branch, not merged; affects results: no

## 2026-11-04 · claude/matchit-matcher
- Add benchmark_matchers.R: synthetic stress-test suite comparing matchers (Expression / MatchIt / optional matchRanges) across confound scenarios, emitting a covariate-balance scorecard + per-scenario matchQC; matchQC.R gains density + ECDF views. Experimental branch, not merged; affects results: no

## 2026-11-03 · claude/matchit-matcher
- Prototype: MatchIt-based covariate matcher (joint expression+length nearest-neighbour) as an alternative to decile expression matching, with QC (balance table + distribution plots); adds matchitMatched.R, matchQC.R, main_matchit.R; experimental branch, not merged; affects results: no
