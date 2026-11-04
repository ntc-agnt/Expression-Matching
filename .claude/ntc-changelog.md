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

## 2026-11-04 · claude/matchit-matcher
- Add matchit_match.R: reusable flagged MatchIt matching workflow with a --covars spec (default expression; add length / any numeric GGA column / --covar-file column, with log2:/log10: transforms), sourcing covariates from makeHeatmap.rds + GGA WORKING_FILE + optional --covar-file, and emitting the full match-quality QC. Supersedes and removes the positional prototype main_matchit.R. Experimental branch, not merged; affects results: no

## 2026-11-04 · claude/matchit-matcher
- Add benchmark_matchers.R: synthetic stress-test suite comparing matchers (Expression / MatchIt / optional matchRanges) across confound scenarios, emitting a covariate-balance scorecard + per-scenario matchQC; matchQC.R gains density + ECDF views. Experimental branch, not merged; affects results: no

## 2026-11-03 · claude/matchit-matcher
- Prototype: MatchIt-based covariate matcher (joint expression+length nearest-neighbour) as an alternative to decile expression matching, with QC (balance table + distribution plots); adds matchitMatched.R, matchQC.R, main_matchit.R; experimental branch, not merged; affects results: no
