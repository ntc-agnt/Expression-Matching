# qc.R
#
# Quality-control for expression-matched gene sets.
#
# Adapted from the Adelman lab Expression-Matching README workflow
# (https://github.com/AdelmanLab/Expression-Matching), which recommends QC via:
#   - summary statistics for the control / target / expression-matched sets
#   - Wilcoxon rank-sum (Mann-Whitney) tests
#   - a violin plot of log2 expression values
#
# expressionMatchQC() operates on the same length-normalized expression values
# used for matching, writing a text report and a violin-plot PNG per call.

expressionMatchQC <- function(geneset, ctrlgeneset, matchedset, dat,
                              out_prefix, plot_title = NULL) {
  geneIDs <- as.character(dat[, 1])
  expr    <- as.numeric(dat[, 2])

  target_values  <- expr[geneIDs %in% as.character(geneset)]
  control_values <- expr[geneIDs %in% as.character(ctrlgeneset)]
  matched_values <- expr[geneIDs %in% as.character(matchedset)]

  ## ---- text report: summary stats + Wilcoxon tests ----
  safe_test <- function(a, b) {
    tryCatch(
      capture.output(print(suppressWarnings(wilcox.test(a, b)))),
      error = function(e) paste("wilcox.test failed:", conditionMessage(e))
    )
  }

  report_path <- paste0(out_prefix, "_QC.txt")
  report <- c(
    paste0("Expression-matching QC report: ", basename(out_prefix)),
    if (!is.null(plot_title)) paste0("Label: ", plot_title),
    paste0("Generated: ", format(Sys.time())),
    "",
    paste0("Set sizes: target=", length(target_values),
           ", control=", length(control_values),
           ", expression-matched=", length(matched_values)),
    "",
    "## Summary statistics (length-normalized expression, counts/kb)",
    "",
    "Control set (all candidate unchanged genes):",
    capture.output(print(summary(control_values))),
    "",
    "Target set (genes of interest):",
    capture.output(print(summary(target_values))),
    "",
    "Expression-matched set:",
    capture.output(print(summary(matched_values))),
    "",
    "## Wilcoxon rank-sum (Mann-Whitney) tests",
    "",
    "Control vs. target (expected: significant, since the full control",
    "set is not expression-matched):",
    safe_test(control_values, target_values),
    "",
    "Expression-matched vs. target (expected: NOT significant if the",
    "matching is good):",
    safe_test(matched_values, target_values)
  )
  writeLines(report, report_path)

  ## ---- violin plot of log2 expression ----
  # Guard against non-finite values (log2 of zero or negative counts).
  safe_log2 <- function(x) {
    x <- log2(x)
    x[is.finite(x)]
  }
  c_l <- safe_log2(control_values)
  m_l <- safe_log2(matched_values)
  t_l <- safe_log2(target_values)

  plot_path  <- paste0(out_prefix, "_violin.png")
  main_title <- if (is.null(plot_title)) "Expression matching QC" else plot_title

  png(plot_path, width = 900, height = 900, res = 150)
  on.exit(dev.off(), add = TRUE)
  tryCatch({
    if (requireNamespace("vioplot", quietly = TRUE)) {
      vioplot::vioplot(
        c_l, m_l, t_l,
        names   = c("Control", "Exp-matched", "Target"),
        col     = "white",
        border  = c("gray", "dodgerblue3", "red2"),
        rectCol = "white",
        lineCol = c("gray", "dodgerblue3", "red2"),
        colMed  = c("gray", "dodgerblue3", "red2"),
        plotCentre = "line", las = 1
      )
      title(main = main_title, ylab = "log2 expression (counts/kb)")
    } else {
      message("Package 'vioplot' not installed; falling back to boxplot for QC.")
      boxplot(
        list(Control = c_l, `Exp-matched` = m_l, Target = t_l),
        col    = "white",
        border = c("gray", "dodgerblue3", "red2"),
        main   = main_title, ylab = "log2 expression (counts/kb)", las = 1
      )
    }
  }, error = function(e) {
    message("QC plot failed: ", conditionMessage(e))
  })

  invisible(list(
    report = report_path,
    plot   = plot_path,
    control_values = control_values,
    target_values  = target_values,
    matched_values = matched_values
  ))
}
