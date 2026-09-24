#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
outfile = if (length(args)) args[[1L]] else "replication-sessionInfo.txt"
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

package_names = c(
  "xplaineff", "data.table", "ggplot2", "scales", "ranger", "reticulate",
  "pkgload", "knitr", "rmarkdown", "testthat"
)
installed = utils::installed.packages()
available = intersect(package_names, rownames(installed))
package_table = if (length(available)) {
  installed[available, c("Package", "Version"), drop = FALSE]
} else {
  data.frame(Package = character(), Version = character())
}

lines = c(
  "xplaineff replication environment",
  paste0("Recorded at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Selected R package versions:",
  capture.output(print(package_table, row.names = FALSE)),
  "",
  capture.output(sessionInfo())
)
writeLines(lines, outfile)
message("Wrote ", outfile)
