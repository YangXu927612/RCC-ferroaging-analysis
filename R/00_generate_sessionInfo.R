# ==============================================================================
# 00_generate_sessionInfo.R
# Capture the full R session information (R version, platform, OS, attached
# packages and their versions) and write it to sessionInfo.txt in the current
# working directory.
#
# Recommended usage:
#   1. Set the working directory to the repository root.
#   2. Run scripts 01 -> 19.
#   3. After all analyses complete, run this script LAST to record the R
#      environment that produced every figure / table.
#   4. Commit the generated sessionInfo.txt to GitHub.
#
# This satisfies the reproducibility requirement of journals such as
# PLOS ONE under their Code Sharing policy.
# ==============================================================================

# ---- 1. Capture full session info ----
session_text <- capture.output(sessionInfo())

# Append timestamp (so reviewers know when the file was generated)
stamp <- format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S UTC")
session_text <- c(
  paste0("# Generated: ", stamp),
  paste0("# R script: 00_generate_sessionInfo.R"),
  "# ---------------------------------------------------------------------",
  session_text
)

# ---- 2. Write to sessionInfo.txt in the working directory ----
out_path <- file.path(getwd(), "sessionInfo.txt")
writeLines(session_text, out_path)

# ---- 3. Also list CRAN / Bioconductor versions of key packages ----
# (Optional but useful for reviewers who want to know exactly which package
#  versions were used; sessionInfo() already reports attached packages,
#  but some loading scripts install packages without attaching them.)
key_pkgs <- c(
  # CRAN
  "ggplot2", "pheatmap", "ggpubr", "ggrepel", "ggsci", "ggsignif",
  "ggvenn", "ComplexHeatmap", "circlize", "patchwork", "RColorBrewer",
  "dplyr", "tidyverse", "reshape2", "readxl", "openxlsx", "DT",
  "e1071", "randomForest", "rms", "regplot", "rmda", "pROC",
  "MCPcounter", "Seurat", "SingleR", "celldex", "DoubletFinder",
  "celda", "WGCNA", "clusterProfiler", "org.Hs.eg.db", "GOplot",
  "scRNAtoolVis", "plot1cell", "Matrix", "parallel", "preprocessCore",
  "viridis", "scales", "broom", "corrplot", "cowplot", "tools",
  "limma", "data.table",
  # Bioconductor
  "limma", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db",
  "WGCNA", "SingleR", "celldex", "MCPcounter"
)
key_pkgs <- unique(key_pkgs)
pkg_versions <- sapply(key_pkgs, function(p) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  paste0(p, "_", v)
})

extra_block <- c(
  "",
  "# ---------------------------------------------------------------------",
  "# Key package versions (NOT INSTALLED means the package was not detected)",
  "# ---------------------------------------------------------------------",
  pkg_versions
)

writeLines(c(session_text, extra_block), out_path)

message("Wrote session information to: ", out_path)
message("\nFirst 20 lines:\n")
cat(paste(readLines(out_path, n = 20), collapse = "\n"), "\n")
