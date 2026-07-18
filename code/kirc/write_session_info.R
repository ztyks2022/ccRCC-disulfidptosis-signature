#!/usr/bin/env Rscript

# Write a reproducibility snapshot for the analysis repository.
# The script scans R files for package calls and also includes packages named
# in the manuscript Methods but not necessarily invoked in the public scripts.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  file.path(getwd(), "code", "write_session_info.R")
}

repo_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_dir, "code"))) {
  repo_dir <- normalizePath(getwd(), mustWork = TRUE)
}

code_dir <- file.path(repo_dir, "code")
out_file <- file.path(repo_dir, "sessionInfo.txt")
r_files <- list.files(code_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
contents <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
text <- paste(contents, collapse = "\n")

extract_all <- function(pattern, x) {
  m <- gregexpr(pattern, x, perl = TRUE)
  hits <- regmatches(x, m)[[1]]
  if (length(hits) == 1 && hits[1] == "") character() else hits
}

lib_hits <- extract_all("\\b(?:library|require)\\s*\\(\\s*['\"]?([A-Za-z][A-Za-z0-9_.]*)['\"]?", text)
lib_pkgs <- sub(".*\\(\\s*['\"]?([A-Za-z][A-Za-z0-9_.]*).*", "\\1", lib_hits)

ns_hits <- extract_all("\\brequireNamespace\\s*\\(\\s*['\"]([A-Za-z][A-Za-z0-9_.]*)['\"]", text)
ns_pkgs <- sub(".*\\(\\s*['\"]([A-Za-z][A-Za-z0-9_.]*)['\"].*", "\\1", ns_hits)

colon_hits <- extract_all("\\b([A-Za-z][A-Za-z0-9_.]*)::[A-Za-z_.][A-Za-z0-9_.]*", text)
colon_pkgs <- sub("::.*", "", colon_hits)

# Dependencies named in Methods/results or required by figure-generation paths
# even when the current public scripts do not call them directly.
manuscript_declared_pkgs <- c(
  "CellChat",
  "ConsensusClusterPlus",
  "TCGAbiolinks",
  "maftools",
  "xgboost",
  "ranger",
  "gbm",
  "SHAPforxgboost",
  "ggpubr",
  "pROC",
  "SeuratObject"
)

pkgs <- sort(unique(c(lib_pkgs, ns_pkgs, colon_pkgs, manuscript_declared_pkgs)))
pkgs <- pkgs[pkgs != ""]
pkgs <- setdiff(pkgs, c("pkg"))

installed <- as.data.frame(
  installed.packages()[, c("Package", "Version", "LibPath")],
  stringsAsFactors = FALSE
)

pkg_table <- data.frame(
  Package = pkgs,
  Version = ifelse(
    pkgs %in% installed$Package,
    installed$Version[match(pkgs, installed$Package)],
    NA_character_
  ),
  Installed = pkgs %in% installed$Package,
  stringsAsFactors = FALSE
)

# Keep only real installed packages in the main table, and separately report any
# referenced packages missing from the current R library.
pkg_table <- pkg_table[order(tolower(pkg_table$Package)), ]
installed_table <- pkg_table[pkg_table$Installed, ]
missing_table <- pkg_table[!pkg_table$Installed, ]

con <- file(out_file, open = "wt")
on.exit(close(con), add = TRUE)

writeLines("Reproducibility snapshot for KIRC disulfidptosis analysis", con)
writeLines(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), con)
writeLines(paste0("Repository: ", normalizePath(repo_dir, mustWork = FALSE)), con)
writeLines("", con)

writeLines("R sessionInfo()", con)
writeLines("===============", con)
writeLines(paste0("R version ", getRversion(), " (", R.version$date, ")"), con)
writeLines(paste0("Platform: ", R.version$platform), con)
writeLines(paste0("Running under: ", paste(Sys.info()[c("sysname", "release")], collapse = " ")), con)
writeLines("", con)
writeLines(paste0("time zone: ", Sys.timezone()), con)
writeLines("", con)

writeLines(
  paste0(
    "Packages used by code/*.R and manuscript-declared analyses (library, ::, requireNamespace calls; ",
    nrow(installed_table),
    " installed packages)"
  ),
  con
)
writeLines("========================================================================", con)
write.table(installed_table, file = con, quote = FALSE, sep = "\t", row.names = FALSE)
writeLines("", con)

if (nrow(missing_table)) {
  writeLines("Packages referenced but not installed in this R library", con)
  writeLines("=======================================================", con)
  write.table(missing_table, file = con, quote = FALSE, sep = "\t", row.names = FALSE)
  writeLines("", con)
}

writeLines("Note", con)
writeLines("====", con)
writeLines(paste0("Package versions were read from the current R environment (R version ", getRversion(), ")."), con)
writeLines("The list covers packages invoked in code/*.R via library(), require(), requireNamespace(), or pkg::fun(), plus manuscript-declared analysis packages that may be called from local/private scripts.", con)
writeLines("Regenerate with:", con)
writeLines("  Rscript code/kirc/write_session_info.R", con)

message("Wrote ", out_file)
