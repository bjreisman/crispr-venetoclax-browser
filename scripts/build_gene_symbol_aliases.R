#!/usr/bin/env Rscript
# Build data/metadata/gene_symbol_aliases.csv: withdrawn CRISPR-library gene symbols
# mapped to their current HGNC symbols.
#
# Why this exists. The screens in this browser were published between 2019 and 2026 and
# use the library symbols of their day (MARCH5, TAZ, GCN1L1, ATP5B, SEPT9, ...). DepMap
# 26Q1 and other current reference data use the 2020+ HGNC symbols (MARCHF5, TAFAZZIN,
# GCN1, ATP5F1B, SEPTIN9). Without a mapping the two never join, so the Gene Explorer
# showed no baseline dependency row for genes as central as MARCH5.
#
# What is kept. Only HGNC "previous symbol" entries that are safe to apply blindly:
#   1. the old symbol maps to exactly one approved symbol,
#   2. the old symbol is not itself an approved symbol for some other gene,
#   3. the new symbol does not already appear as a separate gene in the screen data
#      (otherwise two library entries would silently merge into one row), and
#   4. no two old symbols in the screen data map onto the same new symbol.
# Rules 3 and 4 are checked against the published gene symbols in
# data/processed/gene_scores.rds, so build that cache first, re-run this script after
# adding screens, and review the diff.
#
# Usage:
#   Rscript scripts/build_gene_symbol_aliases.R [path/to/hgnc_complete_set.txt]
# With no argument the current HGNC release is downloaded to data/external/hgnc/.

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

HGNC_URL <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
out_path <- file.path("data", "metadata", "gene_symbol_aliases.csv")
scores_path <- file.path("data", "processed", "gene_scores.rds")

args <- commandArgs(trailingOnly = TRUE)
hgnc_path <- if (length(args)) args[[1]] else file.path("data", "external", "hgnc", "hgnc_complete_set.txt")
if (!file.exists(hgnc_path)) {
  dir.create(dirname(hgnc_path), recursive = TRUE, showWarnings = FALSE)
  message("Downloading HGNC complete set to ", hgnc_path)
  utils::download.file(HGNC_URL, hgnc_path, quiet = TRUE)
}
if (!file.exists(scores_path)) {
  stop("Build data/processed/gene_scores.rds first; the safety checks compare against it.")
}

normalize_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

hgnc <- read_tsv(hgnc_path, col_types = cols(.default = "c"), progress = FALSE) |>
  filter(.data$status == "Approved") |>
  select("symbol", "prev_symbol")

approved <- unique(normalize_key(hgnc$symbol))

renames <- hgnc |>
  filter(!is.na(.data$prev_symbol)) |>
  separate_rows("prev_symbol", sep = "[|]") |>
  transmute(from = normalize_key(.data$prev_symbol), to = normalize_key(.data$symbol)) |>
  filter(.data$from != "", .data$from != .data$to) |>
  distinct()

ambiguous <- renames |> count(.data$from) |> filter(.data$n > 1) |> pull("from")
renames <- renames |> filter(!.data$from %in% ambiguous, !.data$from %in% approved)

# Compare against the symbols as published, not gene_key, so that re-running this script
# after the aliases are already in force still sees the original library symbols.
screen_keys <- unique(normalize_key(readRDS(scores_path)$gene_symbol))
renames <- renames |> filter(.data$from %in% screen_keys, !.data$to %in% screen_keys)

collapsing <- renames |> count(.data$to) |> filter(.data$n > 1) |> pull("to")
renames <- renames |>
  filter(!.data$to %in% collapsing) |>
  arrange(.data$from)

write_csv(renames, out_path)
message(sprintf(
  "Wrote %s: %d withdrawn symbols mapped (dropped %d ambiguous, %d many-to-one).",
  out_path, nrow(renames), length(ambiguous), length(collapsing)
))
