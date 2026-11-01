library(readr)

# This catalog is metadata only: do not bind it into build_gene_scores().
catalog <- read_csv(
  "data/metadata/additional_screen_candidates.csv",
  col_types = cols(.default = col_character())
)
stopifnot(!anyDuplicated(catalog$candidate_id), all(nzchar(catalog$source_url)))
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(catalog, "data/processed/additional_screen_candidates.rds")
message("Additional-screen catalog: ", nrow(catalog), " study entries.")
