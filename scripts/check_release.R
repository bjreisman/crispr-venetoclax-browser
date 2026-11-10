# Check prepared public files without rebuilding or requiring raw supplements.
source("R/load_data.R")

cache_manifest <- readr::read_csv("data/metadata/cache_manifest.csv", show_col_types = FALSE)
expected_caches <- file.path("data/processed", c(
  "gene_scores.rds", "depmap_crispr.rds", "external_context.rds", "additional_screen_candidates.rds"
))
stopifnot(setequal(cache_manifest$path, expected_caches), !anyDuplicated(cache_manifest$path))
for (i in seq_len(nrow(cache_manifest))) {
  path <- cache_manifest$path[i]
  stopifnot(file.exists(path), file.info(path)$size == cache_manifest$bytes[i])
  if (!identical(digest::digest(file = path, algo = "sha256"), cache_manifest$sha256[i])) {
    stop("Cache checksum mismatch: ", path, ". Regenerate the manifest after an intentional rebuild.")
  }
}

scores <- load_browser_data()
stopifnot(all(standard_columns %in% names(scores)), nrow(scores) > 0L)
stopifnot(dplyr::n_distinct(scores$paper) == 16L)
stopifnot(all(c("BAX", "TP53", "PMAIP1") %in% scores$gene_key))
depmap <- load_depmap_crispr()
stopifnot(nrow(depmap) > 0L, all(c("MOLM-13", "MOLM-14", "MV4-11", "OCI-AML2", "THP-1", "OCI-LY-19") %in% depmap$cell_line))
context <- load_external_context()
stopifnot(nrow(context$drug_response) > 0L, ncol(context$beataml_expression$expression) > 0L)
catalog <- readRDS("data/processed/additional_screen_candidates.rds")
catalog_source <- readr::read_csv("data/metadata/additional_screen_candidates.csv", col_types = readr::cols(.default = "c"))
stopifnot(nrow(catalog) == 29L, !anyDuplicated(catalog$candidate_id))
stopifnot(identical(as.data.frame(catalog), as.data.frame(catalog_source)))

# A source-code timestamp newer than a prepared cache must not trigger a rebuild.
# Missing caches must report an actionable error, not silently create partial data.
check_cache_only_loading <- function() {
  fixture <- tempfile("cache-loader-")
  dir.create(file.path(fixture, "data", "processed"), recursive = TRUE)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  saved <- list(marker = "prepared snapshot")
  for (name in c("gene_scores.rds", "external_context.rds")) {
    file <- file.path(fixture, "data", "processed", name)
    saveRDS(saved, file)
    Sys.setFileTime(file, as.POSIXct("2000-01-01", tz = "UTC"))
  }
  env <- new.env(parent = globalenv())
  env$load_browser_data <- load_browser_data
  env$load_external_context <- load_external_context
  environment(env$load_browser_data) <- env
  environment(env$load_external_context) <- env
  env$write_processed_data <- function(...) stop("Unexpected source rebuild")
  env$write_external_context <- function(...) stop("Unexpected context rebuild")
  stopifnot(identical(env$load_browser_data(fixture), saved))
  stopifnot(identical(env$load_external_context(fixture), saved))
  unlink(file.path(fixture, "data", "processed", "gene_scores.rds"))
  unlink(file.path(fixture, "data", "processed", "external_context.rds"))
  for (loader in list(env$load_browser_data, env$load_external_context)) {
    error <- tryCatch({loader(fixture); NULL}, error = identity)
    stopifnot(inherits(error, "error"), grepl("Restore the prepared data", conditionMessage(error)))
  }
  stopifnot(length(list.files(file.path(fixture, "data", "processed"))) == 0L)
}
check_cache_only_loading()

app_env <- new.env(parent = globalenv())
app <- source("app.R", local = app_env)$value
stopifnot(inherits(app, "shiny.appobj"))
stopifnot(all(c("ui", "server") %in% ls(app_env)))

shiny::testServer(app_env$server, {
  session$setInputs(
    gene_select = c("BAX", "TP53"), gene_papers = unique(app_env$scores$paper),
    gene_significant_only = FALSE, effect_axis_type = "raw",
    depmap_cell_lines = c("MOLM-13", "THP-1", "OCI-LY-19"),
    overview_papers = unique(app_env$scores$paper), overview_direction = "All",
    context_datasets = "BeatAML2", context_drugs = "Venetoclax",
    context_expression_mode = "gene", context_expression_gene = "BCL2"
  )
  stopifnot(nrow(gene_view()) > 0L, all(gene_view()$gene_key %in% c("BAX", "TP53")))
  stopifnot(nrow(gene_effect_display()) > 0L, nrow(depmap_view()) > 0L)
  stopifnot(nrow(overview_hits()) > 0L, nrow(expression_auc_view()) > 0L)
  session$setInputs(gene_view_mode = "screens", screen_extra_groups = baseline_screen_group)
  stopifnot(nrow(gene_screen_view()) > 0L, all(gene_screen_view()$screen_group %in% c(venetoclax_screen_groups, baseline_screen_group)))
  stopifnot(!anyDuplicated(gene_screen_view()[c("screen_id", "gene_key")]))
  # Baseline dependency rows are on by default and sit apart from the screens.
  stopifnot(baseline_screen_group %in% gene_screen_view()$screen_group)
  session$setInputs(gene_select = c("MCL1", "BAX"))
  base <- gene_screen_view()[gene_screen_view()$baseline, ]
  stopifnot(all(base$percentile[base$gene_key == "MCL1" & base$paper == "DepMap"] < -0.9), all(base$supported[base$gene_key == "MCL1"]))
  session$setInputs(gene_select = c("BAX", "TP53"))
  session$setInputs(screen_extra_groups = "Other BH3 mimetics")
  stopifnot("Other BH3 mimetics" %in% gene_screen_view()$screen_group, !any(gene_screen_view()$baseline))
  for (metric in names(screen_metrics)) {
    session$setInputs(screen_metric = metric)
    stopifnot(nzchar(output$gene_screen_plot))
  }
  # Expanding a paper adds its timepoints/replicates as uniquely keyed rows.
  collapsed_rows <- nrow(gene_screen_view())
  session$setInputs(screen_metric = "percentile", screen_expand_papers = c("Sharon2019", "Chen2019", "Guièze2019"))
  expanded <- gene_screen_view()
  stopifnot(nrow(expanded) > collapsed_rows, any(expanded$row_level == 1L))
  stopifnot(all(c("↳ day 8", "↳ day 29") %in% expanded$short_label[expanded$paper == "Sharon2019"]))
  stopifnot(!anyDuplicated(screen_panel_rows(expanded, app_env$expandable_papers, input$screen_expand_papers, app_env$screen_detail_parents)$key))
  stopifnot(nzchar(output$gene_screen_plot))
  # ▸/▾ screen labels carry their paper for click-to-expand, and the plot
  # binds the label click handler; children and unexpandable screens do not.
  rows <- screen_panel_rows(gene_screen_view(), app_env$expandable_papers, input$screen_expand_papers, app_env$screen_detail_parents)
  stopifnot(identical(unique(stats::na.omit(rows$toggle_paper[grepl("Sharon 2019", rows$label)])), "Sharon2019"))
  headers <- startsWith(rows$key, "header :: ")
  stopifnot(all(is.na(rows$toggle_paper[!headers & !grepl("^[▸▾] ", rows$label)])), all(rows$toggle_paper[!headers & grepl("^[▸▾] ", rows$label)] %in% app_env$expandable_papers))
  stopifnot(all(rows$toggle_paper[headers] == paste0("group::", sub("^header :: ", "", rows$key[headers]))))
  # Folding a group keeps its heading, drops its screens, and still plots.
  folded <- screen_panel_rows(gene_screen_view(), app_env$expandable_papers, input$screen_expand_papers, app_env$screen_detail_parents, "Venetoclax · AML")
  stopifnot("header :: Venetoclax · AML" %in% folded$key, nrow(folded) < nrow(rows), grepl("screens\\)", folded$label[folded$key == "header :: Venetoclax · AML"]))
  stopifnot(!any(c("chen_molm13_ven", "sharon_molm13r1_ven") %in% folded$key))
  session$setInputs(screen_label_click = "group::Venetoclax · AML")
  stopifnot(identical(collapsed_screen_groups(), "Venetoclax · AML"), nzchar(output$gene_screen_plot))
  session$setInputs(screen_label_click = "group::Venetoclax · AML")
  stopifnot(length(collapsed_screen_groups()) == 0L)
  stopifnot(grepl("screen_label_click", output$gene_screen_plot, fixed = TRUE))
  session$setInputs(screen_label_click = "Sharon2019")
})
# Paper Explorer: each paper opens on a Gene Explorer representative result
# where one exists, every screen button and timepoint chip selects a result,
# and each result renders its plot, summary, and top-hit tables.
shiny::testServer(app_env$server, {
  session$setInputs(paper_metric = "robust_z", paper_top_n = 10, paper_highlight_genes = "BCL2")
  for (paper in unique(app_env$scores$paper)) {
    session$setInputs(paper_select = paper)
    rows <- app_env$paper_results[app_env$paper_results$paper == paper, ]
    stopifnot(!any(rows$main) || rows$main[match(paper_result_selected(), rows$result_id)])
    stopifnot(nzchar(output$paper_screen_buttons$html))
    buttons <- paper_screen_buttons(app_env$paper_results, paper)
    ids <- unlist(lapply(unlist(buttons, recursive = FALSE), function(screen) screen$results$result_id))
    stopifnot(setequal(ids, unique(rows$result_id)))
    for (result in ids) {
      session$setInputs(paper_result_click = result)
      stopifnot(identical(paper_result_selected(), result), nrow(paper_view()) > 0L, !anyDuplicated(paper_view()$gene_key))
      stopifnot(nzchar(output$paper_rank_plot), nzchar(output$paper_hits_resistance), nzchar(output$paper_summary$html))
    }
  }
  # A result from another paper is ignored.
  session$setInputs(paper_select = "Chen2019")
  session$setInputs(paper_result_click = "Chen2019 :: Venetoclax vs DMSO day 16")
  session$setInputs(paper_result_click = "Sharon2019 :: MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day 8")
  stopifnot(identical(paper_result_selected(), "Chen2019 :: Venetoclax vs DMSO day 16"))
  top <- paper_top_hits()$data
  stopifnot(all(c("BAX", "PMAIP1") %in% top$gene_key[top$side == "resistance"]))
  stopifnot(top$n_resistance[top$gene_key == "BAX"] >= 8L)
  for (metric in names(app_env$paper_metric_choices)) {
    session$setInputs(paper_metric = app_env$paper_metric_choices[[metric]])
    stopifnot(nzchar(output$paper_rank_plot))
  }
})
# Every screen-view registry row must match a displayed experiment, and the
# canonical resistance genes must lead the venetoclax screens they appear in.
stopifnot(nrow(dplyr::anti_join(screen_view_registry, dplyr::distinct(scores, paper, experiment), by = c("paper", "experiment"))) == 0L)
screen_support <- summarise_screen_support(app_env$screen_scores, c("BAX", "PMAIP1"))
stopifnot(all(screen_support$n_resistance >= 8L), all(screen_support$n_sensitizing == 0L))
mcl1_support <- summarise_screen_support(app_env$screen_scores, "BCL2L1")
stopifnot(mcl1_support$n_sensitizing >= 3L, mcl1_support$n_sensitizing > mcl1_support$n_resistance)

lock <- jsonlite::read_json("renv.lock", simplifyVector = FALSE)
missing <- names(lock$Packages)[!vapply(names(lock$Packages), requireNamespace, logical(1), quietly = TRUE)]
stopifnot(length(missing) == 0L)
wrong_version <- names(lock$Packages)[vapply(names(lock$Packages), function(package) {
  utils::packageVersion(package) != package_version(lock$Packages[[package]]$Version)
}, logical(1))]
if (length(wrong_version)) stop("Installed versions differ from renv.lock: ", paste(wrong_version, collapse = ", "))
message("Release checks passed: hashes, data schema, catalog, cache-only loading, app construction/server queries (Gene and Paper Explorer), and locked package versions.")
message("Scores: ", nrow(scores), "; studies: ", dplyr::n_distinct(scores$paper), "; displayed experiments: ", dplyr::n_distinct(scores$experiment_id))
