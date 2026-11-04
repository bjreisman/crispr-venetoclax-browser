# Rebuild BrimEdwards2026 gene-level results from the GEO sgRNA counts.
#
# The paper's gene-level AZD4320 table ("Supplemental Data 1, CRISPR Screen
# Results.xlsx") is not in the published supplement, so this script reproduces
# the authors' pipeline (Haematologica Methods "Genome-wide CRISPR/Cas9
# screening"; the same OHSU pipeline documented in the Stewart2026 Table S4
# ReadMe):
#   * drop sgRNAs with 0 plasmid reads;
#   * drop sgRNAs at <= 100 reads (as CPM at the mean contrast library size)
#     in more than half of the contrast samples;
#   * TMM normalization, edgeR estimateDisp + exactTest on the contrast samples;
#   * gene = the sgRNA with the middle p-value (ceiling(n/2) in p order);
#     its logFC, p-value and BH FDR are mid_lfc / mid_pval / mid_fdr;
#   * concordance tiers as described in the Stewart2026 ReadMe, with sgRNA
#     p < 0.005 as "significant" (reproduces the Figure 2A tier counts).
# Checked against Stewart2026 PalboVsDMSO (published sgRNA and gene tabs), see
# validate.R. Requires edgeR (Bioconductor); the browser build itself only
# reads the TSV written here.
#
# Usage: Rscript --vanilla build_brimedwards2026.R <repo root>

be_read_counts <- function(path) {
  raw <- readr::read_tsv(path, show_col_types = FALSE)
  m <- as.matrix(raw[, setdiff(names(raw), c("sgRNA", "Gene"))])
  rownames(m) <- raw$sgRNA
  list(counts = m, gene = raw$Gene)
}

be_guide_test <- function(counts, gene, plasmid_col, control_cols, treatment_cols) {
  x <- counts[, c(control_cols, treatment_cols), drop = FALSE]
  group <- factor(rep(c("control", "treatment"), c(length(control_cols), length(treatment_cols))),
                  levels = c("control", "treatment"))
  all_zero <- rowSums(x) == 0
  # plasmid_col = NULL skips the plasmid filter (count files without a plasmid).
  pre_filter <- if (is.null(plasmid_col)) rep(FALSE, nrow(counts)) else counts[, plasmid_col] == 0
  cpm_threshold <- 100 / (mean(colSums(x)) / 1e6)
  cont_filter <- rowSums(edgeR::cpm(x) <= cpm_threshold) > ncol(x) / 2
  tested <- !all_zero & !pre_filter & !cont_filter
  y <- edgeR::DGEList(x[tested, , drop = FALSE], group = group)
  y <- edgeR::normLibSizes(y, method = "TMM")
  y <- edgeR::estimateDisp(y, stats::model.matrix(~group))
  et <- edgeR::exactTest(y, pair = c("control", "treatment"))$table
  out <- data.frame(sgRNA = rownames(counts), Gene = gene, tested = tested,
                    logFC = NA_real_, PValue = NA_real_, fdr = NA_real_)
  out$logFC[tested] <- et$logFC
  out$PValue[tested] <- et$PValue
  out$fdr[tested] <- stats::p.adjust(et$PValue, "BH")
  out
}

# Tiers per the Stewart2026 ReadMe: significant = sgRNA FDR < 0.05 (default),
# or sgRNA p < sig_p when given. Tier1:
# >= 2 significant guides, all |LFC| >= 2, one direction, and >= 75% of
# filtered guides significant; Tier2: as Tier1 without the 75% rule; Tier3: as
# Tier2 with |LFC| >= 1; Singleton: one significant guide with |LFC| >= 2.
be_gene_summary <- function(guides, sig_p = NULL) {
  guides$sig <- guides$tested & if (is.null(sig_p)) guides$fdr < 0.05 else guides$PValue < sig_p
  guides$sig[is.na(guides$sig)] <- FALSE
  guides |>
    dplyr::group_by(.data$Gene) |>
    dplyr::summarise(
      total = dplyr::n(),
      not_filtered = sum(.data$tested),
      num_sig = sum(.data$sig),
      min_fdr = if (any(.data$tested)) min(.data$fdr, na.rm = TRUE) else NA_real_,
      sig_lfc_min = if (num_sig > 0) min(.data$logFC[.data$sig]) else NA_real_,
      sig_lfc_max = if (num_sig > 0) max(.data$logFC[.data$sig]) else NA_real_,
      sig_min_abslfc = if (num_sig > 0) min(abs(.data$logFC[.data$sig])) else NA_real_,
      num_up = sum(.data$sig & .data$logFC > 0, na.rm = TRUE),
      num_down = sum(.data$sig & .data$logFC < 0, na.rm = TRUE),
      lfc_mean = if (any(.data$tested)) mean(.data$logFC, na.rm = TRUE) else NA_real_,
      mid_index = ceiling(sum(.data$tested) / 2),
      mid_lfc = if (any(.data$tested)) .data$logFC[.data$tested][order(.data$PValue[.data$tested])][mid_index] else NA_real_,
      mid_pval = if (any(.data$tested)) .data$PValue[.data$tested][order(.data$PValue[.data$tested])][mid_index] else NA_real_,
      mid_fdr = if (any(.data$tested)) .data$fdr[.data$tested][order(.data$PValue[.data$tested])][mid_index] else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      one_direction = .data$num_sig > 0 & (.data$num_up == 0 | .data$num_down == 0),
      tiers = dplyr::case_when(
        .data$num_sig >= 2 & .data$one_direction & .data$sig_min_abslfc >= 2 &
          .data$num_sig / .data$not_filtered >= 0.75 ~ "Tier1",
        .data$num_sig >= 2 & .data$one_direction & .data$sig_min_abslfc >= 2 ~ "Tier2",
        .data$num_sig >= 2 & .data$one_direction & .data$sig_min_abslfc >= 1 ~ "Tier3",
        .data$num_sig == 1 & .data$sig_min_abslfc >= 2 ~ "Singleton",
        TRUE ~ "Unassigned"
      )
    ) |>
    dplyr::select(-"mid_index", -"one_direction")
}

# The AZD4320 Figure 2A legend reports Tier1 (9), Tier2 (103), Tier3 (11).
# Guide FDR < 0.05 gives only 3/9/0 on the GEO counts; guide p < 0.005 gives
# 9/104/10 and recovers all nine labeled Tier1 genes, so that cutoff is used.
be_tier_sig_p <- 0.005

be_build <- function(root = ".") {
  count_file <- file.path(root, "BrimEdwards2026/GSE294240_azd4320_crispr_count.txt.gz")
  out_file <- file.path(root, "BrimEdwards2026/brimedwards2026_gene_results.tsv")
  cts <- be_read_counts(count_file)
  cols <- colnames(cts$counts)
  guides <- be_guide_test(
    cts$counts, cts$gene,
    plasmid_col = grep("plasmid", cols, value = TRUE),
    control_cols = grep("_d14_DMSO_", cols, value = TRUE),
    treatment_cols = grep("_d14_4320_", cols, value = TRUE)
  )
  genes <- be_gene_summary(guides, sig_p = be_tier_sig_p)
  readr::write_tsv(genes, out_file)
  invisible(genes)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  genes <- be_build(if (length(args)) args[1] else ".")
  cat("Wrote", nrow(genes), "genes;", "tiers:\n"); print(table(genes$tiers))
}
