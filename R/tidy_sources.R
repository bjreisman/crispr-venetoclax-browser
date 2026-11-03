library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tibble)

standard_columns <- c(
  "paper", "year", "model_system", "cell_line_or_context", "experiment",
  "gene_symbol", "gene_key", "species", "effect_metric", "effect_value",
  "p_value", "fdr", "rank", "direction", "source_level", "source_file"
)

paper_metadata <- tibble::tribble(
  ~paper, ~year, ~model_system, ~cell_line_or_context, ~species, ~rationale,
  "Lin2019", 2019, "Metabolic CRISPR/Cas9 loss-of-function screen", "AML; THP-1 and MOLM-13", "Human",
  "Metabolism-focused CRISPR screens mapping metabolic modifiers of apoptotic priming and drug sensitivity, nominating heme biosynthesis and mitochondrial metabolism.",
  "Chen2019", 2019, "Genome-wide CRISPR/Cas9 loss-of-function screen", "MOLM-13 AML", "Human",
  "Genome-wide screen for knockouts that alter venetoclax response, highlighting mitochondrial organization and CLPB-linked mitochondrial structure.",
  "Nechiporuk2019", 2019, "Genome-wide CRISPR/Cas9 loss-of-function screens", "MOLM-13 AML", "Human",
  "Venetoclax resistance screens identifying the TP53 apoptotic network, including TP53, BAX, and PMAIP1/NOXA.",
  "Guièze2019", 2019, "CRISPR/Cas9 loss-of-function and ORF gain-of-function screens", "B-lymphoid malignancy models", "Human",
  "Study of venetoclax resistance in lymphoid malignancies using parallel knockout and ORF overexpression screens, implicating BCL2-family buffering, lymphoid transcriptional regulators, AMPK/PKA signaling, and oxidative metabolism.",
  "Glytsou2023", 2023, "Genome-wide CRISPR/Cas9 loss-of-function screens", "MOLM-13 AML", "Human",
  "Parallel BH3-mimetic screens showing mitophagy and mitochondrial dynamics contribute to resistance or sensitivity to venetoclax-based combinations.",
  "Deng2022", 2022, "CRISPRa-associated RNA-seq differential expression", "Eµ-Myc/dCas9a-SAM lymphoma and B-ALL", "Mouse",
  "CRISPR activation mouse model used to study aggressive lymphoma and venetoclax resistance; included here as contextual expression output rather than a directly comparable knockout screen.",
  "Portelinha2024", 2024, "Genome-wide CRISPR/Cas9 loss-of-function positive selection screens", "OCI-LY19 B-cell lymphoma", "Human",
  "B-cell lymphoma screens and functional work connecting SETD1B loss to apoptosis resistance and BCL2 independence.",
  "Eide2023", 2023, "Genome-wide CRISPR/Cas9 loss-of-function screen", "OCI-AML2 AML; venetoclax and venetoclax+ruxolitinib", "Human",
  "Venetoclax-combination study pairing primary AML ex vivo drug response with an OCI-AML2 CRISPR/Cas9 screen for venetoclax and venetoclax+ruxolitinib resistance/sensitization mechanisms.",
  "He2026", 2026, "BH3-priming CRISPR/Cas9 knockout screen", "BH3-mimetic-resistant primary AML PDX cells P1-R and P2-R", "Human",
  "Dynamic BH3-priming screen in resistant AML PDX cells nominating ADSS2 and de novo AMP biosynthesis as sensitizers to venetoclax and MCL1 inhibition.",
  "Nakao2023", 2023, "Genome-wide CRISPR/Cas9 loss-of-function screens", "Mouse MLL::AF9 AML mMA9Cas cells; March5 WT and knockout", "Mouse",
  "Genome-wide mouse AML screens identifying the MARCH5/UBE2J2/UBE2K mitochondrial E3 ubiquitin-ligase complex as a venetoclax-sensitization node and NOXA/Pmaip1 as a March5-regulated effector. Browser scores are reproducible DrugZ re-analyses of the published sgRNA counts; MAGeCK MLE fields are retained for audit and figure concordance."
)

normalize_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- dplyr::na_if(x, "")
  x
}

gene_symbol_key_aliases <- c(
  ADSS = "ADSS2",
  ADSS1 = "ADSSL1",
  TOMM70A = "TOMM70"
)

gene_key <- function(x) {
  key <- normalize_gene_symbol(x) |>
    stringr::str_to_upper() |>
    stringr::str_replace_all("[^A-Z0-9]", "")
  alias <- unname(gene_symbol_key_aliases[key])
  dplyr::if_else(!is.na(alias), alias, key)
}

classify_direction <- function(value, positive_label = "resistance/enriched", negative_label = "sensitization/depleted") {
  dplyr::case_when(
    is.na(value) ~ "unknown",
    value > 0 ~ positive_label,
    value < 0 ~ negative_label,
    TRUE ~ "neutral"
  )
}

coerce_number <- function(x) {
  x <- as.character(x)
  x[x %in% c("#DIV/0!", "NA", "NaN", "Inf", "-Inf", "")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

complete_standard <- function(df) {
  missing <- setdiff(standard_columns, names(df))
  for (col in missing) df[[col]] <- NA
  df |>
    mutate(
      gene_symbol = normalize_gene_symbol(.data$gene_symbol),
      gene_key = if_else(is.na(.data$gene_key), gene_key(.data$gene_symbol), as.character(.data$gene_key)),
      effect_value = coerce_number(.data$effect_value),
      p_value = coerce_number(.data$p_value),
      fdr = coerce_number(.data$fdr),
      rank = coerce_number(.data$rank)
    ) |>
    filter(!is.na(.data$gene_symbol), !is.na(.data$gene_key)) |>
    select(all_of(standard_columns), everything())
}

attach_metadata <- function(df, paper) {
  meta <- paper_metadata |> filter(.data$paper == !!paper) |> slice(1)
  df |>
    mutate(
      paper = paper,
      year = meta$year,
      model_system = meta$model_system,
      cell_line_or_context = meta$cell_line_or_context,
      species = meta$species
    )
}

tidy_lin2019 <- function(root = ".") {
  source_file <- file.path(root, "Lin2019/1-s2.0-S1550413119300117-mmc3.csv")
  raw <- readr::read_csv(source_file, skip = 1, show_col_types = FALSE)
  raw |>
    rename(gene_symbol = Genes, gene_description = `Descriptive Gene Name`) |>
    pivot_longer(
      cols = c(THP_ess, MOLM_ess, Avg_ess, THP_sens, MOLM_sens, Avg_sens),
      names_to = "experiment",
      values_to = "effect_value"
    ) |>
    mutate(
      experiment = recode(
        .data$experiment,
        THP_ess = "THP-1 essentiality",
        MOLM_ess = "MOLM-13 essentiality",
        Avg_ess = "Average essentiality",
        THP_sens = "THP-1 drug sensitivity",
        MOLM_sens = "MOLM-13 drug sensitivity",
        Avg_sens = "Average drug sensitivity"
      ),
      effect_value = coerce_number(.data$effect_value),
      effect_metric = "Published three-score",
      p_value = NA_real_,
      fdr = NA_real_,
      rank = abs(.data$effect_value),
      direction = classify_direction(.data$effect_value, "higher score", "lower score"),
      source_level = "gene",
      source_file = source_file
    ) |>
    attach_metadata("Lin2019") |>
    complete_standard()
}

tidy_chen2019 <- function(root = ".") {
  source_file <- file.path(root, "Chen2019/21598290cd190117-sup-216243_2_supp_5484296_pqktn1.xlsx")
  raw <- readxl::read_excel(source_file, sheet = "CRISPR score and p-value")
  raw |>
    rename(gene_symbol = `Target Gene Symbol`) |>
    transmute(
      gene_symbol,
      `Venetoclax day 8 delta CRISPR score` = `delta CS_day8`,
      `Venetoclax day 8 p-value` = `p_value_day8`,
      `Venetoclax day 16 delta CRISPR score` = `delta CS_day16`,
      `Venetoclax day 16 p-value` = `p_value_day16`
    ) |>
    pivot_longer(-gene_symbol, names_to = "measure", values_to = "value") |>
    mutate(
      experiment = if_else(str_detect(.data$measure, "day 8"), "Venetoclax vs DMSO day 8", "Venetoclax vs DMSO day 16"),
      field = if_else(str_detect(.data$measure, "p-value"), "p_value", "effect_value")
    ) |>
    select(-measure) |>
    pivot_wider(names_from = field, values_from = value) |>
    mutate(
      effect_value = coerce_number(.data$effect_value),
      p_value = coerce_number(.data$p_value),
      effect_metric = "Delta CRISPR score",
      fdr = NA_real_,
      rank = rank(-abs(coerce_number(.data$effect_value)), ties.method = "first"),
      direction = classify_direction(coerce_number(.data$effect_value)),
      source_level = "gene",
      source_file = source_file
    ) |>
    attach_metadata("Chen2019") |>
    complete_standard()
}

tidy_mageck_gene_sheet <- function(source_file, sheet, skip, paper, experiment) {
  raw <- readxl::read_excel(source_file, sheet = sheet, skip = skip)
  if (!"id" %in% names(raw)) return(tibble())
  base <- raw |>
    rename(gene_symbol = id) |>
    mutate(across(everything(), as.character))

  # MAGeCK reports one gene-level LFC shared by both arms; keep its sign so
  # enriched genes are not shown as depleted in the negative-selection arm.
  neg <- if ("neg|lfc" %in% names(base)) {
    base |>
      mutate(mageck_raw_score = coerce_number(`neg|lfc`)) |>
      transmute(
        gene_symbol,
        experiment = paste0(experiment, " | negative selection"),
        effect_metric = "MAGeCK negative selection LFC",
        effect_value = .data$mageck_raw_score,
        p_value = coerce_number(`neg|p-value`),
        fdr = coerce_number(`neg|fdr`),
        rank = coerce_number(`neg|rank`),
        direction = classify_direction(.data$mageck_raw_score),
        source_level = "gene",
        source_file = source_file,
        raw_effect_value = .data$mageck_raw_score
      )
  } else {
    tibble()
  }

  pos <- if ("pos|lfc" %in% names(base)) {
    base |>
      mutate(mageck_raw_score = coerce_number(`pos|lfc`)) |>
      transmute(
        gene_symbol,
        experiment = paste0(experiment, " | positive selection"),
        effect_metric = "MAGeCK positive selection LFC",
        effect_value = .data$mageck_raw_score,
        p_value = coerce_number(`pos|p-value`),
        fdr = coerce_number(`pos|fdr`),
        rank = coerce_number(`pos|rank`),
        direction = classify_direction(.data$mageck_raw_score),
        source_level = "gene",
        source_file = source_file,
        raw_effect_value = .data$mageck_raw_score
      )
  } else {
    tibble()
  }

  bind_rows(neg, pos) |>
    attach_metadata(paper) |>
    complete_standard()
}

tidy_nechiporuk2019 <- function(root = ".") {
  source_file <- file.path(root, "Nechiporuk2019/21598290cd190125-sup-216310_2_table_5485035_pqlflj.xlsx")
  bind_rows(
    tidy_mageck_gene_sheet(source_file, "Table S1. Tab 1 of 4", 1, "Nechiporuk2019", "Y. Kosuke library venetoclax vs DMSO day 7"),
    tidy_mageck_gene_sheet(source_file, "Table S1. Tab 3 of 4", 1, "Nechiporuk2019", "Y. Kosuke library venetoclax vs DMSO day 14"),
    tidy_mageck_gene_sheet(source_file, "Table S2. Tab 1 of 2", 1, "Nechiporuk2019", "Brunello library venetoclax vs DMSO day 14"),
    tidy_mageck_gene_sheet(source_file, "Table S3. Tab 1 of 2", 1, "Nechiporuk2019", "Top Y. Kosuke hits")
  )
}

tidy_glytsou2023 <- function(root = ".") {
  source_file <- file.path(root, "Glytsou2023/cd-22-0601_supplementary_table_s1_suppst1.xlsx")
  sheets <- c("MCL1i+V", "MCL1i", "V+A")
  purrr::map_dfr(sheets, \(sheet) {
    tidy_mageck_gene_sheet(source_file, sheet, 0, "Glytsou2023", sheet)
  })
}

tidy_guieze2019 <- function(root = ".") {
  source_file <- file.path(root, "Guieze2019/1-s2.0-S1535610819303733-mmc3.xlsx")
  raw <- readxl::read_excel(source_file, sheet = "STARS output-score (gene)", skip = 2)
  raw |>
    rename(gene_symbol = `Gene Symbol`) |>
    transmute(
      gene_symbol,
      experiment = paste0("STARS replicate ", Replicate),
      effect_metric = "STARS score",
      effect_value = abs(coerce_number(`STARS Score`)),
      p_value = `p-value`,
      fdr = FDR,
      rank = rank(-abs(coerce_number(`STARS Score`)), ties.method = "first"),
      direction = "resistance/enriched",
      source_level = "gene",
      source_file = source_file,
      raw_effect_value = coerce_number(`STARS Score`)
    ) |>
    attach_metadata("Guièze2019") |>
    complete_standard()
}

tidy_guieze2019_orf <- function(root = ".") {
  source_file <- file.path(root, "Guieze2019/1-s2.0-S1535610819303733-mmc4.xlsx")
  raw <- readxl::read_excel(source_file, sheet = "DRP-ORF_OCILY-2017-03-09-LFC", skip = 1)

  raw |>
    mutate(
      gene_symbol = coalesce(normalize_gene_symbol(.data$`gene symbol`), normalize_gene_symbol(.data$`Best gene match`)),
      lfc_venetoclax = coerce_number(.data$`LFC AVERAGE Veneto`),
      lfc_dmso = coerce_number(.data$`LFC-AVERAGE DMSO`),
      z_ave_lfc_veneto = coerce_number(.data$`Z_AVE LFC VENETO`),
      delta_lfc = coerce_number(.data$`Delta LFC`)
    ) |>
    filter(!is.na(.data$gene_symbol), !is.na(.data$lfc_venetoclax)) |>
    add_count(.data$gene_symbol, name = "n_orf_constructs") |>
    group_by(.data$gene_symbol) |>
    slice_max(.data$lfc_venetoclax, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      experiment = "ORF overexpression screen | venetoclax 100 nM day 10",
      effect_metric = "best ORF LFC average venetoclax",
      effect_value = .data$lfc_venetoclax,
      p_value = NA_real_,
      fdr = NA_real_,
      rank = rank(-.data$lfc_venetoclax, ties.method = "first"),
      direction = classify_direction(.data$lfc_venetoclax, "ORF overexpression resistance/enriched", "ORF overexpression depleted"),
      source_level = "gene summary from ORF constructs",
      source_file = source_file,
      best_construct_id = .data$`Construct IDs`,
      best_construct_barcode = .data$`Construct Barcode`,
      best_gene_match = .data$`Best gene match`,
      clones_per_gene = .data$`clones/gene`
    ) |>
    select(
      gene_symbol, experiment, effect_metric, effect_value, p_value, fdr, rank,
      direction, source_level, source_file, n_orf_constructs, best_construct_id,
      best_construct_barcode, best_gene_match, clones_per_gene, lfc_dmso,
      z_ave_lfc_veneto, delta_lfc
    ) |>
    attach_metadata("Guièze2019") |>
    complete_standard()
}

tidy_portelinha_file <- function(root, file_name, experiment) {
  source_file <- file.path(root, "Portelinha2024", file_name)
  raw <- readxl::read_excel(source_file, sheet = "All positively selected genes")
  raw |>
    group_by(gene_symbol = gene) |>
    summarise(
      raw_effect_value = mean(coerce_number(edger_sgrna_level_lfc), na.rm = TRUE),
      p_value = NA_real_,
      fdr = min(coerce_number(edger_gene_level_fdr), na.rm = TRUE),
      n_sgrna = n(),
      .groups = "drop"
    ) |>
    mutate(
      raw_effect_value = if_else(is.infinite(.data$raw_effect_value), NA_real_, .data$raw_effect_value),
      effect_value = abs(.data$raw_effect_value),
      fdr = if_else(is.infinite(.data$fdr), NA_real_, .data$fdr),
      rank = rank(-.data$effect_value, ties.method = "first", na.last = "keep"),
      experiment = experiment,
      effect_metric = "edgeR sgRNA mean logFC",
      direction = "resistance/enriched",
      source_level = "gene summary from sgRNA",
      source_file = source_file
    ) |>
    attach_metadata("Portelinha2024") |>
    complete_standard()
}

tidy_portelinha2024 <- function(root = ".") {
  bind_rows(
    tidy_portelinha_file(root, "jem_20231143_tables1.xlsx", "Positive selection screen table 1"),
    tidy_portelinha_file(root, "jem_20231143_tables2.xlsx", "Positive selection screen table 2")
  )
}

tidy_deng2022 <- function(root = ".") {
  source_file <- file.path(root, "Deng2022/tables/2022-06-29-DE_gene_list_groupEµ-Myc_dCas9a-SAM_vs_B-ALL_GGiner.csv")
  raw <- readr::read_csv(source_file, show_col_types = FALSE)
  gene_col <- if ("Symbol" %in% names(raw)) "Symbol" else names(raw)[2]
  raw |>
    transmute(
      gene_symbol = .data[[gene_col]],
      experiment = "Eµ-Myc/dCas9a-SAM vs B-ALL expression",
      effect_metric = "limma voom logFC",
      effect_value = logFC,
      p_value = P.Value,
      fdr = adj.P.Val,
      rank = rank(-abs(logFC), ties.method = "first"),
      direction = classify_direction(logFC, "higher in Eµ-Myc/dCas9a-SAM", "higher in B-ALL"),
      source_level = "contextual RNA-seq differential expression",
      source_file = source_file
    ) |>
    attach_metadata("Deng2022") |>
    complete_standard()
}

tidy_eide2023_crispr <- function(root = ".") {
  source_file <- file.path(root, "data/external/eide2023/GSE216087_crispr_counts_for_geo.txt.gz")
  if (!file.exists(source_file)) return(tibble())

  raw <- readr::read_tsv(source_file, show_col_types = FALSE)
  # Scale each sample to the median library size before computing fold-changes;
  # raw-count ratios otherwise carry a sequencing-depth offset.
  count_cols <- setdiff(names(raw), c("sgRNA", "Gene"))
  library_sizes <- colSums(raw[count_cols], na.rm = TRUE)
  size_factors <- library_sizes / stats::median(library_sizes)
  raw[count_cols] <- purrr::map2(raw[count_cols], size_factors, `/`)
  condition_specs <- tibble::tribble(
    ~experiment, ~control_pattern, ~treatment_pattern,
    "OCI-AML2 Yusa library day 14 venetoclax vs DMSO", "D14_DMSO", "D14_VEN[0-9]",
    "OCI-AML2 Yusa library day 14 venetoclax+ruxolitinib vs DMSO", "D14_DMSO", "D14_VENRUX",
    "OCI-AML2 Yusa library day 21 venetoclax vs DMSO", "D21_DMSO", "D21_VEN[0-9]",
    "OCI-AML2 Yusa library day 21 venetoclax+ruxolitinib vs DMSO", "D21_DMSO", "D21_VENRUX"
  )

  purrr::pmap_dfr(condition_specs, function(experiment, control_pattern, treatment_pattern) {
    control_cols <- names(raw)[stringr::str_detect(names(raw), control_pattern)]
    treatment_cols <- names(raw)[stringr::str_detect(names(raw), treatment_pattern)]
    if (!length(control_cols) || !length(treatment_cols)) return(tibble())

    raw |>
      transmute(
        gene_symbol = .data$Gene,
        sgrna_lfc = log2((rowMeans(pick(all_of(treatment_cols)), na.rm = TRUE) + 1) /
          (rowMeans(pick(all_of(control_cols)), na.rm = TRUE) + 1))
      ) |>
      group_by(.data$gene_symbol) |>
      summarise(
        effect_value = stats::median(.data$sgrna_lfc, na.rm = TRUE),
        n_sgrna = n(),
        .groups = "drop"
      ) |>
      mutate(
        # No non-targeting guides in this library; center on the median gene.
        effect_value = .data$effect_value - stats::median(.data$effect_value, na.rm = TRUE),
        experiment = experiment,
        effect_metric = "Median-centered sgRNA log2 fold-change from normalized GEO counts",
        p_value = NA_real_,
        fdr = NA_real_,
        rank = rank(-abs(.data$effect_value), ties.method = "first"),
        direction = classify_direction(.data$effect_value),
        source_level = "computed gene summary from sgRNA counts",
        source_file = source_file
      )
  }) |>
    attach_metadata("Eide2023") |>
    complete_standard()
}

tidy_he2026 <- function(root = ".") {
  source_file <- file.path(root, "He2026/43018_2026_1184_MOESM2_ESM.xlsx")
  mageck_file <- file.path(root, "He2026/43018_2026_1184_MOESM4_ESM.xlsx")
  mageck_summary <- if (file.exists(mageck_file)) {
    mageck <- readxl::read_excel(mageck_file, sheet = "Fig.1", range = "U3:W734") |>
      rename(
        gene_symbol = gene,
        mageck_lfc = lfc,
        log_p_value = `log.p.value`
      ) |>
      mutate(
        gene_symbol = normalize_gene_symbol(.data$gene_symbol),
        mageck_lfc = coerce_number(.data$mageck_lfc),
        log_p_value = coerce_number(.data$log_p_value)
      ) |>
      filter(!is.na(.data$gene_symbol), !is.na(.data$mageck_lfc)) |>
      group_by(.data$gene_symbol) |>
      arrange(desc(.data$log_p_value), desc(abs(.data$mageck_lfc)), .by_group = TRUE) |>
      slice(1) |>
      ungroup()

    crispr_scores <- readxl::read_excel(mageck_file, sheet = "Fig.1", range = "Y3:Z734") |>
      rename(
        gene_symbol = Gene,
        crispr_score = CRISPR_score
      ) |>
      mutate(
        gene_symbol = normalize_gene_symbol(.data$gene_symbol),
        crispr_score = coerce_number(.data$crispr_score)
      ) |>
      filter(!is.na(.data$gene_symbol)) |>
      group_by(.data$gene_symbol) |>
      arrange(desc(abs(.data$crispr_score)), .by_group = TRUE) |>
      slice(1) |>
      ungroup()

    mageck |>
      left_join(crispr_scores, by = "gene_symbol") |>
      transmute(
        gene_symbol,
        experiment = "P1-R/P2-R combined BH3 priming screen CytoC-high vs CytoC-low",
        effect_metric = "MAGeCK LFC (CytoC-high/CytoC-low)",
        effect_value = .data$mageck_lfc,
        p_value = 10^(-.data$log_p_value),
        fdr = NA_real_,
        rank = rank(-abs(.data$mageck_lfc), ties.method = "first"),
        direction = classify_direction(.data$mageck_lfc),
        source_level = "MAGeCK gene summary from figure source data",
        source_file = mageck_file,
        raw_effect_value = .data$mageck_lfc,
        log_p_value = .data$log_p_value,
        crispr_score = .data$crispr_score,
        source_sample = "P1-R/P2-R"
      )
  } else {
    tibble()
  }

  if (!file.exists(source_file)) {
    return(mageck_summary |> attach_metadata("He2026") |> complete_standard())
  }

  raw <- readxl::read_excel(source_file, sheet = "supplement table3", skip = 2) |>
    rename(
      sgrna_id = GeneName,
      p1_baseline_rpm = `P1-R aseline`,
      p1_cytoc_high_rpm = `P1-R CytoC-high`,
      p1_cytoc_low_rpm = `P1-R CytoC-low`,
      p2_baseline_rpm = `P2-R Baseline`,
      p2_cytoc_high_rpm = `P2-R CytoC-high`,
      p2_cytoc_low_rpm = `P2-R CytoC-low`
    ) |>
    mutate(
      gene_symbol = str_replace(.data$sgrna_id, "-[0-9]+$", ""),
      across(ends_with("_rpm"), coerce_number)
    )

  condition_specs <- tibble::tribble(
    ~experiment, ~sample_id, ~baseline_col, ~high_col, ~low_col,
    "P1-R BH3 priming screen CytoC-high vs CytoC-low", "P1-R", "p1_baseline_rpm", "p1_cytoc_high_rpm", "p1_cytoc_low_rpm",
    "P2-R BH3 priming screen CytoC-high vs CytoC-low", "P2-R", "p2_baseline_rpm", "p2_cytoc_high_rpm", "p2_cytoc_low_rpm"
  )

  sample_summaries <- purrr::pmap_dfr(condition_specs, function(experiment, sample_id, baseline_col, high_col, low_col) {
    raw |>
      transmute(
        gene_symbol,
        sgrna_id,
        sgrna_lfc_high_vs_low = log2((.data[[high_col]] + 1) / (.data[[low_col]] + 1)),
        baseline_rpm = .data[[baseline_col]],
        cytoc_high_rpm = .data[[high_col]],
        cytoc_low_rpm = .data[[low_col]]
      ) |>
      filter(!is.na(.data$sgrna_lfc_high_vs_low)) |>
      group_by(.data$gene_symbol) |>
      summarise(
        effect_value = stats::median(.data$sgrna_lfc_high_vs_low, na.rm = TRUE),
        n_sgrna = n_distinct(.data$sgrna_id),
        n_screen_observations = n(),
        median_baseline_rpm = stats::median(.data$baseline_rpm, na.rm = TRUE),
        median_cytoc_high_rpm = stats::median(.data$cytoc_high_rpm, na.rm = TRUE),
        median_cytoc_low_rpm = stats::median(.data$cytoc_low_rpm, na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        experiment = experiment,
        source_sample = sample_id
      )
  })

  rpm_summaries <- sample_summaries |>
    mutate(
      effect_metric = "Median sgRNA log2 fold-change (CytoC-high/CytoC-low RPM)",
      p_value = NA_real_,
      fdr = NA_real_,
      rank = rank(-abs(.data$effect_value), ties.method = "first"),
      direction = classify_direction(.data$effect_value),
      source_level = "computed gene summary from sgRNA RPM",
      source_file = source_file
    )

  bind_rows(mageck_summary, rpm_summaries) |>
    attach_metadata("He2026") |>
    complete_standard()
}

tidy_nakao2023 <- function(root = ".") {
  source_file <- file.path(root, "Nakao2023/41375_2023_1879_MOESM2_ESM.xlsx")
  result_file <- file.path(root, "data/processed/nakao2023/gene_results.tsv")
  manifest_file <- file.path(root, "data/processed/nakao2023/manifest.json")
  if (!file.exists(result_file) || !file.exists(manifest_file)) {
    stop(
      "Nakao2023 derived results are missing. Create .venv-nakao2023, install scripts/requirements-nakao2023.txt, then run .venv-nakao2023/bin/python scripts/build_nakao2023_analysis.py --root ."
    )
  }
  if (file.exists(source_file) && file.info(result_file)$mtime < file.info(source_file)$mtime) {
    stop("Nakao2023 derived results are older than the source workbook. Re-run scripts/build_nakao2023_analysis.py before rebuilding browser data.")
  }

  readr::read_tsv(result_file, show_col_types = FALSE) |>
    rename(gene_symbol = Gene) |>
    transmute(
      gene_symbol,
      experiment,
      effect_metric = "DrugZ normalized Z (drug vs DMSO)",
      effect_value,
      p_value,
      fdr,
      rank,
      direction,
      source_level = "re-analyzed gene summary from published sgRNA counts",
      source_file = source_file,
      analysis_id,
      source_sheet,
      n_sgrna,
      drugz_num_obs,
      drugz_p_value_synth,
      drugz_fdr_synth,
      drugz_rank_synth,
      drugz_p_value_supp,
      drugz_fdr_supp,
      drugz_rank_supp,
      mageck_dmso_beta,
      mageck_treatment_beta,
      mageck_beta_delta,
      mageck_dmso_p_value,
      mageck_dmso_fdr,
      mageck_dmso_wald_p_value,
      mageck_dmso_wald_fdr,
      mageck_treatment_p_value,
      mageck_treatment_fdr,
      mageck_treatment_wald_p_value,
      mageck_treatment_wald_fdr
    ) |>
    attach_metadata("Nakao2023") |>
    complete_standard()
}

build_gene_scores <- function(root = ".") {
  bind_rows(
    tidy_lin2019(root),
    tidy_chen2019(root),
    tidy_nechiporuk2019(root),
    tidy_guieze2019(root),
    tidy_guieze2019_orf(root),
    tidy_glytsou2023(root),
    tidy_deng2022(root),
    tidy_portelinha2024(root),
    tidy_eide2023_crispr(root),
    tidy_he2026(root),
    tidy_nakao2023(root)
  ) |>
    left_join(paper_metadata |> select(paper, rationale), by = "paper") |>
    mutate(
      experiment_id = paste(.data$paper, .data$experiment, sep = " :: "),
      significance = case_when(
        .data$paper == "Lin2019" &
          .data$experiment == "Average drug sensitivity" &
          !is.na(.data$effect_value) &
          .data$effect_value < -0.914 ~ "paper threshold: Avg_sens < -0.914",
        .data$paper == "Guièze2019" &
          grepl("ORF overexpression", .data$experiment) &
          !is.na(.data$effect_value) &
          .data$effect_value >= 2 ~ "paper threshold: ORF LFC >= 2",
        !is.na(.data$fdr) & .data$fdr <= 0.05 ~ "FDR <= 0.05",
        !is.na(.data$p_value) & .data$p_value <= 0.05 ~ "p <= 0.05",
        TRUE ~ "not annotated/significant"
      ),
      signed_log10_fdr = case_when(
        !is.na(.data$fdr) & .data$fdr > 0 ~ sign(.data$effect_value) * -log10(.data$fdr),
        !is.na(.data$fdr) & .data$fdr == 0 ~ sign(.data$effect_value) * 300,
        TRUE ~ NA_real_
      )
    )
}

load_detail_for_gene <- function(gene_key_value, paper, root = ".") {
  paper <- as.character(paper)[1]
  if (is.na(paper) || is.na(gene_key_value)) return(tibble())
  missing_detail_note <- tibble(note = "Source-detail table not bundled for deployment. Use the local app to inspect raw/source rows.")

  if (paper == "Chen2019") {
    source_file <- file.path(root, "Chen2019/21598290cd190117-sup-216243_2_supp_5484296_pqktn1.xlsx")
    if (!file.exists(source_file)) return(missing_detail_note)
    raw_counts <- readxl::read_excel(source_file, sheet = "CRISPR screen raw counts")
    scores <- readxl::read_excel(source_file, sheet = "CRISPR score and p-value")
    bind_rows(
      raw_counts |> rename(gene_symbol = `Target Gene Symbol`) |> mutate(source_sheet = "CRISPR screen raw counts"),
      scores |> rename(gene_symbol = `Target Gene Symbol`) |> mutate(source_sheet = "CRISPR score and p-value")
    ) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Nechiporuk2019") {
    source_file <- file.path(root, "Nechiporuk2019/21598290cd190125-sup-216310_2_table_5485035_pqlflj.xlsx")
    if (!file.exists(source_file)) return(missing_detail_note)
    sheets <- c("Table S1. Tab 1 of 4", "Table S1. Tab 2 of 4", "Table S1. Tab 3 of 4", "Table S1. Tab 4 of 4", "Table S2. Tab 1 of 2", "Table S2. Tab 2 of 2", "Table S3. Tab 1 of 2", "Table S3. Tab 2 of 2")
    map_dfr(sheets, \(sheet) {
      x <- readxl::read_excel(source_file, sheet = sheet, skip = 1)
      gene_col <- intersect(c("id", "Gene", "gene"), names(x))[1]
      if (is.na(gene_col)) return(tibble())
      x |> mutate(gene_symbol = .data[[gene_col]], source_sheet = sheet)
    }) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Glytsou2023") {
    source_file <- file.path(root, "Glytsou2023/cd-22-0601_supplementary_table_s1_suppst1.xlsx")
    if (!file.exists(source_file)) return(missing_detail_note)
    map_dfr(c("MCL1i+V", "MCL1i", "V+A"), \(sheet) {
      readxl::read_excel(source_file, sheet = sheet) |>
        rename(gene_symbol = id) |>
        mutate(source_sheet = sheet)
    }) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Guièze2019") {
    lof_file <- file.path(root, "Guieze2019/1-s2.0-S1535610819303733-mmc3.xlsx")
    orf_file <- file.path(root, "Guieze2019/1-s2.0-S1535610819303733-mmc4.xlsx")
    if (!file.exists(lof_file) && !file.exists(orf_file)) return(missing_detail_note)
    bind_rows(
      if (file.exists(lof_file)) {
        readxl::read_excel(lof_file, sheet = "STARS output-score (gene)", skip = 2) |>
          rename(gene_symbol = `Gene Symbol`) |>
          mutate(source_sheet = "STARS output-score (gene)")
      } else {
        tibble()
      },
      if (file.exists(orf_file)) {
        readxl::read_excel(orf_file, sheet = "DRP-ORF_OCILY-2017-03-09-LFC", skip = 1) |>
          mutate(
            gene_symbol = coalesce(normalize_gene_symbol(.data$`gene symbol`), normalize_gene_symbol(.data$`Best gene match`)),
            source_sheet = "ORF overexpression screen"
          )
      } else {
        tibble()
      }
    ) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Lin2019") {
    source_file <- file.path(root, "Lin2019/1-s2.0-S1550413119300117-mmc3.csv")
    if (!file.exists(source_file)) return(missing_detail_note)
    readr::read_csv(source_file, skip = 1, show_col_types = FALSE) |>
      rename(gene_symbol = Genes) |>
      mutate(source_sheet = "Supplementary Table 2", gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Portelinha2024") {
    if (!file.exists(file.path(root, "Portelinha2024/jem_20231143_tables1.xlsx")) ||
        !file.exists(file.path(root, "Portelinha2024/jem_20231143_tables2.xlsx"))) {
      return(missing_detail_note)
    }
    map_dfr(c("jem_20231143_tables1.xlsx", "jem_20231143_tables2.xlsx"), \(file_name) {
      source_file <- file.path(root, "Portelinha2024", file_name)
      readxl::read_excel(source_file, sheet = "All positively selected genes") |>
        mutate(source_sheet = paste(file_name, "All positively selected genes"))
    }) |>
      mutate(gene_symbol = gene, gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Deng2022") {
    source_file <- file.path(root, "Deng2022/tables/2022-06-29-DE_gene_list_groupEµ-Myc_dCas9a-SAM_vs_B-ALL_GGiner.csv")
    if (!file.exists(source_file)) return(missing_detail_note)
    readr::read_csv(source_file, show_col_types = FALSE) |>
      mutate(gene_symbol = Symbol, source_sheet = "Differential expression table", gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "He2026") {
    source_file <- file.path(root, "He2026/43018_2026_1184_MOESM2_ESM.xlsx")
    mageck_file <- file.path(root, "He2026/43018_2026_1184_MOESM4_ESM.xlsx")
    if (!file.exists(source_file) && !file.exists(mageck_file)) return(missing_detail_note)
    bind_rows(
      if (file.exists(mageck_file)) {
        readxl::read_excel(mageck_file, sheet = "Fig.1", range = "U3:W734") |>
          rename(gene_symbol = gene, mageck_lfc = lfc, log_p_value = `log.p.value`) |>
          mutate(
            gene_symbol = normalize_gene_symbol(.data$gene_symbol),
            mageck_lfc = coerce_number(.data$mageck_lfc),
            log_p_value = coerce_number(.data$log_p_value)
          ) |>
          filter(!is.na(.data$gene_symbol), !is.na(.data$mageck_lfc)) |>
          group_by(.data$gene_symbol) |>
          arrange(desc(.data$log_p_value), desc(abs(.data$mageck_lfc)), .by_group = TRUE) |>
          slice(1) |>
          ungroup() |>
          left_join(
            readxl::read_excel(mageck_file, sheet = "Fig.1", range = "Y3:Z734") |>
              rename(gene_symbol = Gene, crispr_score = CRISPR_score) |>
              mutate(
                gene_symbol = normalize_gene_symbol(.data$gene_symbol),
                crispr_score = coerce_number(.data$crispr_score)
              ) |>
              filter(!is.na(.data$gene_symbol)) |>
              group_by(.data$gene_symbol) |>
              arrange(desc(abs(.data$crispr_score)), .by_group = TRUE) |>
              slice(1) |>
              ungroup(),
            by = "gene_symbol"
          ) |>
          mutate(source_sheet = "Figure 1 source data - MAGeCK analysis")
      } else {
        tibble()
      },
      if (file.exists(source_file)) {
        readxl::read_excel(source_file, sheet = "supplement table3", skip = 2) |>
          rename(sgrna_id = GeneName) |>
          mutate(
            gene_symbol = str_replace(.data$sgrna_id, "-[0-9]+$", ""),
            source_sheet = "Supplementary Table 3 sgRNA RPM"
          )
      } else {
        tibble()
      }
    ) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else if (paper == "Nakao2023") {
    source_file <- file.path(root, "Nakao2023/41375_2023_1879_MOESM2_ESM.xlsx")
    if (!file.exists(source_file)) return(missing_detail_note)
    map_dfr(c("Supplemental Table 1", "Supplemental Table 2"), \(sheet) {
      readxl::read_excel(source_file, sheet = sheet, skip = 1) |>
        rename(sgrna_id = sgRNA, gene_symbol = Gene) |>
        mutate(source_sheet = sheet)
    }) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      filter(.data$gene_key == gene_key_value)
  } else {
    tibble()
  }
}
