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
  "Deng2022", 2022, "Genome-wide CRISPR activation (CRISPRa) screens", "Eµ-Myc/dCas9a-SAM/sgBcl-2 mouse lymphoma lines #214 and #216", "Mouse",
  "CRISPR activation mouse model of aggressive MYC/BCL-2 lymphoma. Genome-wide CRISPRa screens found that activating A1 (Bcl2a1; human BFL-1), BCL-XL, or MCL-1 confers venetoclax resistance. The browser includes the published enriched-sgRNA tables as gain-of-function hit lists.",
  "Portelinha2024", 2024, "Genome-wide CRISPR/Cas9 loss-of-function positive selection screens", "OCI-LY19 B-cell lymphoma", "Human",
  "B-cell lymphoma screens and functional work connecting SETD1B loss to apoptosis resistance and BCL2 independence.",
  "Eide2023", 2023, "Genome-wide CRISPR/Cas9 loss-of-function screen", "OCI-AML2 AML; venetoclax and venetoclax+ruxolitinib", "Human",
  "Venetoclax-combination study pairing primary AML ex vivo drug response with an OCI-AML2 CRISPR/Cas9 screen for venetoclax and venetoclax+ruxolitinib resistance/sensitization mechanisms.",
  "He2026", 2026, "BH3-priming CRISPR/Cas9 knockout screen", "BH3-mimetic-resistant primary AML PDX cells P1-R and P2-R", "Human",
  "Dynamic BH3-priming screen in resistant AML PDX cells nominating ADSS2 and de novo AMP biosynthesis as sensitizers to venetoclax and MCL1 inhibition.",
  "Nakao2023", 2023, "Genome-wide CRISPR/Cas9 loss-of-function screens", "Mouse MLL::AF9 AML mMA9Cas cells; March5 WT and knockout", "Mouse",
  "Genome-wide mouse AML screens identifying the MARCH5/UBE2J2/UBE2K mitochondrial E3 ubiquitin-ligase complex as a venetoclax-sensitization node and NOXA/Pmaip1 as a March5-regulated effector. Browser scores are reproducible DrugZ re-analyses of the published sgRNA counts; MAGeCK MLE fields are retained for audit and figure concordance.",
  "Wang2023", 2023, "Genome-wide and RBP domain-focused CRISPR/Cas9 loss-of-function screens", "MOLM-13 AML", "Human",
  "Parallel AML drug screens in MOLM-13 (Brunello genome-wide and a 492-gene RNA-binding-protein domain library) showing that loss of RNA splicing factors such as RBM10 and SR proteins sensitizes to venetoclax, motivating combination with the splicing-kinase inhibitor SM09419. Browser scores are the published venetoclax CRISPR scores (mean sgRNA log2 fold change vs DMSO); no p-values or FDRs were published.",
  "Sharon2019", 2019, "Genome-wide CRISPR/Cas9 loss-of-function screen", "MOLM13-R1 acquired venetoclax-resistant AML (MOLM-13 derived)", "Human",
  "Genome-wide TKO v1 knockout screen in acquired venetoclax-resistant MOLM13-R1 cells (venetoclax 400 nM vs DMSO, days 8-29) showing that loss of mitochondrial translation genes (DAP3, MRPL54, MRPL17, RBFA) restores venetoclax sensitivity; BAX loss was the top enriched hit.",
  "Yokoyama2026", 2026, "Genome-wide CRISPR/Cas9 loss-of-function screens", "Peripheral T-cell lymphoma lines ED40515(-) and DEL; venetoclax, navitoclax and ABBV-467", "Human",
  "Genome-wide knockout screens in two PTCL lines under venetoclax, navitoclax (BCL2/BCL-XL) and the MCL1 inhibitor ABBV-467, nominating MYLIP loss (via LDLR and cholesterol uptake) and the upstream RXRA/RXRB receptors as BH3-mimetic sensitizers. Browser rows are the published MAGeCK gene summaries (drug vs DMSO); the table gives p-values but no FDRs.",
  "Stewart2026", 2026, "Genome-wide CRISPR/Cas9 loss-of-function screens", "OCI-AML2 AML; venetoclax, palbociclib, and venetoclax+palbociclib", "Human",
  "OCI-AML2 Yusa-library knockout screens (day 21 vs DMSO) showing that loss of translation genes sensitizes to venetoclax but not venetoclax+palbociclib, and that the combination overcomes BAX-loss (venetoclax) and RB1-loss (palbociclib) resistance. Venetoclax and combination contrasts were not count-normalized by the authors, so raw LFCs are shifted negative; compare within-screen percentiles.",
  "Diepstraten2026", 2026, "Genome-wide CRISPR activation (CRISPRa) screens, Partita high-density library", "Eµ-Myc/dCas9a-SAM/sgBcl2 mouse lymphoma lines DHL214 and DHL216; venetoclax IC50/IC80", "Mouse",
  "High-density Partita mouse CRISPRa library (10 sgRNAs per TSS, five sublibraries) screened in DHL214/DHL216 lymphoma under venetoclax. Enriched sgRNAs targeted Bcl2, Bcl2l1 and Mcl1, and validation nominated Irx5 (via BCL-XL induction), Rac1, Trip6 and Dock8 as resistance factors. Selection was heavily bottlenecked, so browser rows are a hit-list reanalysis of the GEO counts (sgRNAs at >=5% read frequency, the authors' Fig. 3B threshold); no p-values or FDRs.",
  "BrimEdwards2026", 2026, "Genome-wide CRISPR/Cas9 loss-of-function screen", "OCI-AML2 AML; AZD4320 (dual BCL2/BCL-XL inhibitor)", "Human",
  "OCI-AML2 Yusa-library knockout screen (AZD4320 30 nM vs DMSO, day 14) showing that loss of the m6A writer adaptor RBM15 or the NuRD interactor ZMYND8 confers resistance to AZD4320 and venetoclax, alongside BAX, PMAIP1 and TP53. Browser scores re-derive the authors' edgeR median-p sgRNA summaries from the GEO counts because the gene-level table was not published; the paper's venetoclax comparison uses an earlier OCI-AML2 screen, not a new one."
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
  TOMM70A = "TOMM70",
  # Mouse A1 paralogs share one human ortholog, BCL2A1 (BFL-1).
  BCL2A1A = "BCL2A1",
  BCL2A1B = "BCL2A1",
  BCL2A1C = "BCL2A1",
  BCL2A1D = "BCL2A1"
)

# Withdrawn HGNC symbols still used by the published screen libraries (MARCH5, TAZ,
# ATP5B, SEPT9, ...), mapped to the current symbols that DepMap and other reference data
# use. Regenerate with scripts/build_gene_symbol_aliases.R.
gene_symbol_alias_file <- file.path("data", "metadata", "gene_symbol_aliases.csv")
gene_symbol_alias_cache <- new.env(parent = emptyenv())

hgnc_symbol_aliases <- function() {
  if (is.null(gene_symbol_alias_cache$map)) {
    map <- character()
    if (file.exists(gene_symbol_alias_file)) {
      table <- readr::read_csv(gene_symbol_alias_file, col_types = readr::cols(.default = "c"))
      map <- stats::setNames(table$to, table$from)
    }
    gene_symbol_alias_cache$map <- map
  }
  gene_symbol_alias_cache$map
}

gene_key <- function(x) {
  key <- normalize_gene_symbol(x) |>
    stringr::str_to_upper() |>
    stringr::str_replace_all("[^A-Z0-9]", "")
  # Curated aliases win over the generated HGNC table.
  alias <- unname(gene_symbol_key_aliases[key])
  renamed <- unname(hgnc_symbol_aliases()[key])
  key <- dplyr::if_else(!is.na(renamed), renamed, key)
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

# Published enriched-sgRNA tables (Supplementary Data 2): edgeR/voom, IC80
# venetoclax vs DMSO. Only significant guides are listed, so these are hit lists.
# Each gene key keeps its strongest guide (lowest FDR, then largest logFC);
# enriched_symbols records mouse paralogs merged under one human key.
tidy_deng2022_crispra <- function(root = ".") {
  source_file <- file.path(root, "data/external/deng2022/41467_2022_32485_MOESM5_ESM.xlsx")
  if (!file.exists(source_file)) return(tibble())
  sheets <- c(
    "214 Enriched sgRNAs" = "CRISPRa line #214 venetoclax IC80 (30 nM) vs DMSO | enriched sgRNAs",
    "216 Enriched sgRNAs" = "CRISPRa line #216 venetoclax IC80 (15 nM) vs DMSO | enriched sgRNAs"
  )
  purrr::imap_dfr(sheets, function(experiment, sheet) {
    readxl::read_excel(source_file, sheet = sheet) |>
      transmute(
        gene_symbol = normalize_gene_symbol(.data$Symbol),
        sgrna = .data$sgRNA,
        effect_value = coerce_number(.data$logFC),
        p_value = coerce_number(.data$P.Value),
        fdr = coerce_number(.data$FDR)
      ) |>
      filter(!is.na(.data$gene_symbol)) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
      group_by(.data$gene_key) |>
      arrange(.data$fdr, desc(.data$effect_value), .by_group = TRUE) |>
      mutate(
        n_enriched_sgrna = n(),
        enriched_symbols = paste(sort(unique(.data$gene_symbol)), collapse = "; ")
      ) |>
      slice(1) |>
      ungroup() |>
      mutate(
        experiment = experiment,
        effect_metric = "edgeR sgRNA logFC (strongest enriched guide)",
        rank = rank(.data$fdr, ties.method = "first"),
        direction = classify_direction(.data$effect_value),
        source_level = "published enriched-sgRNA table (hit list)",
        source_file = source_file
      )
  }) |>
    attach_metadata("Deng2022") |>
    complete_standard()
}

# Wang2023 (Cancer Cell 2023; doi 10.1016/j.ccell.2022.12.002; PMC9839614).
# Table S1 gene-level CRISPR scores from MOLM-13 Cas9 knockout screens. Only the
# Venetoclax columns are imported. Each value is the paper's "CRISPR score":
# mean log2 fold change (CPM-normalized drug / DMSO) across a gene's sgRNAs;
# positive = enriched under venetoclax (resistance), negative = depleted
# (sensitization). The table reports no p-values or FDRs.
tidy_wang2023 <- function(root = ".") {
  source_file <- file.path(root, "Wang2023/1-s2.0-S1535610822005888-mmc2.xlsx")
  sheets <- list(
    list(
      sheet = "Genome-wide CRISPR screen",
      skip = 1, # title row above the header
      experiment = "Brunello genome-wide library venetoclax 25 nM vs DMSO day 20"
    ),
    list(
      sheet = "RBP CRISPR screen",
      skip = 0,
      experiment = "RBP domain-focused library venetoclax vs DMSO"
    )
  )
  purrr::map_dfr(sheets, function(spec) {
    readxl::read_excel(source_file, sheet = spec$sheet, skip = spec$skip) |>
      transmute(
        gene_symbol = normalize_gene_symbol(.data$`Target gene`),
        effect_value = coerce_number(.data$Venetoclax)
      ) |>
      filter(!is.na(.data$gene_symbol), !is.na(.data$effect_value)) |>
      mutate(
        experiment = spec$experiment,
        source_sheet = spec$sheet,
        effect_metric = "CRISPR score (mean sgRNA log2 fold change, venetoclax / DMSO)",
        p_value = NA_real_,
        fdr = NA_real_,
        rank = rank(-abs(.data$effect_value), ties.method = "first"),
        direction = classify_direction(.data$effect_value),
        source_level = "gene",
        source_file = source_file
      )
  }) |>
    attach_metadata("Wang2023") |>
    complete_standard()
}

# Sharon2019: author MAGeCK RRA workbook linked from BioGRID ORCS Dataset 406.
# MOLM13-R1 (acquired venetoclax-resistant MOLM-13 clone), TKO v1, venetoclax
# 400 nM vs DMSO; days 8, 16, 22 and 29 are repeat timepoints of one screen.
# Days 16/22/29 have per-day sheets; day 8 exists only in the combined sheet,
# whose day 16/22/29 columns are identical to the per-day sheets. Its "score"
# columns are MAGeCK RRA scores, so effects come only from the "lfc" columns.
# TKO control targets (EGFP, LacZ, luciferase, chr10 promoter/random) are dropped.
sharon2019_control_ids <- c("chr10Promi", "chr10Rand", "EGFP", "LacZ", "luciferase")

sharon2019_experiment <- function(day) {
  paste0("MOLM13-R1 acquired-resistant venetoclax 400 nM vs DMSO day ", day)
}

tidy_sharon2019_day8 <- function(source_file) {
  raw <- readxl::read_excel(source_file, sheet = "All days including D8")
  arm <- function(prefix, label) {
    raw |>
      transmute(
        gene_symbol = .data$Gene,
        experiment = paste0(sharon2019_experiment(8), " | ", label, " selection"),
        effect_metric = paste0("MAGeCK ", label, " selection LFC"),
        effect_value = coerce_number(.data[[paste0("Day_8_", prefix, "|lfc")]]),
        p_value = coerce_number(.data[[paste0("Day_8_", prefix, "|p-value")]]),
        fdr = coerce_number(.data[[paste0("Day_8_", prefix, "|fdr")]]),
        rank = coerce_number(.data[[paste0("Day_8_", prefix, "|rank")]]),
        direction = classify_direction(.data$effect_value),
        source_level = "gene",
        source_file = source_file,
        raw_effect_value = .data$effect_value
      )
  }
  bind_rows(arm("neg", "negative"), arm("pos", "positive")) |>
    attach_metadata("Sharon2019") |>
    complete_standard()
}

tidy_sharon2019 <- function(root = ".") {
  source_file <- file.path(root, "Sharon2019/MaGECK_D8-16-22-29.xlsx")
  per_day <- purrr::map_dfr(c(16, 22, 29), \(day) {
    tidy_mageck_gene_sheet(source_file, paste0("Day_", day), 0, "Sharon2019", sharon2019_experiment(day))
  })
  bind_rows(tidy_sharon2019_day8(source_file), per_day) |>
    filter(!.data$gene_symbol %in% sharon2019_control_ids) |>
    attach_metadata("Sharon2019") |>
    complete_standard()
}

# Excel turned MARCHn/SEPTn/SEP15/DEC1 into dates in both supplementary tables.
# The sgRNA IDs in Supplementary Table 2 keep the RefSeq accession
# (e.g. "43164_NM_017824.4_..." for "5-Mar"). Symbols below are the pre-Excel
# library symbols for those accessions (NCBI nuccore esummary, 2026-11-03), so
# gene keys match the other screens; current HGNC symbols are kept for audit.
# "1-Mar" (MARCH1 + MARC1) and "2-Mar" (MARCH2 + MARC2) each merge two genes in
# the gene table and cannot be resolved, so they are dropped.
yokoyama2026_refseq_symbols <- tibble::tribble(
  ~refseq, ~restored_symbol, ~current_symbol,
  "NM_178450", "MARCH3", "MARCHF3",
  "NM_020814", "MARCH4", "MARCHF4",
  "NM_017824", "MARCH5", "MARCHF5",
  "NM_005885", "MARCH6", "MARCHF6",
  "NM_022826", "MARCH7", "MARCHF7",
  "NM_001002266", "MARCH8", "MARCHF8",
  "NM_138396", "MARCH9", "MARCHF9",
  "NM_001100875", "MARCH10", "MARCHF10",
  "NM_001102562", "MARCH11", "MARCHF11",
  "NM_052838", "SEPT1", "SEPTIN1",
  "NM_006155", "SEPT2", "SEPTIN2",
  "NM_145733", "SEPT3", "SEPTIN3",
  "NM_001256782", "SEPT4", "SEPTIN4",
  "NM_002688", "SEPT5", "SEPTIN5",
  "NM_145802", "SEPT6", "SEPTIN6",
  "NM_001011553", "SEPT7", "SEPTIN7",
  "NM_001300799", "SEPT8", "SEPTIN8",
  "NM_001113491", "SEPT9", "SEPTIN9",
  "NM_178584", "SEPT10", "SEPTIN10",
  "NM_018243", "SEPT11", "SEPTIN11",
  "NM_144605", "SEPT12", "SEPTIN12",
  "NM_207366", "SEPT14", "SEPTIN14",
  "NM_004261", "SEP15", "SELENOF",
  "NM_017418", "DEC1", "DEC1"
)

yokoyama2026_arms <- c(
  ED_venetoclax = "ED40515(-) venetoclax vs DMSO",
  ED_navitoclax = "ED40515(-) navitoclax vs DMSO",
  ED_ABBV467 = "ED40515(-) ABBV-467 vs DMSO",
  DEL_venetoclax = "DEL venetoclax vs DMSO",
  DEL_navitoclax = "DEL navitoclax vs DMSO",
  DEL_ABBV467 = "DEL ABBV-467 vs DMSO"
)

# Supplementary Table 1: MAGeCK gene summaries for six model-drug groups vs
# DMSO. Like tidy_mageck_gene_sheet(), both arms keep the signed gene LFC
# (median sgRNA LFC). The table has no FDR column, so fdr stays NA.
tidy_yokoyama2026 <- function(root = ".") {
  source_file <- file.path(root, "Yokoyama2026/41388_2026_3970_MOESM2_ESM.xlsx")
  if (!file.exists(source_file)) return(tibble())
  date_symbol <- "^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$"

  header <- readxl::read_excel(source_file, sheet = "Supplemental_Table_1. gene", col_names = FALSE, n_max = 3, .name_repair = "minimal")
  group <- tidyr::fill(tibble(g = unlist(header[2, ])), "g")$g
  col_names <- c("id", paste(group[-1], unlist(header[3, -1]), sep = "::"))
  genes <- readxl::read_excel(source_file, sheet = "Supplemental_Table_1. gene", skip = 3, col_names = col_names, col_types = "text")

  sgrna <- readxl::read_excel(source_file, sheet = "Supplementary_Table_2. sgRNA", skip = 3, col_types = "text",
    col_names = c("Gene", "sgrna", names(yokoyama2026_arms)))

  # Restore a date-mangled row only when its sgRNAs carry one RefSeq accession
  # and their median LFC reproduces the gene LFC in all six groups.
  restored <- sgrna |>
    filter(str_detect(.data$Gene, date_symbol)) |>
    mutate(refseq = str_match(.data$sgrna, "_([NX][MR]_[0-9]+)[.][0-9]+_")[, 2]) |>
    group_by(.data$Gene) |>
    filter(n_distinct(.data$refseq) == 1) |>
    summarise(refseq = first(.data$refseq), across(all_of(names(yokoyama2026_arms)), \(x) stats::median(coerce_number(x))), .groups = "drop") |>
    tidyr::pivot_longer(all_of(names(yokoyama2026_arms)), names_to = "group", values_to = "sgrna_median") |>
    inner_join(
      genes |>
        select("id", all_of(paste0(names(yokoyama2026_arms), "::neg|lfc"))) |>
        tidyr::pivot_longer(-"id", names_to = "group", values_to = "gene_lfc") |>
        mutate(group = str_remove(.data$group, "::neg\\|lfc$"), gene_lfc = coerce_number(.data$gene_lfc)),
      by = c(Gene = "id", "group")
    ) |>
    group_by(.data$Gene, .data$refseq) |>
    summarise(max_diff = max(abs(.data$sgrna_median - .data$gene_lfc)), .groups = "drop") |>
    filter(.data$max_diff < 1e-4) |>
    inner_join(yokoyama2026_refseq_symbols, by = "refseq")

  base <- genes |>
    filter(!str_detect(.data$id, "^Control_")) |>
    left_join(restored, by = c(id = "Gene")) |>
    mutate(
      source_symbol = .data$id,
      gene_symbol = if_else(str_detect(.data$id, date_symbol), .data$restored_symbol, .data$id)
    ) |>
    filter(!is.na(.data$gene_symbol))

  purrr::imap_dfr(yokoyama2026_arms, function(experiment, group) {
    col <- \(field) coerce_number(base[[paste0(group, "::", field)]])
    lfc <- col("neg|lfc")
    arm <- function(side, label) {
      tibble(
        gene_symbol = base$gene_symbol,
        experiment = paste0(experiment, " | ", label, " selection"),
        effect_metric = paste("MAGeCK", label, "selection LFC"),
        effect_value = lfc,
        p_value = col(paste0(side, "|p_value")),
        fdr = NA_real_,
        rank = col(paste0(side, "|rank")),
        direction = classify_direction(lfc),
        source_level = "gene",
        source_file = source_file,
        raw_effect_value = lfc,
        source_sheet = "Supplemental_Table_1. gene",
        source_symbol = base$source_symbol,
        restored_refseq = base$refseq,
        restored_current_symbol = base$current_symbol
      )
    }
    bind_rows(arm("neg", "negative"), arm("pos", "positive"))
  }) |>
    attach_metadata("Yokoyama2026") |>
    complete_standard()
}

# Published gene-level edgeR summaries (Table S4, mmc2.xls): OCI-AML2 Yusa v1
# knockout library, day-21 drug vs DMSO (the workbook ReadMe says "day 23").
# effect_value is the authors' mid_lfc: the LFC of the median-p sgRNA, as in the
# Figure 3A volcano; positive = enriched under drug. The authors could not
# normalize the venetoclax and combination counts, so their mid_lfc values are
# shifted to a median near -6 and their FDRs flag most genes. For those two arms
# effect_value is median-centered and p/FDR are not used for significance; the
# published values are kept in published_effect_value/p_value/fdr.
tidy_stewart2026 <- function(root = ".") {
  source_file <- file.path(root, "Stewart2026/mmc2.xls")
  if (!file.exists(source_file)) return(tibble())
  unnormalized <- c("VenVsDMSO_gene", "ComboVsDMSO_gene")
  sheets <- c(
    "VenVsDMSO_gene" = "OCI-AML2 Yusa library day 21 venetoclax 200 nM vs DMSO",
    "ComboVsDMSO_gene" = "OCI-AML2 Yusa library day 21 venetoclax 200 nM + palbociclib 1 µM vs DMSO",
    "PalboVsDMSO_gene" = "OCI-AML2 Yusa library day 21 palbociclib 1 µM vs DMSO (not a BH3 mimetic)"
  )
  purrr::imap_dfr(sheets, function(experiment, sheet) {
    readxl::read_excel(source_file, sheet = sheet) |>
      transmute(
        gene_symbol = normalize_gene_symbol(.data$Gene),
        effect_value = coerce_number(.data$mid_lfc),
        p_value = coerce_number(.data$mid_pval),
        fdr = coerce_number(.data$mid_fdr),
        n_sgrna = coerce_number(.data$total),
        n_sgrna_tested = coerce_number(.data$not_filtered),
        n_sgrna_significant = coerce_number(.data$num_sig),
        n_sgrna_up = coerce_number(.data$num_up),
        n_sgrna_down = coerce_number(.data$num_down),
        sgrna_lfc_mean = coerce_number(.data$lfc_mean),
        sgrna_min_fdr = coerce_number(.data$min_fdr),
        tier = .data$tiers,
        camera_dir = .data$camera_dir,
        camera_p_value = coerce_number(.data$camera_pval),
        camera_fdr = coerce_number(.data$camera_fdr)
      ) |>
      mutate(
        median_centered_mid_lfc = .data$effect_value - stats::median(.data$effect_value, na.rm = TRUE),
        published_effect_value = .data$effect_value,
        published_p_value = .data$p_value,
        published_fdr = .data$fdr,
        rank = rank(.data$p_value, ties.method = "first", na.last = "keep")
      ) |>
      mutate(
        effect_value = if (sheet %in% unnormalized) .data$median_centered_mid_lfc else .data$effect_value,
        p_value = if (sheet %in% unnormalized) NA_real_ else .data$p_value,
        fdr = if (sheet %in% unnormalized) NA_real_ else .data$fdr,
        experiment = experiment,
        effect_metric = if (sheet %in% unnormalized) {
          "edgeR median-p sgRNA log2 fold-change (mid_lfc), median-centered; published counts unnormalized"
        } else {
          "edgeR median-p sgRNA log2 fold-change (mid_lfc)"
        },
        direction = classify_direction(.data$effect_value),
        source_level = "gene",
        source_file = source_file
      )
  }) |>
    attach_metadata("Stewart2026") |>
    complete_standard()
}

# Diepstraten2026 (Sci Adv 2026; doi 10.1126/sciadv.aec0722; PMC13588164).
# Partita mouse CRISPRa library (10 sgRNAs per TSS, five sublibraries each
# transduced separately) in Eu-Myc/dCas9a-SAM/sgBcl2 DHL214 and DHL216 lymphoma
# lines; venetoclax IC50 or IC80 vs DMSO, 2-3 replicates per sublibrary.
# Source: GEO GSE296642_DHL_Partita1-5_ven_counts.xlsx (raw sgRNA counts).
# The authors ran no statistical test for the venetoclax screen; they plotted
# per-sample sgRNA proportions and labelled sgRNAs with frequency >= 5%
# (Fig. 3B legend). Post-selection samples are severely bottlenecked
# (DMSO 0.3-5% of sgRNAs detected, venetoclax 0.1-1%; top 10 sgRNAs hold
# 68-100% of venetoclax reads), so no genome-wide gene score is computed.
# This is a hit-list REANALYSIS of the counts: an sgRNA is a hit when it
# reaches >= 5% of reads in at least one venetoclax replicate (authors'
# threshold) and its mean CPM is higher under venetoclax than DMSO.
# Sublibraries are pooled within a line x dose (library parts are not separate
# screens). Each gene key keeps its strongest sgRNA (highest venetoclax
# frequency); `enriched_symbols` records merged paralogs. Sign: positive =
# activating the gene enriched under venetoclax (resistance). No p/FDR.
tidy_diepstraten2026 <- function(root = ".") {
  source_file <- file.path(root, "Diepstraten2026/GSE296642_DHL_Partita1-5_ven_counts.xlsx")
  if (!file.exists(source_file)) return(tibble())
  hit_frequency <- 0.05
  experiments <- c(
    "214_IC50" = "CRISPRa Partita DHL214 venetoclax IC50 (10 nM) vs DMSO | sgRNAs >=5% frequency (hit list)",
    "214_IC80" = "CRISPRa Partita DHL214 venetoclax IC80 (30 nM) vs DMSO | sgRNAs >=5% frequency (hit list)",
    "216_IC50" = "CRISPRa Partita DHL216 venetoclax IC50 (5 nM) vs DMSO | sgRNAs >=5% frequency (hit list)",
    "216_IC80" = "CRISPRa Partita DHL216 venetoclax IC80 (15 nM) vs DMSO | sgRNAs >=5% frequency (hit list)"
  )
  counts <- purrr::map_dfr(readxl::excel_sheets(source_file), function(sheet) {
    d <- readxl::read_excel(source_file, sheet = sheet)
    names(d)[tolower(names(d)) == "gene"] <- "gene"
    d |>
      select(-any_of("water")) |>
      tidyr::pivot_longer(-c("sgRNA", "gene"), names_to = "sample", values_to = "count") |>
      mutate(sublibrary = sheet)
  }) |>
    mutate(
      count = dplyr::coalesce(as.numeric(.data$count), 0),
      line = stringr::str_extract(.data$sample, "^\\d{3}"),
      condition = stringr::str_match(.data$sample, "^\\d{3}_(.*)_Lib")[, 2] |>
        stringr::str_remove("^s_") |>
        stringr::str_to_upper()
    ) |>
    filter(.data$condition %in% c("DMSO", "IC50", "IC80")) |>
    group_by(.data$sample, .data$sublibrary) |>
    mutate(cpm = .data$count / sum(.data$count) * 1e6) |>
    ungroup()

  dmso <- counts |>
    filter(.data$condition == "DMSO") |>
    group_by(.data$sublibrary, .data$line, .data$sgRNA) |>
    summarise(dmso_cpm = mean(.data$cpm), .groups = "drop")

  counts |>
    filter(.data$condition != "DMSO") |>
    group_by(.data$sublibrary, .data$line, .data$condition, .data$sgRNA, .data$gene) |>
    summarise(
      ven_cpm = mean(.data$cpm),
      max_sgrna_frequency = max(.data$cpm) / 1e6,
      n_replicates_ge5pct = sum(.data$cpm / 1e6 >= hit_frequency),
      n_replicates = n(),
      .groups = "drop"
    ) |>
    left_join(dmso, by = c("sublibrary", "line", "sgRNA")) |>
    mutate(effect_value = log2((.data$ven_cpm + 1) / (.data$dmso_cpm + 1))) |>
    filter(
      .data$max_sgrna_frequency >= hit_frequency,
      .data$effect_value > 0,
      !grepl("^(NTC|Non|non)", .data$gene)
    ) |>
    mutate(
      gene_symbol = normalize_gene_symbol(.data$gene),
      gene_key = gene_key(.data$gene_symbol),
      experiment = unname(experiments[paste(.data$line, .data$condition, sep = "_")])
    ) |>
    group_by(.data$experiment, .data$gene_key) |>
    arrange(desc(.data$max_sgrna_frequency), .by_group = TRUE) |>
    mutate(
      n_enriched_sgrna = n(),
      enriched_symbols = paste(sort(unique(.data$gene_symbol)), collapse = "; ")
    ) |>
    slice(1) |>
    group_by(.data$experiment) |>
    mutate(rank = rank(-.data$max_sgrna_frequency, ties.method = "first")) |>
    ungroup() |>
    transmute(
      experiment = .data$experiment,
      gene_symbol = .data$gene_symbol,
      gene_key = .data$gene_key,
      effect_metric = "log2((mean venetoclax CPM + 1) / (mean DMSO CPM + 1)), strongest sgRNA reaching >=5% read frequency (reanalysis of GEO counts)",
      effect_value = .data$effect_value,
      p_value = NA_real_,
      fdr = NA_real_,
      rank = .data$rank,
      direction = classify_direction(.data$effect_value),
      source_level = "reanalysed sgRNA counts; authors' >=5% frequency threshold (hit list)",
      source_file = source_file,
      sgrna = .data$sgRNA,
      sublibrary = .data$sublibrary,
      max_sgrna_frequency = .data$max_sgrna_frequency,
      n_replicates_ge5pct = .data$n_replicates_ge5pct,
      n_replicates = .data$n_replicates,
      n_enriched_sgrna = .data$n_enriched_sgrna,
      enriched_symbols = .data$enriched_symbols
    ) |>
    attach_metadata("Diepstraten2026") |>
    complete_standard()
}

# OCI-AML2 Yusa v1 knockout screen, AZD4320 30 nM vs DMSO for 14 days
# (GSE294240). The paper's gene-level table (Supplemental Data 1, "CRISPR
# Screen Results.xlsx") is not in the published supplement, so gene results are
# re-derived from the GEO sgRNA counts by scripts/build_brimedwards2026.R, which
# follows the authors' edgeR pipeline (plasmid and CPM filters, TMM, guide
# exactTest, median-p sgRNA per gene, concordance tiers). effect_value is that
# guide's log2 fold-change (mid_lfc), as in the Figure 2A volcano; positive =
# enriched under AZD4320 = knockout confers resistance. p_value/fdr are the
# median guide's p-value and BH FDR, as in the Stewart2026 gene tabs. The
# GEO series has no venetoclax arm; the paper's venetoclax comparison uses an
# earlier OCI-AML2 screen that is not part of this deposit.
tidy_brimedwards2026 <- function(root = ".") {
  source_file <- file.path(root, "BrimEdwards2026/GSE294240_azd4320_crispr_count.txt.gz")
  result_file <- file.path(root, "BrimEdwards2026/brimedwards2026_gene_results.tsv")
  if (!file.exists(source_file)) return(tibble())
  if (!file.exists(result_file) || file.info(result_file)$mtime < file.info(source_file)$mtime) {
    stop("BrimEdwards2026 derived results are missing or older than the GEO counts. Run Rscript --vanilla scripts/build_brimedwards2026.R . (requires edgeR) before rebuilding browser data.")
  }

  readr::read_tsv(result_file, show_col_types = FALSE) |>
    transmute(
      gene_symbol = normalize_gene_symbol(.data$Gene),
      experiment = "OCI-AML2 Yusa library day 14 AZD4320 30 nM vs DMSO",
      effect_metric = "edgeR median-p sgRNA log2 fold-change (mid_lfc), re-derived from GEO counts",
      effect_value = coerce_number(.data$mid_lfc),
      p_value = coerce_number(.data$mid_pval),
      fdr = coerce_number(.data$mid_fdr),
      n_sgrna = coerce_number(.data$total),
      n_sgrna_tested = coerce_number(.data$not_filtered),
      n_sgrna_significant = coerce_number(.data$num_sig),
      n_sgrna_up = coerce_number(.data$num_up),
      n_sgrna_down = coerce_number(.data$num_down),
      sgrna_lfc_mean = coerce_number(.data$lfc_mean),
      sgrna_min_fdr = coerce_number(.data$min_fdr),
      tier = .data$tiers
    ) |>
    # Genes with no sgRNA passing the filters have no score.
    filter(!is.na(.data$effect_value)) |>
    mutate(
      rank = rank(.data$p_value, ties.method = "first", na.last = "keep"),
      direction = classify_direction(.data$effect_value),
      source_level = "re-analyzed gene summary from published sgRNA counts",
      source_file = source_file
    ) |>
    attach_metadata("BrimEdwards2026") |>
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
    tidy_deng2022_crispra(root),
    tidy_portelinha2024(root),
    tidy_eide2023_crispr(root),
    tidy_he2026(root),
    tidy_nakao2023(root),
    tidy_wang2023(root),
    tidy_sharon2019(root),
    tidy_yokoyama2026(root),
    tidy_stewart2026(root),
    tidy_diepstraten2026(root),
    tidy_brimedwards2026(root)
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
    source_file <- file.path(root, "data/external/deng2022/41467_2022_32485_MOESM5_ESM.xlsx")
    if (!file.exists(source_file)) return(missing_detail_note)
    map_dfr(c("214 Enriched sgRNAs", "216 Enriched sgRNAs"), \(sheet) {
      readxl::read_excel(source_file, sheet = sheet) |>
        mutate(gene_symbol = .data$Symbol, source_sheet = paste("Supplementary Data 2", sheet))
    }) |>
      mutate(gene_key = gene_key(.data$gene_symbol)) |>
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
