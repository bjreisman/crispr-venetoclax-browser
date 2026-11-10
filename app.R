library(shiny)
library(bslib)
library(dplyr)
library(DT)
library(plotly)
library(readr)

source("R/load_data.R")
source("R/plots.R")
source("R/screen_views.R")

scores <- slim_browser_scores(load_browser_data())
depmap_crispr <- load_depmap_crispr()
external_context <- load_external_context()
metadata <- paper_metadata
additional_screen_candidates <- readRDS("data/processed/additional_screen_candidates.rds")
default_gene_keys <- intersect(c("BAX", "TP53", "PMAIP1"), unique(scores$gene_key))
max_gene_plot_rows <- 5000
max_paper_table_rows <- 10000
drug_response <- external_context$drug_response
beataml_expression <- external_context$beataml_expression
beataml_gene_sets <- external_context$beataml_gene_sets
eide_crispr_context <- external_context$eide_crispr
zhang_library <- external_context$zhang_library

screen_scores <- build_screen_scores(scores) |>
  mutate(row_level = 0L, row_key = .data$screen_id, baseline = FALSE)
screen_detail_scores <- build_screen_scores(scores, build_screen_detail_registry(scores)) |>
  mutate(row_level = 1L, row_key = .data$screen_id, baseline = FALSE)
expandable_papers <- sort(unique(screen_detail_scores$paper))
screen_detail_parents <- unique(screen_detail_scores$parent_screen_id)
paper_results <- build_paper_result_registry(scores, build_screen_detail_registry(scores))
# Percentile is the x-axis of the Paper Explorer plot, so it is not a y option.
paper_metric_choices <- screen_metric_choices[screen_metric_choices != "percentile"]
screen_extra_group_choices <- c(
  "Baseline dependency (DepMap, Lin2019)" = baseline_screen_group,
  "Other BH3 mimetics" = "Other BH3 mimetics",
  "CRISPRa, ORF, and BH3-priming assays" = "Other assays (direction not comparable)"
)
baseline_scores <- build_baseline_scores(scores, depmap_crispr) |>
  mutate(row_level = 0L, row_key = .data$screen_id)

experiment_stats <- scores |>
  group_by(.data$experiment_id) |>
  summarise(
    effect_mean = mean(.data$effect_value, na.rm = TRUE),
    effect_sd   = sd(.data$effect_value, na.rm = TRUE),
    .groups = "drop"
  )

collapse_selection_pairs <- function(df) {
  df |>
    mutate(
      experiment_base = sub(" \\| (positive|negative) selection$", "", .data$experiment),
      selection_pair = grepl(" \\| (positive|negative) selection$", .data$experiment),
      evidence_value = dplyr::case_when(
        !is.na(.data$fdr) ~ .data$fdr,
        !is.na(.data$p_value) ~ .data$p_value,
        TRUE ~ Inf
      ),
      evidence_rank = dplyr::case_when(
        is.finite(.data$evidence_value) ~ .data$evidence_value,
        !is.na(.data$rank) ~ .data$rank / 1e6,
        TRUE ~ Inf
      )
    ) |>
    group_by(.data$paper, .data$experiment_base, .data$gene_key) |>
    arrange(.data$evidence_rank, desc(abs(.data$effect_value)), .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    select(-any_of(c("experiment_base", "selection_pair", "evidence_value", "evidence_rank")))
}

ui <- page_navbar(
  title = "CRISPR Venetoclax Resistance Browser",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  fillable = FALSE,
  header = tags$head(
    tags$style(HTML("
      .sidebar-card { border: 1px solid #d9e2e7; border-radius: 8px; padding: 14px; background: #fbfcfd; }
      .metric-strip { display: grid; grid-template-columns: repeat(4, minmax(120px, 1fr)); gap: 10px; margin-bottom: 12px; }
      .metric { border: 1px solid #d9e2e7; border-radius: 8px; padding: 10px 12px; background: white; }
      .metric .value { font-size: 1.35rem; font-weight: 700; }
      .small-muted { color: #58656d; font-size: 0.9rem; }
      .rationale-box { padding: 4px 4px 0 4px; }
      .paper-hit-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 16px; }
      .paper-hit-grid h3 { font-size: 1.02rem; margin: 0 0 6px 0; }
      .paper-screens .grp { color: #58656d; font-size: 0.8rem; margin: 10px 0 4px 0; }
      .paper-screens .grp:first-child { margin-top: 2px; }
      .paper-screens .screen-btn { display: block; width: 100%; text-align: left; border: 1px solid #c9d3d9; border-radius: 6px; background: white; padding: 6px 9px; margin-bottom: 6px; font-size: 0.88rem; line-height: 1.3; cursor: pointer; }
      .paper-screens .screen-btn:hover { border-color: #8aa0ad; }
      .paper-screens .screen-btn.on { border-color: #2a6fb0; background: #eaf2fa; color: #1d4f80; }
      .paper-screens .chips { display: flex; flex-wrap: wrap; gap: 4px; margin-top: 5px; }
      .paper-screens .chip { font-size: 0.76rem; padding: 1px 7px; border: 1px solid #c9d3d9; border-radius: 5px; background: white; color: #58656d; cursor: pointer; }
      .paper-screens .chip.on { background: #2a6fb0; border-color: #2a6fb0; color: white; }
      .paper-screens .note { color: #58656d; font-size: 0.78rem; margin: 2px 0 10px 0; }
      .metric .detail { color: #58656d; font-size: 0.85rem; line-height: 1.35; }
      .rationale-box p { font-size: 1.02rem; line-height: 1.55; max-width: 980px; }
      .rationale-box strong { font-size: 1.08rem; }
      .baseline-note { color: #58656d; margin: 0 0 14px 0; max-width: 1000px; }
      .baseline-grid { display: grid; grid-template-columns: 1fr; gap: 16px; align-items: start; }
      .baseline-panel { border: 1px solid #d9e2e7; border-radius: 8px; padding: 12px; background: #ffffff; }
      .baseline-panel h3 { font-size: 1.05rem; margin: 0 0 4px 0; }
      .baseline-panel .small-muted { margin-bottom: 8px; }
      .baseline-empty { color: #58656d; margin: 12px 0 0 0; }
      .gene-summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 10px; margin: 4px 0 12px 0; }
      .gene-summary { border: 1px solid #d9e2e7; border-radius: 8px; padding: 9px 12px; background: #fbfcfd; }
      .gene-summary .gene { font-weight: 700; }
      .gene-summary .detail { color: #58656d; font-size: 0.9rem; line-height: 1.4; }
      @media (max-width: 800px) { .metric-strip { grid-template-columns: repeat(2, minmax(120px, 1fr)); } }
    ")),
    tags$script(HTML("
      $(document).on('click', '.paper-screens [data-result]', function(e) {
        e.stopPropagation();
        Shiny.setInputValue('paper_result_click', $(this).attr('data-result'), {priority: 'event'});
      });
    "))
  ),
  nav_panel(
    "Gene Explorer",
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        div(
          class = "sidebar-card",
          selectizeInput(
            "gene_select",
            "Genes",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Search gene symbols",
              maxOptions = 250,
              plugins = list("remove_button")
            )
          ),
          selectInput("gene_papers", "Papers", choices = sort(unique(scores$paper)), selected = sort(unique(scores$paper)), multiple = TRUE),
          selectInput("depmap_cell_lines", "DepMap cell lines", choices = sort(unique(stats::na.omit(depmap_crispr$cell_line))), selected = sort(unique(stats::na.omit(depmap_crispr$cell_line))), multiple = TRUE),
          radioButtons(
            "gene_view_mode",
            "View",
            choices = c("One row per screen" = "screens", "All results" = "all"),
            selected = "screens"
          ),
          conditionalPanel(
            condition = "input.gene_view_mode == 'screens'",
            radioButtons(
              "screen_metric",
              "Metric",
              choices = screen_metric_choices,
              selected = "percentile"
            ),
            checkboxGroupInput(
              "screen_extra_groups",
              "Also show",
              choices = screen_extra_group_choices,
              selected = baseline_screen_group
            ),
            selectizeInput(
              "screen_expand_papers",
              "Expand papers",
              choices = expandable_papers,
              selected = character(),
              multiple = TRUE,
              options = list(placeholder = "Click a dot, or pick papers", plugins = list("remove_button"))
            )
          ),
          conditionalPanel(
            condition = "input.gene_view_mode == 'all'",
            selectInput(
              "effect_axis_type",
              "Cross-paper axis",
              choices = c(
                "Signed −log₁₀(FDR / p-value)" = "signed_log10",
                "Per-experiment z-score"                      = "zscore",
                "Raw effect value (paper units)"              = "raw"
              ),
              selected = "signed_log10"
            ),
            checkboxInput("gene_significant_only", "Show significant/annotated hits only", FALSE)
          ),
          downloadButton("download_gene_view", "Download view")
        )
      ),
      card(
        card_header("Cross-paper gene effects"),
        uiOutput("gene_screen_summary"),
        conditionalPanel(
          condition = "input.gene_view_mode == 'screens'",
          uiOutput("gene_screen_plot_ui"),
          tags$div(
            class = "small-muted",
            style = "display: flex; gap: 18px; flex-wrap: wrap; margin: 0 0 6px 0;",
            tags$span(HTML("<span style='color:#2a6fb0'>&#9679;</span> FDR ≤ 0.05 or paper threshold")),
            tags$span(HTML("<span style='color:#2a6fb0'>&#9675;</span> Not significant or no statistics")),
            tags$span(HTML("<span style='color:#2a6fb0'>&#9670;</span> In published hit list")),
            tags$span(HTML("<span style='color:#6b6b6b'>&#9670;</span> Baseline dependency (DepMap gene effect ≤ −0.5)")),
            tags$span(HTML("<span style='color:#6b6b6b'>&#9671;</span> Baseline: no dependency, or no threshold (Lin2019)")),
            uiOutput("gene_screen_band_note", inline = TRUE),
            tags$span("Blank: gene not measured or no value on this metric")
          ),
          uiOutput("gene_screen_hidden_note"),
          tags$p(
            class = "small-muted",
            "One representative result per independent screen: the latest timepoint, and for MAGeCK the selection arm matching the gene's direction. ",
            "Percentiles rank the gene among all genes in that screen. Hit-list tables are not genome-wide, so listed genes are shown without a percentile. ",
            "Papers marked ▸ have more results (timepoints, replicates, doses, derived tables): click one of their dots, or pick them under Expand papers, to show each result as an indented row. ",
            "Switch to All results for every arm and table in paper-native units."
          )
        ),
        conditionalPanel(
          condition = "input.gene_view_mode == 'all'",
          plotlyOutput("gene_effect_plot", height = "1520px"),
          uiOutput("gene_effect_na_note")
        )
      ),
      card(
        card_header("Evidence table"),
        DTOutput("gene_table")
      ),
      accordion(
        open = FALSE,
        accordion_panel(
          "Baseline Dependency Context",
          tags$p(
            class = "baseline-note",
            "Baseline CRISPR KO and Lin2019 essentiality scores for the selected genes. These contextualize general growth dependency and are separate from venetoclax-selection effects."
          ),
          div(
            class = "baseline-grid",
            div(
              class = "baseline-panel",
              tags$h3("DepMap CRISPR KO"),
              tags$div(class = "small-muted", "More negative gene effect means stronger baseline dependency."),
              uiOutput("depmap_crispr_plot_ui")
            ),
            div(
              class = "baseline-panel",
              tags$h3("Lin2019 Essentiality"),
              tags$div(class = "small-muted", "Published THP-1, MOLM-13, and average essentiality three-scores."),
              uiOutput("lin2019_essentiality_plot_ui")
            )
          )
        ),
        accordion_panel(
          "Source detail for selected gene and paper",
          fluidRow(
            column(6, selectizeInput("detail_gene", "Gene", choices = NULL, multiple = FALSE, options = list(placeholder = "Search one gene", maxOptions = 250))),
            column(6, selectInput("detail_paper", "Paper", choices = sort(unique(scores$paper))))
          ),
          DTOutput("detail_table")
        )
      )
    )
  ),
  nav_panel(
    "Paper Explorer",
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        div(
          class = "sidebar-card",
          selectInput("paper_select", "Paper", choices = sort(unique(scores$paper))),
          tags$label(class = "control-label", "Screen"),
          uiOutput("paper_screen_buttons"),
          selectInput("paper_metric", "Metric", choices = paper_metric_choices, selected = "robust_z"),
          sliderInput("paper_top_n", "Top hits labelled each way", min = 0, max = 25, value = 10, step = 1),
          selectizeInput(
            "paper_highlight_genes",
            "Highlight/add genes",
            choices = NULL,
            multiple = TRUE,
            options = list(
              placeholder = "Search gene symbols",
              maxOptions = 250,
              plugins = list("remove_button")
            )
          ),
          downloadButton("download_paper_view", "Download view")
        )
      ),
      card(
        card_header("Study rationale"),
        uiOutput("paper_rationale")
      ),
      uiOutput("paper_summary"),
      card(
        card_header("Whole-screen view"),
        tags$p(class = "small-muted", "Every gene in the result, ordered by within-screen percentile. Shaded bands mark the metric's signal threshold; the strongest hits each way and BAX, PMAIP1, and TP53 are labelled. Hover a point for details."),
        plotlyOutput("paper_rank_plot", height = "620px")
      ),
      card(
        card_header("Top hits"),
        tags$p(class = "small-muted", "“Other venetoclax screens” counts the independent venetoclax screens from other papers where the gene is in the same tail (top or bottom 10%, or a published hit list), with how many of those reach FDR ≤ 0.05."),
        div(
          class = "paper-hit-grid",
          div(tags$h3(style = "color:#b73e3e", "Resistance (enriched)"), DTOutput("paper_hits_resistance")),
          div(tags$h3(style = "color:#247a6b", "Sensitizing (depleted)"), DTOutput("paper_hits_sensitizing"))
        )
      ),
      accordion(
        open = FALSE,
        accordion_panel("Full table as published (every arm and column)", DTOutput("paper_table"))
      )
    )
  ),
  nav_panel(
    "Cross-study Overview",
    tags$p(class = "baseline-note", "Exploratory overview: recurrence counts dataset presence, not independently replicated significant sensitization. Assay types and score scales differ; see About Data for interpretation limits."),
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        div(
          class = "sidebar-card",
          selectInput("overview_papers", "Papers", choices = sort(unique(scores$paper)), selected = sort(unique(scores$paper)), multiple = TRUE),
          selectInput("overview_direction", "Direction", choices = c("All", sort(unique(scores$direction))), selected = "All"),
          sliderInput("overview_top_n", "Genes in heatmap", min = 10, max = 80, value = 40, step = 5),
          downloadButton("download_overview", "Download recurrent hits")
        )
      ),
      div(
        class = "metric-strip",
        div(class = "metric", div(class = "small-muted", "Rows"), div(class = "value", textOutput("metric_rows", inline = TRUE))),
        div(class = "metric", div(class = "small-muted", "Genes"), div(class = "value", textOutput("metric_genes", inline = TRUE))),
        div(class = "metric", div(class = "small-muted", "Papers"), div(class = "value", textOutput("metric_papers", inline = TRUE))),
        div(class = "metric", div(class = "small-muted", "Experiments"), div(class = "value", textOutput("metric_experiments", inline = TRUE)))
      ),
      card(
        card_header("Recurrent gene matrix"),
        plotlyOutput("overview_heatmap", height = "720px")
      ),
      card(
        card_header("Top recurrent hits"),
        DTOutput("overview_table")
      )
    )
  ),
  nav_panel(
    "External Context",
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        div(
          class = "sidebar-card",
          selectInput("context_datasets", "Drug-response datasets", choices = sort(unique(drug_response$dataset)), selected = sort(unique(drug_response$dataset)), multiple = TRUE),
          selectizeInput(
            "context_drugs",
            "Venetoclax-containing treatments",
            choices = NULL,
            selected = character(),
            multiple = TRUE,
            options = list(placeholder = "Search treatments", maxOptions = 250, plugins = list("remove_button"))
          ),
          radioButtons(
            "context_expression_mode",
            "Expression feature",
            choices = c("Gene" = "gene", "Gene set" = "gene_set"),
            selected = "gene",
            inline = TRUE
          ),
          conditionalPanel(
            condition = "input.context_expression_mode == 'gene'",
            selectizeInput(
              "context_expression_gene",
              "Expression gene",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = list(placeholder = "Search BeatAML expression genes", maxOptions = 250)
            )
          ),
          conditionalPanel(
            condition = "input.context_expression_mode == 'gene_set'",
            selectizeInput(
              "context_gene_set",
              "Gene set",
              choices = NULL,
              selected = NULL,
              multiple = FALSE,
              options = list(placeholder = "Search curated gene sets", maxOptions = 100)
            )
          ),
          selectInput("context_crispr_experiment", "Eide2023 CRISPR experiment", choices = sort(unique(eide_crispr_context$experiment))),
          sliderInput("context_crispr_top_n", "CRISPR genes shown", min = 10, max = 80, value = 30, step = 5),
          selectizeInput(
            "zhang_gene",
            "Zhang2024 library gene",
            choices = NULL,
            multiple = FALSE,
            options = list(placeholder = "Search library genes", maxOptions = 250)
          ),
          downloadButton("download_context_drug_response", "Download drug response")
        )
      ),
      card(
        card_header("BeatAML2 expression/gene-set score vs venetoclax response"),
        uiOutput("expression_auc_summary"),
        plotlyOutput("expression_auc_plot", height = "780px"),
        DTOutput("expression_auc_table")
      ),
      card(
        card_header("BeatAML2 and Eide2023 ex vivo drug response distributions"),
        tags$p(class = "baseline-note", "Primary-sample drug response context. Lower AUC indicates greater ex vivo sensitivity; these are not CRISPR screen effects. BeatAML2 expression-vs-AUC currently uses BeatAML2 venetoclax response because matched expression is public in that dataset."),
        plotlyOutput("context_drug_plot", height = "860px"),
        DTOutput("context_drug_table")
      ),
      card(
        card_header("Eide2023 GEO CRISPR count summary"),
        tags$p(class = "baseline-note", "Gene-level summaries computed from depth-normalized public sgRNA counts as median sgRNA log2 fold-change versus DMSO, centered on the experiment median. No p-value or FDR is inferred here."),
        plotlyOutput("context_crispr_plot", height = "960px"),
        DTOutput("context_crispr_table")
      ),
      card(
        card_header("Zhang2024 CRISPR library coverage"),
        tags$p(class = "baseline-note", "The public Zhang2024 XLSX contains the screen library design rather than a scored hit table; this panel records which genes were targeted and how many guides were listed."),
        DTOutput("zhang_library_table")
      )
    )
  ),
  nav_panel(
    "About Data",
    card(
      card_header("Work in progress: interpretation limits"),
      tags$p("This browser supports exploration of published datasets. It does not perform a validated pooled effect-size meta-analysis or establish therapeutic efficacy."),
      tags$ul(
        tags$li("Recurrence reflects gene presence across selected datasets, rather than significant, directionally consistent replication across independent studies."),
        tags$li("The significant/annotated flag includes nominal p-values even when a reported FDR exceeds 0.05. Missing FDR values can pass the FDR filter."),
        tags$li("Knockout, ORF overexpression, expression, baseline dependency, and BH3-priming results have different interpretations. A centered z-score need not retain biological direction."),
        tags$li("Eide count summaries are depth-normalized and median-centered but have not undergone a replicate-aware reanalysis. Patient correlations can include repeated specimens or measurements. Symbol matching is not formal human-mouse orthology mapping.")
      )
    ),
    card(
      card_header("Included studies"),
      DTOutput("metadata_table")
    ),
    card(
      card_header("Additional screens to consider"),
      tags$p(class = "baseline-note", "Literature audit checked 31 October 2026, starting with the EndNote library and extending to PubMed, article supplements, and public data repositories. Candidates include direct venetoclax screens, related BH3-mimetic screens, supporting context, and reanalyses. These entries are not included in the gene-score comparisons or screen counts. Data availability and remaining checks are recorded separately."),
      downloadButton("download_additional_screens", "Download full candidate catalog"),
      DTOutput("additional_screen_table")
    ),
    card(
      card_header("Data notes"),
      tags$ul(
        tags$li("Published gene-level scores are the primary cross-study layer."),
        tags$li("Raw count and sgRNA-level tables are exposed in source detail where local files support gene-level lookup."),
        tags$li("DepMap CRISPR KO gene-effect scores are from DepMap Public 26Q1 CRISPRGeneEffect.csv; more negative values indicate stronger baseline dependency."),
        tags$li("DepMap 26Q1 direct matches are available for MOLM-13, MOLM-14, MV4-11, OCI-AML2, THP-1, and OCI-LY19. OCI-Ly1 was not found in the 26Q1 model metadata."),
        tags$li("BeatAML2 and Eide2023 ex vivo drug-response data are contextual patient/sample sensitivity measurements; lower AUC means greater sensitivity and these are not knockout effect sizes."),
        tags$li("BeatAML2 curated gene-set scores are contextual RNA-seq signatures computed as the mean of per-gene z-scored expression across member genes present in the public BeatAML2 expression matrix."),
        tags$li("Eide2023 CRISPR rows are computed from public GEO sgRNA count tables: counts are scaled to the median library size, sgRNA log2 fold-changes versus DMSO are summarized as the gene median, and gene effects are centered on the experiment median (the library has no non-targeting guides); no p-value or FDR is inferred from those counts in this browser."),
        tags$li("He2026 includes the Figure 1 source-data combined MAGeCK gene table with p-values, plus separate P1-R and P2-R sgRNA RPM-derived LFC summaries without p-values/FDRs. LFC is represented as CytoC-high/CytoC-low, so negative values indicate CytoC-low enrichment/sensitizing knockout candidates such as ADSS2."),
        tags$li("Nakao2023 uses a reproducible DrugZ re-analysis of author-supplied sgRNA counts for March5-WT venetoclax, March5-WT S63845, and March5-KO venetoclax screens. Displayed effect is DrugZ normalized Z for drug vs DMSO; positive values indicate resistance/enrichment and negative values sensitization/depletion. Direction-aware DrugZ p-values/FDRs are displayed, while MAGeCK MLE beta statistics remain audit fields in downloads. Non-targeting-control and miRNA pseudo-genes are excluded from browser rows."),
        tags$li("Wang2023 contributes two MOLM-13 venetoclax knockout screens (Brunello genome-wide, and an RNA-binding-protein domain library) as published CRISPR scores, mean sgRNA log2 fold change vs DMSO, without p-values or FDRs."),
        tags$li("Sharon2019 is a genome-wide knockout screen in MOLM13-R1, a MOLM-13 clone with acquired venetoclax resistance, from the authors' MAGeCK workbook (days 8 to 29; the screen view uses day 29). It identifies knockouts that restore sensitivity, such as mitochondrial translation genes."),
        tags$li("Yokoyama2026 contributes peripheral T-cell lymphoma screens (ED40515(−), DEL) under venetoclax, navitoclax, and the MCL1 inhibitor ABBV-467, as MAGeCK LFC and p-values without FDR. Excel-date-corrupted gene IDs were restored from RefSeq accessions where unambiguous."),
        tags$li("Stewart2026 contributes OCI-AML2 screens for venetoclax, venetoclax + palbociclib, and palbociclib alone (edgeR, day 21 vs DMSO). The authors could not normalize the venetoclax and combination contrasts (raw values shifted to a median near −6; FDRs flag most genes), so the browser median-centers those two arms and does not count their FDRs as significant. The published values appear under the Raw effect (as published) and Signed −log₁₀ FDR metrics and in hover text."),
        tags$li("BrimEdwards2026 contributes an OCI-AML2 knockout screen with the dual BCL2/BCL-XL inhibitor AZD4320, re-derived from GEO counts with the authors' edgeR pipeline because the published gene table is not in the supplement."),
        tags$li("Diepstraten2026 contributes Partita CRISPRa venetoclax screens in DHL214 and DHL216 mouse lymphoma as hit lists: sgRNAs reaching at least 5% read frequency (the authors' Figure 3B threshold), reanalysed from GEO counts without p-values or FDRs. Many hits rest on a single replicate."),
        tags$li("Zhang2024 is included as CRISPR library coverage context because the public XLSX contains guide design rather than a scored genome-screen hit table."),
        tags$li("Lin2019 significance is annotated from the paper's published threshold for the average drug-sensitivity score: Avg_sens < -0.914. The supplementary table does not include per-gene p-values or FDRs."),
        tags$li("Deng2022 CRISPR activation (CRISPRa) screens are the published enriched-sgRNA tables (Supplementary Data 2: IC80 venetoclax vs DMSO in mouse lymphoma lines #214 and #216), shown as gain-of-function hit lists with each gene's strongest guide. The paper's lymphoma vs B-ALL RNA-seq comparison is not included because it compares tumour types rather than measuring venetoclax response."),
        tags$li("Gene matching uses normalized symbols and a small curated alias map: ADSS is grouped with ADSS2, ADSS1 with ADSSL1, TOMM70A with TOMM70, and mouse Bcl2a1a/b/c/d (A1) with human BCL2A1. ADSSL1 remains separate from ADSS2. No formal human-mouse ortholog mapping is applied."),
        tags$li(
          tags$strong("Guièze2019 ORF screen:"),
          " Hits are defined as genes with at least one ORF construct exceeding log₂ fold-change ≥ 2 in venetoclax vs. DMSO (replicating the paper's own threshold: ",
          tags$em("\"71 ORFs arising from 70 genes had an LFC greater than 2\""),
          ", Guièze et al. 2019). The library contained 17,255 barcoded ORFs covering 12,952 genes (~1.3 constructs/gene on average), so most genes were represented by a single construct. No cross-construct aggregation was performed in the original study or here. The ",
          tags$code("n_orf_constructs"),
          " field in the source detail view shows how many constructs were present per gene; genes supported by only one construct have no internal replication."
        )
      )
    )
  ),
  nav_panel(
    "About",
    card(
      card_header("About this browser"),
      tags$p("An interactive browser for CRISPR screen data from studies of venetoclax resistance and sensitivity across hematologic malignancy models. Gene-level scores from published supplementary tables are harmonized into a common format to support cross-study comparison."),
      tags$p(tags$strong("Contact: "), "Benjamin Reisman, Johns Hopkins University.")
    ),
    card(
      card_header("Changelog"),
        tags$dl(
        tags$dt("2026-11-10 (v0.1.3)"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("Baseline dependency in the Gene Explorer"), " — ", "The cross-paper chart now ends with a Baseline dependency group, shown by default under Also show: DepMap CRISPR gene effect in MOLM-13, MOLM-14, MV4-11, OCI-AML2, THP-1, and OCI-LY-19, and Lin2019 MOLM-13 and THP-1 essentiality. Each is scored within its own cell line, so a gene that cells need to grow sits on the sensitizing side, where it can look like a venetoclax sensitizer. Grey diamonds mark these rows; they are not counted in the screen summaries.")
          )
        ),
        tags$dt("2026-11-03"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("Paper Explorer redesign"), " — ", "Each paper's screens are buttons in the sidebar, grouped by screen type, with small buttons for timepoints, doses, and replicates; the result shown in the Gene Explorer opens by default. The plot shows every gene by within-screen percentile on a choice of robust z, signed −log₁₀ FDR, or raw effect, labelling the top hits and BAX, PMAIP1, and TP53. A summary strip, top-hit tables with recurrence in other papers' venetoclax screens, and a collapsible full table replace the 11-column table."),
            tags$li(tags$strong("Expandable screens in the Gene Explorer"), " — ", "Screens with several results are marked ▸. Click one of a paper's dots, or pick it under Expand papers, to show each timepoint, replicate, dose, or derived table as an indented row; click again to collapse."),
            tags$li(tags$strong("BrimEdwards2026 and Diepstraten2026"), " — ", "Added the BrimEdwards2026 OCI-AML2 AZD4320 knockout screen (other BH3 mimetics) and the Diepstraten2026 Partita CRISPRa venetoclax hit lists. Updated data-availability notes for nine further candidates in the additional-screens catalog."),
            tags$li(tags$strong("Stewart2026 unnormalized arms"), " — ", "Stewart2026 venetoclax and venetoclax + palbociclib values are now median-centered, and their published FDRs no longer count as significant; this removes about 32,000 artifactual significant rows. The published values remain available under the Raw effect (as published) and Signed −log₁₀ FDR metrics."),
            tags$li(tags$strong("Wang2023, Sharon2019, Yokoyama2026, and Stewart2026"), " — ", "Added four studies from the additional-screens catalog: Wang2023 (MOLM-13 genome-wide and RNA-binding-protein libraries), Sharon2019 (acquired-resistant MOLM13-R1), Yokoyama2026 (peripheral T-cell lymphoma; venetoclax, navitoclax, ABBV-467), and Stewart2026 (OCI-AML2 venetoclax, palbociclib, and their combination). The Gene Explorer screen view now covers 17 venetoclax screens."),
            tags$li(tags$strong("Deng2022 CRISPRa venetoclax screens"), " — ", "Added the Deng2022 genome-wide CRISPR activation screens in mouse lymphoma lines #214 and #216 as gain-of-function hit lists from the published enriched-sgRNA tables (IC80 venetoclax vs DMSO). Mouse A1 paralogs (Bcl2a1a/b/c/d) now map to human BCL2A1, which also merges those Nakao2023 rows. The Deng2022 lymphoma vs B-ALL RNA-seq comparison was removed: it compared tumour types rather than venetoclax response, and its 6,895 FDR-significant genes inflated cross-study recurrence."),
            tags$li(tags$strong("Gene Explorer screen summary view"), " — ", "The Gene Explorer now opens on one row per independent screen, grouped by venetoclax context, with one panel per gene on a shared within-screen percentile axis and a per-gene summary of top/bottom-10% screens. Other BH3 mimetics and CRISPRa, ORF, and BH3-priming assays can be added with toggles; the previous all-results plot remains available under View."),
            tags$li(
              tags$strong("Effect sign, normalization, and rank fixes (v0.1.2)"),
              tags$ul(
                tags$li("MAGeCK negative-selection arms (Nechiporuk2019, Glytsou2023) now keep the signed gene LFC, and positive-selection arms no longer force it positive. Previously enriched resistance genes such as BAX and TP53 appeared as depleted in negative-selection arms."),
                tags$li("Eide2023 count-derived fold-changes are now normalized to the median library size and centered on the median gene, removing a −1.3 to −3.2 log₂ depth offset."),
                tags$li("Portelinha2024 ranks are now computed across genes; previously every gene had rank 1.")
              )
            )
          )
        ),
        tags$dt("2026-08-25"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("Nakao2023 MARCH5 CRISPR screens"), " — ", "Added three mouse mMA9Cas drug-vs-DMSO screen re-analyses from published sgRNA counts. DrugZ normalized Z and directional FDR drive browser rows; MAGeCK MLE results are retained for provenance and paper-figure concordance.")
          )
        ),
        tags$dt("2026-07-03"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("He2026 BH3-priming CRISPR screen"), " — ", "Added He2026 Figure 1 source-data combined MAGeCK gene-level LFC and p-values for the P1-R/P2-R BH3-priming screen, plus separate P1-R and P2-R sgRNA RPM-derived LFC summaries without p-values/FDRs.")
          )
        ),
        tags$dt("2026-05-25"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("BeatAML2 gene-set sensitivity context"), " — ", "Added curated gene-set scoring to the External Context expression-vs-AUC view, allowing pathway-like RNA signatures to be compared against BeatAML2 venetoclax response."),
            tags$li(tags$strong("Lin2019 paper-threshold hits"), " — ", "Annotated Lin2019 average drug-sensitivity hits using the paper's cutoff, Avg_sens < -0.914, while keeping p-value and FDR fields blank because they are not present in the source table."),
            tags$li(tags$strong("External venetoclax sensitivity context"), " — ", "Added Eide2023 CRISPR count summaries to the core gene browser, plus an External Context tab for BeatAML2/Eide2023 ex vivo venetoclax drug response and Zhang2024 CRISPR library coverage.")
          )
        ),
        tags$dt("2026-05-23"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("MAGeCK LFC fallback removed"), " — ", "MAGeCK experiments without a reported LFC column now omit that selection arm entirely, rather than substituting −log₁₀(p-value) as a proxy effect size. Affected: Nechiporuk2019 Y. Kosuke library day 14 negative selection arm, which the authors did not report in the supplementary table."),
            tags$li(tags$strong("Guièze2019 ORF screen documented"), " — ", "Added a data note explaining that ORF hits are defined per the paper's own threshold (any construct with LFC ≥ 2), that the library averaged ~1.3 constructs per gene, and that most genes lack multi-construct validation."),
            tags$li(
              tags$strong("Cross-paper axis normalization"), " — ",
              "Added a ", tags$em("Cross-paper axis"), " dropdown to the Gene Explorer sidebar with three options:",
              tags$ul(
                tags$li(tags$strong("Signed −log₁₀(FDR/p-value)"), " — default; places all studies on a comparable significance scale regardless of the underlying effect metric (STARS score, MAGeCK LFC, delta CRISPR score, etc.)."),
                tags$li(tags$strong("Per-experiment z-score"), " — normalizes each experiment to its own mean and SD across all genes."),
                tags$li(tags$strong("Raw effect value"), " — paper-native units; scales are not directly comparable across studies.")
              )
            ),
            tags$li(tags$strong("Layout improvements"), " — ", "Baseline Dependency Context and Source Detail sections converted to collapsed accordions so the main effects plot is visible on load. Page set to scroll rather than fit content to viewport.")
          )
        ),
        tags$dt("2026-04-27"),
        tags$dd(
          tags$ul(
            tags$li(tags$strong("Lin2019 baseline dependency context"), " — ", "Added Lin2019 essentiality three-scores as a baseline panel alongside DepMap CRISPR KO in the Gene Explorer."),
            tags$li(tags$strong("DepMap dependency plot"), " — ", "Added DepMap Public 26Q1 CRISPR gene-effect scores for MOLM-13, THP-1, and OCI-LY19 as a baseline context panel in the Gene Explorer."),
            tags$li(tags$strong("Initial release"), " — ", "Gene Explorer, Paper Explorer, Cross-study Overview, and About Data tabs covering Chen2019, Deng2022, Glytsou2023, Guièze2019, Lin2019, Nechiporuk2019, and Portelinha2024.")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  all_gene_choices <- gene_choices(scores)
  updateSelectizeInput(session, "gene_select", choices = all_gene_choices, selected = default_gene_keys, server = TRUE)
  updateSelectizeInput(session, "detail_gene", choices = all_gene_choices, selected = default_gene_keys[1], server = TRUE)
  updateSelectizeInput(session, "paper_highlight_genes", choices = all_gene_choices, selected = character(), server = TRUE)
  context_drug_choices <- sort(unique(stats::na.omit(drug_response$drug)))
  updateSelectizeInput(
    session,
    "context_drugs",
    choices = context_drug_choices,
    selected = intersect(c("Venetoclax", "VEN", "VEN + Azacytidine", "VEN + Ruxolitinib"), context_drug_choices),
    server = TRUE
  )
  zhang_choices <- if (nrow(zhang_library)) stats::setNames(zhang_library$gene_key, paste0(zhang_library$gene_symbol, " [", zhang_library$n_sgrna, " sgRNAs]")) else character()
  zhang_default <- if ("ZNF740" %in% zhang_choices) "ZNF740" else character()
  updateSelectizeInput(session, "zhang_gene", choices = zhang_choices, selected = zhang_default, server = TRUE)
  expression_choices <- if (nrow(beataml_expression$genes)) {
    stats::setNames(beataml_expression$genes$gene_key, paste0(beataml_expression$genes$gene_symbol, " [", beataml_expression$genes$gene_key, "]"))
  } else {
    character()
  }
  updateSelectizeInput(session, "context_expression_gene", choices = expression_choices, selected = intersect(c("BCL2", "MCL1", "BCL2L11"), expression_choices)[1], server = TRUE)
  gene_set_choices <- if (!is.null(beataml_gene_sets$coverage) && nrow(beataml_gene_sets$coverage)) {
    gene_set_labels <- beataml_gene_sets$coverage |>
      filter(.data$n_genes_scored > 0) |>
      arrange(.data$category, .data$gene_set) |>
      mutate(label = paste0(.data$gene_set, " [", .data$n_genes_scored, "/", .data$n_genes_defined, " genes]"))
    stats::setNames(gene_set_labels$gene_set, gene_set_labels$label)
  } else {
    character()
  }
  updateSelectizeInput(session, "context_gene_set", choices = gene_set_choices, selected = intersect(c("BCL2 family anti-apoptotic", "Venetoclax resistance candidates"), gene_set_choices)[1], server = TRUE)

  observeEvent(input$gene_select, {
    selected <- input$gene_select
    if (length(selected)) {
      selected_choices <- all_gene_choices[all_gene_choices %in% selected]
      updateSelectizeInput(session, "detail_gene", choices = selected_choices, selected = selected[1], server = TRUE)
    } else {
      updateSelectizeInput(session, "detail_gene", choices = all_gene_choices, selected = character(), server = TRUE)
    }
  })

  # The selected Paper Explorer result; the paper's first screen by default.
  paper_result_selected <- reactiveVal(NULL)
  observeEvent(input$paper_result_click, {
    if (input$paper_result_click %in% paper_results$result_id[paper_results$paper == input$paper_select]) {
      paper_result_selected(input$paper_result_click)
    }
  })

  output$paper_screen_buttons <- renderUI({
    req(input$paper_select)
    selected <- paper_result_selected()
    groups <- paper_screen_buttons(paper_results, input$paper_select)
    multi <- any(vapply(unlist(groups, recursive = FALSE), function(screen) nrow(screen$results) > 1, logical(1)))
    tags$div(
      class = "paper-screens",
      lapply(names(groups), function(group) tagList(
        tags$div(class = "grp", group),
        lapply(groups[[group]], function(screen) {
          ids <- screen$results$result_id
          tags$div(
            class = paste("screen-btn", if (any(ids %in% selected)) "on"), `data-result` = ids[1],
            screen$label,
            if (length(ids) > 1) tags$div(class = "chips", lapply(seq_along(ids), function(i) {
              tags$span(class = paste("chip", if (identical(ids[i], selected)) "on"), `data-result` = ids[i], screen$results$chip_label[i])
            }))
          )
        })
      )),
      if (multi) tags$div(class = "note", "Small buttons are other timepoints, doses, replicates, or tables of that screen; the first is the one shown in the Gene Explorer.")
    )
  })

  observeEvent(input$paper_select, {
    paper_result_selected(paper_default_result(paper_results, input$paper_select))
  }, ignoreInit = FALSE)

  selected_gene_keys <- reactive({
    selected <- input$gene_select
    if (is.null(selected)) selected <- default_gene_keys
    selected
  })

  gene_view <- reactive({
    req(length(selected_gene_keys()) > 0)
    df <- scores |>
      filter(.data$paper %in% input$gene_papers)
    df <- df |> filter(.data$gene_key %in% selected_gene_keys())
    if (isTRUE(input$gene_significant_only)) df <- df |> filter(.data$significance != "not annotated/significant")
    df
  })

  gene_effect_view <- reactive({
    gene_view() |>
      filter(!(.data$paper == "Lin2019" & grepl("essentiality", .data$experiment, ignore.case = TRUE)))
  })

  depmap_view <- reactive({
    if (!nrow(depmap_crispr) || !length(selected_gene_keys())) return(tibble())
    df <- depmap_crispr |>
      filter(.data$gene_key %in% selected_gene_keys(), .data$cell_line %in% input$depmap_cell_lines)
    df
  })

  lin2019_essentiality_view <- reactive({
    if (!length(selected_gene_keys())) return(tibble())
    scores |>
      filter(
        .data$paper == "Lin2019",
        grepl("essentiality", .data$experiment, ignore.case = TRUE),
        .data$gene_key %in% selected_gene_keys()
      )
  })

  gene_effect_display <- reactive({
    df <- collapse_selection_pairs(gene_effect_view())
    switch(
      input$effect_axis_type,
      signed_log10 = df |> mutate(
        display_value = dplyr::case_when(
          !is.na(.data$fdr) & .data$fdr > 0        ~ sign(.data$effect_value) * -log10(.data$fdr),
          !is.na(.data$fdr) & .data$fdr == 0       ~ sign(.data$effect_value) * 300,
          !is.na(.data$p_value) & .data$p_value > 0 ~ sign(.data$effect_value) * -log10(.data$p_value),
          !is.na(.data$p_value) & .data$p_value == 0 ~ sign(.data$effect_value) * 300,
          TRUE ~ NA_real_
        )
      ),
      zscore = df |>
        left_join(experiment_stats, by = "experiment_id") |>
        mutate(
          display_value = (.data$effect_value - .data$effect_mean) /
            dplyr::if_else(!is.na(.data$effect_sd) & .data$effect_sd > 0, .data$effect_sd, 1)
        ) |>
        select(-any_of(c("effect_mean", "effect_sd"))),
      raw = df |> mutate(display_value = .data$effect_value)
    )
  })

  # Every row of the selected result as published (both MAGeCK arms).
  paper_raw_view <- reactive({
    req(input$paper_select, paper_result_selected())
    scores |>
      filter(.data$paper == input$paper_select, .data$experiment %in% paper_results$experiment[paper_results$result_id == paper_result_selected()]) |>
      arrange(desc(abs(.data$effect_value)))
  })

  # One row per gene with within-screen percentile, robust z, and signed FDR.
  paper_view <- reactive({
    req(paper_result_selected() %in% paper_results$result_id)
    paper_result_scores(scores, paper_results, paper_result_selected())
  })

  overview_view <- reactive({
    df <- scores |> filter(.data$paper %in% input$overview_papers)
    if (!identical(input$overview_direction, "All")) df <- df |> filter(.data$direction == input$overview_direction)
    df
  })

  overview_hits <- reactive({
    overview_view() |>
      group_by(.data$gene_symbol, .data$gene_key) |>
      summarise(
        papers = n_distinct(.data$paper),
        experiments = n_distinct(.data$experiment_id),
        best_abs_effect = max(abs(.data$effect_value), na.rm = TRUE),
        best_fdr = suppressWarnings(min(.data$fdr, na.rm = TRUE)),
        directions = paste(sort(unique(.data$direction)), collapse = "; "),
        .groups = "drop"
      ) |>
      mutate(
        best_abs_effect = if_else(is.infinite(.data$best_abs_effect), NA_real_, .data$best_abs_effect),
        best_fdr = if_else(is.infinite(.data$best_fdr), NA_real_, .data$best_fdr)
      ) |>
      arrange(desc(.data$papers), desc(.data$experiments), desc(.data$best_abs_effect))
  })

  context_drug_view <- reactive({
    df <- drug_response
    if (length(input$context_datasets)) df <- df |> filter(.data$dataset %in% input$context_datasets)
    if (length(input$context_drugs)) df <- df |> filter(.data$drug %in% input$context_drugs)
    df |> arrange(.data$dataset, .data$drug, .data$response_value)
  })

  expression_auc_view <- reactive({
    mode <- input$context_expression_mode
    if (is.null(mode)) mode <- "gene"
    if (identical(mode, "gene_set")) {
      req(input$context_gene_set)
      score_matrix <- beataml_gene_sets$scores
      if (!length(score_matrix) || !input$context_gene_set %in% rownames(score_matrix)) return(tibble())
      expr <- tibble(
        rnaseq_sample = colnames(score_matrix),
        expression_value = as.numeric(score_matrix[input$context_gene_set, ])
      )
    } else {
      req(input$context_expression_gene)
      expression_matrix <- beataml_expression$expression
      if (!length(expression_matrix) || !input$context_expression_gene %in% rownames(expression_matrix)) return(tibble())
      expr <- tibble(
        rnaseq_sample = colnames(expression_matrix),
        expression_value = as.numeric(expression_matrix[input$context_expression_gene, ])
      )
    }
    drug_response |>
      filter(.data$dataset == "BeatAML2", .data$drug == "Venetoclax") |>
      inner_join(expr, by = "rnaseq_sample") |>
      filter(!is.na(.data$expression_value), !is.na(.data$response_value)) |>
      arrange(.data$expression_value)
  })

  expression_feature_label <- reactive({
    if (identical(input$context_expression_mode, "gene_set")) {
      return(as.character(input$context_gene_set))
    }
    if (!nrow(beataml_expression$genes) || is.null(input$context_expression_gene)) return(as.character(input$context_expression_gene))
    beataml_expression$genes |>
      filter(.data$gene_key == input$context_expression_gene) |>
      pull(.data$gene_symbol) |>
      dplyr::first(default = as.character(input$context_expression_gene))
  })

  expression_feature_axis_label <- reactive({
    if (identical(input$context_expression_mode, "gene_set")) {
      paste(expression_feature_label(), "BeatAML2 mean z-scored expression")
    } else {
      paste(expression_feature_label(), "BeatAML2 normalized expression")
    }
  })

  expression_feature_note <- reactive({
    if (!identical(input$context_expression_mode, "gene_set")) return(NULL)
    coverage <- beataml_gene_sets$coverage
    if (is.null(coverage) || !nrow(coverage) || is.null(input$context_gene_set)) return(NULL)
    row <- coverage |> filter(.data$gene_set == input$context_gene_set) |> slice(1)
    if (!nrow(row)) return(NULL)
    paste0(
      " Gene-set score is the mean z-scored RNA expression of ",
      row$n_genes_scored, " of ", row$n_genes_defined, " curated genes",
      if (nzchar(row$genes_missing)) paste0("; missing from BeatAML2 expression: ", row$genes_missing) else "",
      "."
    )
  })

  context_crispr_view <- reactive({
    req(input$context_crispr_experiment)
    eide_crispr_context |>
      filter(.data$experiment == input$context_crispr_experiment) |>
      arrange(desc(abs(.data$effect_value)))
  })

  zhang_library_view <- reactive({
    if (!nrow(zhang_library)) return(tibble())
    if (!is.null(input$zhang_gene) && nzchar(input$zhang_gene)) {
      zhang_library |> filter(.data$gene_key == input$zhang_gene)
    } else {
      zhang_library
    }
  })

  gene_screen_view <- reactive({
    req(length(selected_gene_keys()) > 0)
    groups <- c(venetoclax_screen_groups, input$screen_extra_groups)
    baseline <- if (baseline_screen_group %in% groups) {
      baseline_scores |>
        filter(
          .data$gene_key %in% selected_gene_keys(),
          (.data$paper == "DepMap" & .data$cell_line %in% input$depmap_cell_lines) |
            (.data$paper != "DepMap" & .data$paper %in% input$gene_papers)
        )
    }
    bind_rows(
      screen_scores,
      screen_detail_scores |> filter(.data$paper %in% input$screen_expand_papers)
    ) |>
      filter(
        .data$gene_key %in% selected_gene_keys(),
        .data$paper %in% input$gene_papers,
        .data$screen_group %in% groups
      ) |>
      bind_rows(baseline)
  })

  # Clicking a dot expands or collapses that paper's individual results.
  observeEvent(event_data("plotly_click", source = "screen_panels"), {
    paper <- event_data("plotly_click", source = "screen_panels")$customdata
    if (is.null(paper) || !length(paper) || !paper[[1]] %in% expandable_papers) return()
    paper <- paper[[1]]
    current <- input$screen_expand_papers
    expanded <- if (paper %in% current) setdiff(current, paper) else c(current, paper)
    updateSelectizeInput(session, "screen_expand_papers", selected = expanded)
  })

  screen_metric <- reactive({
    metric <- input$screen_metric
    if (is.null(metric) || !metric %in% names(screen_metrics)) "percentile" else metric
  })

  output$gene_screen_band_note <- renderUI({
    tags$span(screen_metrics[[screen_metric()]]$band_note)
  })

  output$gene_screen_hidden_note <- renderUI({
    spec <- screen_metrics[[screen_metric()]]
    hidden <- gene_screen_view() |>
      filter(.data$row_level == 0, is.na(.data[[spec$column]]), !(spec$hit_list_column & .data$hit_list_only)) |>
      distinct(.data$short_label) |>
      pull(.data$short_label)
    if (!length(hidden)) return(NULL)
    tags$p(
      class = "small-muted",
      tags$strong(paste0("Not shown on this metric (no ", if (identical(screen_metric(), "signed_fdr")) "FDR" else "value", "): ")),
      paste(hidden, collapse = "; ")
    )
  })

  output$gene_screen_summary <- renderUI({
    req(length(selected_gene_keys()) > 0)
    support <- summarise_screen_support(screen_scores, selected_gene_keys(), papers = input$gene_papers)
    symbols <- gene_screen_view() |>
      group_by(.data$gene_key) |>
      summarise(gene_display = format_gene_symbols(.data$gene_symbol), .groups = "drop")
    support <- support |>
      left_join(symbols, by = "gene_key") |>
      mutate(gene_display = dplyr::coalesce(.data$gene_display, .data$gene_symbol))
    chips <- lapply(seq_len(nrow(support)), function(i) {
      row <- support[i, ]
      detail <- if (row$n_tested == 0) {
        "Not measured in the selected venetoclax screens."
      } else if (row$n_resistance == 0 && row$n_sensitizing == 0) {
        paste0("No top- or bottom-10% signal in ", row$n_tested, " venetoclax screens.")
      } else {
        parts <- c(
          if (row$n_resistance > 0) paste0("Resistance (top 10%) in ", row$n_resistance, " of ", row$n_tested, " venetoclax screens; FDR ≤ 0.05 in ", row$n_resistance_fdr, "."),
          if (row$n_sensitizing > 0) paste0("Sensitizing (bottom 10%) in ", row$n_sensitizing, " of ", row$n_tested, " venetoclax screens; FDR ≤ 0.05 in ", row$n_sensitizing_fdr, ".")
        )
        paste(parts, collapse = " ")
      }
      tags$div(class = "gene-summary", tags$div(class = "gene", row$gene_display), tags$div(class = "detail", detail))
    })
    tags$div(class = "gene-summary-grid", chips)
  })

  output$gene_screen_plot_ui <- renderUI({
    df <- gene_screen_view()
    if (!nrow(df)) {
      return(tags$p(class = "baseline-empty", "The selected genes were not measured in the selected screens."))
    }
    plotlyOutput("gene_screen_plot", height = paste0(screen_panel_height(df, expandable_papers, input$screen_expand_papers, screen_detail_parents), "px"))
  })

  output$gene_screen_plot <- renderPlotly({
    df <- gene_screen_view()
    validate(need(nrow(df) > 0, "Search for one or more genes to populate this view."))
    plot_screen_panels(
      df,
      gene_order = selected_gene_keys(),
      metric = screen_metric(),
      expandable_papers = expandable_papers,
      expanded_papers = input$screen_expand_papers,
      screen_detail_parents = screen_detail_parents
    )
  })

  output$gene_effect_plot <- renderPlotly({
    df <- gene_effect_display()
    validate(
      need(nrow(df) > 0, "Search for one or more genes to populate this view."),
      need(nrow(df) <= max_gene_plot_rows, paste("Too many rows to plot interactively. Select fewer genes or papers; current rows:", nrow(df)))
    )
    x_labels <- c(
      signed_log10 = "Signed −log₁₀(FDR; †p-value where FDR unavailable; NA if neither present)",
      zscore       = "Per-experiment z-score (relative to all genes in that experiment)",
      raw          = "Effect value (paper-native units — not directly comparable across papers)"
    )
    plot_gene_effects(df, x_label = unname(x_labels[input$effect_axis_type]))
  })

  output$gene_effect_na_note <- renderUI({
    df <- gene_effect_display()
    hidden <- df |> filter(is.na(.data$display_value), !is.na(.data$effect_value))
    if (!nrow(hidden)) return(NULL)
    sources <- hidden |>
      group_by(.data$paper, .data$experiment) |>
      summarise(n = n(), .groups = "drop") |>
      mutate(label = paste0(.data$paper, " — ", .data$experiment, " (", .data$n, " gene", if_else(.data$n == 1, "", "s"), ")")) |>
      pull(label)
    tags$div(
      style = "margin: 6px 12px 10px; padding: 8px 12px; background: #fff8e1; border-left: 3px solid #f9a825; border-radius: 4px; font-size: 0.88rem; color: #5d4200;",
      tags$strong(paste0(nrow(hidden), " row", if (nrow(hidden) > 1) "s" else "", " not shown on this axis")),
      " — no FDR or p-value available (threshold-defined hits). Switch to ",
      tags$strong("Raw effect value"), " to include them.",
      tags$ul(style = "margin: 4px 0 0 0; padding-left: 18px;", lapply(sources, tags$li))
    )
  })

  output$gene_table <- renderDT({
    gene_view() |>
      select(paper, experiment, gene_symbol, species, effect_metric, effect_value, p_value, fdr, rank, direction, source_level) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv"))) |>
      formatSignif(columns = c("effect_value", "p_value"), digits = 3)
  }, server = TRUE)

  output$depmap_crispr_plot_ui <- renderUI({
    if (!nrow(depmap_view())) {
      return(tags$p(class = "baseline-empty", "No DepMap CRISPR KO rows are available for the selected genes and cell lines."))
    }
    n_genes <- depmap_view() |>
      pull(.data$gene_key) |>
      unique() |>
      length()
    plot_height <- max(440, min(900, 260 + 76 * n_genes))
    plotlyOutput("depmap_crispr_plot", height = paste0(plot_height, "px"))
  })

  output$lin2019_essentiality_plot_ui <- renderUI({
    if (!nrow(lin2019_essentiality_view())) {
      return(tags$p(class = "baseline-empty", "No Lin2019 essentiality rows are available for the selected genes."))
    }
    n_genes <- lin2019_essentiality_view() |>
      pull(.data$gene_key) |>
      unique() |>
      length()
    plot_height <- max(440, min(900, 260 + 76 * n_genes))
    plotlyOutput("lin2019_essentiality_plot", height = paste0(plot_height, "px"))
  })

  output$depmap_crispr_plot <- renderPlotly({
    validate(
      need(nrow(depmap_crispr) > 0, "DepMap CRISPR data have not been processed yet."),
      need(nrow(depmap_view()) > 0, "Select one or more genes and matched DepMap cell lines.")
    )
    plot_depmap_crispr(depmap_view())
  })

  output$lin2019_essentiality_plot <- renderPlotly({
    validate(
      need(nrow(lin2019_essentiality_view()) > 0, "No Lin2019 essentiality rows are available for the selected genes.")
    )
    plot_lin2019_essentiality(lin2019_essentiality_view())
  })

  output$detail_table <- renderDT({
    req(input$detail_gene, input$detail_paper)
    detail <- load_detail_for_gene(input$detail_gene, input$detail_paper)
    datatable(detail, filter = "top", options = list(pageLength = 10, scrollX = TRUE))
  }, server = TRUE)

  output$paper_rationale <- renderUI({
    req(input$paper_select)
    row <- metadata |> filter(.data$paper == input$paper_select) |> slice(1)
    tags$div(
      class = "rationale-box",
      tags$strong(paste0(row$paper, " (", row$year, ")")),
      tags$p(row$rationale),
      tags$p(class = "small-muted", paste(row$model_system, row$cell_line_or_context, row$species, sep = " | "))
    )
  })

  # Strongest hits each way, with how often each recurs in the same direction
  # in independent venetoclax screens from other papers.
  paper_top_hits <- reactive({
    df <- paper_view()
    req(nrow(df))
    spec <- paper_metric_spec(df, input$paper_metric)
    value <- df[[spec$column]]
    hits <- select_paper_hits(df, value, max(10, input$paper_top_n))
    other_papers <- setdiff(unique(screen_view_registry$paper), input$paper_select)
    support <- summarise_screen_support(screen_scores, unique(unlist(hits)), papers = other_papers)
    out <- df |>
      mutate(value = value) |>
      filter(.data$gene_key %in% unlist(hits)) |>
      left_join(support |> select(-"gene_symbol"), by = "gene_key") |>
      mutate(
        side = if_else(.data$gene_key %in% hits$resistance, "resistance", "sensitizing"),
        same = if_else(.data$side == "resistance", .data$n_resistance, .data$n_sensitizing),
        same_fdr = if_else(.data$side == "resistance", .data$n_resistance_fdr, .data$n_sensitizing_fdr),
        recurrence = if_else(.data$n_tested == 0, "not measured", paste0(.data$same, " of ", .data$n_tested, if_else(.data$same_fdr > 0, paste0(" · ", .data$same_fdr, " at FDR ≤ 0.05"), "")))
      )
    list(data = out, spec = spec)
  })

  paper_hit_table <- function(which_side) {
    hits <- paper_top_hits()
    df <- hits$data[hits$data$side == which_side, , drop = FALSE]
    df <- df[order(if (which_side == "resistance") -df$value else df$value), , drop = FALSE]
    hit_list <- all(is.na(df$percentile))
    tbl <- tibble::tibble(
      Gene = df$gene_symbol,
      !!hits$spec$label := signif(df$value, 3),
      Percentile = if (hit_list) NULL else sprintf("%+.3f", df$percentile),
      FDR = case_when(is.na(df$fdr) ~ "—", df$fdr <= 0 ~ "0", TRUE ~ format(signif(df$fdr, 2), scientific = TRUE)),
      `Other venetoclax screens` = df$recurrence
    )
    datatable(tbl, rownames = FALSE, options = list(dom = "t", pageLength = nrow(tbl), ordering = FALSE, columnDefs = list(list(className = "dt-left", targets = "_all"))))
  }
  output$paper_hits_resistance <- renderDT(paper_hit_table("resistance"))
  output$paper_hits_sensitizing <- renderDT(paper_hit_table("sensitizing"))

  output$paper_summary <- renderUI({
    df <- paper_view()
    req(nrow(df))
    hit_list <- all(is.na(df$percentile))
    has_stats <- any(!is.na(df$fdr)) || any(grepl("^paper threshold", df$significance))
    count_tile <- function(label, n) {
      div(class = "metric", div(class = "small-muted", label), div(class = "value", if (has_stats) format(n, big.mark = ",") else "—"),
          if (!has_stats) div(class = "detail", "No FDR or threshold in the browser for this result"))
    }
    refs <- df |> filter(.data$gene_key %in% paper_reference_genes)
    ref_text <- if (hit_list) {
      listed <- intersect(paper_reference_genes, refs$gene_key)
      if (length(listed)) paste("In hit list:", paste(listed, collapse = ", ")) else "None in the hit list"
    } else if (nrow(refs)) {
      paste(paste0(refs$gene_symbol, " ", sub("-", "−", sprintf("%+.2f", refs$percentile))), collapse = " · ")
    } else "Not measured"
    div(
      class = "metric-strip",
      div(class = "metric", div(class = "small-muted", if (hit_list) "Genes in published table" else "Genes scored"), div(class = "value", format(nrow(df), big.mark = ","))),
      count_tile("Resistance hits (FDR ≤ 0.05 or paper cutoff)", sum(df$supported & df$effect_value > 0)),
      count_tile("Sensitizing hits (FDR ≤ 0.05 or paper cutoff)", sum(df$supported & df$effect_value < 0)),
      div(class = "metric", div(class = "small-muted", if (hit_list) "Reference genes" else "Reference genes (percentile)"), div(class = "detail", style = "font-size:1rem; font-weight:600; color:#2a3a44", ref_text))
    )
  })

  output$paper_rank_plot <- renderPlotly({
    plot_paper_screen(
      paper_view(),
      metric = input$paper_metric,
      top_n = input$paper_top_n,
      highlight_gene_keys = input$paper_highlight_genes
    )
  })

  output$paper_table <- renderDT({
    paper_raw_view() |>
      slice_head(n = max_paper_table_rows) |>
      select(paper, experiment, gene_symbol, species, effect_metric, effect_value, p_value, fdr, rank, direction, source_level) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")))
  }, server = TRUE)

  output$overview_heatmap <- renderPlotly({
    plot_cross_study_heatmap(overview_view(), top_n = input$overview_top_n)
  })

  output$overview_table <- renderDT({
    datatable(overview_hits(), filter = "top", extensions = "Buttons", options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")))
  }, server = TRUE)

  output$context_drug_plot <- renderPlotly({
    validate(need(nrow(context_drug_view()) > 0, "Select one or more datasets and venetoclax-containing treatments."))
    plot_drug_response_distribution(context_drug_view())
  })

  output$expression_auc_summary <- renderUI({
    df <- expression_auc_view()
    validate(need(nrow(df) >= 3, "Select a gene or gene set with matched BeatAML2 expression and venetoclax AUC data."))
    pearson <- suppressWarnings(stats::cor.test(df$expression_value, df$response_value, method = "pearson"))
    spearman <- suppressWarnings(stats::cor.test(df$expression_value, df$response_value, method = "spearman", exact = FALSE))
    linear_fit <- stats::lm(response_value ~ expression_value, data = df)
    r_squared <- summary(linear_fit)$r.squared
    tags$p(
      class = "baseline-note",
      tags$strong(expression_feature_label()),
      " across ", nrow(df), " BeatAML2 samples with matched RNA-seq and venetoclax response. ",
      "Pearson r = ", signif(unname(pearson$estimate), 3), " (p = ", signif(pearson$p.value, 3), "); ",
      "Spearman rho = ", signif(unname(spearman$estimate), 3), " (p = ", signif(spearman$p.value, 3), "); ",
      "linear R² = ", signif(r_squared, 3), ". ",
      expression_feature_note(),
      "Lower AUC indicates greater sensitivity."
    )
  })

  output$expression_auc_plot <- renderPlotly({
    validate(need(nrow(expression_auc_view()) >= 3, "Select a gene or gene set with matched BeatAML2 expression and venetoclax AUC data."))
    plot_expression_auc(
      expression_auc_view(),
      gene_label = expression_feature_label(),
      x_axis_label = expression_feature_axis_label()
    )
  })

  output$expression_auc_table <- renderDT({
    expression_auc_view() |>
      select(sample_id, rnaseq_sample, drug, treatment_context, expression_value, response_value, response_direction) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv"))) |>
      formatSignif(columns = c("expression_value", "response_value"), digits = 3)
  }, server = TRUE)

  output$context_drug_table <- renderDT({
    context_drug_view() |>
      select(dataset, sample_id, drug, treatment_context, response_metric, response_value, response_direction, disease_context, any_of(c("disease_stage", "specimen_type", "eln_2022", "variant_summary"))) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv"))) |>
      formatSignif(columns = "response_value", digits = 3)
  }, server = TRUE)

  output$context_crispr_plot <- renderPlotly({
    validate(need(nrow(context_crispr_view()) > 0, "No CRISPR rows are available for this experiment."))
    plot_external_crispr_rank(context_crispr_view(), top_n = input$context_crispr_top_n)
  })

  output$context_crispr_table <- renderDT({
    context_crispr_view() |>
      select(dataset, experiment, gene_symbol, effect_metric, effect_value, rank, direction, source_file) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv"))) |>
      formatSignif(columns = "effect_value", digits = 3)
  }, server = TRUE)

  output$zhang_library_table <- renderDT({
    zhang_library_view() |>
      select(dataset, gene_symbol, n_sgrna, example_sgrnas, source_file) |>
      datatable(filter = "top", extensions = "Buttons", options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")))
  }, server = TRUE)

  output$metric_rows <- renderText(format(nrow(overview_view()), big.mark = ","))
  output$metric_genes <- renderText(format(n_distinct(overview_view()$gene_key), big.mark = ","))
  output$metric_papers <- renderText(n_distinct(overview_view()$paper))
  output$metric_experiments <- renderText(n_distinct(overview_view()$experiment_id))

  output$metadata_table <- renderDT({
    datatable(metadata, options = list(pageLength = 7, scrollX = TRUE), rownames = FALSE)
  }, server = TRUE)

  output$additional_screen_table <- renderDT({
    additional_screen_candidates |>
      select(
        Study = study, Priority = priority, Role = inclusion_role,
        `Screen type` = perturbation, Model = model, Treatment = treatment,
        `Data available` = processed_data, `Next step` = browser_action,
        `Evidence checked` = evidence_status, `EndNote records` = endnote_record_ids,
        Source = source_url
      ) |>
      datatable(
        filter = "top", rownames = FALSE,
        options = list(pageLength = 10, scrollX = TRUE)
      )
  }, server = TRUE)

  output$download_additional_screens <- downloadHandler(
    filename = function() "additional_screen_candidates.csv",
    content = function(file) readr::write_csv(additional_screen_candidates, file)
  )

  output$download_gene_view <- downloadHandler(
    filename = function() "gene_explorer_view.csv",
    content = function(file) readr::write_csv(gene_view(), file)
  )
  output$download_paper_view <- downloadHandler(
    filename = function() paste0(input$paper_select, "_", make.names(paper_result_selected()), ".csv"),
    content = function(file) readr::write_csv(paper_raw_view(), file)
  )
  output$download_overview <- downloadHandler(
    filename = function() "cross_study_recurrent_hits.csv",
    content = function(file) readr::write_csv(overview_hits(), file)
  )
  output$download_context_drug_response <- downloadHandler(
    filename = function() "external_context_drug_response.csv",
    content = function(file) readr::write_csv(context_drug_view(), file)
  )
}

shinyApp(ui, server)
