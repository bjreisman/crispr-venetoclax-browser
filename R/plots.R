library(dplyr)
library(ggplot2)
library(plotly)
library(forcats)

effect_color_values <- c(
  "resistance/enriched" = "#b73e3e",
  "sensitization/depleted" = "#247a6b",
  "ORF overexpression resistance/enriched" = "#b73e3e",
  "ORF overexpression depleted" = "#247a6b",
  "higher score" = "#7a5c1f",
  "lower score" = "#326b9a",
  "neutral" = "#767676",
  "unknown" = "#9a9a9a"
)

format_gene_symbols <- function(x) {
  symbols <- sort(unique(stats::na.omit(as.character(x))))
  symbols <- symbols[order(!grepl("^[A-Z0-9.-]+$", symbols), toupper(symbols), symbols)]
  paste(symbols, collapse = " / ")
}

plot_gene_effects <- function(df, x_label = "Effect value") {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  plot_df <- df |>
    group_by(.data$gene_key) |>
    mutate(gene_display = format_gene_symbols(.data$gene_symbol)) |>
    ungroup() |>
    mutate(
      significant = .data$significance != "not annotated/significant",
      label = paste0(
        .data$paper,
        "<br>", .data$experiment,
        "<br>", .data$gene_display,
        "<br>Source symbol: ", .data$gene_symbol,
        "<br>Displayed: ", dplyr::if_else(is.na(.data$display_value), "NA", as.character(signif(.data$display_value, 3))),
        "<br>Raw (", .data$effect_metric, "): ", signif(.data$effect_value, 3),
        "<br>", .data$significance
      ),
      y_label = forcats::fct_reorder(paste(.data$paper, .data$experiment, sep = " | "), .data$year)
    )
  x_range <- range(plot_df$display_value, na.rm = TRUE)
  x_span <- diff(x_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  x_left <- x_range[1] - 0.03 * x_span
  x_right <- x_range[2] + 0.03 * x_span

  p <- ggplot(plot_df, aes(
    x = display_value,
    y = y_label,
    color = gene_display,
    text = label
  )) +
    geom_vline(xintercept = 0, color = "#b9b9b9", linewidth = 0.35) +
    geom_point(data = filter(plot_df, .data$significant), shape = 16, size = 5, alpha = 0.9) +
    geom_point(data = filter(plot_df, !.data$significant), shape = 1, size = 5, stroke = 1.2, alpha = 0.9) +
    annotate("text", x = x_left, y = Inf, label = "Sensitization / depletion", hjust = 0, vjust = 1.4, size = 3.6, color = "#247a6b") +
    annotate("text", x = x_right, y = Inf, label = "Resistance / enrichment", hjust = 1, vjust = 1.4, size = 3.6, color = "#b73e3e") +
    coord_cartesian(xlim = c(x_left, x_right), clip = "off") +
    labs(x = x_label, y = NULL, color = "Gene") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 10),
      plot.margin = margin(26, 16, 20, 8)
    )

  ggplotly(p, tooltip = "text") |>
    layout(legend = list(orientation = "h", y = -0.25))
}

plot_ranked_experiment <- function(df, enriched_n = 8, depleted_n = 8, highlight_gene_keys = character()) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  if (is.null(enriched_n) || is.na(enriched_n)) enriched_n <- 8
  if (is.null(depleted_n) || is.na(depleted_n)) depleted_n <- 8
  enriched_n <- max(0, enriched_n)
  depleted_n <- max(0, depleted_n)
  if (is.null(highlight_gene_keys)) highlight_gene_keys <- character()

  enriched <- df |>
    filter(!is.na(.data$effect_value), .data$effect_value > 0) |>
    arrange(desc(.data$effect_value), .data$rank) |>
    slice_head(n = enriched_n)

  depleted <- df |>
    filter(!is.na(.data$effect_value), .data$effect_value < 0) |>
    arrange(.data$effect_value, .data$rank) |>
    slice_head(n = depleted_n)

  highlighted <- df |>
    filter(.data$gene_key %in% highlight_gene_keys, !is.na(.data$effect_value))

  plot_df <- bind_rows(depleted, enriched, highlighted) |>
    arrange(.data$rank, desc(abs(.data$effect_value))) |>
    distinct(.data$gene_key, .keep_all = TRUE) |>
    arrange(.data$effect_value) |>
    mutate(
      highlighted = .data$gene_key %in% highlight_gene_keys,
      gene_label_text = if_else(.data$highlighted, paste0("<b>", .data$gene_symbol, "</b>"), .data$gene_symbol),
      gene_label = forcats::fct_reorder(.data$gene_label_text, .data$effect_value),
      label = paste0(
        .data$gene_symbol,
        if_else(.data$highlighted, " (highlighted)", ""),
        "<br>", .data$direction,
        "<br>effect: ", signif(.data$effect_value, 3),
        "<br>FDR: ", signif(.data$fdr, 3)
      )
    )

  if (!nrow(plot_df)) return(plotly_empty(type = "scatter", mode = "markers"))

  p <- ggplot(plot_df, aes(x = effect_value, y = gene_label, color = direction, text = label)) +
    geom_vline(xintercept = 0, color = "#b9b9b9", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = effect_value, yend = gene_label), linewidth = 0.5, alpha = 0.65) +
    geom_point(aes(size = highlighted), alpha = 0.85) +
    scale_size_manual(values = c(`FALSE` = 2.4, `TRUE` = 4.2), guide = "none") +
    scale_color_manual(values = effect_color_values, na.value = "#9a9a9a") +
    labs(x = "Effect value", y = NULL, color = "Direction") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggplotly(p, tooltip = "text") |>
    layout(legend = list(orientation = "h", y = -0.18))
}

plot_cross_study_heatmap <- function(df, top_n = 40) {
  if (!nrow(df)) return(plotly_empty(type = "heatmap"))

  top_genes <- df |>
    filter(.data$significance != "not annotated/significant" | abs(.data$effect_value) > 2) |>
    count(.data$gene_symbol, .data$gene_key, sort = TRUE) |>
    slice_head(n = top_n)

  if (!nrow(top_genes)) {
    top_genes <- df |> count(.data$gene_symbol, .data$gene_key, sort = TRUE) |> slice_head(n = top_n)
  }

  matrix_df <- df |>
    semi_join(top_genes, by = c("gene_symbol", "gene_key")) |>
    group_by(.data$gene_symbol, .data$experiment_id) |>
    arrange(desc(abs(.data$effect_value)), .by_group = TRUE) |>
    slice(1) |>
    ungroup() |>
    mutate(score = if_else(is.infinite(.data$effect_value), NA_real_, .data$effect_value))

  max_abs_score <- max(abs(matrix_df$score), na.rm = TRUE)
  if (!is.finite(max_abs_score) || max_abs_score == 0) max_abs_score <- 1

  p <- ggplot(matrix_df, aes(
    x = gene_symbol,
    y = forcats::fct_rev(experiment_id),
    fill = score,
    text = paste0(
      gene_symbol,
      "<br>", experiment_id,
      "<br>", direction,
      "<br>Effect: ", signif(score, 3)
    )
  )) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#247a6b",
      mid = "#f4f4f1",
      high = "#b73e3e",
      midpoint = 0,
      limits = c(-max_abs_score, max_abs_score),
      na.value = "#f3f3f3"
    ) +
    labs(x = NULL, y = NULL, fill = "Effect") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank()
    )

  ggplotly(p, tooltip = "text")
}

plot_depmap_crispr <- function(df) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  plot_df <- df |>
    filter(!is.na(.data$gene_effect)) |>
    group_by(.data$gene_key) |>
    mutate(gene_display = format_gene_symbols(.data$gene_symbol)) |>
    ungroup() |>
    mutate(
      label = paste0(
        .data$gene_display,
        "<br>", .data$cell_line,
        "<br>Gene effect: ", signif(.data$gene_effect, 3),
        "<br>", .data$dependency_strength,
        "<br>", .data$depmap_release
      ),
      gene_display = forcats::fct_reorder(.data$gene_display, .data$gene_effect, .fun = median, na.rm = TRUE)
    )

  if (!nrow(plot_df)) return(plotly_empty(type = "scatter", mode = "markers"))

  y_levels <- levels(plot_df$gene_display)
  x_range <- range(plot_df$gene_effect, na.rm = TRUE)
  x_span <- diff(x_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1

  plot_ly(
    plot_df,
    x = ~gene_effect,
    y = ~gene_display,
    color = ~cell_line,
    colors = c("MOLM-13" = "#f8766d", "OCI-LY-19" = "#2fbf4a", "THP-1" = "#619cff"),
    type = "scatter",
    mode = "markers",
    marker = list(size = 10, opacity = 0.9),
    text = ~label,
    hoverinfo = "text"
  ) |>
    layout(
      xaxis = list(
        title = "DepMap CRISPR gene effect",
        zeroline = FALSE,
        range = c(x_range[1] - 0.08 * x_span, x_range[2] + 0.08 * x_span)
      ),
      yaxis = list(title = "", categoryorder = "array", categoryarray = y_levels, automargin = TRUE),
      legend = list(orientation = "h", x = 0, y = -0.22),
      margin = list(l = 120, r = 24, t = 18, b = 92),
      shapes = list(
        list(type = "line", x0 = 0, x1 = 0, y0 = 0, y1 = 1, xref = "x", yref = "paper", line = list(color = "#b9b9b9", width = 1)),
        list(type = "line", x0 = -0.5, x1 = -0.5, y0 = 0, y1 = 1, xref = "x", yref = "paper", line = list(color = "#555555", width = 1, dash = "dash"))
      ),
      annotations = list(
        list(
          x = -0.5, y = 1.02, xref = "x", yref = "paper",
          text = "dependency threshold", showarrow = FALSE,
          xanchor = "right", font = list(size = 11, color = "#555555")
        )
      )
    )
}

plot_lin2019_essentiality <- function(df) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  plot_df <- df |>
    filter(!is.na(.data$effect_value)) |>
    group_by(.data$gene_key) |>
    mutate(gene_display = format_gene_symbols(.data$gene_symbol)) |>
    ungroup() |>
    mutate(
      context = sub(" essentiality$", "", .data$experiment),
      label = paste0(
        .data$gene_display,
        "<br>", .data$context,
        "<br>Essentiality score: ", signif(.data$effect_value, 3),
        "<br>", .data$effect_metric,
        "<br>Source: Lin2019"
      ),
      gene_display = forcats::fct_reorder(.data$gene_display, .data$effect_value, .fun = median, na.rm = TRUE)
    )

  if (!nrow(plot_df)) return(plotly_empty(type = "scatter", mode = "markers"))

  y_levels <- levels(plot_df$gene_display)
  x_range <- range(plot_df$effect_value, na.rm = TRUE)
  x_span <- diff(x_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1

  plot_ly(
    plot_df,
    x = ~effect_value,
    y = ~gene_display,
    color = ~context,
    colors = c("THP-1" = "#619cff", "MOLM-13" = "#f8766d", "Average" = "#7a7a7a"),
    type = "scatter",
    mode = "markers",
    marker = list(size = 10, opacity = 0.9),
    text = ~label,
    hoverinfo = "text"
  ) |>
    layout(
      xaxis = list(
        title = "Lin2019 essentiality three-score",
        zeroline = FALSE,
        range = c(max(0, x_range[1] - 0.08 * x_span), x_range[2] + 0.08 * x_span)
      ),
      yaxis = list(title = "", categoryorder = "array", categoryarray = y_levels, automargin = TRUE),
      legend = list(orientation = "h", x = 0, y = -0.22),
      margin = list(l = 120, r = 24, t = 18, b = 92),
      shapes = list(
        list(type = "line", x0 = 1, x1 = 1, y0 = 0, y1 = 1, xref = "x", yref = "paper", line = list(color = "#b9b9b9", width = 1, dash = "dot"))
      )
    )
}

plot_drug_response_distribution <- function(df) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  plot_df <- df |>
    filter(!is.na(.data$response_value)) |>
    mutate(
      drug = forcats::fct_reorder(.data$drug, .data$response_value, .fun = median, na.rm = TRUE),
      label = paste0(
        .data$dataset,
        "<br>", .data$drug,
        "<br>Sample: ", .data$sample_id,
        "<br>AUC: ", signif(.data$response_value, 3),
        "<br>", .data$response_direction
      )
    )

  p <- ggplot(plot_df, aes(x = response_value, y = drug, color = dataset, text = label)) +
    geom_boxplot(aes(group = interaction(dataset, drug)), outlier.shape = NA, alpha = 0.25, linewidth = 0.35) +
    geom_jitter(height = 0.18, width = 0, alpha = 0.35, size = 1.6) +
    labs(x = "Drug response AUC (lower = more sensitive)", y = NULL, color = "Dataset") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggplotly(p, tooltip = "text") |>
    layout(legend = list(orientation = "h", y = -0.18))
}

plot_expression_auc <- function(df, gene_label = "Gene", x_axis_label = NULL) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))
  if (is.null(x_axis_label)) x_axis_label <- paste(gene_label, "BeatAML2 normalized expression")

  plot_df <- df |>
    filter(!is.na(.data$expression_value), !is.na(.data$response_value)) |>
    mutate(
      label = paste0(
        "Sample: ", .data$sample_id,
        "<br>RNA-seq sample: ", .data$rnaseq_sample,
        "<br>", gene_label, ": ", signif(.data$expression_value, 3),
        "<br>AUC: ", signif(.data$response_value, 3),
        "<br>Treatment: ", .data$drug
      )
    )

  if (!nrow(plot_df)) return(plotly_empty(type = "scatter", mode = "markers"))

  fit <- stats::lm(response_value ~ expression_value, data = plot_df)
  r_squared <- summary(fit)$r.squared
  x_range <- range(plot_df$expression_value, na.rm = TRUE)
  y_range <- range(plot_df$response_value, na.rm = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span == 0) x_span <- 1
  if (!is.finite(y_span) || y_span == 0) y_span <- 1

  line_df <- tibble(
    expression_value = seq(x_range[1], x_range[2], length.out = 100)
  )
  line_df$response_value <- as.numeric(stats::predict(fit, newdata = line_df))

  plot_ly() |>
    add_markers(
      data = plot_df,
      x = ~expression_value,
      y = ~response_value,
      text = ~label,
      hoverinfo = "text",
      marker = list(color = "#33658a", size = 8, opacity = 0.75),
      name = "Samples"
    ) |>
    add_lines(
      data = line_df,
      x = ~expression_value,
      y = ~response_value,
      line = list(color = "#b73e3e", width = 3),
      hoverinfo = "skip",
      name = "Linear fit"
    ) |>
    layout(
      xaxis = list(
        title = x_axis_label,
        range = c(x_range[1] - 0.05 * x_span, x_range[2] + 0.05 * x_span)
      ),
      yaxis = list(
        title = "Venetoclax AUC (lower = more sensitive)",
        range = c(y_range[1] - 0.08 * y_span, y_range[2] + 0.12 * y_span)
      ),
      annotations = list(
        list(
          x = x_range[1],
          y = y_range[2] + 0.06 * y_span,
          text = paste0("Linear R² = ", signif(r_squared, 3)),
          showarrow = FALSE,
          xanchor = "left",
          yanchor = "top",
          bgcolor = "white",
          bordercolor = "#d9e2e7",
          borderwidth = 1,
          font = list(color = "#2f3b43", size = 13)
        )
      ),
      legend = list(orientation = "h", y = -0.2),
      margin = list(l = 72, r = 24, t = 42, b = 92)
    )
}

plot_external_crispr_rank <- function(df, top_n = 20) {
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))
  if (is.null(top_n) || is.na(top_n)) top_n <- 20

  plot_df <- df |>
    filter(!is.na(.data$effect_value)) |>
    arrange(desc(abs(.data$effect_value))) |>
    slice_head(n = top_n) |>
    arrange(.data$effect_value) |>
    mutate(
      gene_symbol = forcats::fct_reorder(.data$gene_symbol, .data$effect_value),
      label = paste0(
        .data$gene_symbol,
        "<br>", .data$direction,
        "<br>LFC: ", signif(.data$effect_value, 3),
        "<br>Rank: ", .data$rank
      )
    )

  p <- ggplot(plot_df, aes(x = effect_value, y = gene_symbol, color = direction, text = label)) +
    geom_vline(xintercept = 0, color = "#b9b9b9", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = effect_value, yend = gene_symbol), linewidth = 0.5, alpha = 0.65) +
    geom_point(size = 2.8, alpha = 0.9) +
    scale_color_manual(values = effect_color_values, na.value = "#9a9a9a") +
    labs(x = "Median sgRNA log2 fold-change", y = NULL, color = "Direction") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  ggplotly(p, tooltip = "text") |>
    layout(legend = list(orientation = "h", y = -0.18))
}

screen_panel_levels <- function(df) {
  groups <- screen_view_groups[screen_view_groups %in% as.character(df$screen_group)]
  unlist(lapply(groups, function(group) {
    labels <- unique(screen_view_registry$short_label[screen_view_registry$screen_group == group])
    c(paste0("<b>", group, "</b>"), labels[labels %in% df$short_label])
  }))
}

screen_panel_height <- function(df) {
  n_genes <- dplyr::n_distinct(df$gene_key)
  n_rows <- length(screen_panel_levels(df))
  panel_rows <- ceiling(n_genes / 4)
  max(360, 110 + panel_rows * (22 * n_rows + 40))
}

plot_screen_panels <- function(df, gene_order = unique(df$gene_key), metric = "percentile") {
  spec <- screen_metrics[[metric]]
  df <- df |>
    mutate(value = .data[[spec$column]]) |>
    filter(!is.na(.data$value) | (spec$hit_list_column & .data$hit_list_only))
  if (!nrow(df)) return(plotly_empty(type = "scatter", mode = "markers"))

  y_levels <- screen_panel_levels(df)
  finite_values <- df$value[is.finite(df$value)]
  value_range <- range(c(0, finite_values, if (is.finite(spec$band)) c(-spec$band, spec$band)))
  if (identical(metric, "percentile")) value_range <- c(-1, 1)
  span <- diff(value_range)
  hit_list_x <- value_range[2] + 0.17 * span
  show_hit_list <- spec$hit_list_column && any(df$hit_list_only)
  x_limits <- c(value_range[1] - 0.04 * span, if (show_hit_list) hit_list_x + 0.06 * span else value_range[2] + 0.04 * span)

  plot_df <- df |>
    group_by(.data$gene_key) |>
    mutate(gene_display = format_gene_symbols(.data$gene_symbol)) |>
    ungroup() |>
    mutate(
      x = dplyr::if_else(spec$hit_list_column & .data$hit_list_only, hit_list_x, .data$value),
      marker = dplyr::case_when(
        spec$hit_list_column & .data$hit_list_only ~ "In published hit list",
        .data$supported ~ "FDR ≤ 0.05 or paper threshold",
        TRUE ~ "Not significant or no statistics"
      ),
      y = factor(.data$short_label, levels = rev(y_levels)),
      label = paste0(
        .data$gene_display, " · ", .data$short_label,
        "<br>", .data$paper, " — ", .data$experiment,
        dplyr::if_else(
          .data$hit_list_only,
          "<br>Listed in a published hit table (not genome-wide; no percentile or z-score)",
          paste0(
            "<br>Within-screen percentile: ", sprintf("%+.2f", .data$percentile), " of ", format(.data$screen_genes, big.mark = ","), " genes",
            "<br>Robust z-score: ", sprintf("%+.1f", .data$robust_z)
          )
        ),
        "<br>", .data$effect_metric, ": ", signif(.data$effect_value, 3),
        "<br>FDR: ", dplyr::case_when(
          is.na(.data$fdr) ~ "not reported",
          .data$fdr <= 0 ~ "reported as 0",
          TRUE ~ as.character(signif(.data$fdr, 3))
        ),
        "<br>", .data$significance
      )
    )
  plot_df$gene_display <- factor(plot_df$gene_display, levels = unique(plot_df$gene_display[order(match(plot_df$gene_key, gene_order))]))

  marker_shapes <- c(
    "FDR ≤ 0.05 or paper threshold" = 16,
    "Not significant or no statistics" = 1,
    "In published hit list" = 18
  )

  breaks <- if (identical(metric, "percentile")) c(-1, 0, 1) else pretty(value_range, n = 4)
  breaks <- breaks[breaks >= value_range[1] & breaks <= value_range[2]]
  break_labels <- sub("^-", "−", format(breaks, trim = TRUE))
  if (identical(metric, "percentile")) break_labels <- c("−1", "0", "+1")
  if (show_hit_list) {
    breaks <- c(breaks, hit_list_x)
    break_labels <- c(break_labels, "hit list")
  }

  p <- ggplot(plot_df, aes(x = x, y = y, shape = marker, text = label))
  if (is.finite(spec$band)) {
    p <- p +
      annotate("rect", xmin = spec$band, xmax = value_range[2] + 0.02 * span, ymin = 0.5, ymax = length(y_levels) + 0.5, fill = "#b73e3e", alpha = 0.1) +
      annotate("rect", xmin = value_range[1] - 0.02 * span, xmax = -spec$band, ymin = 0.5, ymax = length(y_levels) + 0.5, fill = "#247a6b", alpha = 0.1)
  }
  p <- p +
    geom_vline(xintercept = 0, color = "#b9b9b9", linewidth = 0.35) +
    geom_point(size = 3, color = "#2a6fb0", stroke = 1) +
    scale_shape_manual(values = marker_shapes, name = NULL, drop = FALSE) +
    scale_y_discrete(limits = rev(y_levels), drop = FALSE) +
    scale_x_continuous(limits = x_limits, breaks = breaks, labels = break_labels) +
    facet_wrap(~gene_display, ncol = min(4, dplyr::n_distinct(plot_df$gene_key))) +
    labs(x = spec$axis, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#f0f0ee", linewidth = 0.3),
      panel.border = element_rect(color = "#aab4bb", fill = NA, linewidth = 0.6),
      panel.spacing.x = grid::unit(14, "pt"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      axis.text.y = element_text(size = 10.5)
    )

  ggplotly(p, tooltip = "text") |>
    layout(showlegend = FALSE, margin = list(t = 50, b = 60))
}
