#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure S6A-C: IOS subgroup-level microbiome metrics
# Shannon diversity, genus richness, and Av. paragallinarum relative abundance
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript FigS6_ABC_IOS_subgroup_microbiome_metrics.R
#   Rscript FigS6_ABC_IOS_subgroup_microbiome_metrics.R data/combined_otu_G.xlsx data/GroupID.xlsx results/FigS6_ABC_IOS_subgroup_microbiome_metrics.pdf
#
# Notes:
#   Abundance table: first column = feature names; remaining columns = sample IDs
#   Metadata table: Sample, Organ, Group2
#   Groups retained: HC, EH, MLD, MOD, SEV
#   The target feature for panel C is searched using the pattern "Avibacterium paragallinarum|paragallinarum".
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

abundance_file <- ifelse(length(args) >= 1, args[1], "data/combined_otu_G.xlsx")
metadata_file  <- ifelse(length(args) >= 2, args[2], "data/GroupID.xlsx")
output_pdf     <- ifelse(length(args) >= 3, args[3], "results/FigS6_ABC_IOS_subgroup_microbiome_metrics.pdf")
output_png     <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(abundance_file)) {
  stop("Abundance file not found: ", abundance_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

target_organ <- "IOS"
group_levels <- c("HC", "EH", "MLD", "MOD", "SEV")
target_feature_pattern <- "Avibacterium paragallinarum|paragallinarum"

# ----------------------------
# 1. Read and validate data
# ----------------------------
abundance_raw <- read_excel(abundance_file, sheet = 1)

if (ncol(abundance_raw) < 2) {
  stop("The abundance table should contain one feature column and at least one sample column.")
}

colnames(abundance_raw)[1] <- "Feature"

abundance_df <- abundance_raw %>%
  mutate(Feature = as.character(Feature)) %>%
  filter(!is.na(Feature), Feature != "") %>%
  mutate(across(-Feature, ~ suppressWarnings(as.numeric(.x))))

metadata_raw <- read_excel(metadata_file) %>%
  mutate(across(everything(), as.character))

required_meta_cols <- c("Sample", "Organ", "Group2")
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata_raw))

if (length(missing_meta_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

metadata_df <- metadata_raw %>%
  filter(Organ == target_organ, Group2 %in% group_levels) %>%
  transmute(
    Sample = Sample,
    Organ = Organ,
    Group = Group2
  )

sample_cols <- setdiff(colnames(abundance_df), "Feature")
matched_samples <- intersect(sample_cols, metadata_df$Sample)

if (length(matched_samples) < 2) {
  stop("Fewer than two matched samples were found between abundance table and metadata.")
}

abundance_df <- abundance_df %>%
  select(Feature, all_of(matched_samples))

# ----------------------------
# 2. Calculate Shannon, richness, and Av. paragallinarum abundance
# ----------------------------
calc_shannon <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x <- x[x > 0]
  if (length(x) == 0 || sum(x) == 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log(p))
}

mat <- abundance_df %>%
  column_to_rownames("Feature") %>%
  as.data.frame()

mat[] <- lapply(mat, function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
})

shannon_values <- apply(mat, 2, calc_shannon)
richness_values <- apply(mat, 2, function(x) sum(x > 0, na.rm = TRUE))

target_features <- grep(
  target_feature_pattern,
  rownames(mat),
  ignore.case = TRUE,
  value = TRUE
)

if (length(target_features) == 0) {
  warning("No feature matched the target pattern: ", target_feature_pattern)
  target_abundance <- rep(NA_real_, ncol(mat))
  names(target_abundance) <- colnames(mat)
} else {
  target_raw <- colSums(mat[target_features, , drop = FALSE], na.rm = TRUE)
  total_raw <- colSums(mat, na.rm = TRUE)
  target_abundance <- ifelse(total_raw > 0, target_raw / total_raw * 100, NA_real_)
}

metric_df <- tibble(
  Sample = names(shannon_values),
  Shannon = as.numeric(shannon_values),
  Richness = as.numeric(richness_values),
  Apg_relative_abundance = as.numeric(target_abundance[Sample])
) %>%
  inner_join(metadata_df, by = "Sample") %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Organ = factor(Organ, levels = target_organ)
  ) %>%
  arrange(Group, Sample)

if (nrow(metric_df) == 0) {
  stop("No valid samples were retained for plotting.")
}

if (any(table(metric_df$Group) == 0)) {
  warning("At least one target group has no matched samples.")
}

plot_long <- metric_df %>%
  pivot_longer(
    cols = c(Shannon, Richness, Apg_relative_abundance),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("Shannon", "Richness", "Apg_relative_abundance"),
      labels = c(
        "A  Shannon diversity",
        "B  Genus richness",
        "C  Av. paragallinarum relative abundance"
      )
    )
  )

# ----------------------------
# 3. Statistical analysis
# ----------------------------
kruskal_results <- plot_long %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Group) %>%
  ungroup()

dunn_results <- plot_long %>%
  group_by(Metric) %>%
  dunn_test(Value ~ Group, p.adjust.method = "BH") %>%
  ungroup() %>%
  mutate(
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

sig_results <- dunn_results %>%
  filter(p.signif != "ns") %>%
  group_by(Metric) %>%
  arrange(p.adj, .by_group = TRUE) %>%
  mutate(
    y.position = max(plot_long$Value[plot_long$Metric == first(Metric)], na.rm = TRUE) +
      seq_len(n()) * 0.06 * diff(range(plot_long$Value[plot_long$Metric == first(Metric)], na.rm = TRUE))
  ) %>%
  ungroup()

sig_results$y.position[!is.finite(sig_results$y.position)] <- NA_real_

# ----------------------------
# 4. Plotting function
# ----------------------------
group_colors <- c(
  "HC"  = "#A2E8EA",
  "EH"  = "#9CB9CA",
  "MLD" = "#F8BF99",
  "MOD" = "#BE6C58",
  "SEV" = "#9C7190"
)

plot_metric <- function(data, metric_name, y_label) {
  d <- data %>%
    filter(Metric == metric_name)

  stat_d <- sig_results %>%
    filter(Metric == metric_name, is.finite(y.position))

  p <- ggplot(d, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.65,
      color = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      shape = 21,
      width = 0.15,
      size = 1.8,
      color = "black",
      stroke = 0.3,
      alpha = 0.85
    ) +
    scale_fill_manual(values = group_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      x = NULL,
      y = y_label,
      title = as.character(metric_name)
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0),
      panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
    )

  if (nrow(stat_d) > 0) {
    p <- p +
      stat_pvalue_manual(
        data = stat_d,
        label = "p.signif",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 3.2,
        inherit.aes = FALSE
      )
  }

  p
}

# ----------------------------
# 5. Generate panels
# ----------------------------
p_shannon <- plot_metric(
  data = plot_long,
  metric_name = "A  Shannon diversity",
  y_label = "Genus Shannon index"
)

p_richness <- plot_metric(
  data = plot_long,
  metric_name = "B  Genus richness",
  y_label = "Genus richness"
)

p_apg <- plot_metric(
  data = plot_long,
  metric_name = "C  Av. paragallinarum relative abundance",
  y_label = "Relative abundance (%)"
)

final_plot <- p_shannon + p_richness + p_apg +
  plot_layout(ncol = 3)

# ----------------------------
# 6. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 11.0,
  height = 4.2,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 11.0,
  height = 4.2,
  units = "in",
  dpi = 300
)

write.csv(
  metric_df,
  file = file.path(output_dir, "FigS6_ABC_IOS_microbiome_metrics.csv"),
  row.names = FALSE
)

write.csv(
  kruskal_results,
  file = file.path(output_dir, "FigS6_ABC_IOS_kruskal_results.csv"),
  row.names = FALSE
)

write.csv(
  dunn_results,
  file = file.path(output_dir, "FigS6_ABC_IOS_dunn_BH_results.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Metrics and statistical results saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_FigS6_ABC_IOS.txt")
)
