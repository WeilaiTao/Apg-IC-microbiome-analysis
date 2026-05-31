#!/usr/bin/env Rscript
# ============================================================
# Figure 6B-C: ARG diversity and richness across habitats
# ARG Shannon index and ARG subtype richness
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript Fig6BC_ARG_diversity_richness.R
#   Rscript Fig6BC_ARG_diversity_richness.R data/ALL.arg_summary.cleaned_v3.xlsx data/all_group.xlsx results/Fig6BC_ARG_diversity_richness.pdf
#
# Notes:
#   ARG long table should contain:
#   Sample, Unified_ARG, ARG_Abundance
#
#   Metadata table should contain:
#   Sample and Group2
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(rstatix)
  library(patchwork)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

arg_file    <- ifelse(length(args) >= 1, args[1], "data/ALL.arg_summary.cleaned_v3.xlsx")
group_file  <- ifelse(length(args) >= 2, args[2], "data/all_group.xlsx")
output_pdf  <- ifelse(length(args) >= 3, args[3], "results/Fig6BC_ARG_diversity_richness.pdf")
output_png  <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(arg_file)) {
  stop("ARG file not found: ", arg_file)
}

if (!file.exists(group_file)) {
  stop("Group file not found: ", group_file)
}

group_col <- "Group2"

target_groups <- c(
  "AISLE", "ICin", "ICinlet", "HEin", "HEinlet",
  "CE_HE", "CE_IC", "IOS_HE", "IOS_IC"
)

# ----------------------------
# 1. Read and validate data
# ----------------------------
arg_raw <- read_excel(arg_file, sheet = 1)
required_arg_cols <- c("Sample", "Unified_ARG", "ARG_Abundance")

missing_arg_cols <- setdiff(required_arg_cols, colnames(arg_raw))
if (length(missing_arg_cols) > 0) {
  stop("Missing required ARG columns: ", paste(missing_arg_cols, collapse = ", "))
}

arg_long <- arg_raw %>%
  transmute(
    Sample = as.character(Sample) %>% str_squish(),
    Unified_ARG = as.character(Unified_ARG) %>% str_squish(),
    Abundance = suppressWarnings(as.numeric(ARG_Abundance))
  ) %>%
  filter(!is.na(Sample), Sample != "",
         !is.na(Unified_ARG), Unified_ARG != "") %>%
  mutate(Abundance = ifelse(is.na(Abundance) | Abundance < 0, 0, Abundance))

group_raw <- read_excel(group_file)
names(group_raw) <- trimws(names(group_raw))

if (!"Sample" %in% names(group_raw)) {
  names(group_raw)[1] <- "Sample"
}

if (!group_col %in% names(group_raw)) {
  stop("Missing required group column: ", group_col)
}

group_df <- group_raw %>%
  transmute(
    Sample = as.character(Sample) %>% str_squish(),
    Group = as.character(.data[[group_col]]) %>% str_squish()
  ) %>%
  filter(!is.na(Sample), Sample != "",
         !is.na(Group), Group != "") %>%
  filter(Group %in% target_groups) %>%
  distinct(Sample, .keep_all = TRUE)

arg_long <- arg_long %>%
  filter(Sample %in% group_df$Sample)

if (nrow(arg_long) == 0) {
  stop("No ARG records were retained after matching with metadata.")
}

# ----------------------------
# 2. Build ARG abundance matrix
# ----------------------------
arg_sum <- arg_long %>%
  group_by(Sample, Unified_ARG) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

arg_wide <- arg_sum %>%
  pivot_wider(
    names_from = Unified_ARG,
    values_from = Abundance,
    values_fill = 0
  )

arg_wide <- arg_wide %>%
  filter(Sample %in% group_df$Sample) %>%
  arrange(match(Sample, group_df$Sample))

arg_matrix <- arg_wide %>%
  select(-Sample) %>%
  as.data.frame()

arg_matrix[] <- lapply(arg_matrix, function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x) | x < 0] <- 0
  x
})

# ----------------------------
# 3. Calculate ARG Shannon index and richness
# ----------------------------
calc_shannon <- function(v) {
  v <- as.numeric(v)
  v[is.na(v) | v < 0] <- 0
  s <- sum(v)
  if (s <= 0) return(0)
  p <- v / s
  p <- p[p > 0]
  -sum(p * log(p))
}

metric_df <- tibble(
  Sample = arg_wide$Sample,
  ARG_Shannon = apply(arg_matrix, 1, calc_shannon),
  ARG_Richness = apply(arg_matrix, 1, function(x) sum(x > 0, na.rm = TRUE))
) %>%
  inner_join(group_df, by = "Sample") %>%
  mutate(Group = factor(Group, levels = target_groups))

if (nrow(metric_df) == 0) {
  stop("No valid samples were retained for plotting.")
}

# ----------------------------
# 4. Statistical analysis
# ----------------------------
metric_long <- metric_df %>%
  pivot_longer(
    cols = c(ARG_Shannon, ARG_Richness),
    names_to = "Metric",
    values_to = "Value"
  )

kruskal_results <- metric_long %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Group) %>%
  ungroup()

wilcox_results <- metric_long %>%
  group_by(Metric) %>%
  pairwise_wilcox_test(Value ~ Group, p.adjust.method = "BH") %>%
  ungroup() %>%
  mutate(
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# ----------------------------
# 5. Plotting parameters
# ----------------------------
group_colors <- c(
  "AISLE"   = "#9EC5DC",
  "ICin"    = "#FFAAAA",
  "ICinlet" = "#C6AFE9",
  "HEin"    = "#AFE9AF",
  "HEinlet" = "#AFE9DD",
  "CE_HE"   = "#F7A4C2",
  "CE_IC"   = "#DD588E",
  "IOS_HE"  = "#96BACC",
  "IOS_IC"  = "#C9A69B"
)

plot_box <- function(data, y_var, y_label, panel_title) {
  ggplot(data, aes(x = Group, y = .data[[y_var]], fill = Group)) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.5,
      alpha = 0.55
    ) +
    geom_jitter(
      shape = 21,
      width = 0.15,
      size = 1.8,
      color = "black",
      stroke = 0.3,
      alpha = 0.8
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      x = NULL,
      y = y_label,
      title = panel_title
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
    )
}

# ----------------------------
# 6. Generate panels
# ----------------------------
p_shannon <- plot_box(
  data = metric_df,
  y_var = "ARG_Shannon",
  y_label = "ARG Shannon index",
  panel_title = "B  ARG Shannon index"
)

p_richness <- plot_box(
  data = metric_df,
  y_var = "ARG_Richness",
  y_label = "ARG subtype richness",
  panel_title = "C  ARG subtype richness"
)

final_plot <- p_shannon + p_richness +
  plot_layout(ncol = 2)

# ----------------------------
# 7. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 10.0,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 10.0,
  height = 4.8,
  units = "in",
  dpi = 300
)

write.csv(
  metric_df,
  file = file.path(output_dir, "Fig6BC_ARG_diversity_richness_values.csv"),
  row.names = FALSE
)

write.csv(
  kruskal_results,
  file = file.path(output_dir, "Fig6BC_ARG_diversity_richness_kruskal.csv"),
  row.names = FALSE
)

write.csv(
  wilcox_results,
  file = file.path(output_dir, "Fig6BC_ARG_diversity_richness_wilcox_BH.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Values and statistical results saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig6BC_ARG_diversity_richness.txt")
)











#!/usr/bin/env Rscript
# ============================================================
# Figure 6G: Mobile ARG ratio across habitats
# Ratio of ARGs co-localized with structural MGEs within 10 kb
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript Fig6D_mobile_ARG_ratio_by_habitat.R
#   Rscript Fig6D_mobile_ARG_ratio_by_habitat.R data/summary_by_sample.tsv data/all_group.xlsx results/Fig6D_mobile_ARG_ratio_by_habitat.pdf
#
# Notes:
#   Ratio table should contain:
#   Sample, ARG_Total, ARG_with_NEARBY_structural_MGE_10kb, Ratio
#
#   Metadata table should contain:
#   Sample, Group
#
#   Groups retained: ENV, GUT, IOS
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(scales)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

ratio_file <- ifelse(length(args) >= 1, args[1], "data/summary_by_sample.tsv")
group_file <- ifelse(length(args) >= 2, args[2], "data/all_group.xlsx")
output_pdf <- ifelse(length(args) >= 3, args[3], "results/Fig6D_mobile_ARG_ratio_by_habitat.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(ratio_file)) {
  stop("Ratio file not found: ", ratio_file)
}

if (!file.exists(group_file)) {
  stop("Group file not found: ", group_file)
}

group_levels <- c("ENV", "GUT", "IOS")

# ----------------------------
# 1. Read and validate data
# ----------------------------
ratio_raw <- read_tsv(ratio_file, show_col_types = FALSE, progress = FALSE)

required_ratio_cols <- c(
  "Sample",
  "ARG_Total",
  "ARG_with_NEARBY_structural_MGE_10kb",
  "Ratio"
)

missing_ratio_cols <- setdiff(required_ratio_cols, colnames(ratio_raw))
if (length(missing_ratio_cols) > 0) {
  stop("Missing required ratio columns: ", paste(missing_ratio_cols, collapse = ", "))
}

ratio_df <- ratio_raw %>%
  transmute(
    Sample = str_squish(as.character(Sample)),
    ARG_Total = suppressWarnings(as.numeric(ARG_Total)),
    ARG_with_NEARBY_structural_MGE_10kb = suppressWarnings(as.numeric(ARG_with_NEARBY_structural_MGE_10kb)),
    Ratio = suppressWarnings(as.numeric(Ratio))
  ) %>%
  filter(!is.na(Sample), Sample != "") %>%
  mutate(
    ARG_Total = ifelse(is.na(ARG_Total) | ARG_Total < 0, 0, ARG_Total),
    ARG_with_NEARBY_structural_MGE_10kb = ifelse(
      is.na(ARG_with_NEARBY_structural_MGE_10kb) | ARG_with_NEARBY_structural_MGE_10kb < 0,
      0,
      ARG_with_NEARBY_structural_MGE_10kb
    ),
    Ratio = ifelse(is.na(Ratio), ARG_with_NEARBY_structural_MGE_10kb / ARG_Total, Ratio),
    Ratio = pmax(pmin(Ratio, 1), 0)
  )

group_raw <- read_excel(group_file, .name_repair = "unique")
names(group_raw) <- trimws(names(group_raw))

if (!"Sample" %in% names(group_raw)) {
  names(group_raw)[1] <- "Sample"
}

if (!"Group" %in% names(group_raw)) {
  stop("Missing required metadata column: Group")
}

group_df <- group_raw %>%
  transmute(
    Sample = str_squish(as.character(Sample)),
    Group = str_squish(as.character(Group))
  ) %>%
  filter(!is.na(Sample), Sample != "",
         !is.na(Group), Group != "") %>%
  filter(Group %in% group_levels) %>%
  distinct(Sample, .keep_all = TRUE)

plot_df <- ratio_df %>%
  inner_join(group_df, by = "Sample") %>%
  mutate(Group = factor(Group, levels = group_levels)) %>%
  filter(!is.na(Group), is.finite(Ratio))

if (nrow(plot_df) < 3) {
  stop("Too few samples were retained after merging ratio table and metadata.")
}

# ----------------------------
# 2. Statistical analysis
# ----------------------------
kruskal_result <- kruskal_test(plot_df, Ratio ~ Group)

stat_df <- pairwise_wilcox_test(
  plot_df,
  Ratio ~ Group,
  p.adjust.method = "BH"
) %>%
  mutate(
    p.adj.signif = case_when(
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01  ~ "**",
      p.adj <= 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

cmp_order <- c("ENV-GUT", "ENV-IOS", "GUT-IOS")
y_max <- max(plot_df$Ratio, na.rm = TRUE)
step <- max(0.05, y_max * 0.08)
base <- y_max + step * 0.6

stat_df <- stat_df %>%
  mutate(
    cmp = paste0(group1, "-", group2),
    cmp_rev = paste0(group2, "-", group1),
    cmp_key = ifelse(cmp %in% cmp_order, cmp, cmp_rev),
    cmp_key = factor(cmp_key, levels = cmp_order)
  ) %>%
  arrange(cmp_key) %>%
  mutate(y.position = base + (row_number() - 1) * step)

n_df <- plot_df %>%
  count(Group, name = "n") %>%
  mutate(label = paste0("n = ", n))

n_y <- base + nrow(stat_df) * step
y_upper <- max(1.02, n_y + step, max(stat_df$y.position, na.rm = TRUE) + step)

# ----------------------------
# 3. Plotting parameters
# ----------------------------
group_fill <- c(
  "ENV" = "#D7BADB",
  "GUT" = "#C5E0A4",
  "IOS" = "#CEBFB6"
)

group_point <- c(
  "ENV" = "#9E6FA8",
  "GUT" = "#6E9E49",
  "IOS" = "#7E736B"
)

# ----------------------------
# 4. Generate plot
# ----------------------------
set.seed(123)

p <- ggplot(plot_df, aes(x = Group, y = Ratio, fill = Group)) +
  geom_violin(
    width = 0.90,
    alpha = 0.55,
    color = NA,
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.25,
    outlier.shape = NA,
    alpha = 0.65,
    color = "grey20",
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = Group),
    position = position_jitter(width = 0.12, height = 0),
    size = 2.0,
    alpha = 0.70
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 1.8,
    color = "black"
  ) +
  geom_text(
    data = n_df,
    aes(x = Group, y = n_y, label = label),
    inherit.aes = FALSE,
    size = 3.0
  ) +
  stat_pvalue_manual(
    data = stat_df,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3.5,
    bracket.size = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_fill, drop = FALSE) +
  scale_color_manual(values = group_point, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, y_upper),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Mobile ARG ratio"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
  )

# ----------------------------
# 5. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p,
  width = 4.2,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p,
  width = 4.2,
  height = 4.8,
  units = "in",
  dpi = 300
)

write.csv(
  plot_df,
  file = file.path(output_dir, "Fig6D_mobile_ARG_ratio_values.csv"),
  row.names = FALSE
)

write.csv(
  kruskal_result,
  file = file.path(output_dir, "Fig6D_mobile_ARG_ratio_kruskal.csv"),
  row.names = FALSE
)

write.csv(
  stat_df,
  file = file.path(output_dir, "Fig6D_mobile_ARG_ratio_pairwise_wilcox_BH.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Values and statistical results saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig6D_mobile_ARG_ratio.txt")
)