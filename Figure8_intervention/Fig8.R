#!/usr/bin/env Rscript
# ============================================================
# Figure 8C: Onset interval after baffle intervention
# Days to outbreak in adjacent houses under natural transmission and baffle intervention
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript Fig8C_onset_interval_baffle_intervention.R
#   Rscript Fig8C_onset_interval_baffle_intervention.R data/Fig8C_onset_interval.csv results/Fig8C_onset_interval_baffle_intervention.pdf
#
# Notes:
#   Optional input CSV should contain:
#   Group, Days
#
#   If no input file is provided or the default input file is not found,
#   the embedded onset-interval data are used.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(length(args) >= 1, args[1], "data/Fig8C_onset_interval.csv")
output_pdf <- ifelse(length(args) >= 2, args[2], "results/Fig8C_onset_interval_baffle_intervention.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

group_levels <- c("Natural transmission", "Baffle intervention")

# ----------------------------
# 1. Read data
# ----------------------------
if (file.exists(input_file)) {
  raw_df <- read_csv(input_file, show_col_types = FALSE)
} else {
  natural <- c(0, 6, 1, 1, 3, 4, 5, 1, 5)
  baffle <- c(7, 27, 20, 57)

  raw_df <- data.frame(
    Group = c(
      rep("Natural transmission", length(natural)),
      rep("Baffle intervention", length(baffle))
    ),
    Days = c(natural, baffle)
  )
}

required_cols <- c("Group", "Days")
missing_cols <- setdiff(required_cols, colnames(raw_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

plot_df <- raw_df %>%
  transmute(
    Group = as.character(Group),
    Days = suppressWarnings(as.numeric(Days))
  ) %>%
  filter(Group %in% group_levels, is.finite(Days)) %>%
  mutate(Group = factor(Group, levels = group_levels))

if (nrow(plot_df) == 0) {
  stop("No valid rows were retained after reading the input data.")
}

if (any(table(plot_df$Group) == 0)) {
  stop("Both groups must contain at least one observation.")
}

# ----------------------------
# 2. Statistical analysis
# ----------------------------
wilcox_result <- wilcox.test(
  Days ~ Group,
  data = plot_df,
  exact = FALSE
)

stat_df <- data.frame(
  comparison = "Natural transmission vs Baffle intervention",
  p_value = wilcox_result$p.value
)

# ----------------------------
# 3. Plotting parameters
# ----------------------------
group_colors <- c(
  "Natural transmission" = "#D9534F",
  "Baffle intervention" = "#428BCA"
)

# ----------------------------
# 4. Generate plot
# ----------------------------
set.seed(123)

p <- ggplot(plot_df, aes(x = Group, y = Days, fill = Group)) +
  geom_boxplot(
    width = 0.50,
    outlier.shape = NA,
    alpha = 0.45,
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_jitter(
    shape = 21,
    width = 0.08,
    size = 2.2,
    color = "black",
    stroke = 0.3,
    alpha = 0.85
  ) +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(exact = FALSE),
    label = "p.format",
    label.y = max(plot_df$Days, na.rm = TRUE) * 1.08,
    size = 3.5
  ) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(
    breaks = seq(0, 60, 10),
    limits = c(0, max(65, max(plot_df$Days, na.rm = TRUE) * 1.15)),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = NULL,
    y = "Days to outbreak"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
  )

# ----------------------------
# 5. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p,
  width = 4.6,
  height = 4.2,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p,
  width = 4.6,
  height = 4.2,
  units = "in",
  dpi = 300
)

write.csv(
  plot_df,
  file = file.path(output_dir, "Fig8C_onset_interval_values.csv"),
  row.names = FALSE
)

write.csv(
  stat_df,
  file = file.path(output_dir, "Fig8C_onset_interval_wilcox_result.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Values and statistical results saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig8C_onset_interval.txt")
)