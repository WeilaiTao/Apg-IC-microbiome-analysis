#!/usr/bin/env Rscript
# ============================================================
# Supplementary Figure S5B: Serum immunoglobulins and cytokines
# ELISA-based comparison of IgA, IgG, IgM, IFN-gamma, and IL-1beta
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript FigS5B_serum_ELISA_markers.R
#   Rscript FigS5B_serum_ELISA_markers.R data/elisa.xlsx results/FigS5B_serum_ELISA_markers.pdf
#
# Notes:
#   The input Excel file should contain:
#   Group, IgA, IgG, IgM, IFN-gamma, IL-1beta
#
#   Groups retained: EH, MLD, MOD, SEV
#   Statistical comparisons are performed against the EH group.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(length(args) >= 1, args[1], "data/elisa.xlsx")
output_pdf <- ifelse(length(args) >= 2, args[2], "results/FigS5B_serum_ELISA_markers.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# ----------------------------
# 1. Read and validate data
# ----------------------------
df_raw <- read_excel(input_file)
names(df_raw) <- trimws(names(df_raw))

if (!"Group" %in% colnames(df_raw)) {
  stop("Missing required column: Group")
}

if (ncol(df_raw) < 2) {
  stop("The input file should contain Group and at least one ELISA marker column.")
}

group_levels <- c("EH", "MLD", "MOD", "SEV")

df_long <- df_raw %>%
  mutate(Group = trimws(as.character(Group))) %>%
  pivot_longer(
    cols = -Group,
    names_to = "Marker_raw",
    values_to = "Value"
  ) %>%
  mutate(
    Value = suppressWarnings(as.numeric(Value)),
    Marker_raw = trimws(as.character(Marker_raw))
  ) %>%
  filter(Group %in% group_levels, is.finite(Value)) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    Marker = case_when(
      grepl("IgA", Marker_raw, ignore.case = TRUE) ~ "IgA (ug/mL)",
      grepl("IgG", Marker_raw, ignore.case = TRUE) ~ "IgG (ug/mL)",
      grepl("IgM", Marker_raw, ignore.case = TRUE) ~ "IgM (ug/mL)",
      grepl("IFN", Marker_raw, ignore.case = TRUE) ~ "IFN-gamma (pg/mL)",
      grepl("IL-1|IL1", Marker_raw, ignore.case = TRUE) ~ "IL-1beta (pg/mL)",
      TRUE ~ Marker_raw
    ),
    Marker = factor(
      Marker,
      levels = c(
        "IgA (ug/mL)",
        "IgG (ug/mL)",
        "IgM (ug/mL)",
        "IFN-gamma (pg/mL)",
        "IL-1beta (pg/mL)"
      )
    )
  ) %>%
  filter(!is.na(Marker))

if (nrow(df_long) == 0) {
  stop("No valid rows were retained after filtering.")
}

# ----------------------------
# 2. Statistical analysis
# ----------------------------
comparisons <- list(
  c("EH", "MLD"),
  c("EH", "MOD"),
  c("EH", "SEV")
)

stat_tbl <- df_long %>%
  group_by(Marker) %>%
  pairwise_wilcox_test(
    Value ~ Group,
    comparisons = comparisons,
    p.adjust.method = "BH"
  ) %>%
  ungroup()

# ----------------------------
# 3. Plotting parameters
# ----------------------------
group_colors <- c(
  "EH"  = "#9CB9CA",
  "MLD" = "#F8BF99",
  "MOD" = "#BE6C58",
  "SEV" = "#9C7190"
)

# ----------------------------
# 4. Generate plot
# ----------------------------
p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.7,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.85
  ) +
  facet_wrap(~ Marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    tip.length = 0.015,
    step.increase = 0.10,
    size = 3.5,
    bracket.size = 0.5
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = NULL,
    y = "Concentration",
    caption = "Statistical significance vs. EH group, Wilcoxon rank-sum test"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", color = "grey50", linewidth = 0.5),
    plot.caption = element_text(hjust = 0, color = "grey30")
  )

# ----------------------------
# 5. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p,
  width = 9.5,
  height = 6.2,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p,
  width = 9.5,
  height = 6.2,
  units = "in",
  dpi = 300
)

write.csv(
  df_long,
  file = file.path(output_dir, "FigS5B_serum_ELISA_values.csv"),
  row.names = FALSE
)

write.csv(
  stat_tbl,
  file = file.path(output_dir, "FigS5B_serum_ELISA_wilcox_BH_results.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Statistical results saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_FigS5B_serum_ELISA.txt")
)