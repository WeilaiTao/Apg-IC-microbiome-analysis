```r
#!/usr/bin/env Rscript

# ============================================================
# Supplementary Figure S2A-C: Phase-specific production deviations
# Feed intake, HDEP, and mortality deviations relative to production targets
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation
# in laying hens
#
# Usage:
#   Rscript FigS2_ABC_phase_deviation.R
#   Rscript FigS2_ABC_phase_deviation.R data/IC_hen.xlsx results/FigS2_ABC_phase_deviation.pdf
#
# Required input columns:
#   Day
#   Feed intake
#   Feed-intake target
#   Mortality
#   Mortality target
#   HDEP
#   HDEP target
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(length(args) >= 1, args[1], "data/IC_hen.xlsx")
output_pdf <- ifelse(length(args) >= 2, args[2], "results/FigS2_ABC_phase_deviation.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Please provide the Excel file path as the first argument, for example:\n",
    "Rscript FigS2_ABC_phase_deviation.R data/IC_hen.xlsx"
  )
}

# ----------------------------
# 1. Read and validate data
# ----------------------------
required_cols <- c(
  "Day",
  "Feed intake",
  "Feed-intake target",
  "Mortality",
  "Mortality target",
  "HDEP",
  "HDEP target"
)

raw_df <- read_excel(input_file)

missing_cols <- setdiff(required_cols, colnames(raw_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df_clean <- raw_df %>%
  mutate(
    Day = suppressWarnings(as.numeric(`Day`)),
    `Feed intake` = suppressWarnings(as.numeric(`Feed intake`)),
    `Feed-intake target` = suppressWarnings(as.numeric(`Feed-intake target`)),
    Mortality = suppressWarnings(as.numeric(`Mortality`)),
    `Mortality target` = suppressWarnings(as.numeric(`Mortality target`)),
    HDEP = suppressWarnings(as.numeric(`HDEP`)),
    `HDEP target` = suppressWarnings(as.numeric(`HDEP target`))
  ) %>%
  filter(!is.na(Day)) %>%
  mutate(
    Delta_Feed = `Feed intake` - `Feed-intake target`,
    Delta_HDEP = HDEP - `HDEP target`,
    Excess_Mortality = Mortality - `Mortality target`,
    Phase = case_when(
      Day <= 0 ~ "Pre-outbreak",
      Day >= 1 & Day <= 7 ~ "Early outbreak",
      Day >= 8 & Day <= 16 ~ "Peak outbreak",
      Day >= 17 ~ "Recovery phase",
      TRUE ~ NA_character_
    ),
    Phase = factor(
      Phase,
      levels = c(
        "Pre-outbreak",
        "Early outbreak",
        "Peak outbreak",
        "Recovery phase"
      )
    )
  ) %>%
  filter(!is.na(Phase))

if (nrow(df_clean) == 0) {
  stop("No valid rows were retained after phase assignment.")
}

# ----------------------------
# 2. Plotting parameters
# ----------------------------
phase_colors <- c(
  "Pre-outbreak" = "#4DBBD5",
  "Early outbreak" = "#F39B7F",
  "Peak outbreak" = "#DC0000",
  "Recovery phase" = "#00A087"
)

selected_comparisons <- data.frame(
  group1 = c(
    "Pre-outbreak",
    "Pre-outbreak",
    "Pre-outbreak",
    "Early outbreak"
  ),
  group2 = c(
    "Early outbreak",
    "Peak outbreak",
    "Recovery phase",
    "Peak outbreak"
  ),
  stringsAsFactors = FALSE
)

selected_keys <- paste(selected_comparisons$group1, selected_comparisons$group2, sep = "___")

# ----------------------------
# 3. Helper function
# ----------------------------
plot_phase_deviation <- function(data, y_var, y_label, panel_title) {
  formula_obj <- as.formula(paste(y_var, "~ Phase"))

  # Kruskal-Wallis test is calculated for reporting consistency.
  kw_test <- data %>%
    kruskal_test(formula_obj)

  # Dunn's post hoc test with Benjamini-Hochberg adjustment.
  stat_tbl <- data %>%
    dunn_test(formula_obj, p.adjust.method = "BH") %>%
    mutate(pair_key = paste(group1, group2, sep = "___")) %>%
    filter(pair_key %in% selected_keys) %>%
    add_xy_position(x = "Phase", step.increase = 0.10)

  p <- ggplot(data, aes(x = Phase, y = .data[[y_var]])) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "#B2182B",
      linewidth = 0.75,
      alpha = 0.75
    ) +
    geom_boxplot(
      aes(fill = Phase),
      width = 0.55,
      alpha = 0.35,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.75
    ) +
    geom_jitter(
      aes(fill = Phase),
      shape = 21,
      color = "black",
      stroke = 0.45,
      size = 2.8,
      alpha = 0.70,
      width = 0.15
    ) +
    scale_fill_manual(values = phase_colors) +
    labs(
      x = "Outbreak phase",
      y = y_label,
      title = panel_title
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0,
        margin = margin(b = 8)
      ),
      axis.title = element_text(face = "bold", color = "black"),
      axis.text.x = element_text(
        angle = 25,
        hjust = 1,
        color = "black",
        size = 11
      ),
      axis.text.y = element_text(color = "black", size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  if (nrow(stat_tbl) > 0) {
    p <- p +
      stat_pvalue_manual(
        stat_tbl,
        label = "p.adj.signif",
        hide.ns = TRUE,
        tip.length = 0.015,
        bracket.size = 0.65,
        size = 4
      )
  }

  attr(p, "kruskal_test") <- kw_test
  attr(p, "dunn_test") <- stat_tbl

  return(p)
}

# ----------------------------
# 4. Generate panels
# ----------------------------
p_feed <- plot_phase_deviation(
  data = df_clean,
  y_var = "Delta_Feed",
  y_label = expression(Delta ~ "Feed intake (g/hen/day)"),
  panel_title = "A  Feed intake deviation"
)

p_hdep <- plot_phase_deviation(
  data = df_clean,
  y_var = "Delta_HDEP",
  y_label = expression(Delta ~ "HDEP (%)"),
  panel_title = "B  HDEP deviation"
)

p_mortality <- plot_phase_deviation(
  data = df_clean,
  y_var = "Excess_Mortality",
  y_label = "Excess mortality (birds/day)",
  panel_title = "C  Excess mortality"
)

final_plot <- p_feed + p_hdep + p_mortality +
  plot_layout(ncol = 3)

# ----------------------------
# 5. Export figure
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 15.5,
  height = 5.8,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 15.5,
  height = 5.8,
  units = "in",
  dpi = 300
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)

# ----------------------------
# 6. Export statistical tables and session information
# ----------------------------
kruskal_results <- bind_rows(
  Feed_intake = attr(p_feed, "kruskal_test"),
  HDEP = attr(p_hdep, "kruskal_test"),
  Excess_mortality = attr(p_mortality, "kruskal_test"),
  .id = "Metric"
)

dunn_results <- bind_rows(
  Feed_intake = attr(p_feed, "dunn_test"),
  HDEP = attr(p_hdep, "dunn_test"),
  Excess_mortality = attr(p_mortality, "dunn_test"),
  .id = "Metric"
)

write.csv(
  kruskal_results,
  file = file.path(output_dir, "FigS2_ABC_kruskal_results.csv"),
  row.names = FALSE
)

write.csv(
  dunn_results,
  file = file.path(output_dir, "FigS2_ABC_dunn_BH_results.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_FigS2_ABC.txt")
)

message("Statistical results and session information saved to: ", output_dir)
```
