#!/usr/bin/env Rscript

# ============================================================
# Figure 2B-D: Production impacts during infectious coryza outbreak
# Feed intake, HDEP, and mortality relative to production targets
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation
# in laying hens
#
# Usage:
#   Rscript Fig2_BCD_production_impacts.R
#   Rscript Fig2_BCD_production_impacts.R data/IC_hen.xlsx results/Fig2_BCD_production_impacts.pdf
#
# Notes:
#   The input Excel file should contain the following columns:
#   Day, Feed intake, Feed-intake target, Mortality, Mortality target,
#   HDEP, HDEP target
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(length(args) >= 1, args[1], "data/IC_hen.xlsx")
output_pdf <- ifelse(length(args) >= 2, args[2], "results/Fig2_BCD_production_impacts.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Please provide the Excel file path as the first argument, for example:\n",
    "Rscript Fig2_BCD_production_impacts.R data/IC_hen.xlsx"
  )
}

# ----------------------------
# 1. Read and validate data
# ----------------------------
required_cols <- c(
  "Day", "Feed intake", "Feed-intake target",
  "Mortality", "Mortality target",
  "HDEP", "HDEP target"
)

raw_df <- read_excel(input_file)

missing_cols <- setdiff(required_cols, colnames(raw_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- raw_df %>%
  rename(
    day         = `Day`,
    fi          = `Feed intake`,
    fi_target   = `Feed-intake target`,
    mort        = `Mortality`,
    mort_target = `Mortality target`,
    hdep        = `HDEP`,
    hdep_target = `HDEP target`
  ) %>%
  mutate(
    across(
      c(day, fi, fi_target, mort, mort_target, hdep, hdep_target),
      ~ suppressWarnings(as.numeric(.x))
    )
  ) %>%
  filter(!is.na(day)) %>%
  arrange(day)

if (nrow(df) == 0) {
  stop("No valid rows were retained after reading the input file.")
}

# ----------------------------
# 2. Helper functions
# ----------------------------

# Assign a new segment when days are not consecutive.
# This avoids ribbons being drawn across long gaps or sign changes.
add_segment_by_gap <- function(d, gap = 1.01) {
  if (nrow(d) == 0) return(d)
  d <- d[order(d$day), , drop = FALSE]
  d$segment <- cumsum(c(TRUE, diff(d$day) > gap))
  d
}

make_target_line_plot <- function(plot_df, y_label, panel_title, show_legend = FALSE) {
  surplus_df <- plot_df %>% filter(diff > 0) %>% add_segment_by_gap()
  deficit_df <- plot_df %>% filter(diff < 0) %>% add_segment_by_gap()

  p <- ggplot() +
    {
      if (nrow(deficit_df) > 0)
        geom_ribbon(
          data = deficit_df,
          aes(
            x = day, ymin = value, ymax = target,
            group = segment, fill = "Deficit vs target"
          ),
          alpha = 0.28, color = NA
        )
    } +
    {
      if (nrow(surplus_df) > 0)
        geom_ribbon(
          data = surplus_df,
          aes(
            x = day, ymin = target, ymax = value,
            group = segment, fill = "Surplus vs target"
          ),
          alpha = 0.28, color = NA
        )
    } +
    geom_line(
      data = plot_df,
      aes(x = day, y = value, color = "Observed value"),
      linewidth = 1.05
    ) +
    geom_line(
      data = plot_df,
      aes(x = day, y = target, color = "Target"),
      linewidth = 0.9,
      linetype = 2
    ) +
    geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dashed") +
    scale_fill_manual(
      name = NULL,
      values = c(
        "Deficit vs target" = "#E67D70",
        "Surplus vs target" = "#7CE6B1"
      )
    ) +
    scale_color_manual(
      name = NULL,
      values = c(
        "Observed value" = "#30343F",
        "Target" = "#8A8F99"
      )
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = y_label
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = element_text(margin = margin(r = 6)),
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      legend.position = ifelse(show_legend, "top", "none"),
      legend.justification = "left",
      panel.border = element_rect(fill = NA, color = "#D0D0D0", linewidth = 0.6)
    )

  p
}

# ----------------------------
# 3. Prepare plotting tables
# ----------------------------
fi_df <- df %>%
  transmute(
    day,
    value = fi,
    target = fi_target,
    diff = fi - fi_target
  )

hdep_df <- df %>%
  transmute(
    day,
    value = hdep,
    target = hdep_target,
    diff = hdep - hdep_target
  )

mort_df <- df %>%
  transmute(
    day,
    value = mort,
    target = mort_target,
    diff = mort - mort_target,
    status = ifelse(diff > 0, "Above target", "At or below target")
  )

# ----------------------------
# 4. Generate plots
# ----------------------------
p_fi <- make_target_line_plot(
  fi_df,
  y_label = "Feed intake (g/hen/day)",
  panel_title = "B  Feed intake vs target",
  show_legend = TRUE
) +
  annotate(
    "label",
    x = 0,
    y = max(fi_df$value, fi_df$target, na.rm = TRUE),
    label = "Outbreak day 0",
    vjust = -0.2,
    label.size = NA,
    size = 3.5
  ) +
  coord_cartesian(clip = "off")

p_hdep <- make_target_line_plot(
  hdep_df,
  y_label = "HDEP (%)",
  panel_title = "C  HDEP vs target",
  show_legend = FALSE
)

p_mort <- ggplot(mort_df, aes(x = day, y = value, fill = status)) +
  geom_col(width = 0.85) +
  geom_line(
    aes(x = day, y = target),
    color = "#8A8F99",
    linetype = 2,
    linewidth = 0.9
  ) +
  geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dashed") +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Above target" = "#E67D70",
      "At or below target" = "#C9D9C5"
    )
  ) +
  labs(
    title = "D  Daily mortality vs target",
    x = "Day relative to outbreak onset",
    y = "Mortality (birds/day)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 6)),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "#D0D0D0", linewidth = 0.6)
  )

# Combine panels.
p_final <- p_fi / p_hdep / p_mort +
  plot_layout(heights = c(1.2, 1.2, 1), guides = "collect") &
  theme(legend.position = "top")

# ----------------------------
# 5. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p_final,
  width = 7.2,
  height = 8.8,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p_final,
  width = 7.2,
  height = 8.8,
  units = "in",
  dpi = 300
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)

# Optional: print session information for reproducibility.
writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo_Fig2_BCD.txt"))
