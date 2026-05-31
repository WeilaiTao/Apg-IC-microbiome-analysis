#!/usr/bin/env Rscript
# ============================================================
# Figure 5: Strain-sharing route summary
# Total strain-sharing events and unique strain counts across habitat pairs
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript Fig5_strain_sharing_route_summary.R
#   Rscript Fig5_strain_sharing_route_summary.R data/strain_sharing_routes.csv results/Fig5_strain_sharing_route_summary.pdf
#
# Notes:
#   The input CSV file should contain:
#   Transmission_Route, Event_Count, Unique_Strain_Count
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(length(args) >= 1, args[1], "data/strain_sharing_routes.csv")
output_pdf <- ifelse(length(args) >= 2, args[2], "results/Fig5_strain_sharing_route_summary.pdf")
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
required_cols <- c("Transmission_Route", "Event_Count", "Unique_Strain_Count")

df <- read_csv(input_file, show_col_types = FALSE)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    Transmission_Route = as.character(Transmission_Route),
    Event_Count = suppressWarnings(as.numeric(Event_Count)),
    Unique_Strain_Count = suppressWarnings(as.numeric(Unique_Strain_Count))
  ) %>%
  filter(
    !is.na(Transmission_Route),
    is.finite(Event_Count),
    is.finite(Unique_Strain_Count)
  ) %>%
  arrange(desc(Event_Count))

if (nrow(df) == 0) {
  stop("No valid rows were retained after reading the input file.")
}

df <- df %>%
  mutate(
    Transmission_Route = factor(Transmission_Route, levels = Transmission_Route)
  )

# ----------------------------
# 2. Plotting parameters
# ----------------------------
scale_factor <- max(df$Event_Count, na.rm = TRUE) / max(df$Unique_Strain_Count, na.rm = TRUE) * 0.85

route_colors <- c(
  "ENV_vs_ENV" = "#7EA1C4",
  "GUT_vs_GUT" = "#C06C84",
  "ENV_vs_GUT" = "#6C5B7B",
  "IOS_vs_IOS" = "#E3B489",
  "ENV_vs_IOS" = "#88B04B"
)

used_colors <- route_colors[names(route_colors) %in% as.character(df$Transmission_Route)]
line_color <- "#2C3E50"

# ----------------------------
# 3. Generate plot
# ----------------------------
p <- ggplot(df, aes(x = Transmission_Route)) +
  geom_col(
    aes(y = Event_Count, fill = Transmission_Route),
    width = 0.55,
    alpha = 0.85
  ) +
  geom_line(
    aes(
      y = Unique_Strain_Count * scale_factor,
      group = 1,
      color = "Unique strain count"
    ),
    linewidth = 0.9
  ) +
  geom_point(
    aes(
      y = Unique_Strain_Count * scale_factor,
      color = "Unique strain count"
    ),
    size = 2.8
  ) +
  geom_text(
    aes(y = Event_Count, label = comma(Event_Count)),
    vjust = -0.4,
    size = 3.2
  ) +
  geom_text(
    aes(
      y = Unique_Strain_Count * scale_factor,
      label = Unique_Strain_Count
    ),
    vjust = -1.1,
    size = 3.2,
    color = line_color
  ) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.18)),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Unique strain count"
    )
  ) +
  scale_fill_manual(
    name = "Transmission route",
    values = used_colors
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Unique strain count" = line_color)
  ) +
  labs(
    x = NULL,
    y = "Total sharing events"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.title.y.right = element_text(color = line_color, margin = margin(l = 8)),
    axis.text.y.right = element_text(color = line_color),
    axis.ticks.y.right = element_line(color = line_color),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
    legend.position = "top",
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
  )

# ----------------------------
# 4. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p,
  width = 7.0,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p,
  width = 7.0,
  height = 4.8,
  units = "in",
  dpi = 300
)

write.csv(
  df,
  file = file.path(output_dir, "Fig5_strain_sharing_route_summary_data.csv"),
  row.names = FALSE
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)
message("Processed data saved to: ", output_dir)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig5_strain_sharing_route_summary.txt")
)