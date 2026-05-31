```r
#!/usr/bin/env Rscript

# ============================================================
# Figure 3B-D: Environmental diversity and Avibacterium paragallinarum abundance
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation
# in laying hens
#
# Usage:
#   Rscript Fig3_BCD_environmental_statistics.R
#   Rscript Fig3_BCD_environmental_statistics.R data/Fig3_BCD_environmental_metrics.xlsx results/Fig3_BCD_environmental_statistics.pdf
#
# Required input columns:
#   Sample
#   Group
#   Shannon
#   Genus_richness
#   Apg_abundance
#
# Notes:
#   Group names should preferably be:
#   COR, ICin, ICinlet, HEin, HEinlet
#
#   For compatibility with earlier internal files, the following names are
#   automatically converted:
#   HE-IC -> COR
#   IC    -> ICin
#   HE    -> HEin
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggdist)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(scales)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_file <- ifelse(
  length(args) >= 1,
  args[1],
  "data/Fig3_BCD_environmental_metrics.xlsx"
)

output_pdf <- ifelse(
  length(args) >= 2,
  args[2],
  "results/Fig3_BCD_environmental_statistics.pdf"
)

output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Please provide the Excel file path as the first argument, for example:\n",
    "Rscript Fig3_BCD_environmental_statistics.R data/Fig3_BCD_environmental_metrics.xlsx"
  )
}

# ----------------------------
# 1. Read and validate data
# ----------------------------
required_cols <- c(
  "Sample",
  "Group",
  "Shannon",
  "Genus_richness",
  "Apg_abundance"
)

raw_df <- read_excel(input_file)

missing_cols <- setdiff(required_cols, colnames(raw_df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    "\nRequired columns are: ",
    paste(required_cols, collapse = ", ")
  )
}

group_levels <- c("COR", "ICin", "ICinlet", "HEin", "HEinlet")

group_colors_fill <- c(
  "COR"      = "#9EC5DC",
  "ICin"    = "#FFAAAA",
  "ICinlet" = "#C6AFE9",
  "HEin"    = "#AFE9AF",
  "HEinlet" = "#AFE9DD"
)

group_colors_edge <- c(
  "COR"      = "#3A78A1",
  "ICin"    = "#C44F4F",
  "ICinlet" = "#7A5BAE",
  "HEin"    = "#4E9A4E",
  "HEinlet" = "#3B8C87"
)

df <- raw_df %>%
  mutate(
    Sample = as.character(Sample),
    Group = as.character(Group),
    Group = case_when(
      Group == "HE-IC" ~ "COR",
      Group == "IC" ~ "ICin",
      Group == "HE" ~ "HEin",
      TRUE ~ Group
    ),
    Shannon = suppressWarnings(as.numeric(Shannon)),
    Genus_richness = suppressWarnings(as.numeric(Genus_richness)),
    Apg_abundance = suppressWarnings(as.numeric(Apg_abundance))
  ) %>%
  filter(Group %in% group_levels) %>%
  mutate(Group = factor(Group, levels = group_levels))

if (nrow(df) == 0) {
  stop("No valid rows were retained after group filtering.")
}

# ----------------------------
# 2. Convert to long format
# ----------------------------
df_long <- df %>%
  select(Sample, Group, Shannon, Genus_richness, Apg_abundance) %>%
  pivot_longer(
    cols = c(Shannon, Genus_richness, Apg_abundance),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))

# ----------------------------
# 3. Helper function
# ----------------------------
plot_metric <- function(data, metric_name, y_label, panel_title, per_mille_axis = FALSE) {
  plot_df <- data %>%
    filter(Metric == metric_name) %>%
    filter(!is.na(Value), !is.na(Group))

  if (nrow(plot_df) == 0) {
    stop("No data available for metric: ", metric_name)
  }

  # Kruskal-Wallis test.
  kw_tbl <- plot_df %>%
    kruskal_test(Value ~ Group) %>%
    mutate(Metric = metric_name, .before = 1)

  # Dunn's post hoc test with Benjamini-Hochberg adjustment.
  dunn_tbl <- plot_df %>%
    dunn_test(Value ~ Group, p.adjust.method = "BH") %>%
    mutate(Metric = metric_name, .before = 1)

  # Only significant comparisons are shown to avoid overcrowding.
  stat_tbl <- dunn_tbl %>%
    filter(p.adj < 0.05) %>%
    add_xy_position(x = "Group", step.increase = 0.10)

  p <- ggplot(
    plot_df,
    aes(x = Group, y = Value, fill = Group, color = Group)
  ) +
    ggdist::stat_halfeye(
      side = "right",
      width = 0.40,
      adjust = 1.1,
      slab_alpha = 0.65,
      normalize = "groups",
      point_interval = NULL,
      justification = 1.18,
      position = position_nudge(x = 0.10)
    ) +
    geom_boxplot(
      width = 0.22,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.85
    ) +
    geom_jitter(
      width = 0.06,
      size = 2.8,
      alpha = 0.90,
      stroke = 0
    ) +
    scale_fill_manual(values = group_colors_fill, drop = FALSE) +
    scale_color_manual(values = group_colors_edge, drop = FALSE) +
    labs(
      x = NULL,
      y = y_label,
      title = panel_title
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.text.x = element_text(size = 11, face = "bold", angle = 25, hjust = 1),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    ) +
    coord_cartesian(clip = "off")

  if (per_mille_axis) {
    p <- p +
      scale_y_continuous(
        labels = label_number(scale = 1000, accuracy = 0.1, suffix = "\u2030"),
        expand = expansion(mult = c(0.02, 0.18))
      )
  } else {
    p <- p +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.18))
      )
  }

  if (nrow(stat_tbl) > 0) {
    p <- p +
      stat_pvalue_manual(
        data = stat_tbl,
        label = "p.adj.signif",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.65,
        size = 4,
        hide.ns = TRUE,
        inherit.aes = FALSE
      )
  }

  return(list(
    plot = p,
    kruskal = kw_tbl,
    dunn = dunn_tbl
  ))
}

# ----------------------------
# 4. Generate panels
# ----------------------------
res_shannon <- plot_metric(
  data = df_long,
  metric_name = "Shannon",
  y_label = "Shannon diversity",
  panel_title = "B  Shannon diversity",
  per_mille_axis = FALSE
)

res_richness <- plot_metric(
  data = df_long,
  metric_name = "Genus_richness",
  y_label = "Genus richness",
  panel_title = "C  Genus richness",
  per_mille_axis = FALSE
)

res_apg <- plot_metric(
  data = df_long,
  metric_name = "Apg_abundance",
  y_label = expression(italic("Av. paragallinarum") ~ "abundance (" * "\u2030" * ")"),
  panel_title = expression("D  " * italic("Av. paragallinarum") * " abundance"),
  per_mille_axis = TRUE
)

final_plot <- res_shannon$plot + res_richness$plot + res_apg$plot +
  plot_layout(ncol = 3)

# ----------------------------
# 5. Export figure
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = final_plot,
  width = 13.5,
  height = 5.0,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = final_plot,
  width = 13.5,
  height = 5.0,
  units = "in",
  dpi = 300
)

message("Figure saved to: ", output_pdf)
message("Figure saved to: ", output_png)

# ----------------------------
# 6. Export statistics and session information
# ----------------------------
kruskal_results <- bind_rows(
  res_shannon$kruskal,
  res_richness$kruskal,
  res_apg$kruskal
)

dunn_results <- bind_rows(
  res_shannon$dunn,
  res_richness$dunn,
  res_apg$dunn
)

write.csv(
  kruskal_results,
  file = file.path(output_dir, "Fig3_BCD_kruskal_results.csv"),
  row.names = FALSE
)

write.csv(
  dunn_results,
  file = file.path(output_dir, "Fig3_BCD_dunn_BH_results.csv"),
  row.names = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig3_BCD.txt")
)

message("Statistical results and session information saved to: ", output_dir)
```





```r
#!/usr/bin/env Rscript

# ============================================================
# Figure 3E: Environmental microbiome NMDS
# Bray-Curtis NMDS based on Hellinger-transformed abundance data
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation
# in laying hens
#
# Usage:
#   Rscript Fig3_E_environmental_NMDS.R
#   Rscript Fig3_E_environmental_NMDS.R data/env_group.xlsx data/envmicro_otu_S.xlsx results/Fig3_E_environmental_NMDS.pdf
#
# Required input:
#   1. Group file:
#      - The first two columns should be Group and Sample, or should contain
#        columns named "Group" and "Sample".
#
#   2. Abundance file:
#      - The first column should contain taxon/species names.
#      - The remaining columns should be sample abundance columns.
#      - Rows = taxa/species, columns = samples.
#
# Output:
#   - NMDS plot in PDF and PNG formats
#   - NMDS site scores
#   - PERMANOVA result table
#   - sessionInfo record
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(vegan)
  library(ggplot2)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

group_file <- ifelse(length(args) >= 1, args[1], "data/env_group.xlsx")
abund_file <- ifelse(length(args) >= 2, args[2], "data/envmicro_otu_S.xlsx")
output_pdf <- ifelse(length(args) >= 3, args[3], "results/Fig3_E_environmental_NMDS.pdf")
output_png <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(group_file)) {
  stop(
    "Group file not found: ", group_file, "\n",
    "Example:\n",
    "Rscript Fig3_E_environmental_NMDS.R data/env_group.xlsx data/envmicro_otu_S.xlsx"
  )
}

if (!file.exists(abund_file)) {
  stop(
    "Abundance file not found: ", abund_file, "\n",
    "Example:\n",
    "Rscript Fig3_E_environmental_NMDS.R data/env_group.xlsx data/envmicro_otu_S.xlsx"
  )
}

# ----------------------------
# 1. Read abundance table
# ----------------------------
ab_raw <- read_excel(abund_file)

if (ncol(ab_raw) < 3) {
  stop("The abundance table should contain one taxon column and at least two sample columns.")
}

taxon_col <- colnames(ab_raw)[1]

taxa_names <- as.character(ab_raw[[1]])
taxa_names[is.na(taxa_names) | taxa_names == ""] <- paste0("Taxon_", seq_len(sum(is.na(taxa_names) | taxa_names == "")))

ab_tbl <- ab_raw[, -1, drop = FALSE]

ab_num <- as.data.frame(
  lapply(ab_tbl, function(x) suppressWarnings(as.numeric(as.character(x)))),
  check.names = FALSE
)

rownames(ab_num) <- make.unique(taxa_names)
ab_num[is.na(ab_num)] <- 0

# Transpose to samples x taxa.
mat <- t(as.matrix(ab_num))

# Remove all-zero samples and all-zero taxa.
mat <- mat[rowSums(mat) > 0, , drop = FALSE]
mat <- mat[, colSums(mat) > 0, drop = FALSE]

if (nrow(mat) < 3) {
  stop("Fewer than three non-zero samples were retained after filtering.")
}

# ----------------------------
# 2. Read and align group table
# ----------------------------
grp_raw <- read_excel(group_file)

if (!all(c("Group", "Sample") %in% colnames(grp_raw))) {
  if (ncol(grp_raw) < 2) {
    stop("The group file should contain at least two columns: Group and Sample.")
  }
  grp <- grp_raw %>%
    rename(
      Group = 1,
      Sample = 2
    )
} else {
  grp <- grp_raw %>%
    select(Group, Sample)
}

grp <- grp %>%
  mutate(
    Group = as.character(Group),
    Sample = as.character(Sample)
  ) %>%
  filter(!is.na(Group), !is.na(Sample), Group != "", Sample != "")

common_samples <- intersect(rownames(mat), grp$Sample)

if (length(common_samples) < 3) {
  stop(
    "Fewer than three matched samples were found between abundance and group tables.\n",
    "Matched samples: ", paste(common_samples, collapse = ", ")
  )
}

mat <- mat[common_samples, , drop = FALSE]

grp <- grp %>%
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, rownames(mat)))

# Ensure exact row order correspondence between group metadata and abundance matrix.
if (!identical(grp$Sample, rownames(mat))) {
  stop("Sample order mismatch after alignment. Please check sample names.")
}

# User-defined group order used in the manuscript.
preferred_group_order <- c("HE-IC", "IC", "ICinlet", "HE", "HEinlet")

observed_groups <- unique(grp$Group)
final_group_order <- c(
  preferred_group_order[preferred_group_order %in% observed_groups],
  setdiff(observed_groups, preferred_group_order)
)

grp$Group <- factor(grp$Group, levels = final_group_order)

# ----------------------------
# 3. Hellinger transformation, Bray-Curtis distance, NMDS and PERMANOVA
# ----------------------------
mat_hel <- decostand(mat, method = "hellinger")
bray_dist <- vegdist(mat_hel, method = "bray")

set.seed(123)
nmds <- metaMDS(
  bray_dist,
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  wascores = FALSE,
  trace = FALSE
)

adonis_result <- adonis2(
  bray_dist ~ Group,
  data = grp,
  permutations = 999
)

# ----------------------------
# 4. Prepare plotting data
# ----------------------------
site_scores <- as.data.frame(scores(nmds, display = "sites"))
site_scores$Sample <- rownames(site_scores)

plot_df <- site_scores %>%
  left_join(grp, by = "Sample")

if (!all(c("NMDS1", "NMDS2") %in% colnames(plot_df))) {
  colnames(plot_df)[1:2] <- c("NMDS1", "NMDS2")
}

# Convex hull function.
get_hull_closed <- function(df) {
  if (nrow(df) < 3) return(NULL)
  h <- chull(df$NMDS1, df$NMDS2)
  hull <- df[h, c("NMDS1", "NMDS2", "Group")]
  rbind(hull, hull[1, , drop = FALSE])
}

hull_list <- lapply(split(plot_df, plot_df$Group), get_hull_closed)
hull_df <- do.call(rbind, hull_list[!vapply(hull_list, is.null, logical(1))])

# Color palette.
group_colors <- c(
  "HE-IC" = "#9EC5DC",
  "IC" = "#FFAAAA",
  "ICinlet" = "#C6AFE9",
  "HE" = "#AFE9AF",
  "HEinlet" = "#AFE9DD"
)

# If additional groups exist, assign default ggplot colors.
missing_color_groups <- setdiff(levels(grp$Group), names(group_colors))
if (length(missing_color_groups) > 0) {
  additional_colors <- scales::hue_pal()(length(missing_color_groups))
  names(additional_colors) <- missing_color_groups
  group_colors <- c(group_colors, additional_colors)
}

group_colors <- group_colors[levels(grp$Group)]

# Annotation labels.
r2_value <- adonis_result$R2[1]
p_value <- adonis_result$`Pr(>F)`[1]

label_permanova <- paste0(
  "PERMANOVA: R² = ",
  format(round(r2_value, 3), nsmall = 3),
  ", P = ",
  ifelse(p_value < 0.001, "< 0.001", signif(p_value, 3))
)

label_stress <- paste0("Stress = ", round(nmds$stress, 3))

x_range <- range(plot_df$NMDS1, na.rm = TRUE)
y_range <- range(plot_df$NMDS2, na.rm = TRUE)

x_annot <- x_range[1] + 0.02 * diff(x_range)
y_annot1 <- y_range[2] - 0.06 * diff(y_range)
y_annot2 <- y_range[2] - 0.14 * diff(y_range)

# ----------------------------
# 5. Generate NMDS plot
# ----------------------------
p_nmds <- ggplot() +
  {
    if (!is.null(hull_df))
      geom_polygon(
        data = hull_df,
        aes(x = NMDS1, y = NMDS2, group = Group, fill = Group),
        color = NA,
        alpha = 0.18
      )
  } +
  {
    if (!is.null(hull_df))
      geom_polygon(
        data = hull_df,
        aes(x = NMDS1, y = NMDS2, group = Group, color = Group),
        fill = NA,
        linewidth = 0.9,
        linejoin = "mitre"
      )
  } +
  geom_point(
    data = plot_df,
    aes(x = NMDS1, y = NMDS2, color = Group),
    size = 3.2,
    alpha = 0.95
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey60",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    color = "grey60",
    linewidth = 0.4
  ) +
  annotate(
    "text",
    x = x_annot,
    y = y_annot1,
    label = label_permanova,
    hjust = 0,
    size = 4.5,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_annot,
    y = y_annot2,
    label = label_stress,
    hjust = 0,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_color_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Environmental microbiome NMDS",
    x = "NMDS1",
    y = "NMDS2",
    color = "Group",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )

print(p_nmds)

# ----------------------------
# 6. Export figure and tables
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = p_nmds,
  width = 7.0,
  height = 5.6,
  units = "in"
)

ggsave(
  filename = output_png,
  plot = p_nmds,
  width = 7.0,
  height = 5.6,
  units = "in",
  dpi = 300
)

write.csv(
  plot_df,
  file = file.path(output_dir, "Fig3_E_NMDS_site_scores.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(adonis_result),
  file = file.path(output_dir, "Fig3_E_PERMANOVA_result.csv"),
  row.names = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(output_dir, "sessionInfo_Fig3_E_NMDS.txt")
)

message("NMDS figure saved to: ", output_pdf)
message("NMDS figure saved to: ", output_png)
message("NMDS site scores and PERMANOVA results saved to: ", output_dir)
```

