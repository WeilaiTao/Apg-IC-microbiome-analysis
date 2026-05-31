#!/usr/bin/env Rscript

# ============================================================

# Figure 4D: Shannon diversity of infraorbital sinus microbiomes

# Genus-level Shannon index between HE and IC groups

# Manuscript: Ventilation-shaped farm environments link infectious coryza dissemination,

# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens

# Usage:

# Rscript Fig4D_IOS_Shannon_raincloud.R

# Rscript Fig4D_IOS_Shannon_raincloud.R data/combined_otu_G.xlsx data/GroupID.xlsx results/Fig4D_IOS_Shannon_raincloud.pdf

# Notes:

# Genus abundance table: Genus, sample_1, sample_2, ...

# Metadata table: Sample, Organ, Group1

# This script compares HE and IC samples from IOS only.

# ============================================================

suppressPackageStartupMessages({
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggdist)
library(rstatix)
library(ggpubr)
})

# ----------------------------

# 0. Input and output settings

# ----------------------------

args <- commandArgs(trailingOnly = TRUE)

abundance_file <- ifelse(length(args) >= 1, args[1], "data/combined_otu_G.xlsx")
metadata_file  <- ifelse(length(args) >= 2, args[2], "data/GroupID.xlsx")
output_pdf     <- ifelse(length(args) >= 3, args[3], "results/Fig4D_IOS_Shannon_raincloud.pdf")
output_png     <- sub("\.pdf$", ".png", output_pdf)

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

target_organ  <- "IOS"
target_groups <- c("HE", "IC")

# ----------------------------

# 1. Read and validate data

# ----------------------------

abundance_raw <- read_excel(abundance_file)

if (ncol(abundance_raw) < 2) {
stop("The abundance table should contain one genus column and at least one sample column.")
}

colnames(abundance_raw)[1] <- "Genus"

abundance_df <- abundance_raw %>%
mutate(Genus = as.character(Genus)) %>%
filter(!is.na(Genus), Genus != "") %>%
mutate(across(-Genus, ~ suppressWarnings(as.numeric(.x))))

metadata_raw <- read_excel(metadata_file) %>%
mutate(across(everything(), as.character))

required_meta_cols <- c("Sample", "Organ", "Group1")
missing_meta_cols <- setdiff(required_meta_cols, colnames(metadata_raw))

if (length(missing_meta_cols) > 0) {
stop("Missing required metadata columns: ", paste(missing_meta_cols, collapse = ", "))
}

metadata_df <- metadata_raw %>%
filter(Organ == target_organ, Group1 %in% target_groups) %>%
transmute(Sample = Sample, Organ = Organ, Group = Group1)

sample_cols <- setdiff(colnames(abundance_df), "Genus")
matched_samples <- intersect(sample_cols, metadata_df$Sample)

if (length(matched_samples) < 2) {
stop("Fewer than two matched samples were found between abundance table and metadata.")
}

abundance_df <- abundance_df %>%
select(Genus, all_of(matched_samples))

# ----------------------------

# 2. Calculate Shannon index

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
column_to_rownames("Genus") %>%
as.data.frame()

shannon_values <- apply(mat, 2, calc_shannon)

plot_df <- tibble(
Sample = names(shannon_values),
Shannon = as.numeric(shannon_values)
) %>%
inner_join(metadata_df, by = "Sample") %>%
mutate(
Group = factor(Group, levels = target_groups),
Organ = factor(Organ, levels = target_organ)
) %>%
arrange(Group, Sample)

if (nrow(plot_df) == 0) {
stop("No valid samples were retained for plotting.")
}

if (any(table(plot_df$Group) == 0)) {
stop("At least one target group has no matched samples.")
}

# ----------------------------

# 3. Statistical analysis

# ----------------------------

kruskal_result <- kruskal_test(plot_df, Shannon ~ Group)

wilcox_result <- pairwise_wilcox_test(
plot_df,
Shannon ~ Group,
p.adjust.method = "BH"
) %>%
mutate(
p.signif = case_when(
p.adj < 0.001 ~ "***",
p.adj < 0.01  ~ "**",
p.adj < 0.05  ~ "*",
TRUE          ~ "ns"
)
)

y_max <- max(plot_df$Shannon, na.rm = TRUE)
y_min <- min(plot_df$Shannon, na.rm = TRUE)
y_range <- y_max - y_min

if (is.na(y_range) || y_range == 0) {
y_range <- 1
}

if (nrow(wilcox_result) > 0) {
wilcox_result <- wilcox_result %>%
mutate(y.position = y_max + 0.12 * y_range)
}

# ----------------------------

# 4. Generate plot

# ----------------------------

group_colors <- c("HE" = "#AFC6D6", "IC" = "#D9B9AE")
point_colors <- c("HE" = "#5E95B1", "IC" = "#AC7765")

set.seed(123)

p <- ggplot(plot_df, aes(x = Group, y = Shannon, fill = Group, color = Group)) +
stat_halfeye(
adjust = 0.7,
width = 0.35,
justification = -0.2,
.width = 0,
point_colour = NA,
alpha = 0.65
) +
geom_boxplot(
width = 0.12,
outlier.shape = NA,
linewidth = 0.5,
alpha = 0.45
) +
geom_jitter(
aes(x = as.numeric(Group) - 0.25),
width = 0.08,
size = 2.3,
alpha = 0.85
) +
scale_fill_manual(values = group_colors) +
scale_color_manual(values = point_colors) +
labs(x = NULL, y = "Genus Shannon index") +
scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
facet_wrap(~ Organ) +
theme_classic() +
theme(
legend.position = "none",
strip.background = element_rect(fill = "grey90", color = "grey50"),
panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
)

if (nrow(wilcox_result) > 0) {
p <- p +
stat_pvalue_manual(
wilcox_result,
label = "p.signif",
xmin = "group1",
xmax = "group2",
y.position = "y.position",
tip.length = 0.01,
bracket.size = 0.5,
size = 3.5,
inherit.aes = FALSE
)
}

# ----------------------------

# 5. Export

# ----------------------------

ggsave(
filename = output_pdf,
plot = p,
width = 4.0,
height = 4.2,
units = "in"
)








#!/usr/bin/env Rscript
# ============================================================
# Figure 4E: NMDS of infraorbital sinus microbiomes
# Bray-Curtis NMDS with group-wise axis distributions
#
# Manuscript:
# Ventilation-shaped farm environments link infectious coryza dissemination,
# infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens
#
# Usage:
#   Rscript Fig4E_IOS_NMDS_group_distribution.R
#   Rscript Fig4E_IOS_NMDS_group_distribution.R data/SI_Apg_otu_S.xlsx data/GroupID.xlsx results/Fig4E_IOS_NMDS_group_distribution.pdf
#
# Notes:
#   Abundance table: first column = feature/species names; remaining columns = sample IDs
#   Metadata table: Group, Sample
#   Groups retained: Ctrl, EHA, MLD, MOD, SEV
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(vegan)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# 0. Input and output settings
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

matrix_file <- ifelse(length(args) >= 1, args[1], "data/SI_Apg_otu_S.xlsx")
group_file  <- ifelse(length(args) >= 2, args[2], "data/GroupID.xlsx")
output_pdf  <- ifelse(length(args) >= 3, args[3], "results/Fig4E_IOS_NMDS_group_distribution.pdf")
output_png  <- sub("\\.pdf$", ".png", output_pdf)

output_dir <- dirname(output_pdf)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(matrix_file)) {
  stop("Abundance matrix file not found: ", matrix_file)
}

if (!file.exists(group_file)) {
  stop("Group file not found: ", group_file)
}

# ----------------------------
# 1. Read and validate data
# ----------------------------
mat_raw <- read_excel(matrix_file, col_names = TRUE)

if (ncol(mat_raw) < 2) {
  stop("The abundance table should contain one feature column and at least one sample column.")
}

feature_names <- as.character(mat_raw[[1]])
mat_df <- as.data.frame(mat_raw[, -1, drop = FALSE])
colnames(mat_df) <- trimws(colnames(mat_df))
rownames(mat_df) <- make.unique(feature_names)

mat_t <- t(as.matrix(mat_df))
mode(mat_t) <- "numeric"
mat_t[is.na(mat_t)] <- 0

all_samples <- rownames(mat_t)

grp_raw <- read_excel(group_file, col_names = TRUE)

if (ncol(grp_raw) < 2) {
  stop("The metadata table should contain at least two columns: Group and Sample.")
}

colnames(grp_raw)[1:2] <- c("Group", "Sample")

grp <- grp_raw %>%
  transmute(
    Group = as.character(Group),
    Sample = as.character(Sample)
  )

group_levels <- c("Ctrl", "EHA", "MLD", "MOD", "SEV")

grp <- grp %>%
  filter(Group %in% group_levels)

keep_samples <- all_samples[all_samples %in% grp$Sample]

if (length(keep_samples) < 2) {
  stop("Fewer than two matched samples were found between abundance table and metadata.")
}

final_matrix <- mat_t[keep_samples, , drop = FALSE]

final_group_info <- data.frame(
  Sample = keep_samples,
  Group = factor(grp$Group[match(keep_samples, grp$Sample)], levels = group_levels),
  stringsAsFactors = FALSE
)

final_group_info <- final_group_info %>%
  filter(!is.na(Group))

final_matrix <- final_matrix[final_group_info$Sample, , drop = FALSE]

if (nrow(final_matrix) < 2) {
  stop("Fewer than two valid samples were retained after group filtering.")
}

# ----------------------------
# 2. Bray-Curtis distance, NMDS, and PERMANOVA
# ----------------------------
set.seed(123)

bc_dist <- vegdist(final_matrix, method = "bray")
nmds <- metaMDS(
  bc_dist,
  k = 2,
  trymax = 50,
  autotransform = FALSE,
  trace = FALSE
)

nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
nmds_scores$Sample <- rownames(nmds_scores)

plot_df <- nmds_scores %>%
  inner_join(final_group_info, by = "Sample") %>%
  mutate(Group = droplevels(Group))

adonis_result <- adonis2(bc_dist ~ Group, data = plot_df)

# ----------------------------
# 3. Plotting parameters
# ----------------------------
group_colors <- c(
  "Ctrl" = "#1F78B4",
  "EHA"  = "#4CAF50",
  "MLD"  = "#FFC107",
  "MOD"  = "#FF9800",
  "SEV"  = "#F44336"
)

used_colors <- group_colors[names(group_colors) %in% levels(plot_df$Group)]

stat_label <- paste0(
  "R2 = ", round(adonis_result$R2[1], 3),
  ", P = ", formatC(adonis_result$`Pr(>F)`[1], format = "e", digits = 2)
)

stress_label <- paste0("Stress = ", round(nmds$stress, 3))

label_x <- min(plot_df$NMDS1) + 0.05 * diff(range(plot_df$NMDS1))
label_y1 <- max(plot_df$NMDS2) - 0.05 * diff(range(plot_df$NMDS2))
label_y2 <- max(plot_df$NMDS2) - 0.12 * diff(range(plot_df$NMDS2))

# ----------------------------
# 4. Generate plots
# ----------------------------
p_nmds <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Group)) +
  stat_ellipse(
    aes(group = Group, fill = Group),
    type = "norm",
    geom = "polygon",
    alpha = 0.25,
    color = NA
  ) +
  geom_point(size = 2.6, alpha = 0.9) +
  scale_color_manual(values = used_colors, drop = TRUE) +
  scale_fill_manual(values = used_colors, drop = TRUE) +
  labs(
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    color = "Group",
    fill = "Group"
  ) +
  annotate(
    "text",
    x = label_x,
    y = label_y1,
    label = stat_label,
    size = 3.4,
    hjust = 0
  ) +
  annotate(
    "text",
    x = label_x,
    y = label_y2,
    label = stress_label,
    size = 3.4,
    hjust = 0
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey50", linewidth = 0.5)
  )

p_nmds1 <- ggplot(plot_df, aes(x = Group, y = NMDS1, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.28, alpha = 0.65, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85) +
  scale_fill_manual(values = used_colors, drop = TRUE) +
  scale_color_manual(values = used_colors, drop = TRUE) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

p_nmds2 <- ggplot(plot_df, aes(x = Group, y = NMDS2, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.28, alpha = 0.65, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.85) +
  scale_fill_manual(values = used_colors, drop = TRUE) +
  scale_color_manual(values = used_colors, drop = TRUE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

combined_plot <- ((p_nmds1 + plot_spacer()) / (p_nmds + p_nmds2)) +
  plot_layout(widths = c(0.5, 4.5), heights = c(0.6, 4.4))

# ----------------------------
# 5. Export
# ----------------------------
ggsave(
  filename = output_pdf,
  plot = combined_plot,
  width = 7.0,
  height = 6.2,
  units = "in"
)
