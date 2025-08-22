library(data.table)
library(dplyr)
library(colorRamp2)
library(ComplexHeatmap)
library(ggplot2)

source("utils.R")
source("import.R")

################################################################################

dir.create("outputs/heatmaps", recursive = T, showWarnings = F)
dir.create("outputs/qc", recursive = T, showWarnings = F)
dir.create("outputs/networks", recursive = T, showWarnings = F)
dir.create("outputs/hockeys", recursive = T, showWarnings = F)

# note: this shows that a single gene can contain multiple idr's
# of the same cluster, which i'm going to collapse into a single
# node. the reason is that i think we're primarily interested in
# the aggregate effect the idr may have on association, and having
# a single gene driving that effect may be a bit misleading.

# old
# temp
hpa_location_map_old <- fread("~/code/misc_reference_data/proteinatlas.tsv") %>%
  select(
    `Gene`,
    `Subcellular location`
  ) %>%
  rename(
    gene = `Gene`,
    location = `Subcellular location`
  ) %>%
  tidyr::separate_longer_delim(
    location,
    ","
  ) %>%
  filter(
    location != "" & location %in% names(location_colors)
  ) %>%
  filter(gene %in% colnames(cge_cor_ranked))

hpa_location_map <- hpa %>%
  select(
    `Gene name`,
    `Main location`
  ) %>%
  rename(
    gene = `Gene name`,
    location = `Main location`
  ) %>%
  tidyr::separate_longer_delim(
    location,
    ";"
  ) %>%
  filter(
    location != "" & location %in% names(location_colors)
  ) %>%
  filter(gene %in% colnames(cge_cor_ranked))

gin_cluster_map_redundant <- gin %>%
  rename(
    cluster = `Cluster Number`,
    gene = hgnc_symbol
  ) %>%
  select(
    gene, idr, cluster
  ) %>%
  filter(gene %in% colnames(cge_cor_ranked)) %>%
  data.table(key = "gene")

gin_cluster_map <- gin %>%
  rename(
    cluster = `Cluster Number`,
    gene = hgnc_symbol
  ) %>%
  select(
    gene, cluster
  ) %>%
  distinct %>%
  filter(gene %in% colnames(cge_cor_ranked)) %>%
  data.table(key = "gene")

################################################################################

MAX_RANK <- 20
GIN_GENES <- unique(gin_cluster_map$gene)
HPA_GENES <- unique(hpa_location_map$gene)

################################################################################

# note, not using tidyverse just to avoid memory issues.
cge_cor_edges <- cge_cor_ranked[
  HPA_GENES,
  .SD, .SDcols = c("gene", HPA_GENES)
] %>%
  melt(
    id.vars = "gene",
    variable.name = "gene_b",
    value.name = "rank"
  )
setnames(cge_cor_edges, "gene", "gene_a")
cge_cor_edges[, gene_b := as.character(gene_b)]
cge_cor_edges[, rank := rank - 1]
cge_cor_edges <- cge_cor_edges[rank > 0 & rank <= MAX_RANK,]
cge_cor_edges <- merge(
  cge_cor_edges,
  hpa_location_map,
  by.x = "gene_a",
  by.y = "gene",
  all.x = T,
  allow.cartesian = T
)
setnames(cge_cor_edges, "location", "location_a")
cge_cor_edges <- merge(
  cge_cor_edges,
  hpa_location_map,
  by.x = "gene_b",
  by.y = "gene",
  all.x = T,
  allow.cartesian = T
)
setnames(cge_cor_edges, "location", "location_b")

# sanity check.
edges_per_gene <- cge_cor_edges %>%
  group_by(gene_a, location_a) %>%
  summarize(
    n_connections = n(),
  ) %>%
  ungroup %>%
  arrange(desc(n_connections)) %>%
  rename(gene = gene_a, location = location_a)

location_sizes <- hpa_location_map %>%
  group_by(location) %>%
  summarize(size = n()) %>%
  pull(size, location)

# temp
# location_sizes_old <- hpa_location_map_old %>%
#   group_by(location) %>%
#   summarize(size = n()) %>%
#   pull(size, location)
# 
# hpa_location_map_nuc <- hpa_location_map %>%
#   filter(location == "Nucleoli") %>%
#   pull(gene)
# 
# hpa_location_map_nuc_old <- hpa_location_map_old %>%
#   filter(location == "Nucleoli") %>%
#   pull(gene)
# 
# setdiff(hpa_location_map_nuc, hpa_location_map_nuc_old)
# setdiff(hpa_location_map_nuc_old, hpa_location_map_nuc)

# surprisingly log
cge_cor_location_edges <- cge_cor_edges %>%
  group_by(location_a, location_b) %>%
  summarize(
    n_edges = n()
  ) %>%
  ungroup %>%
  mutate(
    edges_normed = ifelse(
      location_a == location_b,
      n_edges / (location_sizes[as.character(location_a)] * (location_sizes[as.character(location_b)] - 1) / 1000 ),
      n_edges / (location_sizes[as.character(location_a)] * location_sizes[as.character(location_b)] / 1000 )
    )
  )

cairo_pdf(
  paste0(
    "outputs/qc/qq_plot_rank_", MAX_RANK, "_edges_hpa.pdf"
  ),
  height = 6,
  width = 6
)
cge_cor_location_edges %>%
  filter(location_a <= location_b) %>% # ignore double counting.
  pull(edges_normed) %>%
  log10 %>%
  {
    qqnorm(.)
    qqline(.)
  }
dev.off()

cge_cor_location_edges_mat <- cge_cor_location_edges %>%
  reshape2::dcast(
    location_a ~ location_b,
    value.var = "edges_normed",
    fill = 0
  ) %>%
  data.table %>%
  as.matrix(rownames = "location_a")

################################################################################

edge_per_thousand_breaks <- cge_cor_location_edges %>%
  filter(location_a <= location_b) %>% # ignore double counting.
  pull(edges_normed) %>%
  log10 %>%
  quantile(
    .,
    pnorm(c(0, 1.8))
  ) %>%
  {
    seq(from = round(2 * .[1] - .[2], 2), to = round(.[2], 2), length.out = 5)
  } %>%
  unname

exp_color_function_bi <- colorRamp2(
  edge_per_thousand_breaks,
  c("royalblue4", "dodgerblue3", "white", "indianred1", "firebrick"),
)

cairo_pdf(
  paste0(
    "outputs/heatmaps/hpa_cluster_rank_",
    MAX_RANK,
    "_edge_matrix_bicolor.pdf"
  ),
  height = 8,
  width = 12
)
Heatmap(
  log10(cge_cor_location_edges_mat),
  name = "log10(Edges per 1000 Gene Pairs)",
  #heatmap_legend_param = list(
  #  at = edge_per_thousand_breaks,
  #  labels = round(10 ^ edge_per_thousand_breaks)
  #),
  col = exp_color_function_bi,
  clustering_distance_columns = "spearman",
  clustering_distance_rows    = "spearman",
  height = unit(4, "inches"),
  width  = unit(4, "inches"),
  column_title_rot = 90,
  show_column_names = F,
  top_annotation = columnAnnotation(
    `Subcellular Location` = ifelse(
      colnames(cge_cor_location_edges_mat) %in% names(location_colors),
      colnames(cge_cor_location_edges_mat),
      "NA"
    ) %>%
      factor(
        levels = names(location_colors)
      ),
    col = list(
      `Subcellular Location` = location_colors
    ),
    gp = gpar(
      lwd = 0.5,
      col = "black"
    )
  ),
  bottom_annotation = columnAnnotation(
    labels = anno_text(
      x = colnames(cge_cor_location_edges_mat),
      rot = 90,
      gp = gpar(
        font = "arial",
        fontsize = 10,
        col = "black"
      ),
    )
  ),
  show_row_names = F,
  right_annotation = rowAnnotation(
    labels = anno_text(
      x =  rownames(cge_cor_location_edges_mat),
      which = "row",
      gp = gpar(
        font = "arial",
        fontsize = 10,
        col = "black"
      ),
    )
  ),
  left_annotation = rowAnnotation(
    `Subcellular Location` = ifelse(
      row.names(cge_cor_location_edges_mat) %in% names(location_colors),
      row.names(cge_cor_location_edges_mat),
      "NA"
    ) %>%
      factor(
        levels = names(location_colors)
      ),
    col = list(
      `Subcellular Location` = location_colors
    ),
    gp = gpar(
      lwd = 0.5,
      col = "black"
    )
  )
) %>%
  print
dev.off()

################################################################################
# rank plot

# 1.917991 edges per thousand gene pairs
# 180222 edges total
expected_edges_per_thousand <- nrow(cge_cor_edges) / (
  nrow(hpa_location_map) * (nrow(hpa_location_map) - 1) / 1000
)

location_location_edges_ranked <- cge_cor_location_edges %>%
  filter(location_a <= location_b) %>%
  mutate(
    rank = row_number(desc(edges_normed)),
    log_edges_normed = log10(edges_normed)
  )

fwrite(
  location_location_edges_ranked,
  "outputs/exports/location_location_edges_ranked.csv"
)

ggplot(
  location_location_edges_ranked,
  aes(
    x = rank,
    y = edges_normed
  )
) +
  geom_hline(
    yintercept = expected_edges_per_thousand,
    linetype = 2,
    linewidth = 1.2,
    color = "grey20"
  ) +
  geom_line() +
  geom_point(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      mutate(color = location_colors[location_a]),
    aes(
      color = color
    ),
    shape = 16,
    size = 4,
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      mutate(color = "black"),
    aes(
      label = paste0(
        rank,
        ". ",
        location_a,
        " - ",
        location_b
      ),
      color = color,
    ),
    max.overlaps = Inf,
    min.segment.length = 0,
    nudge_x = 10,
    nudge_y = 0.1
  ) +
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(1, 7, 1),
    limits = c(1, 7)
  ) +
  scale_x_continuous(
    name = "HPA Location Fitness Similarity Rank",
    breaks = c(
      1, max(location_location_edges_ranked$rank)
    )
  ) +
  coord_fixed(
    ratio = diff(range(location_location_edges_ranked$rank)) /
      diff(c(1, 7))
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/location_fitness_similarity_rank_",
    MAX_RANK,
    "_self_similarity_highlight.pdf"
  ),
  height = 8,
  width  = 8,
  device = cairo_pdf
)

ggplot(
  location_location_edges_ranked,
  aes(
    x = rank,
    y = log_edges_normed
  )
) +
  geom_hline(
    yintercept = log10(expected_edges_per_thousand),
    linetype = 2,
    linewidth = 1.2,
    color = "grey20"
  ) +
  geom_line() +
  geom_point(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      mutate(color = location_colors[location_a]),
    aes(
      color = color
    ),
    shape = 16,
    size = 4,
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      mutate(color = "black"),
    aes(
      label = paste0(
        rank,
        ". ",
        location_a,
        " - ",
        location_b
      ),
      color = color,
    ),
    max.overlaps = Inf,
    min.segment.length = 0,
    nudge_x = 10,
  ) +
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Log Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(0, 0.9, 0.1),
    limits = c(0, 0.9)
  ) +
  scale_x_continuous(
    name = "HPA Location Fitness Similarity Rank",
    breaks = c(
      1, max(location_location_edges_ranked$rank)
    )
  ) +
  coord_fixed(
    ratio = diff(range(location_location_edges_ranked$rank)) /
      diff(c(0, 0.9))
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/location_fitness_similarity_rank_",
    MAX_RANK,
    "_self_similarity_highlight_log.pdf"
  ),
  height = 8,
  width  = 8,
  device = cairo_pdf
)
