library(data.table)
library(dplyr)
library(colorRamp2)
library(ComplexHeatmap)
library(ggplot2)

source("utils.R")
source("import.R")

dir.create("outputs/heatmaps", recursive = T, showWarnings = F)
dir.create("outputs/qc", recursive = T, showWarnings = F)
dir.create("outputs/networks", recursive = T, showWarnings = F)
dir.create("outputs/hockeys", recursive = T, showWarnings = F)

# note: this shows that a single gene can contain multiple idr's
# of the same cluster, which i'm going to collapse into a single
# node. the reason is that i think we're primarily interested in
# the aggregate effect the idr may have on association, and having
# a single gene driving that effect may be a bit misleading.

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

################################################################################

# note, not using tidyverse just to avoid memory issues.
cge_cor_edges <- cge_cor_ranked[
  GIN_GENES,
  .SD, .SDcols = c("gene", GIN_GENES)
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
  gin_cluster_map,
  by.x = "gene_a",
  by.y = "gene",
  all.x = T,
  allow.cartesian = T
)
setnames(cge_cor_edges, "cluster", "cluster_a")
cge_cor_edges <- merge(
  cge_cor_edges,
  gin_cluster_map,
  by.x = "gene_b",
  by.y = "gene",
  all.x = T,
  allow.cartesian = T
)
setnames(cge_cor_edges, "cluster", "cluster_b")

# sanity check.
edges_per_gene <- cge_cor_edges %>%
  group_by(gene_a, cluster_a) %>%
  summarize(
    n_connections = n(),
  ) %>%
  ungroup %>%
  arrange(desc(n_connections)) %>%
  rename(gene = gene_a, cluster = cluster_a)

cluster_sizes <- gin_cluster_map %>%
  group_by(cluster) %>%
  summarize(size = n()) %>%
  pull(size, cluster)

# surprisingly log
cge_cor_cluster_edges <- cge_cor_edges %>%
  group_by(cluster_a, cluster_b) %>%
  summarize(
    n_edges = n()
  ) %>%
  ungroup %>%
  mutate(
    edges_normed = ifelse(
      cluster_a == cluster_b,
      n_edges / (cluster_sizes[as.character(cluster_a)] * (cluster_sizes[as.character(cluster_b)] - 1) / 1000 ),
      n_edges / (cluster_sizes[as.character(cluster_a)] * cluster_sizes[as.character(cluster_b)] / 1000 )
    ),
    log_edges_normed = log10(edges_normed)
  )

################################################################################
# rank plot

expected_edges_per_thousand <- nrow(cge_cor_edges) / (
  nrow(gin_cluster_map) * (nrow(gin_cluster_map) - 1) / 1000
)

cluster_cluster_edges_ranked <- cge_cor_cluster_edges %>%
  filter(cluster_a <= cluster_b) %>%
  mutate(rank = row_number(desc(edges_normed)))

ggplot(
  cluster_cluster_edges_ranked,
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
  geom_point(
    data = cluster_cluster_edges_ranked %>%
      filter(!cluster_a == 11 | !cluster_b == 11) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]),
    color = "black",
    shape = 16,
    size = 2,
    alpha = 0.8
  ) +
  geom_point(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == 11 & cluster_b == 11) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]),
    aes(
      color = color
    ),
    shape = 16,
    size = 4,
    alpha = 0.8
  ) +
  ggrepel::geom_text_repel(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == 11 & cluster_b == 11) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]),
    aes(
      label = paste0(
        rank,
        ". Cluster ",
        cluster_a,
        " - Cluster ",
        cluster_b
      ),
      color = color,
    ),
    max.overlaps = Inf,
    min.segment.length = 0,
    nudge_x = 100,
  ) +
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(-.1, 9, 0.1),
    labels = round(10 ^ seq(-.1, 9, 0.1), 2),
    limits = range(cluster_cluster_edges_ranked$log_edges_normed)
  ) +
  scale_x_continuous(
    name = "Cluster Fitness Similarity Rank",
    breaks = c(
      1, max(cluster_cluster_edges_ranked$rank)
    )
  ) +
  coord_fixed(
    ratio = diff(range(cluster_cluster_edges_ranked$rank)) / 
      diff(range(cluster_cluster_edges_ranked$log_edges_normed)) / 3
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/cluster_fitness_similarity_rank_",
    MAX_RANK,
    "_self_similarity_highlight_log_grant-version.pdf"
  ),
  height = 6,
  width  = 15,
  device = cairo_pdf
)
