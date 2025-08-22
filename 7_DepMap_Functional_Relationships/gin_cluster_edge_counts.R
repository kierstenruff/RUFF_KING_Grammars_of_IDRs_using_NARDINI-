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
    )
  )

cairo_pdf(
  paste0(
    "outputs/qc/qq_plot_rank_", MAX_RANK, "_edges.pdf"
  ),
  height = 6,
  width = 6
)
cge_cor_cluster_edges %>%
  filter(cluster_a <= cluster_b) %>% # ignore double counting.
  pull(edges_normed) %>%
  log10 %>%
  {
    qqnorm(.)
    qqline(.)
  }
dev.off()

cge_cor_cluster_edges_mat <- cge_cor_cluster_edges %>%
  reshape2::dcast(
    cluster_a ~ cluster_b,
    value.var = "edges_normed",
    fill = 0
  ) %>%
  data.table %>%
  as.matrix(rownames = "cluster_a")

################################################################################

edge_per_thousand_breaks <- cge_cor_cluster_edges %>%
  filter(cluster_a <= cluster_b) %>% # ignore double counting.
  pull(edges_normed) %>%
  log10 %>%
  quantile(
    pnorm(c(-2.7, 2.7))
  ) %>%
  {
    seq(from = round(.[1], 1), to = round(.[2], 1), length.out = 5)
  } %>%
  unname

exp_color_function_bi <- colorRamp2(
  edge_per_thousand_breaks,
  c("royalblue4", "dodgerblue3", "white", "indianred1", "firebrick"),
)

cairo_pdf(
  paste0(
    "outputs/heatmaps/gin_cluster_rank_",
    MAX_RANK,
    "_edge_matrix_bicolor.pdf"
  ),
  height = 8,
  width = 12
)
Heatmap(
  log10(cge_cor_cluster_edges_mat),
  name = "log10(Edges per 1000 Gene Pairs)",
  #heatmap_legend_param = list(
  #  at = edge_per_thousand_breaks,
  #  labels = round(10 ^ edge_per_thousand_breaks)
  #),
  col = exp_color_function_bi,
  clustering_distance_columns = "euclidean",
  clustering_distance_rows    = "euclidean",
  height = unit(4, "inches"),
  width  = unit(4, "inches"),
  column_title_rot = 90,
  show_column_names = F,
  top_annotation = columnAnnotation(
    `Subcellular Location` = ifelse(
      colnames(cge_cor_cluster_edges_mat) %in% names(locations),
      locations[colnames(cge_cor_cluster_edges_mat)],
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
      x = paste0(
        "Cluster ",
        colnames(cge_cor_cluster_edges_mat)
      ),
      rot = 90,
      gp = gpar(
        font = "arial",
        fontsize = 10,
        col = COLORS[GROUPS[colnames(cge_cor_cluster_edges_mat)]]
      ),
    )
  ),
  show_row_names = F,
  right_annotation = rowAnnotation(
    labels = anno_text(
      x = paste0(
        "Cluster ",
        rownames(cge_cor_cluster_edges_mat)
      ),
      which = "row",
      gp = gpar(
        font = "arial",
        fontsize = 10,
        col = COLORS[GROUPS[rownames(cge_cor_cluster_edges_mat)]]
      ),
    )
  ),
  left_annotation = rowAnnotation(
    `Subcellular Location` = ifelse(
      row.names(cge_cor_cluster_edges_mat) %in% names(locations),
      locations[row.names(cge_cor_cluster_edges_mat)],
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
# rank plot, version 1 no annotations

# 2.043575 edges per thousand
# 105014 edges
expected_edges_per_thousand <- nrow(cge_cor_edges) / (
  nrow(gin_cluster_map) * (nrow(gin_cluster_map) - 1) / 1000
)

cluster_cluster_edges_ranked <- cge_cor_cluster_edges %>%
  filter(cluster_a <= cluster_b) %>%
  mutate(
    rank = row_number(desc(edges_normed)),
    log_edges_normed = log10(edges_normed)
  )

ggplot(
  cluster_cluster_edges_ranked,
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
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
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
      filter(cluster_a == cluster_b) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]) %>%
      filter(rank >= 250),
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
    nudge_x = -100,
    nudge_y = -0.5
  ) +
  ggrepel::geom_text_repel(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]) %>%
      filter(rank < 250),
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
    nudge_y = 0.5
  ) +
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(0, 8, 2),
    limits = c(0, 8)
  ) +
  scale_x_continuous(
    name = "Cluster Fitness Similarity Rank",
    breaks = c(
      1, max(cluster_cluster_edges_ranked$rank)
    )
  ) +
  coord_fixed(
    ratio = diff(range(cluster_cluster_edges_ranked$rank)) / 8
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/cluster_fitness_similarity_rank_",
    MAX_RANK,
    "_self_similarity_highlight.pdf"
  ),
  height = 8,
  width  = 8,
  device = cairo_pdf
)

################################################################################
# rank plot, version 2 subcellular location annotations

inter_location_edges <- fread("outputs/exports/location_location_edges_ranked.csv") %>%
  filter(location_a == location_b) %>%
  mutate(
    label = paste(
      location_a,
      " - ",
      location_b
    )
  )

ggplot(
  cluster_cluster_edges_ranked,
  aes(
    x = rank,
    y = edges_normed
  )
) +
  geom_segment(
    aes(
      y    = expected_edges_per_thousand,
      yend = expected_edges_per_thousand,
      x = -50,
      xend = max(cluster_cluster_edges_ranked$rank) + 50,
      color = location_colors[location_a]
    ),
    linetype = 2,
    linewidth = 1,
    color = "grey20",
  ) +
  geom_segment(
    data = inter_location_edges,
    aes(
      y = edges_normed,
      yend = edges_normed,
      x = -50,
      xend = max(cluster_cluster_edges_ranked$rank) + 50,
      color = location_colors[location_a]
    ),
    linetype = 2,
    linewidth = 0.8,
  ) +
  geom_line() +
  geom_point(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]),
    aes(
      color = color
    ),
    shape = 16,
    size = 4,
    alpha = 0.8
  ) +
  # start cluster cluster label
  geom_segment(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(
        section = case_when(
          cluster_a == 26 ~ 'a',
          edges_normed > 2.2 ~ 'b',
          .default = 'c'
        ),
        xend = case_when(
          section == 'a' ~ 50,
          section == 'b' ~ 250,
          section == 'c' ~ 150,
        ),
        color = COLORS[GROUPS[as.character(cluster_a)]],
      ) %>%
      group_by(section) %>%
      mutate(
        section_rn = row_number(rank),
        yend = case_when(
          section == 'a' ~ 7.8,
          section == 'b' ~ 7.5 - section_rn * 0.15,
          section == 'c' ~ 1.7 - section_rn * 0.15,
        )
      ),
    aes(
      xend = xend,
      yend = yend,
      color = color
    )
  ) +
  geom_text(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(
        section = case_when(
          cluster_a == 26 ~ 'a',
          edges_normed > 2.2 ~ 'b',
          .default = 'c'
        ),
        x = case_when(
          section == 'a' ~ 50 + 5,
          section == 'b' ~ 250 + 5,
          section == 'c' ~ 150 - 5,
        ),
        color = COLORS[GROUPS[as.character(cluster_a)]],
        label = paste0(
          rank,
          ". Cluster ",
          cluster_a,
          " - Cluster ",
          cluster_b
        )
      ) %>%
      group_by(section) %>%
      mutate(
        section_rn = row_number(rank),
        y = case_when(
          section == 'a' ~ 7.8,
          section == 'b' ~ 7.5 - section_rn * 0.16,
          section == 'c' ~ 1.7 - section_rn * 0.16,
        ),
        hjust = ifelse(
          section %in% c('a', 'b'),
          0,
          1
        )
      ),
    aes(
      label = label,
      x = x,
      y = y,
      hjust = hjust,
      color = color
    ),
    family = "Arial",
    size = 3.5
  ) +
  # end cluster cluster label
  # start location location label
  geom_segment(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      select(
        location_a, edges_normed
      ) %>%
      rbind(
        data.frame(
          location_a = "Overall edges per 1000 gene pairs",
          edges_normed = expected_edges_per_thousand
        )
      ) %>%
      mutate(
        section_rn = row_number(desc(edges_normed)),
        xend = max(cluster_cluster_edges_ranked$rank) + 70,
        x = max(cluster_cluster_edges_ranked$rank) + 50,
        y = edges_normed,
        yend = 4.2 - 0.25 * section_rn,
        color = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          "grey20",
          location_colors[location_a]
        ),
        label = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          paste0(location_a, " (", edges_normed, ")"),
          paste0(location_a, " - ", location_a, " (", edges_normed, ")")
        )
      ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = color
    )
  ) +
  geom_text(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      select(
        location_a, edges_normed
      ) %>%
      rbind(
        data.frame(
          location_a = "Overall edges per 1000 gene pairs",
          edges_normed = expected_edges_per_thousand
        )
      ) %>%
      mutate(
        section_rn = row_number(desc(edges_normed)),
        xend = max(cluster_cluster_edges_ranked$rank) + 70 + 5,
        x = max(cluster_cluster_edges_ranked$rank) + 50,
        y = edges_normed,
        yend = 4.2 - 0.25 * section_rn,
        color = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          "grey20",
          location_colors[location_a]
        ),
        label = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          paste0(location_a, " (", round(edges_normed, 2), ")"),
          paste0(location_a, " (", round(edges_normed, 2), ")")
        )
      ),
    aes(
      x = xend,
      y = yend,
      color = color,
      label = label
    ),
    hjust = 0,
    family = "Arial",
    size = 3.5
  ) +
  # end location location label
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(0, 8, 2),
    limits = c(0, 8)
  ) +
  scale_x_continuous(
    name = "Cluster Fitness Similarity Rank",
    breaks = c(
      1, max(cluster_cluster_edges_ranked$rank)
    ),
    limits = c(
      -50, max(cluster_cluster_edges_ranked$rank) + 400
    ),
    expand = c(0, 0)
  ) +
  coord_fixed(
    ratio = (diff(range(cluster_cluster_edges_ranked$rank)) + 100) / 8
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/cluster_fitness_similarity_rank_",
    MAX_RANK,
    "_self_similarity_highlight_subcellular-location-labelled.pdf"
  ),
  height = 8,
  width  = 14,
  device = cairo_pdf
)

################################################################################
# log transform

ggplot(
  cluster_cluster_edges_ranked,
  aes(
    x = rank,
    y = log_edges_normed
  )
) +
  geom_segment(
    aes(
      y    = log10(expected_edges_per_thousand),
      yend = log10(expected_edges_per_thousand),
      x = -50,
      xend = max(cluster_cluster_edges_ranked$rank) + 50,
      color = location_colors[location_a]
    ),
    linetype = 2,
    linewidth = 1,
    color = "grey20",
  ) +
  geom_segment(
    data = inter_location_edges,
    aes(
      y = log_edges_normed,
      yend = log_edges_normed,
      x = -50,
      xend = max(cluster_cluster_edges_ranked$rank) + 50,
      color = location_colors[location_a]
    ),
    linetype = 2,
    linewidth = 0.8,
  ) +
  geom_line() +
  geom_point(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(color = COLORS[GROUPS[as.character(cluster_a)]]),
    aes(
      color = color
    ),
    shape = 16,
    size = 4,
    alpha = 0.8
  ) +
  # start cluster cluster label
  geom_segment(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(
        section = case_when(
          rank < 15 ~ 'a',
          rank < 151 ~ 'b',
          .default = 'c'
        ),
        xend = case_when(
          section == 'a' ~ 50,
          section == 'b' ~ 300,
          section == 'c' ~ 130,
        ),
        color = COLORS[GROUPS[as.character(cluster_a)]],
      ) %>%
      group_by(section) %>%
      mutate(
        section_rn = row_number(rank),
        yend = case_when(
          section == 'a' ~ 0.91 - section_rn * 0.021,
          section == 'b' ~ 0.90 - section_rn * 0.021,
          section == 'c' ~ 0.20 - section_rn * 0.021,
        )
      ),
    aes(
      xend = xend,
      yend = yend,
      color = color
    )
  ) +
  geom_text(
    data = cluster_cluster_edges_ranked %>%
      filter(cluster_a == cluster_b) %>%
      mutate(
        section = case_when(
          rank < 15 ~ 'a',
          rank < 151 ~ 'b',
          .default = 'c'
        ),
        x = case_when(
          section == 'a' ~ 50 + 5,
          section == 'b' ~ 300 + 5,
          section == 'c' ~ 130 - 5,
        ),
        color = COLORS[GROUPS[as.character(cluster_a)]],
        label = paste0(
          rank,
          ". Cluster ",
          cluster_a,
          " - Cluster ",
          cluster_b
        )
      ) %>%
      group_by(section) %>%
      mutate(
        section_rn = row_number(rank),
        y = case_when(
          section == 'a' ~ 0.91 - section_rn * 0.021,
          section == 'b' ~ 0.90 - section_rn * 0.021,
          section == 'c' ~ 0.20 - section_rn * 0.021,
        ),
        hjust = ifelse(
          section %in% c('a', 'b'),
          0,
          1
        )
      ),
    aes(
      label = label,
      x = x,
      y = y,
      hjust = hjust,
      color = color
    ),
    family = "Arial",
  ) +
  # end cluster cluster label
  # start location location label
  geom_segment(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      select(
        location_a, log_edges_normed
      ) %>%
      rbind(
        data.frame(
          location_a = "Overall edges per 1000 gene pairs",
          log_edges_normed = log10(expected_edges_per_thousand)
        )
      ) %>%
      mutate(
        section_rn = row_number(desc(log_edges_normed)),
        xend = max(cluster_cluster_edges_ranked$rank) + 70,
        x = max(cluster_cluster_edges_ranked$rank) + 50,
        y = log_edges_normed,
        yend = 0.6 - 0.03 * section_rn,
        color = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          "grey20",
          location_colors[location_a]
        ),
        label = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          paste0(location_a, " (", log_edges_normed, ")"),
          paste0(location_a, " - ", location_a, " (", log_edges_normed, ")")
        )
      ),
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = color
    )
  ) +
  geom_text(
    data = location_location_edges_ranked %>%
      filter(location_a == location_b) %>%
      select(
        location_a, log_edges_normed
      ) %>%
      rbind(
        data.frame(
          location_a = "Overall edges per 1000 gene pairs",
          log_edges_normed = log10(expected_edges_per_thousand)
        )
      ) %>%
      mutate(
        section_rn = row_number(desc(log_edges_normed)),
        xend = max(cluster_cluster_edges_ranked$rank) + 70 + 5,
        x = max(cluster_cluster_edges_ranked$rank) + 50,
        y = log_edges_normed,
        yend = 0.6 - 0.03 * section_rn,
        color = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          "grey20",
          location_colors[location_a]
        ),
        label = ifelse(
          location_a == "Overall edges per 1000 gene pairs",
          paste0(location_a, " (", round(log_edges_normed, 2), ")"),
          paste0(location_a, " (", round(log_edges_normed, 2), ")")
        )
      ),
    aes(
      x = xend,
      y = yend,
      color = color,
      label = label
    ),
    hjust = 0,
    family = "Arial",
  ) +
  # end location location label
  scale_color_identity() +
  scale_y_continuous(
    name = paste0(
      "Log10 Rank ", MAX_RANK, " Edges per 1000 Gene Pairs"
    ),
    breaks = seq(-0.1, 0.9, 0.2),
    limits = c(-0.1, 0.9)
  ) +
  scale_x_continuous(
    name = "Cluster Fitness Similarity Rank",
    breaks = c(
      1, max(cluster_cluster_edges_ranked$rank)
    ),
    limits = c(
      -50, max(cluster_cluster_edges_ranked$rank) + 400
    ),
    expand = c(0, 0)
  ) +
  coord_fixed(
    ratio = (diff(range(cluster_cluster_edges_ranked$rank)) + 100) / 1
  ) +
  base_theme()

ggsave(
  paste0(
    "outputs/hockeys/cluster_fitness_similarity_rank_",
    MAX_RANK,
    "_log_self_similarity_highlight_subcellular-location-labelled.pdf"
  ),
  height = 8,
  width  = 14,
  device = cairo_pdf
)

################################################################################
# polr1f partners

cge_cor_ranked %>%
  select(gene, POLR1F) %>%
  merge(
    gin_cluster_map_redundant,
  ) %>%
  arrange(POLR1F) %>%
  filter(
    gene %in% c(
      "",
      "EP300",
      "MED12L",
      "AFF4",
      "AFF1",
      "MED25",
      "CCNK1"
    )
  ) %>%
  mutate(group = GROUPS[cluster])

top_instead

cge_cor_ranked %>%
  select(gene, POLR2) %>%
  merge(
    gin_cluster_map_redundant,
  ) %>%
  arrange(POLR2A) %>%
  filter(
    gene %in% c(
      "CREBBP",
      "EP300",
      "MED12L",
      "AFF4",
      "AFF1",
      "MED25",
      "CCNK1"
    )
  ) %>%
  mutate(group = GROUPS[cluster])











































