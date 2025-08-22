dt_fus %>%
  colnames %>%
  writeLines

# note the methods for AJ's paper is wrong.
# The bait norm is actually not performed.
dt_fus_processed <- dt_fus %>%
  rename(
    gene = `Gene Symbol`
  ) %>%
  select(
    gene,
    `AN3CA No Ligase Ctrl Rep1`,
    `AN3CA 1A-WT-TbID Rep1`,
    `AN3CA 1A-FUS IDR-TbID Rep1`,
    `AN3CA No Ligase Ctrl Rep2`,
    `AN3CA 1A-WT-TbID Rep2`,
    `AN3CA 1A-FUS IDR-TbID Rep2`,
    `AN3CA No Ligase Ctrl Rep3`,
    `AN3CA 1A-WT-TbID Rep3`,
    `AN3CA 1A-FUS IDR-TbID Rep3`,
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_lfc`,
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_logp`
  ) %>%
  mutate(
    wt_1 = `AN3CA 1A-WT-TbID Rep1` /
       `AN3CA No Ligase Ctrl Rep1`,
    wt_2 = `AN3CA 1A-WT-TbID Rep2` /
       `AN3CA No Ligase Ctrl Rep2`,
    wt_3 = `AN3CA 1A-WT-TbID Rep3` /
       `AN3CA No Ligase Ctrl Rep3`,
    fus_1 = `AN3CA 1A-FUS IDR-TbID Rep1` /
             `AN3CA No Ligase Ctrl Rep1`,
    fus_2 = `AN3CA 1A-FUS IDR-TbID Rep2` /
             `AN3CA No Ligase Ctrl Rep2`,
    fus_3 = `AN3CA 1A-FUS IDR-TbID Rep3` /
             `AN3CA No Ligase Ctrl Rep3`,
    wt_1  = log2(wt_1 ),
    wt_2  = log2(wt_2 ),
    wt_3  = log2(wt_3 ),
    fus_1 = log2(fus_1),
    fus_2 = log2(fus_2),
    fus_3 = log2(fus_3),
  ) %>%
  select(
    gene,
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_lfc`,
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_logp`,
    starts_with("wt"),
    starts_with("fus"),
  )
wt_genes <- dt_fus_processed %>%
  filter(wt_1 > 0 & wt_2 > 0 & wt_3 > 0) %>%
  select(gene) %>%
  merge(gin_cluster_map, by = "gene") %>%
  group_by(cluster) %>%
  summarize(
    raw_frequency = n(),
    normalized_frequency = length(unique(gene)) /
      cluster_sizes[as.character(unique(cluster))],
  ) %>%
  ungroup %>%
  mutate(
    raw_proportion = raw_frequency / sum(raw_frequency),
    normalized_proportion = normalized_frequency /
      sum(normalized_frequency),
    group = "WT"
  )
fus_genes <- dt_fus_processed %>%
  filter(fus_1 > 0 & fus_2 > 0 & fus_3 > 0) %>%
  select(gene) %>%
  merge(gin_cluster_map, by = "gene") %>%
  group_by(cluster) %>%
  summarize(
    raw_frequency = n(),
    normalized_frequency = length(unique(gene)) /
      cluster_sizes[as.character(unique(cluster))],
  ) %>%
  ungroup %>%
  mutate(
    raw_proportion = raw_frequency / sum(raw_frequency),
    normalized_proportion = normalized_frequency /
      sum(normalized_frequency),
    group = "FUS"
  )

ggplot(
  rbind(
    wt_genes,
    fus_genes
  ),
  aes(
    x = cluster,
    y = normalized_frequency,
    fill = COLORS[GROUPS[as.character(cluster)]]
  )
) +
  geom_bar(
    aes(group = factor(group, levels = c("WT", "FUS"))),
    stat = "identity",
    color = "black",
    position = "dodge"
  ) +
  geom_text(
    aes(
      label = group,
      x = cluster + ifelse(
        group == "WT",
        -0.25,
        0.25
      ),
      y = 0.05,
    ),
    size = 2.5,
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ) +
  scale_x_continuous(
    name = "Cluster",
    breaks = c(0:29),
  ) +
  scale_y_continuous(
    name = "IDRs / Cluster Size",
    expand = c(0, 0),
    limits = c(0, max(wt_genes$normalized_frequency) * 1.05)
  ) +
  scale_fill_identity() +
  base_theme()
    
ggsave(
  "aj_wt_vs_fus_bars.pdf",
  height = 6,
  width = 11,
  device = cairo_pdf
)


ggplot(
  rbind(
    wt_genes,
    fus_genes
  ),
  aes(
    x = group,
    y = raw_proportion,
    fill = COLORS[GROUPS[as.character(cluster)]],
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
  ) +
  scale_x_discrete(
    name = "ARID1A IDR",
    limits = c("WT", "FUS")
  ) +
  scale_y_continuous(
    name = "Proportion of interacting IDRs",
    expand = c(0, 0),
    limits = c(0, 1.05)
  ) +
  scale_fill_identity() +
  base_theme()

ggsave(
  "aj_wt_vs_fus_bars_prop.pdf",
  height = 7,
  width = 6,
  device = cairo_pdf
)


ggplot(
  rbind(
    wt_genes,
    fus_genes,
    data.frame(
      cluster = c("NA", "NA"),
      raw_frequency = c(NA, NA),
      normalized_frequency = c(NA, NA),
      raw_proportion = c(0.9055161, 0.9000666),
      normalized_proportion = c(NA, NA),
      group = c("WT", "FUS")
    )
  ) %>%
    filter(cluster %in% c(10, 11, 28, "NA")),
  aes(
    x = group,
    y = raw_proportion,
    fill = ifelse(
      as.character(cluster) %in% names(GROUPS),
      COLORS[GROUPS[as.character(cluster)]],
      "#aaaaaa"
    ) %>%
      factor(
        levels = c(
        "#aaaaaa",
          COLORS
        )
      )
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
  ) +
  scale_x_discrete(
    name = "ARID1A IDR",
    limits = c("WT", "FUS")
  ) +
  scale_y_continuous(
    name = "Proportion of interacting IDRs",
    expand = c(0, 0),
    limits = c(0, 1.05)
  ) +
  scale_fill_identity() +
  base_theme()

ggsave(
  "aj_wt_vs_fus_bars_prop_v2.pdf",
  height = 7,
  width = 6,
  device = cairo_pdf
)
















