library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(readxl)

dt_mut <- read_xlsx("data/aj-apms/Table S2-1A-WT.DBDmut.delIDR1.CBR_TbID.xlsx", sheet = "Table S2-1A-WT.mutants_TurboID")
dt_fus <- read_xlsx("data/aj-apms/Table S2-1A-WT.DBDmut.delIDR1.CBR_TbID.xlsx", sheet = "Table S2-1A-WT.FUSIDR_TbID")
dt_idr <- read_xlsx("data/aj-apms/Table S2-1A-WT.DBDmut.delIDR1.CBR_TbID.xlsx", sheet = "Table S2-1A-WT.42YS.AQG_TbID")
names(dt_mut) %>% writeLines()
names(dt_fus) %>% writeLines()
names(dt_idr) %>% writeLines()
coldata <- fread("data/aj-apms/turbo_coldata.tsv")

BAF_SYMBOLS <- c(
  "ACTB",
  "ACTL6A",
  "ACTL6B",
  "ARID1A",
  "ARID1B",
  "ARID2",
  "BICRA",
  "BICRAL",
  "BCL7B",
  "BCL7A",
  "BCL7C",
  "BRD7",
  "BRD9",
  "DPF1",
  "DPF2",
  "DPF3",
  "PBRM1",
  "PHF10",
  "SMARCA2",
  "SMARCA4",
  "SMARCB1",
  "SMARCC1",
  "SMARCC2",
  "SMARCD1",
  "SMARCD2",
  "SMARCD3",
  "SMARCE1",
  "SS18",
  "SS18L1"
)

################################################################################

dt_idr %>%
  filter(
    `1A.42YS_vs_1A.WT_lfc`  < -1 &
    `1A.42YS_vs_1A.WT_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
  ) %>%
  ggplot(
    aes(
      x = 1,
      y = freq,
      fill = color
    )
  ) +
  geom_bar(
    stat = "identity",
    color = "black"
  ) +
  scale_fill_identity() +
  base_theme()

dt_idr %>%
  filter(
    `1A.AQGscram_vs_1A.WT_lfc`  < -1 &
    `1A.AQGscram_vs_1A.WT_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
  ) %>%
  ggplot(
    aes(
      x = 1,
      y = freq,
      fill = color
    )
  ) +
  geom_bar(
    stat = "identity",
    color = "black"
  ) +
  scale_fill_identity() +
  base_theme()



cge_cor_edges %>%
  filter(cluster_a == 21 & rank <= 20) %>%
  rename(cluster = cluster_b) %>%
  group_by(cluster) %>%
  summarize(n = n()) %>%
  ungroup %>%
  mutate(
    freq = n / sum(n),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
  ) %>%
  ggplot(
    aes(
      x = cluster,
      y = freq,
      fill = color
    )
  ) +
  geom_bar(
    stat = "identity",
    color = "black"
  ) +
  scale_fill_identity() +
  base_theme()



dt_1a <- read_xlsx("data/aj-apms/Table S1-IP.MassSpec.xlsx", sheet = "Table S1-1A WT.mutants.IP.MS")
dt_1b <- read_xlsx("data/aj-apms/Table S1-IP.MassSpec.xlsx", sheet = "Table S1-1B WT.mutant_IP.MS")

plot_bar <- function(contrast) {
  dt_1a %>%
    select(
      `Gene Symbol`,
      contains(paste0(contrast, "_lfc")),
      contains(paste0(contrast, "_logp"))
    ) %>%
    filter(
      if_all(
        contains("_logp"),
        ~.x > 0.6
      ) &
      if_all(
        contains("_lfc"),
        ~.x < -1
      )
    ) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
  ) %>%
  ggplot(
    aes(
      x = 1,
      y = freq,
      fill = color
    )
  ) +
  geom_bar(
    stat = "identity",
    color = "black"
  ) +
  ggtitle(contrast) +
  scale_fill_identity() +
  base_theme()
}
################################################################################

rbindlist(list(
dt_fus %>%
  filter(
    `AN3CA 1A-delIDR1-TbID_vs_AN3CA 1A-WT-TbID_lfc`  < -1 &
    `AN3CA 1A-delIDR1-TbID_vs_AN3CA 1A-WT-TbID_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(`Gene Symbol`, cluster) %>%
  summarize(
    freq = n()
  ) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
    contrast = "AN3CA 1A-delIDR1 vs 1A-WT Turbo-ID 1, Lost Proteins"
  ),

dt_mut %>%
  filter(
    `1A-delIDR1-TbID_vs_1A-WT-TbID_lfc`  < -1 &
    `1A-delIDR1-TbID_vs_1A-WT-TbID_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(`Gene Symbol`, cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
    contrast = "AN3CA 1A-delIDR1 vs 1A-WT Turbo-ID 2, Lost Proteins"
  ),

dt_mut %>%
  filter(
    `1A-CBR-TbID_vs_1A-WT-TbID_lfc`  < -1 &
    `1A-CBR-TbID_vs_1A-WT-TbID_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(`Gene Symbol`, cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
    contrast = "AN3CA 1A-CBR vs 1A-WT Turbo-ID 2, Lost Proteins"
  ),

dt_fus %>%
  filter(
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_lfc`  < -1 &
    `AN3CA 1A-FUS IDR-TbID_vs_AN3CA 1A-WT-TbID_logp` > 2
  ) %>%
  select(`Gene Symbol`) %>%
  merge(
    gin_cluster_map,
    by.x = "Gene Symbol",
    by.y = "gene"
  ) %>%
  group_by(`Gene Symbol`, cluster) %>%
  summarize(freq = n()) %>%
  ungroup %>%
  mutate(
    freq = freq / cluster_sizes[as.character(cluster)],
    freq = freq / sum(freq),
    group = GROUPS[as.character(cluster)],
    color = COLORS[group],
    contrast = "AN3CA 1A-FUS IDR vs 1A-WT Turbo-ID 1, Lost Proteins"
  )
)) %>%
  group_by(`Gene Symbol`) %>%
  mutate(n_occ = n()) %>%
  ungroup %>%
  filter(n_occ > 3) %>%
  ggplot(
    aes(
      x = contrast,
      y = freq,
      fill = color
    )
  ) +
  geom_bar(
    stat = "identity",
    color = "black"
  ) +
  scale_fill_identity() +
  base_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )










