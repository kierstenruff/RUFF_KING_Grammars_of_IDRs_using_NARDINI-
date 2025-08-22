library(data.table)
library(ggplot2)
library(dplyr)

dir.create("export/motifnet", recursive = T, showWarnings = F)

source("utils.R")

################################################################################

gin %>%
  group_by(hgnc_symbol) %>%
  summarize(
    max_idr_length = max(`IDR Length`),
    max_micd = max(`Min Inter Clust Dist`),
    top_idr_cluster = `Cluster Number`[
      which(
        `Min Inter Clust Dist`[`IDR Length` == max_idr_length] == max(
          `Min Inter Clust Dist`[`IDR Length` == max_idr_length]
        )
      )
    ],
    cluster_group = GROUPS[top_idr_cluster]
  ) %>%
  group_by(cluster_group) %>%
  group_split() %>%
  lapply(
    function(x) {
      cluster <- unique(x$cluster_group)
      #fwrite(
      #  list(x$hgnc_symbol),
      #  paste0("export/motifnet/nodes_cluster_", cluster, ".txt"),
      #  col.names = F
      #)
      print(nrow(x))
    }
  )

cge_cor_ranked %>%
  rename(gene_a = gene) %>%
  melt(
    id.vars = "gene_a",
    variable.name = "gene_b",
    value.name = "rank"
  ) %>%
  mutate(gene_b = as.character(gene_b)) %>%
  filter(rank < 5) %>%
  filter(gene_b > gene_a) %>%
  select(gene_a, gene_b) %>%
  fwrite("export/motifnet/edges_rank4.txt", sep = '\t', col.names = F)


















