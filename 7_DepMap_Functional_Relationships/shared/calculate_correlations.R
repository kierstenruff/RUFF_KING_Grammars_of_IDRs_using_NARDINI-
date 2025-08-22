library(data.table)
library(dplyr)

version <- "24q4"

dir.create(paste0(version, "/processed"), showWarnings = F)

CRISPRGeneEffect <- paste0(version, "/CRISPRGeneEffect.csv") %>%
  fread()
CRISPRGeneEffect_mat <- CRISPRGeneEffect %>%
  as.matrix(rownames = 1)

# Number of valid genes per model (24q2):
# CRISPRGeneEffect_mat %>% apply(1, function(x) sum(!is.na(x))) %>% table
# 17348 17349 17787 17931 18107 18169 18171 18172 18173 18174 18175 18176 18286
#     2   116   834   189     1     1     1     1     1     2    21   115     1
# 18287 18419 18427 18428 18442 18443 
#     1     1     1    25     1     6 
# Number of valid models per gene:
# CRISPRGeneEffect_mat %>% apply(2, function(x) sum(!is.na(x))) %>% table
# 175   176   177   178   179   316   486  1057  1176  1177  1197  1199  1200
#   1     1     1     7   502    15   129   156     1    95     1     1     1
# 1201  1202  1316  1317  1318  1319  1320 
#    2   421     1     2     7    74 17025 
# from this we only filter out genes in which the number of models with
# a non-NA value is less than 1000
colnames(CRISPRGeneEffect_mat) <- colnames(CRISPRGeneEffect_mat) %>%
  gsub(" \\(\\d+\\)", "", .)
CRISPRGeneEffect_mat <- CRISPRGeneEffect_mat[
  ,
  (
    CRISPRGeneEffect_mat %>% apply(2, function(x) sum(!is.na(x)))
  ) > 1000
]
CRISPRGeneEffect_cor <- CRISPRGeneEffect_mat %>%
  cor(., use = "pairwise.complete.obs")
fwrite(
  data.table(CRISPRGeneEffect_cor, keep.rownames = "gene"),
  paste0(version, "/processed/cge_cor.csv")
)

#CRISPRGeneEffect_cor <- fread("cge_cor.csv") %>%
#  as.matrix(rownames = 1)
col_ranks <- apply(CRISPRGeneEffect_cor, 2, function(x) row_number(desc(x)))
CRISPRGeneEffect_cor_ranked <- pmin(
  col_ranks,
  t(col_ranks)
)
row.names(CRISPRGeneEffect_cor_ranked) <- colnames(CRISPRGeneEffect_cor_ranked)

fwrite(
  data.table(CRISPRGeneEffect_cor_ranked, keep.rownames = "gene"),
  paste0(version, "/processed/cge_cor_ranked.csv")
)
