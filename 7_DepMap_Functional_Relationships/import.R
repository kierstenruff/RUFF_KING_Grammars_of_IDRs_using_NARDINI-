library(data.table)
library(dplyr)

if (file.exists("data/biomart_hgnc_uniprotswissprot.csv")) {
  mart <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  biomart_hgnc_uniprotswissprot <- biomaRt::getBM(
    attributes = c(
      "hgnc_symbol",
      "uniprotswissprot"
    ),
    filters = c("biotype"),
    values = list(biotype="protein_coding"),
    mart = mart
  )
  biomart_hgnc_uniprotswissprot %>%
    fwrite(
      "data/biomart_hgnc_uniprotswissprot.csv"
    )
} else {
  biomart_hgnc_uniprotswissprot <- fread(
    "data/biomart_hgnc_uniprotswissprot.csv"
  )
}

cge_cor_ranked <- fread(
  "~/code/depmap_data/24q4/processed/cge_cor_ranked.csv",
  key = "gene"
)

gin <- fread("data/GIN_IDRome_all_human_gte_70_and_nonlinkers_gte_50_idrs_dist_gte_1.5.csv") %>%
  #filter(`Min Inter Clust Dist` > 1.5) %>%
  tidyr::separate_wider_delim(
    cols = "Gene",
    " IDR",
    names = c(
      "hgnc_manual",
      "idr"
    )
  ) %>%
  merge(
    .,
    biomart_hgnc_uniprotswissprot,
    by.x = "Uniprot",
    by.y = "uniprotswissprot",
    all.x = T
  ) %>%
  mutate(
    hgnc_symbol = ifelse(
      hgnc_symbol %in% colnames(cge_cor_ranked),
      hgnc_symbol,
      hgnc_manual
    )
  ) %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol %in% colnames(cge_cor_ranked)) %>%
  select(-hgnc_manual) %>%
  distinct

gin_features <- fread("data/gin_feature_scores.csv") %>%
  tidyr::separate_wider_delim(
    cols = "Gene",
    " IDR",
    names = c(
      "hgnc_manual",
      "idr"
    )
  ) %>%
  merge(
    .,
    biomart_hgnc_uniprotswissprot,
    by.x = "Uniprot",
    by.y = "uniprotswissprot",
    all.x = T
  ) %>%
  mutate(
    hgnc_symbol = ifelse(
      hgnc_symbol %in% colnames(cge_cor_ranked),
      hgnc_symbol,
      hgnc_manual
    )
  ) %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol %in% colnames(cge_cor_ranked)) %>%
  select(-hgnc_manual) %>%
  distinct
  
hpa <- fread("data/subcellular_location.tsv")












