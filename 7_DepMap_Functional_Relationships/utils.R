DESCRIPTIONS <- c(
  "0"  = "Blocks of glycine and polar residues",
  "1"  = "Blocks of proline and polar residues",
  "2"  = "Blocks of negative, proline, and glycine residues",
  "3"  = "Small negative blocks",
  "4"  = "Weak negative fraction and high N fraction",
  "5"  = "Blocks of positive, proline, and negative residues",
  "6"  = "S-tracts",
  "7"  = "Large negative blocks",
  "8"  = "Blocks of negative and alanine residues",
  "9"  = "Blocks of positive and negative residues",
  "10" = "Weak positive charge, enriched in M and polar residues",
  "11" = "Q-tracts",
  "12" = "Blocks of positive, negative, and polar residues",
  "13" = "Blocks of negative, proline, and polar residues",
  "14" = "Blocks of positive and proline residues",
  "15" = "Well-mixed T-rich sequences",
  "16" = "Blocks of polar and positive",
  "17" = "Weak positive charge",
  "18" = "Large negative blocks with positive blocks",
  "19" = "High negative fraction of specifically Es",
  "20" = "Blocks of alanine and polar residues",
  "21" = "Blocks of proline, glycine, alanine, and polar residues",
  "22" = "Well-mixed proline / glycine",
  "23" = "Strong K blocks",
  "24" = "Weak negative fraction",
  "25" = "Blocks of positive, polar, and hydrophobic residues",
  "26" = "R blocks",
  "27" = "High proline fraction",
  "28" = "High aromatic fraction",
  "29" = "High R fraction "
)

GROUPS <- c(
  "2"  = "Negative",
  "3"  = "Negative",
  "4"  = "Negative",
  "7"  = "Negative",
  "8"  = "Negative",
  "13" = "Negative",
  "19" = "Negative",
  "24" = "Negative",
  "5"  = "Negative & Positive",
  "9"  = "Negative & Positive",
  "12" = "Negative & Positive",
  "18" = "Negative & Positive",
  "23" = "Negative & Positive",
  "10" = "Positive",
  "14" = "Positive",
  "17" = "Positive",
  "25" = "Positive",
  "26" = "Positive",
  "29" = "Positive",
  "1"  = "Proline",
  "21" = "Proline",
  "27" = "Proline",
  "0"  = "Glycine",
  "22" = "Glycine",
  "20" = "Alanine",
  "16" = "Polar",
  "11" = "Polar",
  "15" = "Polar",
  "6"  = "Polar",
  "28" = "Aromatic"
)

COLORS <- c(
  "Negative & Positive" = "black",
  "Negative" = "#e23b36",
  "Positive" = "#1c9ad7",
  "Proline"  = "#7c297f",
  "Glycine"  = "#d865a6",
  "Alanine"  = "#8c8c8c",
  "Polar"    = "#1f8241",
  "Aromatic" = "#f58220"
)  
   
FEATURES <- c(
  "pol-pol",
  "pol-hyd",
  "pol-pos",
  "pol-neg",
  "pol-aro",
  "pol-ala",
  "pol-pro",
  "pol-gly",
  "hyd-hyd",
  "hyd-pos",
  "hyd-neg",
  "hyd-aro",
  "hyd-ala",
  "hyd-pro",
  "hyd-gly",
  "pos-pos",
  "pos-neg",
  "pos-aro",
  "pos-ala",
  "pos-pro",
  "pos-gly",
  "neg-neg",
  "neg-aro",
  "neg-ala",
  "neg-pro",
  "neg-gly",
  "aro-aro",
  "aro-ala",
  "aro-pro",
  "aro-gly",
  "ala-ala",
  "ala-pro",
  "ala-gly",
  "pro-pro",
  "pro-gly",
  "gly-gly",
  "Frac A",
  "Frac C",
  "Frac D",
  "Frac E",
  "Frac F",
  "Frac G",
  "Frac H",
  "Frac I",
  "Frac K",
  "Frac L",
  "Frac M",
  "Frac N",
  "Frac P",
  "Frac Q",
  "Frac R",
  "Frac S",
  "Frac T",
  "Frac V",
  "Frac W",
  "Frac Y",
  "Frac K+R",
  "Frac D+E",
  "Frac Polar",
  "Frac Aliphatic",
  "Frac Aromatic",
  "R/K Ratio",
  "E/D Ratio",
  "Frac Chain Expanding",
  "FCR",
  "NCPR",
  "Hydrophobicity",
  "Disorder Promoting",
  "Iso point",
  "PPII",
  "A Patch",
  "C Patch",
  "D Patch",
  "E Patch",
  "F Patch",
  "G Patch",
  "H Patch",
  "I Patch",
  "K Patch",
  "L Patch",
  "M Patch",
  "N Patch",
  "P Patch",
  "Q Patch",
  "R Patch",
  "S Patch",
  "T Patch",
  "V Patch",
  "Y Patch",
  "RG Frac"
)

locations <- c(
  "7"  = "Nucleoli",
  "3"  = "Vesicles",
  "19" = "Plasma membrane",
  "24" = "Plasma membrane",
  "18" = "Nuclear speckles",
  "12" = "Nuclear bodies",
  "23" = "Nucleoli",
  "26" = "Nuclear speckles",
  "29" = "Plasma membrane",
  "17" = "Microtubules",
  "10" = "Nucleoplasm",
  "22" = "Vesicles",
  "20" = "Intermediate filaments",
  "16" = "Cell Junctions",
  "11" = "Nucleoplasm",
  "15" = "Vesicles",
  "6"  = "Intermediate filaments",
  "28" = "Cell Junctions"
)
location_colors <- c(
  "Nucleoli"         = "#293C8F",
  "Nuclear speckles" = "#189BD7",
  "Nucleoplasm"      = "#218342",
  "Nuclear bodies"   = "olivedrab4",
  "Intermediate filaments" = "goldenrod1",
  "Microtubules"    = "orange1",
  "Cell Junctions"  = "darkorange",
  "Endoplasmic reticulum" = "orangered2",
  "Plasma membrane" = "#E23C36",
  "Vesicles"        = "firebrick4",
  "NA" = "grey70"
)

#location_colors <- c(
#  "Nucleoli"         = "#639",
#  "Nuclear speckles" = "#36b",
#  "Nucleoplasm"      = "#09c",
#  "Nuclear bodies"   = "#0bc",
#  "Intermediate filaments" = "#4d8",
#  "Microtubules"    = "#9d5",
#  "Cell Junctions"  = "#e94",
#  "Endoplasmic reticulum" = "#c66",
#  "Plasma membrane" = "#a35",
#  "Vesicles"        = "#817",
#  "NA" = "grey70"
#)






IDR_FEATURES <- list(
  "22" = c("gly-gly", "pro-gly"),
  "26" = c("R Patch", "Frac R"),
  "7" = c("neg-neg", "E Patch"),
  "23" = c("K Patch", "pos-pos"),
  "11" = c("Q Patch", "Frac Q"),
  "21" = c("pol-gly", "pro-gly"),
  "15" = c("Frac T", "T Patch"),
  "8" = c("neg-ala", "neg-neg"),
  "20" = c("A Patch", "pol-ala"),
  "2" = c("pro-gly", "neg-pro"),
  "1" = c("pol-pro", "pro-pro"),
  "5" = c("pos-neg", "pos-pro"),
  "18" = c("pos-neg", "neg-neg"),
  "0" = c("pol-gly", "gly-gly"),
  "14" = c("pos-pro", "pro-pro"),
  "28" = c("aro-aro", "Frac Y", "Frac Aromatic"),
  "19" = c("Frac E", "Frac D+E"),
  "13" = c("neg-pro", "pol-neg"),
  "10" = c("hyd-hyd", "Frac M"),
  "27" = c("Frac P", "PPII"),
  "6" = c("Frac S", "S Patch"),
  "3" = c("pol-neg", "neg-neg"),
  "12" = c("pol-pos", "pos-neg"),
  "29" = c("R/K Ratio", "Frac R"),
  "9" = c("pos-neg", "neg-neg"),
  "4" = c("Frac N", "Iso point"),
  "17" = c("Iso point"),
  "16" = c("pol-pol", "pol-pos"),
  "25" = c("pos-pos", "pol-pos"),
  "24" = c("Iso point")
)

################################################################################

base_theme <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.key.spacing.y = unit(4, "pt"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      axis.text    = element_text(color = "black", size = 12, family = "Arial"),
      axis.title   = element_text(color = "black", size = 20, family = "Arial"),
      legend.title = element_text(color = "black", size = 20, family = "Arial"),
      title = element_text(color = "black", size = 14, family = "Arial"),
      text  = element_text(color = "black", size = 14, family = "Arial")
    )
}


















