# RUFF_KING_Grammars_of_IDRs_using_NARDINI+

This repository contains the analysis code associated with the **GIN: Grammars Inferred using NARDINI+** project, led by **Kiersten M. Ruff** and **Matthew R. King**. 

The manuscript is has been published in Cell [here](https://www.cell.com/cell/fulltext/S0092-8674(25)01191-2).

## [NARDINI+ Mutant Generator](https://colab.research.google.com/drive/1VoWfOvSDgZ04ZeGXazLDndD9c2Hm4YCy#scrollTo=xLTTCTcIeP8C)

Our NARDINI+ framework can be used to redesign IDRs to mutate specific non-random grammars and test their effect on IDR function. 
The main process is to identify a non-random feature and mutate that feature while keeping most other grammar features similar to the wild-type sequence.
Thus, one can identify the effect of a specific grammar feature on IDR function. 
We have created a Google Colab notebook that allows users to input an IDR sequence of interest. Then, the notebook identifies the non-random features, allows users to pick which non-random feature they want to mutate, and outputs a mutated version of the input sequence based on that choice.
The alpha version of our NARDINI+ Mutant Generator can be found [here](https://colab.research.google.com/drive/1VoWfOvSDgZ04ZeGXazLDndD9c2Hm4YCy#scrollTo=xLTTCTcIeP8C).


Examples of mutating IDRs based on their non-random features to determine IDR function can be found in the following published manucripts:

1. [A. Patil, A.R. Strom, J.A. Paulo, C.K. Collings, K.M. Ruff, M.K. Shinn, A Sankar, K.S. Cervantes, T. Wauer, J.D. St Laurent, G Xu, L.A. Becker, S.P. Gygi, R.V. Pappu, C.P. Brangwynne, Cigall Kadoch. (2023). A disordered region controls cBAF activity via condensation and partner recruitment. Cell 186 (22), 4936-4955. e26.](https://www.cell.com/cell/fulltext/S0092-8674(23)00965-0)

2. [C. Hoffmann, K.M. Ruff, I.A. Edu, M.K. Shinn, J.V. Tromm, M.R. King, A. Pant, H. Ausserwöger, J.R. Morgan, T.P.J. Knowles, R.V. Pappu, D. Milovanovic. (2025). Synapsin condensation is governed by sequence-encoded molecular grammars. Journal of molecular biology 437 (8), 168987.](https://www.sciencedirect.com/science/article/pii/S0022283625000531)

3. [K. Meyer, K. Yserentant, R. Cheloor-Kovilakam, K.M. Ruff, C. Chung, X. Shu, B. Huang, O.D. Weiner. (2025). YAP charge patterning mediates signal integration through transcriptional co-condensates. Nature Communications 16 (1), 7454.](https://www.nature.com/articles/s41467-025-62157-3)

![](http://drive.google.com/uc?export=view&id=1uE1njkXIShciesx9rspyXIkwaBkC7dme)

![](http://drive.google.com/uc?export=view&id=1P9zmKcgINsXN1B_qLToLgRRRaOPS2FtC)

## GIN Resource

The GIN resource includes molecular grammar analyses of all 24,508 predicted human IDR of length ≥ 30. 
Here, molecular grammar refers to the non-random amino acid composition and the non-random patterning of distinct pairs of amino acid types with respect to one another.
Molecular grammars for each IDR are reported as a 90-feature z-score vector (ZSV). 
These ZSVs were utilized to cluster the human IDRome into 30 distinct grammar clusters which we refer to as GIN clusters. 
Each GIN cluster has a distinct set of grammar features that define it. 
To faciliate the use of GIN as a resource, we created three Google Colab notebooks.

## Google Colab Notebook: [NARDINI+_from_accession](https://colab.research.google.com/drive/15O00GXapuDmD8AijSvkFfRm2fo0xgFPv#scrollTo=DNWtoXwYr71X)

This Google Colab notebook can be found [here](https://colab.research.google.com/drive/15O00GXapuDmD8AijSvkFfRm2fo0xgFPv#scrollTo=DNWtoXwYr71X).

In this notebook, users input a list of proteins as either a list of gene names or Uniprot accessions. The list can be comma separated or uploaded as a file with one gene / accession per line. 
This notebook takes advantage of the fact that GIN was created by analyzing 24,508 IDRs from the human IDRome and thus can extract the IDRs from the user inputed list of genes / accessions. 
The notebook produces two major outputs: (1) a summary of the GIN cluster annotations and IDR information and (2) the ZSVs for all IDRs within the list of genes / accessions. 
Users also have additional output options. These include a schematic of the proteins of interest. Here, IDRs are colored and labeled by their GIN cluster. Domains, downloaded from Uniprot, are shown in yellow 
and labeled by the domain name. To visualize exceptional grammars within the IDRs of interest, users can also plot and download hierarchically clustered heatmaps of the ZSVs.

Users also have access to extracting IDRs and their sequence grammars from eight other species. Here, the compositonal features of the selected IDRome are used as the prior distirbution.
However, if a species besides Homo sapiens is selected the GIN clusters and sequence schematics are not shown given that clustering was only done for the human IDRome. 

![](http://drive.google.com/uc?export=view&id=1puG8zqJhMcVwBRYtCVrJauLuDzbQAF9G)

## Google Colab Notebook: [NARDINI+_from_fasta](https://colab.research.google.com/drive/1Lmb0pm5iFUOC4_ecBnFdmmcT0EOjLvfu#scrollTo=F6reO_2GvIrv)

This Google Colab notebook can be found [here](https://colab.research.google.com/drive/1Lmb0pm5iFUOC4_ecBnFdmmcT0EOjLvfu#scrollTo=F6reO_2GvIrv).

In this notebook, users input an IDR sequence or a list of IDR sequences in FASTA format. This notebook is helpful if users want to specificy their own definitions of IDR sequences. 
However, in this case NARDINI+ must be run for each inputed IDR sequence. Thus, NARDINI+ is broken down into its two components: (1) calculation of compositional ZSVs and (2) calculation of patterning ZSVs. 
Patterning ZSVs are calculated using NARDINI and this involves creating 100,000 sequence scrambles of the IDR of interest to use as a prior distribution. 
Due to this, generating the patterning ZSV for a given IDR can take minutes to run depending on the length of the IDR. 
To facilitate reasonable run times, we have set the maximum number of IDRs to 20 and the maximum length of the IDR to 1000 for the calculation of the patterning ZSVs.
For more IDRs or longer IDRs, users can extract patterning ZSVs using [NARDINI](https://www.sciencedirect.com/science/article/pii/S0022283621006100) within [localCIDER](https://github.com/Pappulab/localCIDER).
The final outputs of this notebook are a summary of the GIN cluster annotations and the ZSVs for all IDRs in the inputed list. 

Users also have access to extracting IDRs and their sequence grammars from eight other species. Here, the compositonal features of the selected IDRome are used as the prior distirbution.
However, if a species besides Homo sapiens is selected the GIN clusters are not shown since human IDRome was utilized as the prior to construct these clusters. 

## Google Colab Notebook: [GIN cluster and NARDINI+ analysis for orthologs](https://colab.research.google.com/drive/1RNsWk-eP1z_DyKr9lI6Q7nBJjR2xY04H#scrollTo=Nc60ua4_k707)

This Google Colab notebook can be found [here](https://colab.research.google.com/drive/1RNsWk-eP1z_DyKr9lI6Q7nBJjR2xY04H#scrollTo=Nc60ua4_k707)

In this notebook, users input a human uniprot accession number and then the associated human GIN clusters are mapped and NARDINI+ z-score vectors are shown for all human and orthologous IDRs. This notebook takes advantage of previously generated orthologous sequences extracted using BLASTp. Additionally, the notebook can rank order specific grammar features for each species and select the IDR with the top rank per species. For the latter, heatmaps are generated for (1) the rank of the top scoring IDRs (2) the percent identity of the top scoring IDRs to the human IDR, and (3) the z-score / grammar value of the top scoring IDRs. The total number of proteins per species refers to the number of unique proteins that are orthologs to human IDR containing proteins. This choice was made for tractability.

![](http://drive.google.com/uc?export=view&id=1flV6XzMRifANlxcL-PhKDzdK2u3l4Isr)
