# RUFF_KING_Grammars_of_IDRs_using_NARDINI+

This repository contains the analysis code associated with the **GIN: Grammars of IDRs using NARDINI+** project, led by **Kiersten M. Ruff** and **Matthew R. King**. 


## GIN Resource

The GIN resource inlcudes molecular grammar analyses of all 24,508 predicted human IDR of length ≥ 30. 
Here, molecular grammar refers to the non-random amino acid composition and the non-random patterning of distinct pairs of amino acid types with respect to one another.
Molecular grammars for each IDR are reported as a 90-feature z-score vector (ZSV). 
These ZSVs were utilized to cluster the human IDRome into 30 distinct grammar clusters which we refer to as GIN clusters. 
Each GIN cluster has a distinct set of grammar features that define it. 
To faciliate the use of GIN as a resource, we created two Google Colab notebooks.


## Google Colab Notebook: [NARDINI+_from_accession](https://colab.research.google.com/drive/15O00GXapuDmD8AijSvkFfRm2fo0xgFPv#scrollTo=DNWtoXwYr71X)

In this notebook, users input a list of proteins as either a list of gene names or Uniprot accessions. The list can be comma separated or uploaded as a file with one gene / accession per line. 
This notebook takes advantage of the fact that GIN was created by analyzing 24,508 IDRs from the human IDRome and thus can extract the IDRs from the user inputed list of genes / accessions. 
The notebook produces two major outputs: (1) a summary of the GIN cluster annotations and IDR information and (2) the ZSVs for all IDRs within the list of genes / accessions. 
Users also have additional output options. These include a schematic of the proteins of interest. Here, IDRs are colored and labeled by their GIN cluster. Domains, downloaded from Uniprot, are shown in yellow 
and labeled by the domain name. To visualize exceptional grammars within the IDRs of interest, users can also plot and download hierarchically clustered heatmaps of the ZSVs.

Users also have access to extracting IDRs and their sequence grammars from eight other species. Here, the compositonal features of selected IDRome are used as the prior distirbution.
However, if a species besides Homo sapiens in selected the GIN clusters and sequence schematics are not shown given that clustering was only done for the human IDRome. 

## Google Colab Notebook: [NARDINI+_from_fasta](https://colab.research.google.com/drive/1Lmb0pm5iFUOC4_ecBnFdmmcT0EOjLvfu#scrollTo=F6reO_2GvIrv)

In this notebook, users input an IDR sequence or a list of IDR sequences in FASTA format. This notebook is helpful if users want to specificy their own definitions of IDR sequences. 
However, in this case NARDINI+ must be run for each inputed IDR sequence. Thus, NARDINI+ is broken down into its two components: (1) calculation of compositional ZSVs and (2) calculation of patterning ZSVs. 
Patterning ZSVs are calculated using NARDINI and this involves creating 100,000 sequence scrambles of the IDR of interest to use as a prior distribution. 
Due to this, generating the patterning ZSV for a given IDR can take minutes to run depending on the length of the IDR. 
To facilitate reasonable run times, we have set the maximum number of IDRs to 20 and the maximum length of the IDR to 1000 for the calculation of the patterning ZSVs.
For more IDRs or longer IDRs, users can extract patterning ZSVs using [NARDINI](https://www.sciencedirect.com/science/article/pii/S0022283621006100) within [localCIDER](https://github.com/Pappulab/localCIDER).
The final outputs of this notebook are a summary of the GIN cluster annotations and the ZSVs for all IDRs in the inputed list. 

Users also have access to extracting IDRs and their sequence grammars from eight other species. Here, the compositonal features of selected IDRome are used as the prior distirbution.
However, if a species besides Homo sapiens in selected the GIN clusters are not shown since human IDRome was utilized as the prior to construct these clusters. 

