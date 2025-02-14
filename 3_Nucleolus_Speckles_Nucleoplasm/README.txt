0_bsp_analyze_IDRs_by_clusters_and_location
- Python notebook that analyzes IDRs by cluster and location
- Creates plots for Figures S4 and S5
- Also gets data for SanKey plots of Figures 4C, 4E, S3B

1_bsp_nucleolus_complex_clusters
- Python notebook that connects ribosomal biogenesis processes to IDR grammars
- Uses graph_env (see below)
- Creates plot for Figure 4C 

2_bsp_nucleolus_complex_ranked_grammars
- Python notebook that determined highest scoring grammars in IDRome
- Creates plots for Figures 5A

3_bsp_nuclear_speckle_complex_clusters
- Python notebook that connects spliceosome processes to IDR grammars
- Uses graph_env (see below)
- Creates plot for Figure 4E

4_bsp_nuclear_speckles_complex_ranked_grammars
- Python notebook that determined highest scoring grammars in IDRome
- Creates plots for Figures 5C

5_bsp_nucleoplasm_complex_clusters
- Python notebook that connects nucleoplasm complexes to IDR grammars
- Uses graph_env (see below)
- Creates plot for Figure S3B

6_bsp_nucleoplasm_complex_ranked_grammars
- Python notebook that determined highest scoring grammars in IDRome
- Creates plots for Figures 5E

7_bsp_analyze_process_complex_idrs
- Python notebook that gets distributions of GIN types by process
- Creates plots for Figures 2D, 2F, S3C

D3Blocks environment
conda create --name graph_env
pip install d3blocks
pip install matplotlib
pip install scipy
pip install pandas
pip install jupyter
pip install localcider
pip install seaborn
conda install salilab::dssp
pip install biopython


