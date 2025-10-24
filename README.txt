__________________________________________________________________________

This repository contains a pipeline for analyzing drug combination synergy 
in A549 cells by propagating drug-target signals through a cancer-specific 
network, using Random Walk with Restart and graph-regularized NMF to learn 
metagene signatures. Subsequent steps cluster drug combinations based on 
the metagene signatures, identify and analyze highly synergistic clusters, 
and pinpoint top genes for enrichment analysis.

__________________________________________________________________________

PIPELINE OVERVIEW

Step 1 : 1_Network Propagation.R
• Language : R
• Purpose : Propagate drug-target signals through an A549-specific
network using Random Walk with Restart
• Main inputs : 
- gene_interaction_cancer_subnetwork.csv
- subnetwork_nodes.csv
- Drug_Combination_A549.xlsx
- DrugCombInationID.csv
• Main output : 
- Network_Propagation_Result_alpha0.5.csv

Step 2 : 2_Propagated scores differentiating between drug combinations with low and high synergy (Fig 2 and 3).R
• Language: R
• Purpose: Statistically test and visualize the ability of network-propagated gene scores to distinguish between drug combinations in the lowest (Q1) and highest (Q4) quartiles of synergy scores (HSA, ZIP, BLISS, and LOEWE).
• Main Inputs:
drug_combination_infomation.csv (contains synergy scores, equivalent to Supplementary Data S1.csv)
Network_Popagation_Result_alpha0.5(-log).csv (contains the negative log10 of the network propagation results with alpha=0.5, equivalent to -log10 of data in Supplementary Data S3.csv)
• Main Outputs:
Fig. 2 of the main text
Fig. 3 of the main text

Step 3 : 3_GNMF.py
• Language : Python
• Purpose : Run graph-regularised NMF (λ = 0.9) to learn metagenes,
select the optimal component count (elbow method)
• Main inputs : Network_Propagation_Result_alpha0.5(-log).csv
gene_interaction_cancer_subnetwork.csv
Drug_Combination_A549.xlsx
• Main output : 
H.csv
W.csv

Step 4 : 4_Identifying highly synergistic clusters (Table 1 and Suppl Table S1).R
• Language: R
• Purpose: Identify highly synergistic clusters of drug combinations by performing hierarchical clustering on metagene signatures and then determining which resulting clusters have an average synergy score (HSA, ZIP, BLISS, or LOEWE) greater than the Q3 of all drug combinations.
• Main Inputs:
drug_combination_infomation.csv (Contains synergy scores, equivalent to Supplementary Data S1.csv)
Network_Popagation_Result_alpha0.5(-log).csv (Negative log10 of network propagation results, equivalent to -log10 of data in Supplementary Data S3.csv)
H.csv (the H matrix from a graph-regularized non-negative matrix factorization (GNMF), representing metagene signatures)
• Main Output:
highly_synergistic_clusters_results.csv: A combined data frame listing identified highly synergistic clusters. The table is structured by the dendrogram height (h) at which the cluster was cut, and includes cluster statistics (size, mean, standard deviation of synergy, main drugs) for each synergy metric (HSA, ZIP, BLISS, LOEWE).

Step 5 : 5_Analyzing 7 highly synergistic clusters (Fig 4, 5, 6).R
• Language: R
• Purpose: Perform a detailed analysis on the seven highly synergistic clusters of drug combinations at a specific dendrogram cut height (h). It generates visualizations (boxplots and heatmaps) of the metagene signatures for these clusters. Then it focuses on analyzing individual drugs (e.g., Dasatinib, Paclitaxel, Quinacrine).
• Main Inputs:
drug_combination_infomation.csv (Contains synergy scores, equivalent to Supplementary Data S1.csv)
Network_Popagation_Result_alpha0.5(-log).csv (Negative log10 of network propagation results, equivalent to -log10 of data in Supplementary Data S3.csv)
H.csv (the H matrix, or metagene signatures, used for clustering)
• Main Outputs:
Fig. 4 of the main text
Fig. 5 of the main text
Fig. 6 of the main text

Step 6 : 6_Analyzing partner drugs (Fig 7).R
• Language: R
• Purpose: Analyze a list of nine highly frequent drug partners (Lapatinib, Vemurafenib, Sunitinib, Axitinib, Vismodegib, Sorafenib, Lenalidomide, Vandetanib, and Imatinib) by visualizing the distribution of their synergy scores and their associated metagene signatures across all drug combinations they are involved in.
• Main Inputs:
drug_combination_infomation.csv (Contains synergy scores, equivalent to Supplementary Data S1.csv)
Network_Popagation_Result_alpha0.5(-log).csv (Negative log10 of network propagation results, equivalent to -log10 of data in Supplementary Data S3.csv)
H.csv (the H matrix, or metagene signatures)
• Main Outputs:
Fig. 7 of the main text

Step 7 : 7_Identifying top genes in Metagene 2 (Suppl Table S10).R
• Language: R
• Purpose: Identify the 200 genes that contribute the most to Metagene 2, based on the W matrix. This list will be used for KEGG enrichment analysis.
• Main Input:
W.csv (The matrix W consists of metagenes (1-26) as columns and genes as rows, indicating the weight of each gene for each metagene.).
• Main Output:
200metagene2.csv: A list of the top 200 genes contributing to Metagene 2.
__________________________________________________________________________


