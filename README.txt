__________________________________________________________________________

README!

This repository contains a three-step analysis pipeline
for classifying drug combinations into clusters with potential
therapeutic benefits to non-small cell lung cancer. Each step
is an independent script that must be executed in the order shown below.
__________________________________________________________________________

PIPELINE OVERVIEW

Step 1 : 1_Network Propagation.R
• Language : R
• Purpose : Propagate drug-target signals through an A549-specific
network using Random Walk with Restart
• Main inputs : 
gene_interaction_cancer_subnetwork.csv
subnetwork_nodes.csv
Drug_Combination_A549.xlsx
DrugCombInationID.csv
• Main output : 
Network_Propagation_Result_alpha0.5.csv

Step 2 : 2_Propagated scores differentiating between drug combinations with low and high synergy (Fig 2 and 3).R
• Language: R
• Purpose: Statistically test and visualize the ability of network-propagated gene scores to distinguish between drug combinations in the lowest (Q1) and highest (Q4) quartiles of synergy scores (HSA, ZIP, BLISS, and LOEWE).
• Main Inputs:
drug_combination_infomation.csv (contains synergy scores, equivalent to Supplementary Data S1.csv)
Network_Popagation_Result_alpha0.5(-log).csv (contains the negative log10 of the network propagation results with alpha=0.5, equivalent to -log10 of data in Supplementary Data S3.csv)
• Main Outputs:
Fig2 of the main text
Fig3 of the main text

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
Network_Popagation_Result_alpha0.5(-log).csv (Negative log10 of network propagation results, equivalent to -log10 of Supplementary Data S3.csv)
H.csv (the H matrix from a graph-regularized non-negative matrix factorization (GNMF), representing metagene signatures)
• Main Output:
highly_synergistic_clusters_results.csv: A combined data frame listing identified highly synergistic clusters. The table is structured by the dendrogram height (h) at which the cluster was cut, and includes cluster statistics (size, mean, standard deviation of synergy, main drugs) for each synergy metric (HSA, ZIP, BLISS, LOEWE).


__________________________________________________________________________

DETAIL OF REQUIRED INPUT FILES

• gene_interaction_cancer_subnetwork.csv
Edge list (source, target) for the 2 100-node A549 protein-interaction sub-network

• subnetwork_nodes.csv
Mapping table: node index → HGNC gene symbol

• Drug_Combination_A549.xlsx
Metadata for 607 drug pairs – IDs, Drug A/B names, target node indices,
and four synergy scores (synergy_zip, synergy_loewe, synergy_hsa, synergy_bliss)

Column Meaning

block_id: numerical drug combination ID from drugcomb.org
DrugBank_ID_row: drug ID of row compound in DrugBank
drug_row: standardized name of row compound refer to drugbank.com
TargetGene_row: the row number of the target genes of row compound which is matching in the nodes of the network. See in subnetwork_nodes.csv
DrugBank_ID_col: drug ID of column compound in DrugBank
drug_col: standardized name of column compound refer to drugbank.com
TargetGene_col: the row number of the target genes of column compound which is matching in the nodes of the network. See in subnetwork_nodes.csv
CombineTargets: sorted and combined the row number of the target genes of both row and column compound
uniqueCombineTargets: sorted and remove duplicated of the row number of the target genes of combinations
uniqueCombineTargetGenes: converted the row number of the target genes of combinations to gene symbols
DNA_targetingByDrugs: 0;Neither drug has a DNA targeting, 1;one of the drugs has DNA targeting, 2;Both drugs target DNA. 
cell_line_name: A549 (human non-small cell lung cancer cell line)
synergy_zip: ZIP score is the % inhibition that is more than the expected effect when the two drugs do not potentiate each other
synergy_loewe: LOEWE synergy score is the % inhibition that is more than the expected effect when the two drugs are the same (drug was combined with itself)
synergy_hsa: HSA synergy score is the % inhibition that is more than expected as the highest single drug effect
synergy_bliss: BLISS synergy score is the % inhibition that is more than the expected effect when the two drugs act in probabilistic independence

• DrugCombInationID.csv
One-column list of the 607 drug-pair IDs (used as column headers)

• Network_Propagation_Result_alpha0.5.csv
Generated by Step 1; propagated profile matrix (genes × drug pairs)

• Network_Propagation_Result_alpha0.5(-log).csv
Log-transformed copy of the file above












































