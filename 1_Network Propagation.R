###############################################################
# Network Propagation via Random Walk with Restart (RWR)
# -------------------------------------------------------------
# Purpose:
#   Propagate drug effect signals from their target genes across a cancer‑specific interaction network
#
# Core equation implemented:
#   xT = α [ I − (1 − α) A D⁻¹ ]⁻¹  x0
#
# where
#   x0 : initial vector (encoded drug targets) (pre-propagate profile)
#   A  : adjacency matrix of undirected cancer sub‑network
#   D  : diagonal out‑degree matrix (used to row‑normalize A)
#   α  : restart probability (trade‑off between restart & walk)
#   I  : identity matrix
#   xT : the final state of the network after the propagation (propagated profile)
#
###############################################################

# ---------------------------#
# 1. Load required packages  #
# ---------------------------#
library(igraph, quietly = TRUE)  # graph objects + adjacency helpers
library(readxl)                  # read Excel (drug‑pair metadata)
library(expm)                    # advanced linear‑algebra functions
                                 # (solve() from base is enough, but expm is loaded
                                 #  in case matrix exponentials are later required)

# -----------------------------------------#
# 2. Build the cancer sub‑network (matrix) #
# -----------------------------------------#

# Load edge list (source, target) for cancer‑specific PPI sub‑network
subnetwork <- read.csv("filesData/gene_interaction_cancer_subnetwork.csv", header = FALSE)

# Convert edge list → igraph object (undirected)
subnetwork <- graph_from_data_frame(subnetwork, directed = FALSE)

# Extract binary adjacency matrix A (|V| × |V|)
A <- get.adjacency(subnetwork)

# -------------------------------------------------------------#
# 3. Load drug‑pair information and list of network node names #
# -------------------------------------------------------------#

df.drugcomb <- read_xlsx(path = "filesData/Drug_Combination_A549.xlsx")  # drug A/B & their targets
DrugCombID <- read.csv("filesData/DrugCombInationID.csv", header = F)
nodes <- read.csv("filesData/subnetwork_nodes.csv", header = TRUE)       # node IDs / gene symbols

# -----------------------------------------------------------------------------#
# 4. Allocate empty results dataframe (rows: nodes, cols: drug combinations)   #
# -----------------------------------------------------------------------------#

df.empty <- data.frame(matrix(
  ncol = nrow(df.drugcomb),  # one column per drug pair
  nrow = nrow(A)             # one row per network node
))

# Use drug‑pair ID (column 1 of df.drugcomb) as column names
colnames(df.empty) <- DrugCombID[, 1]

# Use gene / node identifiers as row names
rownames(df.empty) <- nodes[, 1]

# ------------------------------------------------------------------------#
# 5. Helper function to compute out‑degree for every node in adjacency A  #
# ------------------------------------------------------------------------#
outdegree <- function(A) {
  # Sums each row to count outgoing edges
  apply(A, 1, sum)
}

# Build diagonal matrix D of out‑degrees
D <- diag(outdegree(A), nrow = nrow(A), ncol = nrow(A))

# Identity matrix of same size
I <- diag(nrow(A))

# Set the restart probability (α)
alpha <- 0.5

# ---------------------------------------------------------#
# 6. MAIN LOOP: propagate every drug combination sequentially
# ---------------------------------------------------------#
for (loop in 1:nrow(df.drugcomb)) {

  # -----------------------------------------------------#
  # 6.1 Create initial vector x0 for this drug pair      #
  # -----------------------------------------------------#

  # x0a, x0b will hold partial contributions of drug A and drug B, respectively
  x0a <- matrix(0, nrow = nrow(A), ncol = 1)
  x0b <- matrix(0, nrow = nrow(A), ncol = 1)

  # ---------------- Drug A targets ---------------- #
  # Column 4 of df.drugcomb contains comma‑separated node indices
  CBTDM <- paste(df.drugcomb[loop, 4], sep = ",")
  STDM  <- strsplit(CBTDM, split = ",")
  MTCTDM <- matrix(unlist(STDM))
  nMT <- nrow(MTCTDM)  # number of target nodes for drug A

  # Evenly distribute 0.5 across all targets of drug A
  for (nmt in 1:nMT) {
    x0a[strtoi(MTCTDM[nmt, 1]), 1] <- (0.5 / nMT)
  }

  # ---------------- Drug B targets ---------------- #
  # Column 7 holds the targets for drug B in identical format
  CBTDM <- paste(df.drugcomb[loop, 7], sep = ",")
  STDM  <- strsplit(CBTDM, split = ",")
  MTCTDM <- matrix(unlist(STDM))
  nMT <- nrow(MTCTDM)  # number of target nodes for drug B

  # Evenly distribute the remaining 0.5 across drug B targets
  for (nmt in 1:nMT) {
    x0b[strtoi(MTCTDM[nmt, 1]), 1] <- (0.5 / nMT)
  }

  # Combine both drugs to form full initial vector
  x0 <- x0a + x0b

  # -----------------------------------------------------#
  # 6.2 Random Walk with Restart
  # -----------------------------------------------------#

  # Compute normalized adjacency A·D⁻¹
  normalized_A <- A %*% solve(D)

  # Closed‑form solution for steady‑state vector:
  #   xT = α [ I - (1 - α) normalized_A ]⁻¹  x0
  xT <- (alpha * solve(I - ((1 - alpha) * normalized_A))) %*% x0

  # -----------------------------------------------------#
  # 6.3 Store propagated profile and show progress       #
  # -----------------------------------------------------#
  df.empty[, loop] <- as.vector(xT)

  print(paste("Please wait... ", loop, "/", nrow(df.drugcomb)))
}

# -----------------------------------------------------------#
# 7. Write complete results matrix to disk for downstream use
# -----------------------------------------------------------#

write.csv(df.empty, "Network_Popagation_Result_alpha0.5.csv", row.names=T)

###############################################################
# End of Network Propagation
###############################################################
