# Spearman correlation & permutation testing workflow
# --------------------------------------------------------------------
# Purpose ▸ For each GNMF‑derived cluster ("Group") and each of four
#           synergy score metrics, identify genes whose *propagated*
#           scores strongly correlate (|ρ| ≥ 0.7) with the synergy score
#           inside that cluster, and assess whether the observed
#           correlation could arise by chance using a non‑parametric
#           permutation test (10 000 trials, incl. observed).
# Output  ▸ A data‑frame `df.report` listing every (Group, Synergy_Score,
#           Gene) triplet that meets the threshold, together with its
#           raw p‑value and FDR‑adjusted q‑value.
# --------------------------------------------------------------------

# ── 1. Load GNMF clustering assignments ---------------------------------
# Each row = drug combination; columns include a factor `Groups` (1‒14)
# and synergy scores.  Row‑names are sample IDs that must match the
# propagation matrix later.
df.info <- read.csv('GNMF_Clustering_alpha0.5_lamda0.9.csv',
                   header = TRUE, row.names = 1, check.names = FALSE)

# Prefix sample IDs with "X" so they align with the column names created
# by `read.csv` (R adds "X" when a name starts with a digit).
paste0('X', rownames(df.info)) -> rownames(df.info)

# ── 2. Load and reformat network‑propagated gene scores -------------------
# The file stores genes as rows; we transpose so genes become *columns*
# and samples (drug combinations) become *rows* to match `df.info`.
df.profile <- read.csv('Network_Popagation_Result_alpha0.5(log).csv',
                       header = TRUE, row.names = 1)
as.data.frame(t(df.profile)) -> df.profile

# Re‑order `df.profile` rows so they follow the exact order in `df.info`.
df.profile[rownames(df.info), ] -> df.profile

# ── 3. Define analysis parameters ----------------------------------------
Four_synergy_score <- c('synergy_zip', 'synergy_loewe',
                        'synergy_hsa', 'synergy_bliss')  # metrics to scan
cor.threshold      <- 0.7                                 # |ρ| cut‑off

# Container to collect all significant gene–group–score triples
# Columns: group ID, score name, gene symbol, observed ρ, p‑value, q‑value
# (q‑value will be filled in globally at the end).
df.report <- data.frame(Group               = integer(),
                        Synergy_Score       = character(),
                        Gene                = character(),
                        Observed_Correlation= numeric(),
                        p.value             = numeric(),
                        p.adjust            = numeric(),
                        stringsAsFactors    = FALSE)

# ── 4. Main nested loops --------------------------------------------------
# Outer loop ▸ each synergy score metric
# Inner loop ▸ each GNMF cluster (1‒14)
for (SynScore in 1:4) {
  synergy_score <- Four_synergy_score[SynScore]
  
  for (group in 1:14) {
    count_above_threshold <- 0      # progress counter (optional)
    cor.genes  <- c()              # indices of genes meeting threshold
    coef.genes <- c()              # their observed ρ values
    
    # ── 4a. Scan every gene column for |ρ| ≥ threshold --------------------
    for (i in 1:ncol(df.profile)) {
      coef <- cor(df.profile[which(df.info$Groups == group), i],
                  df.info    [which(df.info$Groups == group), synergy_score],
                  method = "spearman")
      
      if (abs(coef) >= cor.threshold) {
        count_above_threshold <- count_above_threshold + 1
        cor.genes  <- c(cor.genes,  i)
        coef.genes <- c(coef.genes, coef)
        
        # OPTIONAL diagnostic plot: scatter + regression line
        plot(df.profile[which(df.info$Groups == group), i],
             df.info    [which(df.info$Groups == group), synergy_score],
             main = paste(colnames(df.profile)[i],
                          "- Spearman Correlation:", round(coef, 3)),
             xlab = paste("Propagated Score(logE) of",
                           colnames(df.profile)[i],
                           "Gene in Group", group),
             ylab = paste(synergy_score, "Score"))
        fit <- lm(df.info[which(df.info$Groups == group), synergy_score] ~
                   df.profile[which(df.info$Groups == group), i])
        abline(fit, col = "blue", lwd = 2)
      }
    }
    
    # ── 4b. Permutation test for every retained gene ----------------------
    # Only when ≥1 gene passed the correlation threshold
    if (count_above_threshold > 0) {
      for (i in 1:length(cor.genes)) {
        count.greater <- 0          # how many permuted |ρ| ≥ observed |ρ|
        member.num    <- length(which(df.info$Groups == group))  # samples
        
        # Perform 9 999 random label permutations (total 10 000 incl. obs.)
        for (rep in 1:9999) {
          index <- sample(1:nrow(df.profile), size = member.num,
                          replace = FALSE)
          coef  <- cor(df.profile[index, cor.genes[i]],
                       df.info   [index, synergy_score],
                       method = "spearman")
          if (abs(coef) >= abs(coef.genes[i])) count.greater <- count.greater + 1
        }
        
        p.value <- count.greater / 10000  # permutation p‑value
        
        # Append result row to the cumulative report
        df.report <- rbind(df.report,
                           data.frame(Group               = group,
                                      Synergy_Score       = synergy_score,
                                      Gene                = colnames(df.profile)[cor.genes[i]],
                                      Observed_Correlation= coef.genes[i],
                                      p.value             = p.value,
                                      p.adjust            = NA,  # placeholder
                                      stringsAsFactors    = FALSE))
      }
    }
    
    # OPTIONAL console feedback
    print(paste0("Count in Group ", group, " with ", synergy_score,
                 ": ", count_above_threshold))
  }
}

# ── 5. Multiple‑testing correction ---------------------------------------
# Control false‑discovery rate (Benjamini‑Hochberg) across *all* tests
# so q‑values are comparable between groups and scores.
df.report$p.adjust <- p.adjust(df.report$p.value, method = 'fdr')

# ── 6. Inspect non‑significant findings (q ≥ 0.05) -----------------------
insignificant_genes <- which(df.report$p.adjust >= 0.05)
print(df.report[insignificant_genes, ])

# The final object `df.report` remains in the workspace for downstream use
# (e.g. exporting to CSV, filtering q < .05, enrichment analysis).
df.report
