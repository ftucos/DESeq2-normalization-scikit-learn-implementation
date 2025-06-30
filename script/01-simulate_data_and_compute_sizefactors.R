#=====================================================================
#  Simulating an RNA‑seq count matrix and normalising it with DESeq2  
#--------------------------------------------------------------------
#  This script                                                     
#    1. Simulates a toy RNA‑seq experiment (20 samples × 2,000 genes) 
#    2. Splits the data into a **training** and a **test** set         
#    3. Calculates size factors with a *custom* version of           
#       DESeq2::estimateSizeFactorsForMatrix()                       
#    4. Demonstrates how to normalise new (test) samples so that     
#       they are on the same scale as the training data             
#                                                                     
#  The only deviation from the original DESeq2 function is that      
#  the final rescaling step (which forces the geometric mean of      
#  the size factors to 1) has been **disabled**.  Keeping the size    
#  factors on the library‑size scale is helpful when normalising     
#  additional samples later on.                                      
#=====================================================================


# 0. Setup
# ------------------------------
library(DESeq2)

set.seed(17)                                 # ensure reproducibility
source("R/custom_estimate_size_factors.R")   # custom version of estimateSizeFactorsForMatrix()


# ------------------------------
# 1. Simulate a toy dataset
# ------------------------------
# Parameters
n_genes   <- 2000          # number of genes
n_samples <- 20            # number of samples

# 1.1 Simulate sequencing depth for each sample
#     Depths range from 5 million to 50 million reads
library_sizes <- round(runif(n_samples, min = 5e6, max = 50e6))
names(library_sizes) <- paste0("Sample", seq_len(n_samples))

# 1.2 Simulate baseline gene expression
#     Draw gene‑specific means from a Gamma distribution to add heterogeneity
gene_means <- rgamma(n_genes, shape = 2, scale = 1)
names(gene_means) <- paste0("Gene", seq_len(n_genes))

#     Convert means to *relative* expression so that they sum to 1
#     These will act as per‑gene proportions
gene_props <- gene_means / sum(gene_means)

# 1.3 Generate the raw count matrix
#     For each sample, counts_ij ~ Poisson(lambda = gene_prop_i × library_size_j)
count_matrix <- sapply(library_sizes, function(depth) {
  rpois(n_genes, lambda = gene_props * depth)
})
rownames(count_matrix) <- names(gene_means)
colnames(count_matrix) <- names(library_sizes)

# Inspect simulated data -------------------------------------------------
counts_df   <- as.data.frame(count_matrix)

write.csv(counts_df, "data/simulated_counts.csv", row.names = TRUE)

total_counts <- colSums(counts_df)
summary_df   <- data.frame(Sample      = names(library_sizes),
                           LibrarySize = library_sizes,
                           TotalCounts = total_counts)

print(summary_df)   # quick sanity‑check

# ------------------------------
# 2. Train / test split
# ------------------------------
counts_train <- count_matrix[, 1:10]  # training set
counts_test  <- count_matrix[, 11:20] # test set

# ------------------------------
# 3. Normalise the training set (median‑ratio method)
# ------------------------------
# custom_estim_size_factors() is a trimmed‑down version of
# DESeq2::estimateSizeFactorsForMatrix() *without* the final rescaling
# that forces the geometric mean of the size factors to 1.  This keeps
# the totals sizes of new samples closer to the sizes of the training data

# ratio is the default approach that excludes from normalization all genes with 0 count in at least one sample
# poscounts instead handles zeros by ignoring them in the log geometric mean calculation
size_factors_train <- custom_estim_size_factors(counts_train, type = "ratio")
size_factors_train
# export them in results
write.csv(data.frame(Sample = names(size_factors_train),
                      SizeFactor = size_factors_train),
          "data/size_factors_train.csv", row.names = FALSE)

# Geometric mean of size factors (should be close but not necessary one)
exp(mean(log(size_factors_train)))

# 3.1 Apply normalisation: counts_ij / size_factor_j
norm_counts_train <- t(t(counts_train) / size_factors_train)
write.csv(as.data.frame(norm_counts_train),
          "data/norm_counts_train.csv", row.names = TRUE)

colSums(norm_counts_train)  # library sizes now comparable across samples

# Visual check: size factors vs. original library sizes (in millions of reads)
plot(size_factors_train, colSums(counts_train) / 1e6,
     xlab = "Size factor", ylab = "Library size (millions)")

# ------------------------------
# 4. Sanity check – provide an explicit pseudo‑reference sample on the same training data
# ------------------------------
# Build the pseudo‑reference sample as the geometric mean across training samples
# no need to handle zeros because we are using the default ratio method
reference_pseudosample <- exp(rowMeans(log(counts_train)))

size_factors_train_check <- custom_estim_size_factors(counts_train,
                                                      type     = "ratio",
                                                      geoMeans = reference_pseudosample)
# DESeq2::estimateSizeFactorsForMatrix() would output exactly 1 when geoMeans is supplied
exp(mean(log(size_factors_train_check)))

# Confirm that both approaches give the same normalised counts
norm_counts_train_check <- t(t(counts_train) / size_factors_train_check)
all(norm_counts_train_check == norm_counts_train)  # TRUE

# ------------------------------
# 5. Normalise a single sample from the training set (one‑off normalisation)
# ------------------------------
# Extract the first training sample as a one‑column matrix
counts_train_sample_1 <- as.matrix(counts_train[, 1])

size_factor_train_sample_1 <- custom_estim_size_factors(counts_train_sample_1,
                                                        type     = "ratio",
                                                        geoMeans = reference_pseudosample)

counts_normalised_train_sample_1 <- counts_train_sample_1 / size_factor_train_sample_1

# Verify match with previously normalised counts
all(counts_normalised_train_sample_1 == norm_counts_train[, 1])  # TRUE

# ------------------------------
# 7. Normalise the test set so it is on the same scale
# ------------------------------
size_factors_test <- custom_estim_size_factors(counts_test,
                                               type     = "ratio",
                                               geoMeans = reference_pseudosample)
size_factors_test
write.csv(data.frame(Sample = names(size_factors_test),
                     SizeFactor = size_factors_test),
          "data/size_factors_test.csv", row.names = FALSE)

# Geometric mean will again differ from 1 (as expected)
exp(mean(log(size_factors_test)))

# Apply normalisation to the test counts
norm_counts_test <- t(t(counts_test) / size_factors_test)
write.csv(as.data.frame(norm_counts_test),
          "data/norm_counts_test.csv", row.names = TRUE)

cat("Test Samples normalized library sizes:\n")
colSums(norm_counts_test)  # should be comparable to training libraries
cat("Train Samples normalized library sizes:\n")
colSums(norm_counts_train)  

# Size factors should still reflect the underlying library sizes
plot(size_factors_test, colSums(counts_test) / 1e6,
     xlab = "Size factor (test)", ylab = "Library size (millions)")
