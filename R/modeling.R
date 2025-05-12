library(phyloseq)
library(rstan)
library(posterior)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = 4)

cat(Sys.time(), ": Starting script\n")

qc <- readRDS("../results/data__qc.rds")

X <- qc$X
N <- qc$N
sample_ids <- rownames(X)

# Load phyloseq object and apply same sample filter
ps <- readRDS("RISK_CCFA.rds")
ps_qc <- prune_samples(sample_ids, ps)

# Ccollapse to top 100 taxa
top_taxa <- names(sort(taxa_sums(ps_qc), decreasing = TRUE))[1:100]
ps_qc <- prune_taxa(top_taxa, ps_qc)

# Count matrix
Y <- t(as.matrix(otu_table(ps_qc)))

# Build Stan input
stan_data <- list(
  N = nrow(Y),
  K = ncol(Y),
  P = ncol(X),
  Y = Y,
  X = X
)

# Fit
cat(Sys.time(), ": Compiling model\n")
sm <- stan_model("dm_regression.stan")

cat(Sys.time(), ": Sampling...\n")
fit <- sampling(sm, data = stan_data, chains = 4, iter = 1000, warmup = 500, refresh = 100)

# ------------ 3. Save results ------------------------------
cat(Sys.time(), ": Saving fit\n")
saveRDS(fit, "../results/fit_dm_regression_qc.rds")

cat(Sys.time(), ": Done!\n")
