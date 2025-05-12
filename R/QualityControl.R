library(phyloseq)
library(tidyverse)
library(ape)

## QC THRESHOLDS 
MIN_LIBSIZE <- 8000                 # ≥ 8 k reads
KEEP_SITE   <- "Terminal ileum"     # biopsy location to keep
EXCLUDE_ABX <- TRUE                 # drop recent antibiotic users
PREV_FILTER <- 0.005                # 0.5 % OTU prevalence
RESULT_DIR  <- "../results"         # output folder


ps <- readRDS("RISK_CCFA.rds")

## SAMPLE-LEVEL QC 
# minimum library size
ps <- prune_samples(sample_sums(ps) >= MIN_LIBSIZE, ps)

# keep only chosen biopsy site
ps <- subset_samples(ps, biopsy_location == KEEP_SITE)

#  Remove recent antibiotic users
if (EXCLUDE_ABX)
  ps <- subset_samples(ps, antibiotics == "false")

# Build metadata frame & alias
meta <- data.frame(sample_data(ps)) |>
  mutate(sample_site = biopsy_location,
         age = scale(age)[,1])            # centre & scale age

# predictors of interest
predictor_pool <- c("diagnosis","age","sex","sample_site")

# drop predictors that are single-level after QC
surviving <- predictor_pool[
  vapply(meta[predictor_pool],
         \(col) is.numeric(col) || nlevels(as.factor(col)) > 1,
         logical(1))
]

cat("Predictors kept:", paste(surviving, collapse = ", "), "\n")

# drop samples with any NA among surviving predictors
ps   <- prune_samples(complete.cases(meta[ , surviving]), ps)
meta <- meta[rownames(sample_data(ps)), ]        # aligned

cat("Samples after QC :", nsamples(ps), "\n")

## OTU PREVALENCE FILTER
prev <- apply(otu_table(ps), 1, \(x) mean(x > 0))
ps   <- prune_taxa(prev >= PREV_FILTER, ps)
cat("OTUs after filter:", ntaxa(ps), "\n")

##TREE to SPLITS 
tree <- ape::multi2di(phy_tree(ps))
get_children <- \(tr,n) tr$edge[tr$edge[,1]==n,2]
get_tips     <- \(tr,n){ k <- get_children(tr,n)
if(!length(k)) return(n)
unlist(lapply(k,get_tips,tr=tr)) }

n_tip  <- length(tree$tip.label)
raw    <- lapply((1:tree$Nnode)+n_tip, \(node){
  kids <- get_children(tree,node)
  if(length(kids)<2) return(NULL)
  list(left=get_tips(tree,kids[1]),
       right=get_tips(tree,kids[2])) })
splits <- purrr::compact(raw);  V <- length(splits)

otu <- t(as.matrix(otu_table(ps)))
N   <- nrow(otu)
y_left  <- matrix(0L, N, V)
y_total <- matrix(0L, N, V)
for(v in seq_len(V)){
  y_left[,v]  <- rowSums(otu[, splits[[v]]$left , drop=FALSE])
  y_total[,v] <- y_left[,v] +
    rowSums(otu[, splits[[v]]$right, drop=FALSE])
}

## DESIGN MATRIX & BATCH -----------------------------
form <- as.formula(paste("~", paste(surviving, collapse = " + ")))
X    <- model.matrix(form, data = meta)
P    <- ncol(X)

batch <- as.integer(droplevels(meta$run_prefix))
H     <- length(levels(droplevels(meta$run_prefix)))

dir.create(RESULT_DIR, showWarnings = FALSE)

saveRDS(list(N=N, V=V, P=P, y_left=y_left, y_total=y_total, X=X),
        file.path(RESULT_DIR, "data__qc.rds"))

