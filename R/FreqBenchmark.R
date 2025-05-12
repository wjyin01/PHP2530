library(phyloseq)
library(vegan)
library(DESeq2)
library(tidyverse)

# PERMANOVA (Global Differences)

otu_mat <- t(as(otu_table(ps_qc), "matrix"))  # taxa x samples
otu_mat <- t(otu_mat)  # samples x taxa
meta_df <- as(sample_data(ps_qc), "data.frame")

# Bray-Curtis distance
bc_dist <- vegdist(otu_mat, method = "bray")

# PERMANOVA with 999 permutations
adonis_result <- adonis2(bc_dist ~ diagnosis, data = meta_df, permutations = 999)
print(adonis_result)

# Wilcoxon rank-sum (Pairwise CD vs control)

# Filter for CD vs no IBD
keep_samples <- meta_df$diagnosis %in% c("CD", "no")
ps_sub <- prune_samples(keep_samples, ps_qc)
meta_sub <- as(sample_data(ps_sub), "data.frame")
otu_sub <- t(as(otu_table(ps_sub), "matrix"))

# Grouping
group <- meta_sub$diagnosis

# Compute p-values
pvals <- apply(otu_sub, 1, function(x) {
  tryCatch(wilcox.test(x ~ group)$p.value, error = function(e) NA)
})

# Adjust for FDR
padj <- p.adjust(pvals, method = "BH")

wilcox_results <- tibble(
  Taxon = rownames(otu_sub),
  p_value = pvals,
  padj = padj
) %>% arrange(padj)

print(head(wilcox_results, 10))

# DESeq2 Differential Abundance (Multivariable)

dds <- phyloseq_to_deseq2(ps_qc, ~ diagnosis + age + sex)
dds <- DESeq(dds)

# Get results for CD vs no
meta_qc <- as(sample_data(ps_qc), "data.frame")
table(meta_qc$diagnosis)  # check reference level

# Relevel diagnosis if needed
meta_qc$diagnosis <- relevel(factor(meta_qc$diagnosis), ref = "no")
sample_data(ps_qc)$diagnosis <- meta_qc$diagnosis

dds <- DESeq(phyloseq_to_deseq2(ps_qc, ~ diagnosis + age + sex))
res <- results(dds, contrast = c("diagnosis", "CD", "no"))  # CD vs no

# Summarize
res_df <- as.data.frame(res) %>%
  rownames_to_column("Taxon") %>%
  arrange(padj)

print(head(res_df, 10))
