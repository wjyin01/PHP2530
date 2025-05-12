library(phyloseq)
library(vegan)
library(DESeq2)
library(tidyverse)

# PERMANOVA (Global Differences)

otu_mat <- t(as(otu_table(ps_qc), "matrix"))  # taxa x samples
otu_mat <- t(otu_mat)  # samples x taxa
meta_df <- data.frame(sample_data(ps_qc)) %>%
  mutate(diagnosis = factor(diagnosis))


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


# Visualization
library(ggplot2)

# PCoA on Bray-Curtis
ordination <- ordinate(ps_qc, method = "PCoA", distance = "bray")

# Plot
plot_ordination(ps_qc, ordination, color = "diagnosis") +
  geom_point(size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCoA (Bray-Curtis) of Microbiome Composition",
       color = "Diagnosis")

  
# Create effect size: log fold change between medians
group_vals <- sample_data(ps_qc)$diagnosis
otu_df <- as.data.frame(t(otu_table(ps_qc)))

group1 <- "CD"
group2 <- "no"
keep <- group_vals %in% c(group1, group2)

logfc <- apply(otu_df[, keep], 1, function(x) {
  log2(median(x[group_vals[keep] == group1]) + 1) - 
    log2(median(x[group_vals[keep] == group2]) + 1)
})

wilcox_df <- wilcox_results %>%
  filter(group1 == group1, group2 == group2) %>%
  mutate(log2FC = logfc[match(taxon, names(logfc))],
         sig = padj < 0.05)

# Volcano plot
ggplot(wilcox_df, aes(x = log2FC, y = -log10(pval), color = sig)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Wilcoxon Volcano Plot: CD vs Control",
       x = "Log2 Fold Change (Median)",
       y = "-log10(p-value)",
       color = "Significant")
