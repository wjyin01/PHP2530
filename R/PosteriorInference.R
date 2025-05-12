# Summarize posterior draws for diagnosis effect
beta_draws <- posterior$beta[, 2, ]  # diagnosis effect (P=2)
beta_summary <- apply(beta_draws, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>% as.data.frame() %>%
  rename(lower = `2.5%`, median = `50%`, upper = `97.5%`)

# Attach taxon IDs and taxonomy labels
otu_mat <- as(otu_table(ps_qc), "matrix")
if (taxa_are_rows(ps_qc)) otu_mat <- t(otu_mat)
top100_taxa <- names(sort(colSums(otu_mat), decreasing = TRUE))[1:100]
beta_summary$TaxonID <- top100_taxa

tax_map <- tax_table(ps_qc) %>% as.data.frame() %>% rownames_to_column("Taxon")
tax_map$label <- ifelse(!is.na(tax_map$Genus) & tax_map$Genus != "",
                        tax_map$Genus,
                        ifelse(!is.na(tax_map$Family),
                               paste0("Family: ", tax_map$Family),
                               tax_map$Taxon))
beta_summary <- beta_summary %>%
  left_join(tax_map[, c("Taxon", "label")], by = c("TaxonID" = "Taxon")) %>%
  mutate(taxon = paste0(label, " (", TaxonID, ")"))

# Forest plot
ggplot(beta_summary, aes(x = median, y = fct_reorder(taxon, median))) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Diagnosis Effect on Taxa",
       x = "Posterior Median (95% CI)", y = "Taxon")



significant_taxa <- beta_summary %>%
  filter(lower > 0 | upper < 0) %>%
  arrange(desc(abs(median)))
print(significant_taxa[, c("taxon", "median", "lower", "upper")])



theta_hat_mean <- apply(posterior$theta_hat, c(2, 3), mean)  # N x K
colnames(theta_hat_mean) <- top100_taxa
theta_df <- as.data.frame(theta_hat_mean)
theta_df$diagnosis <- sample_data(ps_qc)$diagnosis

# Select top 5 most abundant taxa
top5 <- names(sort(colMeans(theta_df[, top100_taxa]), decreasing = TRUE))[1:5]

# Get readable labels
tax_label_map <- tax_map[top5, "label"]; names(tax_label_map) <- top5

theta_long <- theta_df %>%
  pivot_longer(cols = all_of(top5), names_to = "Taxon", values_to = "Theta") %>%
  mutate(TaxonLabel = tax_label_map[Taxon])

# Boxplot
ggplot(theta_long, aes(x = diagnosis, y = Theta)) +
  geom_boxplot() +
  facet_wrap(~TaxonLabel, scales = "free_y") +
  theme_minimal()
