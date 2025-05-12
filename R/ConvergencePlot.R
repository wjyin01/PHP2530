library(phyloseq)
library(bayesplot)
library(tidyverse)

posterior <- rstan::extract(fit)
posterior_array <- as.array(fit)

# Identify top 3 significant β[2, k]

beta_draws <- posterior$beta[, 2, ]  # β for covariate 2 (diagnosis)
ci <- apply(beta_draws, 2, quantile, c(0.025, 0.975))
signif_idx <- which(ci[1, ] > 0 | ci[2, ] < 0)[1:3]

beta_pars <- paste0("beta[2,", signif_idx, "]")
phi_pars <- "phi[1]"

# Generate readable labels
tax_df <- tax_table(ps_qc) %>% as.data.frame()
tax_df$label <- ifelse(!is.na(tax_df$Genus) & tax_df$Genus != "",
                       tax_df$Genus,
                       ifelse(!is.na(tax_df$Family) & tax_df$Family != "",
                              paste0("Family: ", tax_df$Family),
                              rownames(tax_df)))
tax_label_map <- tax_df$label
names(tax_label_map) <- rownames(tax_df)

param_labels <- c()
for (k in signif_idx) {
  taxon_name <- tax_label_map[k]
  covariate <- colnames_X[2]  # assuming β[2, ] = diagnosis
  stan_name <- paste0("beta[2,", k, "]")
  readable <- paste("Effect of", covariate, "on", taxon_name)
  names(readable) <- stan_name
  param_labels <- c(param_labels, readable)
}
param_labels[phi_pars] <- "Sample 1 Dispersion (phi[1])"


# Trace Plot
mcmc_trace(posterior_array,
           pars = names(param_labels),
           facet_args = list(labeller = labeller(parameter = param_labels))) +
  theme_minimal(base_size = 11) +
  labs(title = "Trace Plots for Selected Parameters",
       y = "Parameter Value", x = "Iteration")


# ACF
mcmc_acf(posterior_array, pars = c("beta[2,1]", "phi[1]"))

