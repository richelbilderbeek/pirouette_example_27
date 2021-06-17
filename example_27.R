#
# Difference from standard:
# - Uses copy_tree for twin tree
# - Uses unlinked node-subst
#
library(pirouette)
library(beautier)
library(testthat)
library(ggplot2)

# Constants
is_testing <- is_on_ci()
example_no <- 27
rng_seed <- 314
folder_name <- paste0("example_", example_no)
crown_age <- 10
n_phylogenies <- 20
if (is_testing) {
  n_phylogenies <- 2
}

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies)) {
  set.seed(314 - 1 + i)
  phylogenies[[i]] <- create_exemplary_dd_tree(n_taxa = n_taxa, crown_age = crown_age)
}
expect_equal(length(phylogenies), n_phylogenies)

# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = length(phylogenies),
  folder_name = folder_name
)
for (i in seq_along(pir_paramses)) {
  pir_paramses[[i]]$alignment_params$sim_tral_fun <-
    get_sim_tral_with_uns_nsm_fun(
      branch_mutation_rate = 0.1,
      node_mutation_rate = 0.1
    )
  pir_paramses[[i]]$twinning_params$sim_twin_tree_fun <-
    create_copy_twtr_from_true_fun()
}
expect_equal(length(pir_paramses), n_phylogenies)
if (is_testing) {
  pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Plot total runs
pir_plots(
  pir_outs = pir_outs
) + ggtitle(paste("Number of replicates: ", n_phylogenies)); ggsave("errors.png", width = 7, height = 7)
  

# Save individual runs
expect_equal(length(pir_paramses), length(pir_outs))
expect_equal(length(pir_paramses), length(phylogenies))
for (i in seq_along(pir_outs)) {
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = dirname(pir_paramses[[i]]$alignment_params$fasta_filename)
  )
}


