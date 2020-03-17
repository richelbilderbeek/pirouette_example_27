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
n_phylogenies <- 5
if (is_testing) {
  n_phylogenies <- 2
}

# Create simulation function
sim_dd_tree_fun <- function(crown_age) {
  extinction_rate <- 0.1
  n_taxa <- 6
  n_0 <- 2 # Initial number of species at stem/crown of tree
  diff <- (log(n_taxa) - log(n_0)) / crown_age
  speciation_rate <- 3.0 * (diff + extinction_rate)
  carrying_capacity <- n_taxa # clade-level
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species
  phylogeny
}
sim_tree_fun <- pryr::partial(
  sim_dd_tree_fun,
  crown_age = crown_age
)

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies)) {
  set.seed(314 - 1 + i)
  phylogenies[[i]] <- sim_tree_fun()
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
) + ggsave(filename = file.path(folder_name, "errors.png"), width = 7, height = 7)

# Save indiidual runs
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


