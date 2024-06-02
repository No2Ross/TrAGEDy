library(tidyverse)
library(dyngen)
library(dyno)
library(unixtools)

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(50000)

setwd("/datastore/Ross/TrAGEDy_V2")

set.seed(2)

backbone <- backbone_linear()

config <-
  initialise_model(
    backbone = backbone,
    num_cells = 1500,
    num_tfs = nrow(backbone$module_info),
    num_targets = 250,
    num_hks = 150,
    simulation_params = simulation_default(
      census_interval = 10, 
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 300),
      total_time = 500
      
    )
  )


model_common <-
  config %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()


model_wt <- model_common %>%
  generate_cells()

plot_gold_mappings(model_wt, do_facet = FALSE)

plot_backbone_modulenet(model_common)

plot_backbone_statenet(model_common)

model_common$module_network

m4_genes <- model_common$feature_info %>% filter(module_id %in% c("Burn1")) %>% pull(feature_id)
b5_genes <- model_common$feature_info %>% filter(module_id %in% c("B5")) %>% pull(feature_id)

model_ko <- model_common
model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 300,
  timepoint = 0, 
  genes = b5_genes,
  num_genes = length(b5_genes),
  multiplier = 0,
)


model_ko <- model_ko %>%
  generate_cells()

plot_backbone_modulenet(model_ko)

plot_gold_mappings(model_ko, do_facet = FALSE)

model_ko <- generate_experiment(model_ko)

plot_gold_mappings(model_ko, do_facet = FALSE)

plot_simulations(model_ko)

model_wt <- generate_experiment(model_wt)

plot_simulations(model_wt)

model_comb <-
  combine_models(list(WT = model_wt, KO = model_ko)) %>% 
  generate_experiment()

plot_gold_mappings(model_comb, do_facet = FALSE)
plot_simulations(model_comb)

plot_heatmap(model_wt, features_oi = 50, trajectory = T)

WT <- as_seurat(model_wt)
KO <- as_seurat(model_ko)

library(dyno)

saveRDS(WT, "simulation/start_shared/start_shared_wt.rds")
saveRDS(KO, "simulation/start_shared/start_shared_ko.rds")


