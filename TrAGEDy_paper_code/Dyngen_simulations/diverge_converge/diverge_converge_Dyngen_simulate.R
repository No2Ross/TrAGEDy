library(tidyverse)
library(dyngen)

library(unixtools)

#set.tempdir("pooling-scRNA/")
ulimit::memory_limit(50000)

setwd("/datastore/Ross/TrAGEDy_V2")


set.seed(5)

backbone <- backbone_bifurcating_converging()

init <- initialise_model(
  backbone = backbone,
  num_cells = 1500,
  num_tfs = nrow(backbone$module_info),
  num_targets = 250,
  num_hks = 150,
  simulation_params = simulation_default(census_interval = 1, ssa_algorithm = ssa_etl(tau = 300 / 3600)),
  verbose = FALSE
)

model_common <-
  init %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()


plot_backbone_modulenet(model_common)

plot_backbone_statenet(model_common)

model_1 <- model_common
model_2 <- model_common

c1_genes <- model_1$feature_info %>% filter(module_id %in% c("C1")) %>% pull(feature_id)
d1_genes <- model_2$feature_info %>% filter(module_id %in% c("D1")) %>% pull(feature_id)
f1_genes <- model_common$feature_info %>% filter(module_id %in% c("F1")) %>% pull(feature_id)

model_1$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100,
  genes = c(c1_genes),
  num_genes = length(c(c1_genes)),
  multiplier = 0,
  timepoint = 0
)

# model_1$simulation_params$experiment_params <- simulation_type_knockdown(
#   num_simulations = 100,
#   genes = f1_genes,
#   num_genes = length(f1_genes),
#   multiplier = 0,
#   timepoint = 0
# )
# 
model_2$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100,
  genes = c(d1_genes),
  num_genes = length(c(d1_genes)),
  multiplier = 0,
  timepoint = 0
)

# model_2$simulation_params$experiment_params <- simulation_type_knockdown(
#   num_simulations = 100,
#   genes = f1_genes,
#   num_genes = length(f1_genes),
#   multiplier = 0,
#   timepoint = 0
# )

model_d <- model_1 %>%
  generate_cells() %>%
  generate_experiment()


model_c <- model_2%>%
  generate_cells()%>%
  generate_experiment()

model_common <- model_common %>%
  generate_cells() %>%
  generate_experiment()

model_d_seurat <- as_seurat(model_d)
model_c_seurat <- as_seurat(model_c)

model_comb <-
  combine_models(list(first = model_d, second = model_c)) %>% 
  generate_experiment()

dyno_comb <- as_dyno(model_common)

dynplot::plot_dimred(dyno_comb)
dynplot::plot_heatmap(dyno_comb)

dyno_d <- as_dyno(model_d)

dynplot::plot_dimred(dyno_d)
dynplot::plot_heatmap(dyno_d)

dyno_c <- as_dyno(model_c)

dynplot::plot_dimred(dyno_c)
dynplot::plot_heatmap(dyno_c)

saveRDS(model_d_seurat, "simulation/div_con/diverge_converge_d.rds")
saveRDS(model_c_seurat, "simulation/div_con/diverge_converge_c.rds")


