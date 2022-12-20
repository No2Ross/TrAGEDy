library(tidyverse)
library(dyngen)


set.seed(3)

backbone <- backbone_bifurcating_converging()

init <- initialise_model(
  backbone = backbone,
  num_cells = 1000,
  num_tfs = nrow(backbone$module_info),
  num_targets = 250,
  num_hks = 250,
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

model_1$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100,
  genes = c(c1_genes),
  num_genes = length(c(c1_genes)),
  multiplier = 0,
  timepoint = 0
)

model_2$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100,
  genes = c(d1_genes),
  num_genes = length(c(d1_genes)),
  multiplier = 0,
  timepoint = 0
)


model_d <- model_1 %>%
  generate_cells() %>%
  generate_experiment()


model_c <- model_2%>%
  generate_cells()%>%
  generate_experiment()


model_d_seurat <- as_seurat(model_d)
model_c_seurat <- as_seurat(model_c)


saveRDS(model_d_seurat, "TrAGEDy_tests/diverge_converge_d.rds")
saveRDS(model_c_seurat, "TrAGEDy_tests/diverge_converge_c.rds")