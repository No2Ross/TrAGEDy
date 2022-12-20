library(tidyverse)
library(dyngen)

set.seed(2)

backbone <- backbone_linear()

config <-
  initialise_model(
    backbone = backbone,
    num_cells = 2000,
    num_tfs = nrow(backbone$module_info),
    num_targets = 300,
    num_hks = 300,
    simulation_params = simulation_default(
      census_interval = 10, 
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 300),
      total_time = 300
      
    )
  )


model_1_config <-
  config %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()

model_2_config <-
  config %>%
  generate_tf_network() %>%
  generate_feature_network() %>% 
  generate_kinetics() %>%
  generate_gold_standard()


plot_backbone_modulenet(model_1_config)

plot_backbone_statenet(model_1_config)


plot_backbone_modulenet(model_2_config)

plot_backbone_statenet(model_2_config)

model_1 <- model_1_config %>%
  generate_cells()%>%
  generate_experiment()


model_2 <- model_2_config %>%
  generate_cells() %>%
  generate_experiment()


model_1_seurat <- as_seurat(model_1)
model_2_seurat <- as_seurat(model_2)

