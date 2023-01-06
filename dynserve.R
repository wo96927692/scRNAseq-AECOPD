#'Title: dynverse 
#'Author: Han Zhifa
#'Date: 2021/12/12
#'Refer: https://dynverse.org/users/2-quick_start/
#'
#### package and dir ####
## packages
#devtools::install_github("dynverse/dyno")
#install.packages("dynplot")
#install.packages("GA")
#install.packages("shinyjs", force = T)
#dynwrap::test_singularity_installation(detailed = TRUE)
#dynwrap::test_docker_installation(detailed = TRUE)
library(dyno)
#dyn.load("/opt/R/4.1.0/lib/R/library/haven/libs/haven.so")
#is.loaded("libiconv")
#install.packages("haven")
library(tidyverse)
#install.packages("Rcpp")
library(Rcpp)
try(library(Seurat, lib.loc = "/home/hanzf/.conda/envs/scs10x/lib/R/library"), 
    silent = T)
require(Seurat, lib.loc = "/opt/R/4.1.0/lib/R/library")
library(SeuratDisk)
library(hdf5r)

## dir
setwd("/ssd/hanzf/research/scs10x/")
results_path <- "results/cluster/dynserve/"
dir.create("results/cluster/dynserve/")
vsnDir <- results_path

#### data ####
combined_dyno <- readRDS(commandArgs(T)[2])# readRDS("results/cluster/combined_celltypes.rds")
## subset
if (T) {
  combined_dyno <- combined_dyno[, combined_dyno$celltypes %in% 
                              c("NKG2C+NK", "NK")]
}

combined_dyno <- 
  combined_dyno[-grep("^RP[SL]",rownames(combined_dyno)), ]
combined_dyno <- 
  combined_dyno[-grep("^MT-",rownames(combined_dyno)), ]

dataset <- wrap_expression(
  counts = t(as.matrix(combined_dyno@assays$RNA@counts)),
  expression = t(as.matrix(combined_dyno@assays$RNA@data))
)
dataset <- add_prior_information(
  dataset,
  start_id = "G10-GTGCTGGTCTAGTCAG-1"
)


guidelines <- guidelines_shiny(dataset, port = 8790, 
                               host = "222.28.169.211")
methods_selected <- guidelines$methods_selected

model <- infer_trajectory(dataset, methods_selected[1])

# chose a method: pca, mds, tsne, ica, lle, landmark_mds, mds_sammon, mds_isomds, mds_smacof, umap, dm_diffusionMap
model <- model %>% 
  add_dimred(dyndimred::dimred(dataset$expression, method = "dm_diffusionMap", 
                               ndim = 3), 
             expression_source = dataset$expression)
save(model, file = paste0(results_path, methods_selected[1], ".rds"))

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = combined_dyno$celltypes, 
  size_cells = 1.5, alpha_cells = 0.6, 
  label_milestones = T
) +# Coloring by grouping 
scale_color_manual(values = allcolour[c(3,10)])

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  color_cells = "feature",
  feature_oi = "KLRK1",
  color_density = "grouping",
  grouping = fibroblast_reprogramming_treutlein$grouping,
  label_milestones = T, 
  alpha_cells = 0.5
) # add background color

plot_dimred(
  model, 
  expression_source = dataset$expression,
  feature_oi = "BIN2", 
  size_cells = 2.5, alpha_cells = 0.7
) # Coloring by expression

model <- model %>% 
  add_root_using_expression(c("KLRC1"), dataset$expression)

model <- label_milestones_markers(
  model,
  markers = list(
    `NK` = c("KLRC1"),
    `  NKG2C+ NK` = c("KLRC2")
  ),
  dataset$expression
)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = combined_dyno$celltypes,
  features_oi = 50
)+# Coloring by grouping 
  scale_color_manual(values = allcolour[c(3,10)])
## No features of interest provided, selecting the top 50 features automatically
## Using dynfeature for selecting the top 50 features
## Coloring by grouping

branching_milestone <- model$milestone_network %>% 
  group_by(from) %>% filter(n() > 1) %>% pull(from) %>% first()

branch_feature_importance <- calculate_branching_point_feature_importance(
  model, expression_source=dataset$expression, 
  milestones_oi = branching_milestone)

branching_point_features <- branch_feature_importance %>% 
  top_n(20, importance) %>% pull(feature_id)

plot_heatmap(
  model,
  expression_source = dataset$expression,
  features_oi = branching_point_features
)
## Coloring by milestone

space <- dyndimred::dimred_mds(dataset$expression)
map(c("KLRC1", "KLRC2"), function(feature_oi) {
  plot_dimred(model, dimred = space, 
              expression_source = dataset$expression, 
              feature_oi = feature_oi, 
              label_milestones = FALSE) +
    theme(legend.position = "none") +
    ggtitle(feature_oi)
}) %>% patchwork::wrap_plots()
