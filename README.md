![](TrAGEDy_example/logo.png)
<p align="center"> Trajectory Alignment of Gene Expression Dynamics </p>

# Introduction to TrAGEDy method
While the study of an organism or biological process by itself is important, understanding the effect of the environment or conditions on an organism or biological process is also critical. Whether this be genetic, nutrient, environmental or disease conditions; understanding what genes define the differences across these conditions can help us better understand the underlying process/organism in general. 

Using TrAGEDy we can find the common alignment between the two conditions, identifying what cells are undergoing a similar biological process between the two conditions.

TrAGEDy can then arrange the cells on a common pseudotime axis in order to pull out genes which are differentially expressed between the conditions at different points in the shared process.

# Installing TrAGEDy

To install TrAGEDy, download the 'TrAGEDy_functions.R' script from this github and source it within your R session
```
source("/path/to/TrAGEDy_functions.R")
```

# Worked Example

To aid in the use of TrAGEDy we have layed out a worked example using a scRNA-seq datasets from Briggs et al, 2020. In this datasets there are cells from the bloodstream form stages of the parasite Tryanosoma brucei under a WT and knockout condition. Under WT conditions, the parasite transitions from a slender to stumpy form, but when the gene ZC3H20 is knocked out, this transition is blocked. 

For simplicity, we will only align the first WT replicate with the ZC3H20 KO dataset. As shown in the TrAGEDy paper, replicates can be aligned with TrAGEDy and merged into one dataset before alignment between the conditions.

We need to library some dependencies before starting TrAGEDy

```
library(slingshot)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(stats)
library(stringr)
library(rgl)
```

## Step 1 - Create trajectories 

TrAGEDy makes no assumptions about what Trajectory Inference package is used, it only requires that the trajectories be linear in nature. For this analysis we used Slingshot from Street et al, 2018 to construct a trajectory based off of reduced dimension embeddings from PHATE (Moon et al, ). Before a trajectory is made, the user first needs to supply a vector of genes which will be used to build the trajectory and calculate gene expression dissimilarity scores between the conditions.

We supply the objects with pseudotime already included, as well as the feature space gene list. While we are setting up our Zenodo, please contact us for the WT and KO files.

```
WT_sce <- readRDS("path/to/WT_sce_merged.rds")
KO_sce <- readRDS("path/to/KO_sce.rds")

features <- read.delim("path/to/WT_ZC3H20_KO_feature_space.csv", sep = ",")[,2]
```

## Step 2 - Create interpolated points

TrAGEDy uses cellAlign's method for creating interpolated points with some minor adjustments. First, the user decides how many interpolated points will be created across the trajectory and how big the window size should be. We chose 50 interpolated points for our analysis and we wanted the window size to mean adjacent interpolated points had mixes of cells being considered.

```
pseudo_end <- min(c(max(KO_sce$slingPseudotime_1, WT_sce$slingPseudotime_1)))
window <- pseudo_end / 45

WT_cell_pseudotime <- matrix(WT_sce$slingPseudotime_1, dimnames =list(WT_sce@colData@rownames))
KO_cell_pseudotime <- matrix(KO_sce$slingPseudotime_1, dimnames =list(KO_sce@colData@rownames))
WT_ID <- data.frame(WT_sce$cell_type, row.names =WT_sce@colData@rownames)
KO_ID <- data.frame(KO_sce$cell_type, row.names =KO_sce@colData@rownames)

#Create Interpolated points across pseudotime 
WT_tree <- nodePseudotime(WT_cell_pseudotime,WT_ID, 50, "WT")
KO_tree <- nodePseudotime(KO_cell_pseudotime,KO_ID, 50, "KO")
```

We then use the window size to get gene expression values for the interpolated points. A larger window size means interpolated point gene expression will be influenced by more cells further away from it and vice versa. As cell density is not uniform across trajectories, it is unwise to apply the same window size to all the interpolated points as this will lead to some interpolated points gene expression only being influenced by a small amount of cells, leading to a unsmooth gene expression profile. TrAGEDy avoids this, by increasing the window size around interpolated points with few cells in it's pseudotime area.


```
KO_cell_pseudo <- data.frame("ID" = KO_sce@colData@rownames, "pseudo" = KO_sce$slingPseudotime_1)
KO_node_pseudo <- data.frame("ID" = row.names(KO_tree$pseudotime), "pseudo" = KO_tree$pseudotime$pseudotime)

WT_cell_pseudo <- data.frame("ID" = WT_sce@colData@rownames, "pseudo" = WT_sce$slingPseudotime_1)
WT_node_pseudo <- data.frame("ID" = row.names(WT_tree$pseudotime), "pseudo" = WT_tree$pseudotime$pseudotime)

KO_node_pseudotime <- matrix(KO_tree$pseudotime$pseudotime , dimnames = list(row.names(KO_tree$pseudotime)), )
WT_node_pseudotime <- matrix(WT_tree$pseudotime$pseudotime , dimnames = list(row.names(WT_tree$pseudotime)), )

#Get gene expression values for the interpolated points
KO_node_exp_mtx <- nodeExpressionEstimate(KO_sce@assays@data@listData$logcounts, KO_node_pseudotime, KO_cell_pseudotime, window, adjust.window = T)
WT_node_exp_mtx <- nodeExpressionEstimate(WT_sce@assays@data@listData$logcounts, WT_node_pseudotime, WT_cell_pseudotime, window, adjust.window = T)
```

We then subset the gene expression matrices if the interpolated points to only include genes which are in our feature space

```
KO_node_exp_mtx  <- KO_node_exp_mtx[ features, ]
WT_node_exp_mtx <- WT_node_exp_mtx[ features, ]
```

## Step 3 - Find Alignment path

We next need to find the optimal path through the data. First we find how much dissimialtity there is between the interpolated points of the conditions. The user chooses what method of calculating dissimilarity is used (Euclidean distance, Pearson, Spearman) but we recommend Spearman's correlation, as it performed best on the simulated datasets. 

A common way to find the path through a process is with Dynamic Time Warping (DTW), but DTW requires every point in the process be matched, but for our process being analysed, the WT has cells (stumpy's) that will not be represented in the ZC3H20 KO, and thus matching them together would be nonsensical. TrAGEDy avoids this issue by finding the lowest dissimilarity points at the start and end of the process and cutting out everything that falls outside of it. This means we can apply DTW fully, returing us the optimal path through the data.

```
penalty_mtx <- dis_mtx_calculator(as.matrix(WT_node_exp_mtx), as.matrix(KO_node_exp_mtx), "spearman")

#Find optimal path through the dissimilarity matrix then cut and matches that have high dissimilarity
path_uncut <- pathfind(penalty_mtx, cut_type = "minimum", method = "mean")
```
We can then cut high dissimilarity matches that occur somewhere between the start and end points

```
path_cut <- cut_deviate(path_uncut, penalty_mtx, method = "mean")
```

We can then visualise how the alignments look through various means.

The path of the data through the matrix

```
PlotAlignment(path_cut, penalty_mtx)
```

![](TrAGEDy_example/PlotAlignment.png)

or with visualisation of the pseudotime of the interpolated points and their matches

```
PlotOutput(WT_tree, KO_tree, path_cut)
```

![](TrAGEDy_example/plotOutput_unaligned.png)


## Step 4 - Adjust pseudotime of interpolated points and cells
Using the alignment path, we can then adjust the pseudotimes of the interpolated points, with matched interpolated points having a similar pseudotime value.

```
adjustedPseudotime <- chunk_node(WT_tree$pseudotime, KO_tree$pseudotime, path_cut)

WT_tree_new <- WT_tree
WT_tree_new$pseudotime <- data.frame(adjustedPseudotime$condition_1$pseudotime, row.names = row.names(WT_tree$pseudotime) )

KO_tree_new <- KO_tree
KO_tree_new$pseudotime <- data.frame(adjustedPseudotime$condition_2$pseudotime, row.names = row.names(KO_tree$pseudotime) )

```

We can then visualise the changes in pseudotime based on the alignment path

```
PlotOutput(WT_tree_new, KO_tree_new , path_cut)
```

![](TrAGEDy_example/plotOutput_aligned.png)


Finally we adjust the pseudotime of the individual cells based on the new pseudotimes of the interpolated points
```
KO_cell_pseudo <- data.frame(KO_sce$slingPseudotime_1, row.names = KO_sce@colData@rownames)
KO_node_pseudo_new <- data.frame(row.names(adjustedPseudotime$condition_2),adjustedPseudotime$condition_2$pseudotime, adjustedPseudotime$condition_2$alignment, row.names =row.names(adjustedPseudotime$condition_2) )
KO_node_pseudo <- data.frame( pseudotime = KO_tree$pseudotime, row.names = row.names(KO_tree$pseudotime))

KO_cell_pseudo_new <- pseudo_cell_align(KO_cell_pseudo , KO_node_pseudo_new, KO_node_pseudo, window)
KO_cell_pseudo_new <- KO_cell_pseudo_new[order(match(row.names(KO_cell_pseudo_new), row.names(KO_cell_pseudo))),]

WT_cell_pseudo <- data.frame(WT_sce$slingPseudotime_1, row.names = WT_sce@colData@rownames)
WT_node_pseudo_new <- data.frame(row.names(adjustedPseudotime$condition_1),adjustedPseudotime$condition_1$pseudotime, adjustedPseudotime$condition_1$align, row.names =row.names(adjustedPseudotime$condition_1) )
WT_node_pseudo <- data.frame(pseudotime = WT_tree$pseudotime, row.names = row.names(WT_tree$pseudotime))

WT_cell_pseudo_new <- pseudo_cell_align(WT_cell_pseudo , WT_node_pseudo_new, WT_node_pseudo, window)
WT_cell_pseudo_new <- WT_cell_pseudo_new[order(match(row.names(WT_cell_pseudo_new), row.names(WT_cell_pseudo))),]
```
We save the results back into the object
```
KO_sce$oldPseudotime <- KO_sce$slingPseudotime_1
KO_sce$newPseudotime <- KO_cell_pseudo_new$pseudotime
WT_sce$oldPseudotime <- WT_sce$slingPseudotime_1
WT_sce$newPseudotime <- WT_cell_pseudo_new$pseudotime
KO_sce$Status <- as.factor(KO_cell_pseudo_new$status)
WT_sce$Status <- as.factor(WT_cell_pseudo_new$status)
```

## Step 5 - Perform TrajDE

Having found the optimal path, we can then identify what genes are differentially expressed (DE) between the two conditions before they diverge from one another. To do this we utilise a sliding window (akin to soft clustering) where the user defines how many windows of comparison will be made, and how many matched interpolated points will be shared between the windows as it slides across the aligned process. The user can also define what statistical test they would like to use (T test, Mann Whitney U test) and what log fold change or minimum percentage thresholds they would like to use.
```
output <- TrajDE(list(WT_sce, KO_sce), list(WT_node_pseudo, KO_node_pseudo), path_cut, n_windows = 4, overlap = 1, p_val = 0.05, min.pct = 0.1, logfc = 0.5, all.genes = F, test_use = "wilcox", correct = T)
```




