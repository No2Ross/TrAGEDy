# TrAGEDy
![](TrAGEDy_example/logo.png)
Trajectory Alignment of Gene Expression Dynamics

# Introduction to TrAGEDy method
While the study of an organism or biological process by itself is important, understanding the effect of the environment or conditions on an organism or biological process is also critical. Whether this be genetic, nutrient, environmental or disease conditions; understanding what genes define the differences across these conditions can help us better understand the underlying process/organism in general. 

Using TrAGEDy we can find the common alignment between the two conditions, identifying what cells are undergoing a similar biological process between the two conditions.

TrAGEDy can then arrange the cells on a common pseudotime axis in order to pull out genes which are differentially expressed between the conditions at different points in the shared process.

# Worked Example

To aid in the use of TrAGEDy we have layed out a worked example using a scRNA-seq datasets from Briggs et al, 2020. In this datasets there are cells from the bloodstream form stages of the parasite Tryanosoma brucei under a WT and knockout condition. Under WT conditions, the parasite transitions from a slender to stumpy form, but when the gene ZC3H20 is knocked out, this transition is blocked. 

For simplicity, we will only align the first WT replicate with the ZC3H20 KO dataset. As shown in the TrAGEDy paper, replicates can be aligned with TrAGEDy and merged into one dataset before alignment between the conditions.

#Step 1 - Create trajectories 

TrAGEDy makes no assumptions about what Trajectory Inference package is used, it only requires that the trajectories be linear in nature. For this analysis we used Slingshot from Street et al, 2018 to construct a trajectory based off of reduced dimension embeddings from PHATE (Moon et al, ). Before a trajectory is made, the user first needs to supply a vector of genes which will be used to build the trajectory and calculate gene expression dissimilarity scores between the conditions.

#Step 2 - Create interpolated points

The next stage uses cellAlign's method for creating interpolated points with some minor adjustments. First, the user decides how many interpolated points will be created across the trajectory ... A larger window size means interpolated point gene expression will be influenced by more cells further away from it and vice versa. 

As cell density is not uniform across trajectories, it is unwise to apply the same window size to all the interpolated points as this will lead to some interpolated points gene expression only being influenced by a small amount of cells, leading to a unsmooth gene expression profile. TrAGEDy avoids this, by increasing the window size around interpolated points with few cells in it's pseudotime area.


#Step 3 - Find Alignment path

We next need to find the optimal path through the data. First we find how much dissimialtity there is between the interpolated points of the conditions. The user chooses what method of calculating dissimilarity is used (Euclidean distance, Pearson, Spearman) but we recommend Spearman's correlation, as it performed best on the simulated datasets. 

A common way to find the path through a process is with Dynamic Time Warping (DTW), but DTW requires every point in the process be matched, but for our process being analysed, the WT has cells (stumpy's) that will not be represented in the ZC3H20 KO, and thus matching them together would be nonsensical. TrAGEDy avoids this issue by finding the lowest dissimilarity points at the start and end of the process and cutting out everything that falls outside of it. This means we can apply DTW fully, returing us the optimal path through the data.

#Step 4 - Perform TrajDE

Having found the optimal path, we can then identify what genes are differentially expressed (DE) between the two conditions before they diverge from one another. To do this we utilise a sliding window (akin to soft clustering) where the user defines how many windows of comparison will be made, and how many matched interpolated points will be shared between the windows as it slides across the aligned process. The user can also define what statistical test they would like to use (T test, Mann Whitney U test) and what log fold change or minimum percentage thresholds they would like to use.



