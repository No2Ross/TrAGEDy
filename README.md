# TrAGEDy
Trajectory Alignment of Gene Expression Dynamics

# Introduction
While the study of an organism or biological process by itself is important, understanding the effect of the environment or conditions on an organism or biological process is also critical. Whether this be genetic, nutrient, environmental or disease conditions; understanding what genes define the differences across these conditions can help us better understand the underlying process/organism in general. 

Using TrAGEDy we can find the common alignment between the two conditions, identifying what cells are undergoing a similar biological process bewteen the two conditions.

TrAGEDy can then arrange the cells on a common pseudotime axis in order to pull out genes which are differentially expressed between the conditions at different points in the shared process.

# Method

The first step in the TrAGEDy process is to perform Trajectory Inference (TI) on the datasets. This can be down with any method which generates pseudotime values for the cells(A). We then sample gene expression at different points of pseudotime on each process by creating pseudoCells across each conditions pseudotime axis. The gene expression values of the pseudoCells are calculated from the surrounding cells, the closer a cell is in pseudotime to the pseudoCell, the more it contributes to its gene expression (B).

We then match pseudoCells between each condition that have similar gene expression profiles, leaving pseudoCells which are not part of the common underlying process unmatched (C). Using these matches, we adjust the pseudotime values of the pseudoCells, creating a common pseudotime axis which reflects the shared biological process between the two conditions (D).


