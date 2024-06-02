# For genes2genes version 0.2.0

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from genes2genes import Main
from genes2genes import ClusterUtils
from genes2genes import TimeSeriesPreprocessor
from genes2genes import VisualUtils
import seaborn as sb
import pickle
import os
from optbinning import ContinuousOptimalBinning
import matplotlib.pyplot as plt

from scipy import sparse



input_dir = "/Users/rosslaidlaw/R/TrAGEDy_V2/simulated_datasets/objects/"

adata_c = ad.read_h5ad(input_dir + 'diverge_converge_c_processed.h5ad') # PAM dataset
adata_d = ad.read_h5ad(input_dir + 'diverge_converge_d_processed.h5ad') # LPS dataset

adata_c.X = adata_c.layers['logcounts']
adata_d.X = adata_d.layers['logcounts']

time_c = (adata_c.obs['slingPseudotime']-np.min(adata_c.obs['slingPseudotime']))/(np.max(adata_c.obs['slingPseudotime'])-np.min(adata_c.obs['slingPseudotime']))
time_d = (adata_d.obs['slingPseudotime']-np.min(adata_d.obs['slingPseudotime']))/(np.max(adata_d.obs['slingPseudotime'])-np.min(adata_d.obs['slingPseudotime']))

adata_c.obs['time'] = time_c
adata_d.obs['time'] = time_d


#Scale the pseudotime between 0 and 1

x = np.asarray(adata_c.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

x = np.asarray(adata_d.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

print(adata_c.obs['seurat_clusters'].value_counts())
print(adata_d.obs['seurat_clusters'].value_counts())

col = np.array(sb.color_palette('colorblind', 7))[range(7)]
joint_cmap_wt={'0':col[0], '1':col[1] , '2':col[2] , '3':col[3], '4':col[4], '5':col[5], '6':col[6]}
joint_cmap_ko={'0':col[0], '1':col[1] , '2':col[2] , '3':col[3], '4':col[4], '5':col[5], '6':col[6]}

# vs = VisualUtils.get_celltype_composition_across_time(adata_wt, adata_ko, n_points=20,
#                                                                   ANNOTATION_COLNAME='seurat_clusters', optimal_binning=False,
#                                                                   ref_cmap=joint_cmap_wt, query_cmap=joint_cmap_ko)

gene_list = np.array(pd.read_csv(input_dir + "diverge_converge_feature_space.csv"))[:,1]

output_dir = "/Users/rosslaidlaw/R/TrAGEDy_V2/simulated_datasets/plots/"

alignment_subset = "MMMMMMMMMMMMMMMMMMMM"

if __name__ == '__main__':
    aligner = Main.RefQueryAligner(adata_c, adata_d, gene_list, 20)  #
    aligner.align_all_pairs()

    aligner.get_aggregate_alignment()
    plt.savefig(output_dir + "DivConv_heatmap_genes2genes_subsetGenes.png")
    plt.show()


# alignment_subset = "MMMMMMMMMMMMMMMMMMMM"
#
# if __name__ == '__main__':
#     aligner = Main.RefQueryAligner(adata_c, adata_d, adata_c.var.index, 20)  #
#     aligner.align_all_pairs()
#
#     aligner.get_aggregate_alignment()
#     plt.savefig(output_dir + "DivConv_heatmap_genes2genes_allGenes.png")
#     plt.show()