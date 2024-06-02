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

adata_ko = ad.read_h5ad(input_dir + 'start_shared_ko_processed.h5ad')
adata_wt = ad.read_h5ad(input_dir + 'start_shared_wt_processed.h5ad')

adata_wt.X = adata_wt.layers['logcounts']
adata_ko.X = adata_ko.layers['logcounts']

time_wt = (adata_wt.obs['slingPseudotime']-np.min(adata_wt.obs['slingPseudotime']))/(np.max(adata_wt.obs['slingPseudotime'])-np.min(adata_wt.obs['slingPseudotime']))
time_ko = (adata_ko.obs['slingPseudotime']-np.min(adata_ko.obs['slingPseudotime']))/(np.max(adata_ko.obs['slingPseudotime'])-np.min(adata_ko.obs['slingPseudotime']))

adata_wt.obs['time'] = time_wt
adata_ko.obs['time'] = time_ko


#Scale the pseudotime between 0 and 1

x = np.asarray(adata_wt.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

x = np.asarray(adata_ko.obs.time)
optb = ContinuousOptimalBinning(name='pseudotime', dtype="numerical")
optb.fit(x, x)
print(len(optb.splits))

print(adata_wt.obs['seurat_clusters'].value_counts())
print(adata_ko.obs['seurat_clusters'].value_counts())

col = np.array(sb.color_palette('colorblind', 7))[range(7)]
joint_cmap_wt={'0':col[0], '1':col[1] , '2':col[2] , '3':col[3], '4':col[4], '5':col[5], '6':col[6]}
joint_cmap_ko={'0':col[0], '1':col[1] , '2':col[2] , '3':col[3], '4':col[4], '5':col[5], '6':col[6]}

# vs = VisualUtils.get_celltype_composition_across_time(adata_wt, adata_ko, n_points=20,
#                                                                   ANNOTATION_COLNAME='seurat_clusters', optimal_binning=False,
#                                                                   ref_cmap=joint_cmap_wt, query_cmap=joint_cmap_ko)

gene_list = np.array(pd.read_csv(input_dir + "start_shared_feature_space.csv"))[:,1]

output_dir = "/Users/rosslaidlaw/R/TrAGEDy_V2/simulated_datasets/plots/"

alignment_subset = "MMMMMMMMMMMMMMMMMMMM"

if __name__ == '__main__':
    aligner = Main.RefQueryAligner(adata_wt, adata_ko, gene_list, 20)  #
    aligner.align_all_pairs()

    aligner.get_aggregate_alignment()
    plt.savefig(output_dir + "startShared_heatmap_genes2genes_subsetGenes.png")
    plt.show()


# alignment_subset = "MMMMMMMMMMMMMMMMMMMM"
#
# if __name__ == '__main__':
#     aligner = Main.RefQueryAligner(adata_wt, adata_ko, adata_wt.var.index, 20)  #
#     aligner.align_all_pairs()
#
#     aligner.get_aggregate_alignment()
#     plt.savefig(output_dir + "startShared_heatmap_genes2genes_allGenes.png")
#     plt.show()