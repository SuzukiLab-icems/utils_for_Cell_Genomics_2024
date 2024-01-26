import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as pl
from matplotlib import rcParams

"""Data Loading"""
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=300, facecolor='white', format='eps')

results_raw = 'result/results_raw.h5ad'
results_normalized = 'result/results_normalized.h5ad'
results_TSNE = 'result/results_TSNE.h5ad'
results_DEGs = 'result/results_DEGs.h5ad'
results_Summary = 'result/results_Summary.h5ad'

adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    var_names='gene_symbols',
    cache=True
)

adata.write(results_raw)

"""preProcessing"""
sc.pl.highest_expr_genes(adata, n_top=20, save=True, show=False)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=True, show=False)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt' ,save=True, show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts' ,save=True, show=False)
sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', save=True, show=False)

adata = adata[adata.obs.n_genes_by_counts < 3000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata, save=True, show=False)

adata.raw = adata

adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

"""PCA Analysis"""
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='Rd3')
sc.pl.pca_variance_ratio(adata, log=True, save=True, show=False)

adata.write(results_normalized)

"""Clustering with tSNE"""
sc.pp.neighbors(adata, n_neighbors=7, n_pcs=30)
sc.tl.tsne(adata)
sc.tl.leiden(adata)
sc.pl.tsne(adata, color=['leiden', 'Rd3', 'Rd3l'], save=True, show=False, edges=False)

adata.write(results_TSNE)

"""DEGs Extraction"""
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, save=True, show=False)

adata.write(results_DEGs)

"""Cell Type Annotation"""
sc.tl.leiden(adata, resolution=1.0)
#Rod: 1,4: Rho
#Cone: 11: Opn1mw
#Bipolar: 2,3,6,9,10: Grm6, Isl1
#Horizontal Cells: 5,: Slc32a1
#AMACRINE CELLS: 15,17,20: Gad1
#MULLER CELLS and Astrocytes: 0,14,21: Slc1a3, Glul, Gfap
#ENDOTHELIAL CELLS: 13: Cldn5
#Pericytes: 8: Acta2
#Fibroblasts: 16: Col1a1
#Microglia: 12: C1qa, P2ry12
#Monocyte: 7, 19: Lgals3, Ms4a7
#T cells: 18: Cd3g, Cd3d, Cd3e
####MArkerList####
marker=[
'Rho',
'Opn1mw',
'Grm6',
'Isl1',
'Slc32a1',
'Gad1',
'Slc1a3',
'Glul',
'Gfap',
'Cldn5',
'Acta2',
'Col1a1',
'C1qa',
'P2ry12',
'Lgals3',
'Ms4a7',
'Cd3g',
'Cd3d',
'Cd3e'
]
#####
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden'], save=False , show=False)

adata.obs['leiden'].cat.categories
adata.obs['leiden_anno'] = adata.obs['leiden']
adata.uns['leiden_anno_colors'] = adata.uns['leiden_colors']

adata.obs['leiden_anno'].cat.categories = [
'Mullar glia/Astrocytes.1',     	#0
'Rods.1',      						#1
'Bipolar cells.1',    				#2
'Bipolar cells.2',      			#3
'Rods.2',      						#4
'Horizontal cells',     			#5
'Bipolar cells.3',   				#6
'Monocytes.1',     					#7
'Pericytes',      					#8
'Bipolar cells.4',      			#9
'Bipolar cells.5',     				#10
'Cones',   							#11
'Microglia',     					#12
'Endothelial cells',      			#13
'Mullar glia/Astrocytes.2',      	#14
'Amacrine cells.1',     			#15
'Fibroblasts',   					#16
'Amacrine cells.2',      			#17
'T cells',     						#18
'Monocytes.2',     					#19
'Amacrine cells.3',     			#20
'Mullar glia/Astrocytes.3',        	#21
]

list = (
'Rods.1',      						#1
'Rods.2',      						#4
'Cones',   							#11
'Amacrine cells.1',     			#15
'Amacrine cells.2',      			#17
'Amacrine cells.3',     			#20
'Bipolar cells.1',    				#2
'Bipolar cells.2',      			#3
'Bipolar cells.3',   				#6
'Bipolar cells.4',      			#9
'Bipolar cells.5',     				#10
'Horizontal cells',     			#5
'Mullar glia/Astrocytes.1',     	#0
'Mullar glia/Astrocytes.2',      	#14
'Mullar glia/Astrocytes.3',        	#21
'Microglia',     					#12
'Monocytes.1',     					#7
'Monocytes.2',     					#19
'T cells',     						#18
'Endothelial cells',      			#13
'Fibroblasts',   					#16
'Pericytes',      					#8
)
adata.obs['leiden_anno_ordered'] = adata.obs['leiden_anno'].cat.reorder_categories(list, inplace=True)

zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['leiden_anno_colors'])
new_colors[[0,1]] = zeileis_colors[[11,11]] #Rods
new_colors[[2]] = zeileis_colors[[23]] #Cones
new_colors[[3,4,5]] = zeileis_colors[[8,8,8]] #Amacrine Cells
new_colors[[6,7,8,9,10]] = zeileis_colors[[6,6,6,6,6]] #Bipolar cells
new_colors[[11]] = zeileis_colors[[0]] #Horizontal cells
new_colors[[12,13,14]] = zeileis_colors[[26,26,26]] #Mullar glia/Astrocytes
new_colors[[15]] = zeileis_colors[[19]] #Microglia
new_colors[[16,17]] = zeileis_colors[[18,18]] #Monocyte
new_colors[[18]] = zeileis_colors[[12]] #T cells
new_colors[[19]] = zeileis_colors[[24]] #Endothelial cells
new_colors[[20]] = zeileis_colors[[27]] #Fibroblasts
new_colors[[21]] = zeileis_colors[[25]] #Pericytes
adata.uns['leiden_anno_colors'] = new_colors
    
sc.pl.tsne(adata, color=['leiden_anno', 'Rd3'], save=True, show=False, edges=False) # for Figure.S4B

adata.write(results_Summary)

"""Generate Expression Matrix with Trajectory Analysis"""
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden_anno']  == 'Rods.1')[0]
sc.tl.dpt(adata)
gene_names = [
'Rd3',
'Rd3l',
'Rho',#Rod: 1,4:
'Opn1mw',#Cone: 11:
'Grm6', 'Isl1',#Bipolar: 2,3,6,9,10:
'Slc32a1',#Horizontal Cells: 5,:
'Gad1',#AMACRINE CELLS: 15,17,20:
'Slc1a3', 'Glul', 'Gfap',#MULLER CELLS and Astrocytes: 0,14,21:
'Cldn5',#ENDOTHELIAL CELLS: 13:
'Acta2',#Pericytes: 8:
'Col1a1',#Fibroblasts: 16:
'P2ry12',#Microglia: 12: C1qa,
'Lgals3', 'Ms4a7',#Monocyte: 7, 19:
'Cd3g', 'Cd3d', 'Cd3e',#T cells: 18:
]
paths = [('retinal cells', [
'Rods.1',      						#1->0
'Rods.2',      						#4->1
'Cones',   							#11->2
'Amacrine cells.1',     			#15->3
'Amacrine cells.2',      			#17->4
'Amacrine cells.3',     			#20->5
'Bipolar cells.1',    				#2->6
'Bipolar cells.2',      			#3->7
'Bipolar cells.3',   				#6->8
'Bipolar cells.4',      			#9->9
'Bipolar cells.5',     				#10->10
'Horizontal cells',     			#5->11
'Mullar glia/Astrocytes.1',     	#0->12
'Mullar glia/Astrocytes.2',      	#14->13
'Mullar glia/Astrocytes.3',        	#21->14
'Microglia',     					#12->15
'Monocytes.1',     					#7->16
'Monocytes.2',     					#19->17
'T cells',     						#18->18
'Endothelial cells',      			#13->19
'Fibroblasts',   					#16->20
'Pericytes',      					#8->21
]
)]

sc.pl.draw_graph(adata, color=['leiden_anno', 'dpt_pseudotime'], legend_loc='on data', save=True)

adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['leiden_anno']
adata.uns['clusters_colors'] = adata.uns['leiden_anno_colors']

_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)

"""Generate Expression Matrix of Marker Genes for Figure.S4C"""
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50, #Smoothing (If you want to obtain raw count data matrix, pls use n_avg=1. But, for generating co-expression network, smoothing process is important.)
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='clusters',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('./figures/Figure.S4B-S4C_retina_marker_gene_exp_matrix.csv'.format(descr))

"""Visualize Expression of Marker Genes for Figure.S4C"""
data = pd.read_csv('./figures/Figure.S4B-S4C_retina_marker_gene_exp_matrix.csv',index_col=0)
data["cluster"] = np.floor(data["groups"])

label = pd.read_csv('./figures/Figure.S4C_Cell_Type_Annotation_Table.csv') #<<Should create by yourself
anno_data = label.merge(data, how="left", on="cluster")

anno_data.to_csv('./figures/Figure.S4C_Cell_Type_Annotated_Exp_Matrix_of_Retinal_Cell_Marker_Genes.csv')

order = (
'Rods',
'Cones',
'Amacrine cells',
'Bipolar cells',
'Horizontal cells',
'Mullar glia/Astrocytes',
'Microglia',
'Monocytes',
'T cells',
'Endothelial cells',
'Fibroblasts',
'Pericytes',
)

marker_order = [
'Rd3',
'Rho',#Rod: 1,4:
'Opn1mw',#Cone: 11:
'Gad1',#AMACRINE CELLS: 15,17,20:
'Grm6', 'Isl1',#Bipolar: 2,3,6,9,10:
'Slc32a1',#Horizontal Cells: 5,:
'Slc1a3', 'Glul', 'Gfap',#MULLER CELLS and Astrocytes: 0,14,21:
'P2ry12',#Microglia: 12: C1qa,
'Lgals3', 'Ms4a7',#Monocyte: 7, 19:
'Cd3g', 'Cd3d', 'Cd3e',#T cells: 18:
'Cldn5',#ENDOTHELIAL CELLS: 13:
'Col1a1',#Fibroblasts: 16:
'Acta2',#Pericytes: 8:
]

anno_data_Z = anno_data.drop(["cluster", "cell_type(detail)", "groups","Rd3l"], axis=1)
anno_data_Z = anno_data_Z.set_index("cell_type").groupby("cell_type").mean()
anno_data_Z = anno_data_Z.apply(stats.zscore, axis=0)
anno_data_Z = anno_data_Z.reindex(order)
anno_data_Z = anno_data_Z.loc[:,marker_order]

anno_data_Z.to_csv('./figures/Figure.S4C_normalized_marker_gene_exp_matrix_in_retina.csv')

sns.set(style="ticks",
    font_scale=2.5,
    font='arial',
    rc = {'figure.figsize':(18,15), 'axes.linewidth':1.5})

sns.heatmap(
	data = anno_data_Z.transpose(),
	cmap="rocket_r",
	center=2, vmin=0, vmax=4
    )

plt.xlabel("Cell Types", fontsize=34, rotation=0)
plt.ylabel("Marker Genes", fontsize=34)
plt.xticks(fontsize=34, rotation=30, ha='right')
plt.yticks(fontsize=34)
plt.tight_layout()
plt.savefig("./figures/Figure.S4C_Cell_Type_Annotation_Zscore_Matrix.eps",format='eps',dpi=150)
plt.savefig("./figures/Figure.S4C_Cell_Type_Annotation_Zscore_Matrix.png",format='png',dpi=150)
plt.close()

