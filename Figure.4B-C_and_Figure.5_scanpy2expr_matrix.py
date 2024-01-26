#!/usr/bin/python

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

adata = sc.read_10x_mtx(
    'filtered_feature_bc_matrix/',
    var_names='gene_symbols',
    cache=True
)

adata.var_names_make_unique()
adata.write(results_raw)

"""preProcessing"""
sc.pl.highest_expr_genes(adata, n_top=20, save=True, show=False)

sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, save=True, show=False)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt' , save=True, show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts' , save=True, show=False)
sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', save=True, show=False)

adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

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
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
#You can detect Sertoli Cells by 'sc.pp.neighbors(adata, n_neighbors=4, n_pcs=50)'.
#As a Sertoli cells marker, Vim is the best to show descrete cell type.
sc.tl.tsne(adata)
sc.pl.tsne(adata, color=['Rd3'])
sc.tl.leiden(adata)
sc.pl.tsne(adata, color=['leiden', 'Rd3', 'Aurka', 'Prm2', 'Ssxb1', 'Dmrt1', 'Piwil1'], save=True, show=False, edges=False)
adata.write(results_TSNE)

"""DEGs Extraction"""
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False, save=True, show=False)

"""Cell Type Annotation"""
sc.tl.leiden(adata, resolution=1.0)
###Cell Type###
#Spermatogonia:12:Dmrt1
#pL/Z:16:Ly6k,Tex101,Sycp1
#P:0,3,13,9,17:Piwil1
#D:1:Pou5f2
#M:6:Aurka,Ccna1
#RS.1:5:Ssxb1
#RS.2:10:Acrv1
#RS.3:14,15:Spaca4
#ES.1:11:Hils1
#ES.2:4,7:Tnp1
#ES.3:2,8:Prm2

marker=[
'Dmrt1',
'Ly6k',
'Tex101',
'Sycp1',
'Piwil1',
'Ybx3',
'Pou5f2',
'Aurka',
'Ccna1',
'Ssxb1',
'Acrv1',
'Spaca4',
'Hils1',
'Tnp1',
'Prm2',
]

sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden'], save=False , show=False)

adata.obs['leiden'].cat.categories
adata.obs['leiden_anno'] = adata.obs['leiden']
adata.uns['leiden_anno_colors'] = adata.uns['leiden_colors']

adata.obs['leiden_anno'].cat.categories = [
'P.1',     	#0
'D',      	#1
'ES.4',     #2
'P.2',      #3
'ES.2',     #4
'RS.1',     #5
'M',   		#6
'ES.3',     #7
'ES.5',     #8
'P.3',      #9
'RS.2',     #10
'ES.1',   	#11
'SG',     	#12
'P.4',      #13
'RS.4',     #14
'RS.5',     #15
'pL/Z',   	#16
'P.5',      #17
'RS.3',     #18
]

list = (
'SG',     	#12>0
'pL/Z',   	#16>1
'P.1',     	#0>2
'P.2',      #3>3
'P.3',      #9>4
'P.4',      #13>5
'P.5',      #17>6
'D',      	#1>7
'M',   		#6>8
'RS.1',     #5>9
'RS.2',     #10>10
'RS.3',     #18>11
'RS.4',     #14>12
'RS.5',     #15>13
'ES.1',   	#11>14
'ES.2',     #4>15
'ES.3',     #7>16
'ES.4',     #2>17
'ES.5',     #8>18
)
adata.obs['leiden_anno_ordered'] = adata.obs['leiden_anno'].cat.reorder_categories(list, inplace=True)

zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['leiden_anno_colors'])
new_colors[[0]] = zeileis_colors[[0]] #SPGs
new_colors[[1]] = zeileis_colors[[6]] #pL/Z
new_colors[[2,3,4,5,6]] = zeileis_colors[[7,7,7,7,7]] #P
new_colors[[7]] = zeileis_colors[[8]] #D
new_colors[[8]] = zeileis_colors[[21]] #M
new_colors[[9,10,11,12,13]] = zeileis_colors[[23,23,23,23,23]] #RS
new_colors[[14,15,16,17,18]] = zeileis_colors[[11,11,11,11,11]] #ES
adata.uns['leiden_anno_colors'] = new_colors

sc.pl.tsne(adata, color=['leiden_anno', 'Rd3'], save=True, show=False, edges=False) # for Figure.4B

"""PAGA"""
sc.tl.paga(adata, groups='leiden_anno')
sc.pl.paga(adata, threshold=0.03, show=False, save=True)
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden_anno', 'Rd3','Rd3l'], legend_loc='on data', save=True, show=False)
sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=False, save=True)

"""Generate Expression Matrix with Trajectory Analysis"""
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden_anno']  == 'SG')[0]
sc.tl.dpt(adata)
gene_names = [
'Rd3',
'Dmrt1', #SG
'Ly6k','Tex101','Sycp1', #pL/Z
'Piwil1','Ybx3', #P
'Pou5f2', #D
'Aurka', 'Ccna1',#M
'Ssxb1', #RS.1
'Acrv1', #RS.2
'Spaca4', #RS.3
'Hils1', #ES.1
'Tnp1', #ES.2
'Prm2', #ES.3
]
paths = [('spermatogenesis', [
'SG',     	#12>0
'pL/Z',   	#16>1
'P.1',     	#0>2
'P.2',      #3>3
'P.3',      #9>4
'P.4',      #13>5
'P.5',      #17>6
'D',      	#1>7
'M',   		#6>8
'RS.1',     #5>9
'RS.2',     #10>10
'RS.3',     #18>11
'RS.4',     #14>12
'RS.5',     #15>13
'ES.1',   	#11>14
'ES.2',     #4>15
'ES.3',     #7>16
'ES.4',     #2>17
'ES.5',     #8>18
]
)]

sc.pl.draw_graph(adata, color=['leiden_anno', 'dpt_pseudotime'], legend_loc='on data', save=True)

adata.obs['distance'] = adata.obs['dpt_pseudotime']
adata.obs['clusters'] = adata.obs['leiden_anno']
adata.uns['clusters_colors'] = adata.uns['leiden_anno_colors']

_, axs = pl.subplots(ncols=3, figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)

"""Generate Expression Matrix of Marker Genes for Figure.4C and S4A"""
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
    data.to_csv('./figures/Figure.4C_marker_expr_table.csv'.format(descr))

"""Generate Expression Matrix of Whole Genes for Figure.5B, 5E, 5F and 5G"""
for ipath, (descr, path) in enumerate(paths):
	_, data = sc.pl.paga_path(
		adata, path, adata.raw.var.index,
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
	data.to_csv('./figures/Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv')

