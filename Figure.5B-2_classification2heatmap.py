import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import scipy.stats as stats
import seaborn as sns
import scanpy as sc
import warnings
warnings.simplefilter('ignore')

'''Data Parsing'''
anno_data = pd.read_csv('IO/Figure.5B-1_cell_type_annotated_expr_matrix.csv', index_col=0)
candidate_list = pd.read_csv('IO/data_from_Figure.5A/Figure.5A-4.sig.csv', index_col=0)
cluster = pd.read_csv('IO/Figure.5B-3_correspondence_table_between_genes_and_cell_type.csv', index_col=0)

'''Data Processing'''
#Labeling Cell Type ID
label = {
  "Broad,"   : "00",
  "SG,"      : "01",
  "SG,SPC,"  : "02",
  "SPC,"     : "03",
  "SG,M,"    : "04",
  "SPC,M,"   : "05",
  "M,"       : "06",
  "SG,RS,"   : "07",
  "SPC,RS,"  : "08",
  "M,RS,"    : "09",
  "RS,"      : "10",
  "RS,ES,"   : "11",
  "SG,RS,ES,": "12",
  "SG,ES,"   : "13",
  "SPC,ES,"  : "14",
  "M,ES,"    : "15",
  "ES,"      : "16"
}
cluster_id = cluster.replace(label).rename(columns = {"cell_type":"cluster_id"})
cluster = cluster.merge(cluster_id, how="left", on="index").rename(columns = {"index":"Gene", "cell_type":"cell_type_1"})

#Processing for Z-Score normalization
candidate_list = candidate_list.loc[:,["m.m.symbol"]]
candidate_list = candidate_list.rename(columns={"m.m.symbol":"index"})
anno_data_Z = anno_data.drop(["cluster", "cell_type_1", "cell_type_2", "groups"], axis=1)
anno_data_Z = anno_data_Z.set_index("cell_type_3").transpose().reset_index()
anno_data_Z = candidate_list.merge(anno_data_Z, how="left", on="index")
anno_data_Z = anno_data_Z.set_index("index").transpose()
anno_data_Z = anno_data_Z.dropna(axis=1)

#Execute Z-Score normalization
anno_data_Z = anno_data_Z.reset_index().groupby("index").mean()
anno_data_Z = anno_data_Z.apply(stats.zscore, axis=0)
anno_data_Z.to_csv('IO/Figure.5B-4.zscore_generalized_expr_matrix_in_each_cell_cluster.csv')
#>Figure.5B-2: zscore in each cell
#>Figure.5B-4: zscore in each cell cluster

'''Generate Summarized Table'''
#Data Processing
anno_data_Z = anno_data_Z.reset_index().rename(columns={"index":"cell_type"})
df_concat = pd.DataFrame(columns=["cell_type", "Z-Score", "Gene"])
for gene in anno_data_Z.columns[1:]:
  df_tmp = anno_data_Z.loc[:, ["cell_type", gene]]
  df_tmp = df_tmp.rename(columns = {gene : "Z-Score"})
  df_tmp["Gene"] = gene
  df_concat = pd.concat([df_concat,df_tmp])
data = df_concat.merge(cluster, how="left", on="Gene").set_index("Gene")
data.to_csv('IO/Figure.5B-5_summary_of_cell-gene-exp.csv')

'''Visualization'''
#preprocessing for Visualization
annotation = {
'SG'	:0,
'pL/Z'	:1,
'P.1'	:2,
'P.2'	:3,
'P.3'	:4,
'P.4'	:5,
'P.5'	:6,
'D'		:7,
'M'		:8,
'RS.1'	:9,
'RS.2'	:10,
'RS.3'	:11,
'RS.4'	:12,
'RS.5'	:13,
'ES.1'	:14,
'ES.2'	:15,
'ES.3'	:16,
'ES.4'	:17,
'ES.5'	:18
}
data = data.replace(annotation)
data['gene_rank'] = data["cluster_id"] + "_" + data.index.astype(str)
data_raw = data
data = data.reset_index().pivot_table(
	values="Z-Score",
	index="gene_rank",
	columns="cell_type"
	)
data.to_csv('IO/Figure.5B-6_normalzied_expr_table_for_heatmap.csv')
data = data.drop(["10_Rd3"], axis=0) #remove Rd3

#Definition of color(for col_colors) from Scanpy library
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
col_colors = np.array(data.columns[0:]).astype("<U7")
col_colors[[0]] = zeileis_colors[[0]] #SPGs
col_colors[[1]] = zeileis_colors[[6]] #pL/Z
col_colors[[2,3,4,5,6]] = zeileis_colors[[7,7,7,7,7]] #P
col_colors[[7]] = zeileis_colors[[8]] #D
col_colors[[8]] = zeileis_colors[[21]] #M
col_colors[[9,10,11,12,13]] = zeileis_colors[[23,23,23,23,23]] #RS
col_colors[[14,15,16,17,18]] = zeileis_colors[[11,11,11,11,11]] #ES

#Definition of color(for raw_colors) from Scanpy library
row_labels = data_raw.sort_values(by="gene_rank", ascending=True).drop_duplicates(subset="gene_rank")
row_colors = np.array(row_labels["cluster_id"]).astype("<U7")
row_colors = np.where(row_colors == "00", zeileis_colors[[15]], row_colors)
row_colors = np.where(row_colors == "01", zeileis_colors[[0]], row_colors)
row_colors = np.where(row_colors == "02", zeileis_colors[[7]], row_colors)
row_colors = np.where(row_colors == "03", zeileis_colors[[7]], row_colors)
row_colors = np.where(row_colors == "04", zeileis_colors[[21]], row_colors)
row_colors = np.where(row_colors == "05", zeileis_colors[[21]], row_colors)
row_colors = np.where(row_colors == "06", zeileis_colors[[21]], row_colors)
row_colors = np.where(row_colors == "07", zeileis_colors[[23]], row_colors)
row_colors = np.where(row_colors == "08", zeileis_colors[[23]], row_colors)
row_colors = np.where(row_colors == "09", zeileis_colors[[23]], row_colors)
row_colors = np.where(row_colors == "10", zeileis_colors[[23]], row_colors)
row_colors = np.where(row_colors == "11", zeileis_colors[[9]], row_colors)
row_colors = np.where(row_colors == "12", zeileis_colors[[9]], row_colors)
row_colors = np.where(row_colors == "13", zeileis_colors[[11]], row_colors)
row_colors = np.where(row_colors == "14", zeileis_colors[[11]], row_colors)
row_colors = np.where(row_colors == "15", zeileis_colors[[11]], row_colors)
row_colors = np.where(row_colors == "16", zeileis_colors[[11]], row_colors)

#Visualization
sns.set(style="ticks",
    font_scale=2.5,
    font='arial',
    rc = {'figure.figsize':(20,20), 'axes.linewidth':1.5})

g = sns.clustermap(data,
      cmap="GnBu",
      center=0,
      vmin=-2, vmax=2,
      col_cluster=False,
      row_cluster=False,
      col_colors=col_colors,
      row_colors=row_colors,
      xticklabels=False,
      yticklabels=False
    )
    
g.fig.subplots_adjust(right=0.7)
g.ax_cbar.set_position((0.8, .2, .03, .4))
g.ax_cbar.set_title('Z-Score', rotation=0)

plt.tight_layout()
plt.savefig('IO/Figure.5B.testis_expression_profile_of_RD3_interactors.eps',format='eps',dpi=300)
plt.savefig('IO/Figure.5B.testis_expression_profile_of_RD3_interactors.png',format='png',dpi=300)

plt.close()
