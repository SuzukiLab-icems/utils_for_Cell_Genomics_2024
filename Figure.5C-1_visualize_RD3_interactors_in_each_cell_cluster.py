import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import scipy.stats as stats
import seaborn as sns
import scanpy as sc
import warnings
warnings.simplefilter('ignore')

'''Data Processing'''
df_table = pd.read_csv('IO/data_from_Figure.5B/Figure.5B-5_summary_of_cell-gene-exp.csv').loc[:,['Gene','cell_type_1','cluster_id']]
df_table = df_table[df_table["Gene"] != "Rd3"] #Remove RD3
anno = df_table.drop(["Gene"], axis=1).drop_duplicates().sort_values("cluster_id", ascending=True)
stats = df_table.drop_duplicates().groupby("cell_type_1").count().reset_index()
stats = stats.loc[:,["cell_type_1","cluster_id"]].rename(columns={"cluster_id":"count"})
df_component = anno.merge(stats, how="left", on="cell_type_1")

'''Visualization'''
#Configuration (color definition according to Scanpy)
sns.set(style="ticks",
    font_scale=2.5,
    font='arial',
    rc = {'figure.figsize':(8,10), 'axes.linewidth':1.5})
zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
pie_colors = np.array(df_component["cell_type_1"]).astype("<U7")
pie_colors = np.where(pie_colors == "Broad,", zeileis_colors[[15]], pie_colors)
pie_colors = np.where(pie_colors == "SG,", zeileis_colors[[0]], pie_colors)
pie_colors = np.where(pie_colors == "SG,SPC,", zeileis_colors[[7]], pie_colors)
pie_colors = np.where(pie_colors == "SPC,", zeileis_colors[[7]], pie_colors)
pie_colors = np.where(pie_colors == "SG,M,", zeileis_colors[[21]], pie_colors)
pie_colors = np.where(pie_colors == "SPC,M,", zeileis_colors[[21]], pie_colors)
pie_colors = np.where(pie_colors == "M,", zeileis_colors[[21]], pie_colors)
pie_colors = np.where(pie_colors == "SG,RS,", zeileis_colors[[23]], pie_colors) #RS
pie_colors = np.where(pie_colors == "M,RS,", zeileis_colors[[23]], pie_colors)
pie_colors = np.where(pie_colors == "RS,", zeileis_colors[[23]], pie_colors)
pie_colors = np.where(pie_colors == "RS,ES,", zeileis_colors[[9]], pie_colors) #RS,ES
pie_colors = np.where(pie_colors == "SG,ES,", zeileis_colors[[11]], pie_colors)
pie_colors = np.where(pie_colors == "ES,", zeileis_colors[[11]], pie_colors) #ES

#Plot
sns.barplot(data = df_component,
            x="count", y="cell_type_1",
            palette=pie_colors,
            linewidth=1.0,
            edgecolor="black"
)
plt.xlabel("Number of Proteins", fontsize=36)
plt.ylabel("Highly Expressing Cell Types", fontsize=36)
plt.xticks(fontsize=36)
plt.yticks(fontsize=36)
plt.xlim([0.0,95.0])
plt.tight_layout()
sns.despine()
plt.savefig('IO/Figure.5C_barplot_for_components_of_RD3_interactors_in_each_cell_cluster.eps',format='eps',dpi=300)
plt.savefig('IO/Figure.5C_barplot_for_components_of_RD3_interactors_in_each_cell_cluster.png',format='png',dpi=300)
plt.close()

'''Stats'''
print(f'Number of Rd3 interactors highly expressing in spermatids: {df_component[df_component["cluster_id"] >= 7]["count"].sum()}')
