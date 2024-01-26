import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as pl
from matplotlib import rcParams
import warnings
warnings.simplefilter('ignore')

'''Data Parsing'''
#testicular cell expression matrix
data = pd.read_csv('IO/data_from_Figure.4B-C/Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv')
data["cluster"] = np.floor(data["groups"])
data = data.drop(["Unnamed: 0"],axis=1).set_index("cluster")
#cell type annotation (Create by yourself according to the result of scRNA-seq data analysis)
label = pd.read_csv('IO/Figure.5B_cluster.csv')
anno_data = label.merge(data, how="left", on="cluster")

'''Generate cell type-annotated expression matrix as Input for Figure.5B-2_classification2heatmap.py'''
anno_data.to_csv('IO/Figure.5B-1_cell_type_annotated_expr_matrix.csv')

'''Data Processing'''
#Normalized Expression Level by Z-score
anno_data_Z = anno_data.drop(["cluster", "cell_type_2", "cell_type_3", "groups"], axis=1)
anno_data_Z = anno_data_Z.set_index("cell_type_1").groupby("cell_type_1").mean()
anno_data_Z = anno_data_Z.apply(stats.zscore, axis=0)
anno_data_Z.to_csv('IO/Figure.5B-2_zscore_generalized_cell_type_annotated_expr_matrix.csv')
anno_data_Z = anno_data_Z.transpose()

#Calculate Cell Type Index score
sigma, p = 0.98, 1
ex_SG = anno_data_Z.loc[:, ["SPC", "M", "RS", "ES"]]
ex_SPC = anno_data_Z.loc[:, ["SG", "M", "RS", "ES"]]
ex_M = anno_data_Z.loc[:, ["SPC", "SG", "RS", "ES"]]
ex_RS = anno_data_Z.loc[:, ["SPC", "M", "SG", "ES"]]
ex_ES = anno_data_Z.loc[:, ["SPC", "M", "RS", "SG"]]
SG = anno_data_Z[anno_data_Z["SG"] - ex_SG.mean(axis=1) >= sigma*p]
SPC = anno_data_Z[anno_data_Z["SPC"] - ex_SPC.mean(axis=1) >= sigma*p]
M = anno_data_Z[anno_data_Z["M"] - ex_M.mean(axis=1) >= sigma*p]
RS = anno_data_Z[anno_data_Z["RS"] - ex_RS.mean(axis=1) >= sigma*p]
ES = anno_data_Z[anno_data_Z["ES"] - ex_ES.mean(axis=1) >= sigma*p]

#Summarization
drop_list = pd.concat([SG,SPC,M,RS,ES])
drop_list = drop_list[~drop_list.index.duplicated()]
Broad = anno_data_Z.drop(index=drop_list.index)
SG["cell_type"] = "SG"
SPC["cell_type"] = "SPC"
M["cell_type"] = "M"
RS["cell_type"] = "RS"
ES["cell_type"] = "ES"
Broad["cell_type"] = "Broad"
summary = pd.concat([SG, SPC, M, RS, ES, Broad])
overlap = summary.reset_index().loc[:,["index", "cell_type"]]
overlap['cell_type'] = overlap['cell_type'].astype(str) + ','
overlap = overlap.groupby("index").sum().reset_index()
overlap.to_csv('IO/Figure.5B-3_correspondence_table_between_genes_and_cell_type.csv')

