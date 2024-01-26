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
df_list = pd.read_csv('IO/data_from_Figure.5B/Figure.5B-5_summary_of_cell-gene-exp.csv').loc[:,['Gene','cell_type_1']].drop_duplicates()
df_sig = pd.read_csv('IO/data_from_Figure.5A/Figure.5A-4.sig.csv',index_col=0)
df_sig = df_list.rename(columns={"Gene":"m.m.symbol"}).merge(df_sig,how="left",on="m.m.symbol")
columns = [
'cell_type_1',
'Gene Symbol',
'Accession',
'Ensembl Gene ID',
'm.m.symbol',
'ensemble',
'KO_rep.1',
'KO_rep.2',
'KO_rep.3',
'RD3_rep.1',
'RD3_rep.2',
'RD3_rep.3',
'statistics',
'p-value',
'-Log10(p-value)',
'Log2FC'
]
df_sig = df_sig.loc[:,columns]
fixed_columns = [
'cell_type',
'hSymbol',
'hAccession',
'hEnsembl_ID',
'mSymbol',
'mEnsembl_ID',
'KO_rep.1',
'KO_rep.2',
'KO_rep.3',
'RD3_rep.1',
'RD3_rep.2',
'RD3_rep.3',
'statistics',
'p-value',
'-Log10(p-value)',
'Log2FC'
]
df_sig.columns = fixed_columns
spermatid = ['ES,','M,RS,','RS,','RS,ES,','SG,ES,','SG,RS,']
df_spermatid = df_sig.set_index("cell_type").loc[spermatid]
df_spermatid = df_spermatid[df_spermatid["mSymbol"] != "Rd3"] #Remove Rd3
df_spermatid.to_csv('IO/Figure.5C_summary_table_of_spermatid_specific_RD3_interactors.csv')
