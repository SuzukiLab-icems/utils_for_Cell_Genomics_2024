import sys
import os
import glob
import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import warnings
warnings.simplefilter('ignore')

def ztest(df):
  #data processing
  df_ztest = df.loc[:,["Accession","KO_rep.1", "KO_rep.2", "KO_rep.3","RD3_rep.1", "RD3_rep.2","RD3_rep.3"]]
  df_ztest_T = df_ztest.set_index("Accession").transpose()
  df_ztest_T.insert(0, "pheno_data", ["KO", "KO", "KO", "RD3", "RD3", "RD3"])
  df_pheno_KO = df_ztest_T[df_ztest_T["pheno_data"] == "KO"].drop("pheno_data",axis=1).transpose()
  df_pheno_RD3 = df_ztest_T[df_ztest_T["pheno_data"] == "RD3"].drop("pheno_data",axis=1).transpose()
  #define values for statistical analysis of RD3 KO
  µ_KO = df_pheno_KO.mean(axis=1)
  sigma_KO = df_pheno_KO.std(axis=1, ddof=0)
  n_KO = df_pheno_KO.count(axis=1)
  #define values for statistical analysis of RD3 OE
  µ_RD3 = df_pheno_RD3.mean(axis=1)
  sigma_RD3 = df_pheno_RD3.std(axis=1, ddof=0)
  n_RD3 = df_pheno_RD3.count(axis=1)
  #perform z-test for RD3 KO
  for rep in df_pheno_KO.columns[0:]:
    if is_numeric_dtype(df_pheno_KO[rep]):
      x = df_pheno_KO[rep]
      Z_KO = abs(x - µ_KO) / (sigma_KO / np.sqrt(n_KO))
      df_pheno_KO["pvalues_" + rep] = stats.norm.sf(Z_KO, loc=0, scale=1)
  #perform z-test for RD3 OE
  for rep in df_pheno_RD3.columns[0:]:
    if is_numeric_dtype(df_pheno_RD3[rep]):
      x = df_pheno_RD3[rep]
      Z_RD3 = (x - µ_RD3) / (sigma_RD3 / np.sqrt(n_RD3))
      df_pheno_RD3["pvalues_" + rep] = stats.norm.cdf(Z_RD3, loc=0, scale=1)
  #data summary
  df_ztest = df_pheno_KO.merge(df_pheno_RD3, how="left", on="Accession")
  df_ztest = df_ztest.reset_index()
  df_ztest = df_ztest.loc[:,
    ["Accession",
    "pvalues_KO_rep.1", "pvalues_KO_rep.2", "pvalues_KO_rep.3",
    "pvalues_RD3_rep.1", "pvalues_RD3_rep.2", "pvalues_RD3_rep.3"]
    ]
  df_ztest.to_csv('./result/Figure.5A-1.ztest_table.csv')

def cleanup(df):
  #data processing
  df_ztest = pd.read_csv('./result/Figure.5A-1.ztest_table.csv')
  df_ttest = df.merge(df_ztest, how="left", on="Accession")
  df_temp = df.loc[:, ["Accession"]]
  df_KO_rep_1 = df_ttest.loc[:, ["Accession", "KO_rep.1", "pvalues_KO_rep.1"]]
  df_KO_rep_2 = df_ttest.loc[:, ["Accession", "KO_rep.2", "pvalues_KO_rep.2"]]
  df_KO_rep_3 = df_ttest.loc[:, ["Accession", "KO_rep.3", "pvalues_KO_rep.3"]]
  df_RD3_rep_1 = df_ttest.loc[:, ["Accession", "RD3_rep.1", "pvalues_RD3_rep.1"]]
  df_RD3_rep_2 = df_ttest.loc[:, ["Accession", "RD3_rep.2", "pvalues_RD3_rep.2"]]
  df_RD3_rep_3 = df_ttest.loc[:, ["Accession", "RD3_rep.3", "pvalues_RD3_rep.3"]]
  #outlier removal of RD3 KO
  df_KO_rep_1_in = df_KO_rep_1[df_KO_rep_1["pvalues_KO_rep.1"] > 0.01]
  df_KO_rep_2_in = df_KO_rep_2[df_KO_rep_2["pvalues_KO_rep.2"] > 0.01]
  df_KO_rep_3_in = df_KO_rep_3[df_KO_rep_3["pvalues_KO_rep.3"] > 0.01]
  df_KO = df_temp.merge(
    df_KO_rep_1_in, how="left", on="Accession"
    ).merge(df_KO_rep_2_in, how="left", on="Accession"
    ).merge(df_KO_rep_3_in, how="left", on="Accession"
  )
  #outlier removal of RD3 OE
  df_RD3_rep_1_in = df_RD3_rep_1[df_RD3_rep_1["pvalues_RD3_rep.1"] > 0.01]
  df_RD3_rep_2_in = df_RD3_rep_2[df_RD3_rep_2["pvalues_RD3_rep.2"] > 0.01]
  df_RD3_rep_3_in = df_RD3_rep_3[df_RD3_rep_3["pvalues_RD3_rep.3"] > 0.01]
  df_RD3 = df_temp.merge(
    df_RD3_rep_1_in, how="left", on="Accession"
    ).merge(df_RD3_rep_2_in, how="left", on="Accession"
    ).merge(df_RD3_rep_3_in, how="left", on="Accession"
  )
  #generate cleaned up data for t-test analysis
  df_ttest = df_temp.merge(
    df_KO, how="left", on="Accession"
    ).merge(df_RD3, how="left", on="Accession"
  )
  df_ttest = df_ttest.loc[:,
    [
    "Accession",
    "KO_rep.1", "KO_rep.2", "KO_rep.3",
    "RD3_rep.1", "RD3_rep.2", "RD3_rep.3"
    ]
  ]
  df_ttest.to_csv('./result/Figure.5A-2.ttest_table.csv')

def ttest():
  #data processing
  df_ttest = pd.read_csv('./result/Figure.5A-2.ttest_table.csv', index_col=0)
  df_T = df_ttest.set_index("Accession").transpose()
  df_T.insert(0, "pheno_data", ["KO", "KO", "KO", "RD3", "RD3", "RD3"])
  df_pheno_KO = df_T[df_T["pheno_data"] == "KO"]
  df_pheno_RD3 = df_T[df_T["pheno_data"] == "RD3"]
  #t-test: In this study, to detect many interactors for performing Hub-Explorer study, we have never applied multiple testing corrections for p-values.
  f=open('./result/Figure.5A-3.statistics.txt', "w", encoding='utf-8', newline='\n')
  for accession in df_T.columns[1:]:
    if is_numeric_dtype(df_T[accession]):
      if (len(df_pheno_KO[accession].dropna()) == 0) or (len(df_pheno_RD3[accession].dropna()) == 0):
        _ = f.write(f'{accession},statistics:,NaN,pvalue:,NaN\n')
      else:
        statistics, pvalue = stats.ttest_ind(df_pheno_KO[accession].dropna(), df_pheno_RD3[accession].dropna(), equal_var=False)
        _ = f.write(f'{accession},statistics:,{statistics},p-value:,{pvalue}\n')
  f.close()

def significant_proteins(df_stat):
  df_annotation = pd.read_csv('annotation_table_for_Figrue.5A.csv').loc[:, ["Input", "Symbol", "Ensembl ID"]] #This annotation file should be prepared manually.
  df_annotation = df_annotation.rename(columns={'Input':'Gene Symbol', 'Symbol':'m.m.symbol', 'Ensembl ID':'ensemble'})
  df_stat = df_annotation.merge(df_stat, how='left', on='Gene Symbol')
  df_sig = df_stat[(df_stat['Log2FC'] >= 1.5) & (df_stat['p-value'] < 0.05)]
  df_unsig = df_stat[(df_stat['Log2FC'] < 1.5) | (df_stat['p-value'] >= 0.05)]
  df_sig.to_csv('./result/Figure.5A-4.sig.csv')
  df_unsig.to_csv('./result/Figure.5A-5.unsig.csv')

def volcano_plot():
  #In our study, we manually put RD3 data as a positive control onto Figure.5A_VolcanoPlot.eps(png) because Foldchange value of RD3 is 'inf'.
  df_sig = pd.read_csv('./result/Figure.5A-4.sig.csv')
  df_unsig = pd.read_csv('./result/Figure.5A-5.unsig.csv')
  sns.set(style="ticks", font_scale=2.5, font='arial',
  rc = {'figure.figsize':(12,8), 'axes.linewidth':1.5})
  sns.scatterplot(data=df_sig, x="Log2FC", y="-Log10(p-value)", color="orchid", s=250, edgecolor='magenta')
  sns.scatterplot(data=df_unsig, x="Log2FC", y="-Log10(p-value)",color="azure", s=50, edgecolor='dimgrey')
  plt.xlabel("Log2FC(RD3 vs. KO)", fontsize=34)
  plt.ylabel("-Log10(p-value)", fontsize=34)
  plt.xticks(fontsize=34)
  plt.yticks(fontsize=34)
  plt.axhline(y=-1.0 * np.log10(0.05), c='grey', lw=1, ls='--', zorder=0)
  plt.axvline(x=1.5, c='grey', lw=1, ls='--', zorder=0)
  plt.xlim([-9.0,14.0])
  plt.ylim([0,7])
  plt.tight_layout()
  sns.despine()
  plt.savefig('./result/Figure.5A_VolcanoPlot.eps',format='eps',dpi=300)
  plt.savefig('./result/Figure.5A_VolcanoPlot.png',format='png',dpi=300)
  plt.close()

def main():
	print('processing..')
	sample_name = {
		"Abundances (Scaled): F5: Sample, KO": "KO_rep.1",
		"Abundances (Scaled): F7: Sample, KO": "KO_rep.2",
		"Abundances (Scaled): F9: Sample, KO": "KO_rep.3",
		"Abundances (Scaled): F6: Sample, Res": "RD3_rep.1",
		"Abundances (Scaled): F8: Sample, Res": "RD3_rep.2",
		"Abundances (Scaled): F10: Sample, Res": "RD3_rep.3"}
	df = pd.read_excel('221215_Noguchi_RD3_SPOT_trap_n3_FusionIT.xlsx', sheet_name=None)["Proteins"]
	df = df.rename(columns=sample_name)
	df = df.loc[:,["Accession","Gene Symbol","Ensembl Gene ID","Entrez Gene ID","KO_rep.1","KO_rep.2","KO_rep.3","RD3_rep.1", "RD3_rep.2","RD3_rep.3"]]
	df_template = df.loc[:, ["Accession", "Gene Symbol", "Ensembl Gene ID", "Entrez Gene ID"]]
	print('outlier removing..')
	ztest(df)
	cleanup(df)
	print('proceeding t-test..')
	ttest()
	print('generating summary..')
	df_ttest = pd.read_csv('./result/Figure.5A-2.ttest_table.csv')
	df_pvalue = pd.read_table('./result/Figure.5A-3.statistics.txt', sep=',', usecols=[0,2,4], header=None, names=['Accession', 'statistics', 'p-value'])
	df_stat = df_template.merge(df_ttest, how="left", on="Accession")
	df_stat = df_stat.merge(df_pvalue, how="left", on="Accession")
	df_stat['-Log10(p-value)'] = -1 * np.log10(df_stat['p-value'])
	df_stat["Log2FC"] = \
		np.log2(df_stat.loc[:, ["RD3_rep.1", "RD3_rep.2", "RD3_rep.3"]].sum(axis=1) / df_stat.loc[:, ["RD3_rep.1", "RD3_rep.2", "RD3_rep.3"]].count(axis=1)) \
		- np.log2(df_stat.loc[:, ["KO_rep.1", "KO_rep.2", "KO_rep.3"]].sum(axis=1) / df_stat.loc[:, ["KO_rep.1", "KO_rep.2", "KO_rep.3"]].count(axis=1))
	significant_proteins(df_stat)
	volcano_plot()
	df_stat.to_csv('./result/Figure.5A-6.result_table.csv')

if __name__ == "__main__":
	main()
	print('done.')










