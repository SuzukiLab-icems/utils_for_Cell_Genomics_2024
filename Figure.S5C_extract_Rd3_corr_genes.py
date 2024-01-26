import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
from scipy.stats import spearmanr
import warnings
warnings.simplefilter('ignore')

'''configuration'''
expr_matrix = 'IO/Figure.4B-C_testicular_gene_expr_table_for_Figure.5.csv'
targeted_gene = 'Rd3'

def generate_Rd3_coexpression_network(expr_matrix, targeted_gene):
    gene_1 = targeted_gene
    expr_matrix = pd.read_csv(expr_matrix, index_col=0)
    opposite_list = expr_matrix.columns.to_list()
    opposite_list.remove(gene_1)
    GCN = []
    for gene_2 in opposite_list:
        correlation, _ = spearmanr(expr_matrix[gene_1], expr_matrix[gene_2])
        GCN.append([gene_1, gene_2, str(correlation)])
    GCN = pd.DataFrame(GCN,columns=["gene_1","gene_2","rho"])
    GCN = GCN[(GCN["rho"] != "nan")&(GCN["gene_2"] != "groups")]
    return GCN

def corr2list(expr_matrix, targeted_gene):
	corr_table = generate_Rd3_coexpression_network(expr_matrix, targeted_gene)
	corr_table["rho"] = corr_table["rho"].astype("float")
	corr_table = corr_table.sort_values("rho", ascending=True)
	high_corr = corr_table[corr_table["rho"] >= 0.8]
	low_corr = corr_table[corr_table["rho"] < 0.8]
	low_corr["rank"] = low_corr["rho"].rank(axis=0, method='first',ascending=True, pct=False)
	high_corr["rank"] = high_corr["rho"].rank(axis=0, method='first',ascending=True, pct=False) + low_corr["rank"].max()
	return high_corr, low_corr
	
def Rank_Plot(high_corr, low_corr):
    sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(14,8), 'axes.linewidth':2.5})
    sns.scatterplot(
      data=high_corr,
      x="rank", y="rho",
      color="orchid", s=200, edgecolor='magenta'
      )
    sns.scatterplot(
      data=low_corr,
      x="rank", y="rho",
      color="lightgrey", s=50, edgecolor='lightgrey'
      )
    plt.xlabel("Rank(Gene)", fontsize=40)
    plt.ylabel("spearman(ρ)", fontsize=40)
    plt.xticks(fontsize=34)
    plt.yticks(fontsize=34)
    plt.axhline(y=0.8, color="black", lw=1.5, ls="--")
    plt.xlim([-1000,25000])
    plt.ylim([-1.1,1.1])
    plt.tight_layout()
    sns.despine()
    plt.savefig('IO/Figure.S5C_Rd3_correlated_genes.eps',format='eps',dpi=300)
    plt.savefig('IO/Figure.S5C_Rd3_correlated_genes.png',format='png',dpi=300)
    plt.close()
    return

if __name__ == "__main__":
	high_corr, low_corr = corr2list(expr_matrix, targeted_gene)
	high_corr.to_csv('IO/Figure.S5C_GCN_for_Rd3.csv')
	Rank_Plot(high_corr, low_corr)
	
