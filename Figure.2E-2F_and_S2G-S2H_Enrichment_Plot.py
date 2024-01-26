"""About EnrichmentScore　deriviation
In the following sgRNA count matrix,

	A(sample)		B(input)
1	A1				B1
2	A2				B2

We calculate ES according to the following processes:
a. modified FoldChange(mFC) = CPM(A1) / (CPM(B1) + 1)
b. Enrichment Score(ES) = log2(1 + mFC)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns

"""Data Processing.."""
df1 = pd.read_table("./PoolA_count_result/PoolA_1st.txt", sep="\t")
df2 = pd.read_table("./PoolA_count_result/PoolA_2nd.txt", sep="\t")
df3 = pd.read_table("./PoolA_count_result/PoolA_3rd.txt", sep="\t")
df4 = pd.read_table("./PoolB_count_result/PoolB_1st.txt", sep="\t")
df5 = pd.read_table("./PoolB_count_result/PoolB_2nd.txt", sep="\t")
df6 = pd.read_table("./PoolB_count_result/PoolB_3rd.txt", sep="\t")

df1 = df1.rename(columns = { df1.columns[2]: '1st' })
df2 = df2.rename(columns = { df2.columns[2]: '2nd' })
df3 = df3.rename(columns = { df3.columns[2]: '3rd' })
df4 = df4.rename(columns = { df4.columns[2]: '1st' })
df5 = df5.rename(columns = { df5.columns[2]: '2nd' })
df6 = df6.rename(columns = { df6.columns[2]: '3rd' })

df1 = pd.concat([df1, df4])
df2 = pd.concat([df2, df5])
df3 = pd.concat([df3, df6])

df1 = df1.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df2 = df2.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df3 = df3.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()

df_2nd_vs_1st = df2.merge(df1.loc[:, ["Gene", "1st"]], on="Gene")
df_3rd_vs_1st = df3.merge(df1.loc[:, ["Gene", "1st"]], on="Gene")

df_2nd_vs_1st["Log2FC"] = \
np.log2(1 + 1000000 * (df_2nd_vs_1st["2nd"] / df_2nd_vs_1st["2nd"].sum())) \
- np.log2(1 + 1000000 * (df_2nd_vs_1st["1st"] / df_2nd_vs_1st["1st"].sum()))

df_3rd_vs_1st["Log2FC"] = \
np.log2(1 + 1000000 * (df_3rd_vs_1st["3rd"] / df_3rd_vs_1st["3rd"].sum())) \
- np.log2(1 + 1000000 * (df_3rd_vs_1st["1st"] / df_3rd_vs_1st["1st"].sum()))

df_2nd_vs_1st["ES"] = \
	np.log2(
		1 + \
		1000000 * (df_2nd_vs_1st["2nd"] / df_2nd_vs_1st["2nd"].sum()) / (1 + 1000000 * (df_2nd_vs_1st["1st"] / df_2nd_vs_1st["1st"].sum()))
		)

df_3rd_vs_1st["ES"] = \
	np.log2(
		1 + \
		1000000 * (df_3rd_vs_1st["3rd"] / df_3rd_vs_1st["3rd"].sum()) / (1 + 1000000 * (df_3rd_vs_1st["1st"] / df_3rd_vs_1st["1st"].sum()))
		)

df_2nd_vs_1st = df_2nd_vs_1st.sort_values(by="ES", ascending=False)
df_2nd_vs_1st["Rank"] = df_2nd_vs_1st["ES"].rank(method='first', ascending=False, na_option='bottom')
df_2nd_vs_1st.to_csv("./2nd_vs_1st_merged.csv")

df_3rd_vs_1st = df_3rd_vs_1st.sort_values(by="ES", ascending=False)
df_3rd_vs_1st["Rank"] = df_3rd_vs_1st["ES"].rank(method='first', ascending=False, na_option='bottom')
df_3rd_vs_1st.to_csv("./3rd_vs_1st_merged.csv")

"""Generate Figure.2E-2F.."""
sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(8,15), 'axes.linewidth':2})
font_size = 34
width = 2
    
fig, (ax1, ax2) = plt.subplots(2, 1)

df_2nd_vs_1st_high_exp = df_2nd_vs_1st[df_2nd_vs_1st["Log2FC"]>=1.5]
df_2nd_vs_1st_low_exp = df_2nd_vs_1st[df_2nd_vs_1st["Log2FC"]<1.5]

sns.scatterplot(data=df_2nd_vs_1st_high_exp, x="Rank", y="ES",
    color="dodgerblue", s=200, edgecolor='royalblue', alpha=0.25, ax=ax1)
sns.scatterplot(data=df_2nd_vs_1st_low_exp, x="Rank", y="ES",
    color="whitesmoke", s=100, edgecolor='silver', alpha=0.1, ax=ax1)

#3rd_vs_1st
df_3rd_vs_1st_high_exp = df_3rd_vs_1st[df_3rd_vs_1st["Log2FC"]>=1.5]
df_3rd_vs_1st_low_exp = df_3rd_vs_1st[df_3rd_vs_1st["Log2FC"]<1.5]

sns.scatterplot(data=df_3rd_vs_1st_high_exp, x="Rank", y="ES",
    color="mediumblue", s=200, edgecolor='midnightblue', alpha=0.25, ax=ax2)
sns.scatterplot(data=df_3rd_vs_1st_low_exp, x="Rank", y="ES",
    color="whitesmoke", s=100, edgecolor='silver', alpha=0.1, ax=ax2)

ax1.set_xlim([-1000,24000])
ax1.set_ylim([-1,15])
ax1.set_ylabel("Enrichment Score", fontsize=font_size)
ax1.set_xlabel("Gene Rank", fontsize=font_size)
ax1.tick_params(labelsize = font_size, width = width )
ax1.axhline(c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)

ax2.set_xlim([-1000,24000])
ax2.set_ylim([-1,15])
ax2.set_ylabel("Enrichment Score", fontsize=font_size)
ax2.set_xlabel("Gene Rank", fontsize=font_size)
ax2.tick_params(labelsize = font_size, width = width )
ax2.axhline(c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)

plt.tight_layout()
sns.despine()

plt.savefig("./Figure.2E-F.eps",format='eps',dpi=300)
plt.savefig("./Figure.2E-F.png",format='png',dpi=300)
plt.close()


"""Generate Figure.S2D-S2E.."""
sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(8,8), 'axes.linewidth':2})
font_size = 34
width = 2
dot_size = 100

df_2nd_vs_1st["detection"] = df_2nd_vs_1st["2nd"].apply(lambda x: "detected" if x >= 1 else "undetected")
df_3rd_vs_1st["detection"] = df_3rd_vs_1st["3rd"].apply(lambda x: "detected" if x >= 1 else "undetected")

sns.scatterplot(
	data = df_2nd_vs_1st,
	x = "Log2FC", y = "ES",
	s = dot_size,
	hue = "detection",
	palette = ["whitesmoke", "dodgerblue"],
	hue_order = ["undetected", "detected"],
	edgecolor = "darkgray",
	lw = width
	)
plt.axvline(x = 1.5, c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)
plt.axhline(y = 2.0, c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.xlabel("Log2FoldChange(2nd vs. 1st)", fontsize=font_size)
plt.ylabel("Enrichment Score(2nd vs. 1st)", fontsize=font_size)
plt.tight_layout()
sns.despine()
plt.savefig("Figure.S2G_1st_vs_INPUT.eps", dpi=300)
plt.savefig("Figure.S2G_1st_vs_INPUT.png", dpi=300)
plt.close()

sns.scatterplot(
	data = df_3rd_vs_1st,
	x = "Log2FC", y = "ES",
	s = dot_size,
	hue = "detection",
	palette = ["whitesmoke", "mediumblue"],
	hue_order = ["undetected", "detected"],
	edgecolor = "darkgray",
	lw = width
	)
plt.axvline(x = 1.5, c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)
plt.axhline(y = 2.0, c='grey', lw=width, ls='--', alpha = 0.6, zorder=1)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
plt.xlabel("Log2FoldChange(3rd vs. 1st)", fontsize=font_size)
plt.ylabel("Enrichment Score(3rd vs. 1st)", fontsize=34)
plt.tight_layout()
sns.despine()
plt.savefig("Figure.S2H_2nd_vs_INPUT.eps", dpi=300)
plt.savefig("Figure.S2H_2nd_vs_INPUT.png", dpi=300)
plt.close()
