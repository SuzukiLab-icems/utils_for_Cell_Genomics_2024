#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

"""Data Processing.."""
df = pd.read_csv("Figure.3B-4A_adult_C57BL6_expr_table.csv")
df_testis = df[df["Anatomical entity name"] == "testis"]
df_testis["log2(TPM+1)"] = np.log2(1 + df_testis["TPM"])

df_Ovol2 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000037279"]
df_Ovol2["gene"] = "Ovol2"
df_Rnf215 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000003581"]
df_Rnf215["gene"] = "Rnf215"
df_Cldn34c4 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000043569"]
df_Cldn34c4["gene"] = "Cldn34c4"
df_Olfr743 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000094285"]
df_Olfr743["gene"] = "Olfr743"
df_Rd3 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000049353"]
df_Rd3["gene"] = "Rd3"
df_Tacr1 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000030043"]
df_Tacr1["gene"] = "Tacr1"
df_Cebpg = df_testis[df_testis["Gene ID"] == "ENSMUSG00000056216"]
df_Cebpg["gene"] = "Cebpg"
df_Rbm26 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000022119"]
df_Rbm26["gene"] = "Rbm26"
df_P2ry2 = df_testis[df_testis["Gene ID"] == "ENSMUSG00000032860"]
df_P2ry2["gene"] = "P2ry2"

df_fig = pd.concat([
    df_Ovol2,
    df_Rnf215,
    df_Cldn34c4,
    df_Olfr743,
    df_Rd3,
    df_Tacr1,
    df_Cebpg,
    df_Rbm26,
    df_P2ry2],
    axis=0
    )

df_fig.to_csv("Figure.3B_expression_table.csv")

"""Plotting.."""
sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(16,8), 'axes.linewidth':1.5})
font_size = 46

sns.barplot(
    data = df_fig,
    x="gene", y="log2(TPM+1)",
    linewidth = 3,
    edgecolor = "black",
	palette="magma",
    saturation=1,
    )


sns.swarmplot(
    data = df_fig,
    x="gene", y="log2(TPM+1)",
    linewidth = 3,
    edgecolor = "black",
    size = 20,
	color="whitesmoke",
    )


plt.tick_params(labelsize = font_size, width = 1.5)
plt.xticks(fontsize=font_size, rotation=90, ha ='right')
plt.yticks(fontsize=font_size, rotation=90)
plt.xlabel("gene", fontsize=font_size)
plt.ylabel("Log2(TPM+1)", fontsize=font_size)
plt.ylim([-1,10])
plt.tight_layout()

plt.savefig("Figure.3B.eps",format='eps',dpi=300)
plt.savefig("Figure.3B.png",format='png',dpi=300)
