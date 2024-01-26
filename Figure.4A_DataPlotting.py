#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

#File Load
df = pd.read_csv("Figure.3B-4A_adult_C57BL6_expr_table.csv")
n_tissues = np.unique(df["Anatomical entity name"].tolist()).size
n_samples = df.groupby(["Gene ID","Anatomical entity name"]).count()["TPM"].sum()
print(f'n_tissues = {n_tissues}\nn_samples = {n_samples}')

#Filtering / Rd3 Extraction
df = df.loc[:,["Gene ID","Anatomical entity name","TPM"]]
df = df[df["Gene ID"] == "ENSMUSG00000049353"]

#Calculation
df["log2(TPM+1)"] = np.log2(1 + df["TPM"])
#Prep
df = df.sort_values(by="log2(TPM+1)", ascending=False)
df.to_csv('Figure.4A_M.M.Bgee_RD3.csv')

#Plot
sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(18,10), 'axes.linewidth':1.5})
color_dict = dict(
{
"Ammon's horn":"lavender",
"adult mammalian kidney":"lavender",
"brain":"lavender",
"cerebellar cortex":"lavender",
"cerebellum":"lavender",
"colon":"lavender",
"cortex of kidney":"lavender",
"duodenum":"lavender",
"esophagus":"lavender",
"granulocyte":"lavender",
"heart":"lavender",
"hindlimb stylopod muscle":"lavender",
"ileum":"lavender",
"jejunum":"lavender",
"lip":"lavender",
"liver":"lavender",
"lung":"lavender",
"muscle tissue":"lavender",
"ovary":"lavender",
"pancreas":"lavender",
"placenta":"lavender",
"primary visual cortex":"lavender",
"proximal tubule":"lavender",
"quadriceps femoris":"lavender",
"retinal neural layer":"darkslateblue",
"right kidney":"lavender",
"skeletal muscle tissue":"lavender",
"spleen":"lavender",
"stomach":"lavender",
"superior frontal gyrus":"lavender",
"testis":"magenta",
"thymus":"lavender",
"uterus":"lavender",
"zone of skin":"lavender"
}
)

sns.violinplot(
    data = df,
    x="Anatomical entity name", y="log2(TPM+1)",
    inner=None,
    scale="width",
    bw=0.25,
    linewidth = 3,
    edgecolor = "black",
    size = 20,
    color="whitesmoke",
    saturation=1
    )

sns.swarmplot(
    data = df,
    x="Anatomical entity name", y="log2(TPM+1)",
    linewidth = 1.5,
    edgecolor = "black",
    size = 10,
    palette=color_dict
    )

plt.tick_params(labelsize = 32, width = 1.5)
plt.xticks(fontsize=32, rotation=90)
plt.yticks([0,2,4,6,8],fontsize=32)
plt.ylabel("Log2(TPM+1)", fontsize=32)
plt.xlabel("Tissue", fontsize=32)
plt.ylim([-1,9])
plt.tight_layout()
sns.despine()

plt.savefig("Figure.4A_Violin_Plot_for_RD3_Expression.eps",format='eps',dpi=300)
plt.savefig("Figure.4A_Violin_Plot_for_RD3_Expression.png",format='png',dpi=300)
plt.close()

	

