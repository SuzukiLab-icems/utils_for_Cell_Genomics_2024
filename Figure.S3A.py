#!/usr/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

df = pd.read_csv("Figure.S3A_raw_data.csv")

sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(10,9), 'axes.linewidth':1.5})
    
sns.violinplot(
    data = df,
    x="gene", y="2**delta_Ct",
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
    x="gene", y="2**delta_Ct",
    linewidth = 3,
    edgecolor = "black",
    size = 20,
    palette="magma"
    )

plt.tick_params(labelsize = 32, width = 1.5)
plt.xticks(fontsize=32, rotation=20, ha ='right')
plt.yticks(fontsize=32)
plt.xlabel("gene", fontsize=32)
plt.ylabel("2-∆Ct(target vs Actb)", fontsize=32)
plt.ylim([0.000000001, 1.0])
plt.yscale('log')
plt.tight_layout()
sns.despine()

plt.savefig("Figure.S3A.eps",format='eps',dpi=300)
plt.savefig("Figure.S3A.png",format='png',dpi=300)

