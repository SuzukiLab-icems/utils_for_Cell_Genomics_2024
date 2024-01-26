#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

df = pd.read_csv("./figures/Figure.4C_marker_expr_table.csv", index_col=0)
df = df.loc[:,["groups","Rd3"]].rename(columns={"groups":"cell types"})
df["cell types"] = np.floor(df["cell types"])
df["Rd3"] = stats.zscore(df["Rd3"])

sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(20,8), 'axes.linewidth':1.5})

sns.regplot(
    data = df,
    x="cell types",y="Rd3",
    order=3,
    x_jitter=.3,
    color='lavender',
    scatter_kws={'edgecolor':'darkslateblue','s':250},
    line_kws={'color':'mediumorchid','linewidth':5,'ls':'-.'},
    ci=None)

plt.tick_params(labelsize = 38, width = 1.5)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],fontsize=38, rotation=0, ha ='center')
plt.yticks([-1,0,1,2,3],fontsize=38)
plt.xlabel("Cell Type", fontsize=38)
plt.ylabel("Rd3 Expression (Z-Score)", fontsize=38)
plt.ylim([1.1*df["Rd3"].min(),1.1*df["Rd3"].max()])
plt.xlim([-1,19])
plt.tight_layout()
sns.despine()

plt.savefig("./figures/Figure.4C_Rd3_Expr_Profiling.eps",format='eps',dpi=300)
plt.savefig("./figures/Figure.4C_Rd3_Expr_Profiling.png",format='png',dpi=300)
plt.close()
