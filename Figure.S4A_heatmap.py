#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

#File Loading >> manipulating
df = pd.read_csv("./figures/Figure.4C_marker_expr_table.csv")
df["groups"] = np.floor(df["groups"])
df = df.groupby('groups', as_index=False).mean()
df = df.drop(["Unnamed: 0","Rd3"], axis=1).set_index("groups")
df = stats.zscore(df)

#Visualization
sns.set(style="ticks",
    font_scale=2.5,
    font='helvetica',
    rc = {'figure.figsize':(20,8), 'axes.linewidth':1.5})

sns.heatmap(data=df.transpose(),
    cmap="GnBu",
    center=0,
    vmin=-2, vmax=2,
    xticklabels=50, yticklabels=1,
    )

plt.xlabel("Cell Types", fontsize=24)
plt.ylabel("Marker", fontsize=24)
plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],fontsize=24,rotation=90)
plt.yticks(fontsize=24, rotation=150, ha="right")
plt.tight_layout()

plt.savefig("./figures/Figure.S4A_Cell_Types_Marker.eps",format='eps',dpi=300)
plt.savefig("./figures/Figure.S4A_Cell_Types_Marker.png",format='png',dpi=300)
