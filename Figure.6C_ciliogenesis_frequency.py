import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

df = pd.read_csv("Figure.6C_summary_table.csv",index_col=0)

for type_A in np.unique(df["RD3"].tolist()):
	for type_B in np.unique(df["RD3"].tolist()):
		statistics, pvalue = stats.mannwhitneyu(
			df[df["RD3"] == type_A]["%freq"],
			df[df["RD3"] == type_B]["%freq"]
			)
		with open("Figure.6C_frequency_stats.txt", "a") as f:
			f.write(
				f'>{type_A}-{type_B}\ntotal cell:{df[df["RD3"] == type_A]["cell_number"].sum()}/"stats:"{statistics}/"pvalue:"{pvalue}\n'
				)
			f.close()

sns.set(style="ticks",
        font_scale=2.5,
        font='arial',
        rc = {'figure.figsize':(9.5,12), 'axes.linewidth':1.5})

sns.violinplot(data=df.reset_index(),
    palette="mako",
    x="RD3",
    y="%freq",
    order=[
    "WT","KO","OE"
    ],
    zorder=0
    )

sns.stripplot(data=df.reset_index(),
    color="white",
    x="RD3",
    y="%freq",
    s=15,
    order=[
	"WT","KO","OE"
    ],
    edgecolor='black',
    linewidth=1.0,
    zorder=1
    )
    
plt.xticks(fontsize=38, rotation=0)
plt.yticks(fontsize=38)

plt.ylim([-9,51])
plt.xlabel("RD3", fontsize=38, ha="right")
plt.ylabel("%ciliated cells", fontsize=38)
plt.tight_layout()
sns.despine()

plt.savefig("Figure.6C_frequency.eps",format='eps',dpi=300)
plt.savefig("Figure.6C_frequency.png",format='png',dpi=300)

plt.close()
