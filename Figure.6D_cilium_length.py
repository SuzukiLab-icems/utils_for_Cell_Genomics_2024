import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats

df = pd.read_csv("Figure.6D_length_summary_table_30-50.csv",index_col=0)
df_count = df.drop_duplicates()

for type_A in np.unique(df["RD3"].tolist()):
	for type_B in np.unique(df["RD3"].tolist()):
		statistics, pvalue = stats.mannwhitneyu(
			df[df["RD3"] == type_A]["Length(µm)"],
			df[df["RD3"] == type_B]["Length(µm)"]
			)
		with open("Figure.6D_length_stats.txt", "a") as f:
			f.write(
				f'>{type_A}-{type_B}\ntotal cell:{df_count[df_count["RD3"] == type_A]["cell_number"].sum()}/total cilium:{df[df["RD3"] == type_A]["RD3"].count()}/"stats:"{statistics}/"pvalue:"{pvalue}\n'
				)
			f.close()

sns.set(style="ticks",
        font_scale=2.5,
        font='arial',
        rc = {'figure.figsize':(9.5,12), 'axes.linewidth':1.5})

sns.ecdfplot(
	data=df.reset_index(),
	stat = "proportion",
    palette="mako",
    hue="RD3",
    lw=3,
    x="Length(µm)",
    hue_order=["WT","KO","OE"],
    zorder=0
    )
    
plt.xticks(fontsize=38, rotation=0)
plt.yticks(fontsize=38)

plt.xlabel("Length(µm)", fontsize=38, ha="right")
plt.ylabel("Cumulative Fraction", fontsize=38)
plt.tight_layout()
sns.despine()

plt.savefig("Figure.6D_length.eps",format='eps',dpi=300)
plt.savefig("Figure.6D_length.png",format='png',dpi=300)

plt.close()
