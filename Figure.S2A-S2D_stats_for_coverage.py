#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns
import glob

sns.set(style="ticks", font_scale=2.5, font='arial',
rc = {'figure.figsize':(12,12), 'axes.linewidth':1.5})
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.0

#import files
files = glob.glob("*/mageck_result/*.count.txt")
for file in files:
	sample=file.split("/")[0]
	df = pd.read_table(file, index_col=0)
	df = df[df["sample1"]>0]
	df_count = df.groupby(by="Gene").count().rename(columns={"sample1":"count"})
	count = df_count.reset_index().groupby(by="count").count()
	print(f'{file}:{count}')
	sns.countplot(data=df_count, x="count", palette="viridis")
	plt.xlabel("sgRNA Counts per Gene" ,fontsize=38)
	plt.ylabel("Number of Genes",fontsize=38)
	plt.xticks(fontsize=38,rotation=0)
	plt.yticks(fontsize=38,rotation=0)
	plt.tight_layout()
	plt.savefig(f'Figure.S2C-S2D_{sample}_coverage.eps', dpi=300)
	plt.savefig(f'Figure.S2C-S2D_{sample}_coverage.png', dpi=300)
	plt.close()
	#distribution
	df_sum = df.groupby(by="Gene").sum().rename(columns={"sample1":"total_count"})
	df_dist = df_sum.merge(df_count.loc[:,["count"]].reset_index(),how="left",on="Gene")
	df_dist["Log2(CPM+1)"] = np.log2(1 + 1000000 * df_dist["total_count"] / df_dist["total_count"].sum())
	sns.histplot(data=df_dist,x="Log2(CPM+1)",hue="count", palette="viridis")
	plt.xlabel("Log2(CPM+1)" ,fontsize=38)
	plt.ylabel("Number of Genes",fontsize=38)
	plt.xticks(fontsize=38,rotation=0)
	plt.yticks(fontsize=38,rotation=0)
	plt.tight_layout()
	plt.axvline(x=2, linewidth=1.5, c='mediumorchid')
	plt.savefig(f'Figure.S2A-S2B_{sample}_count_histgram.eps', dpi=300)
	plt.savefig(f'Figure.S2A-S2B_{sample}_count_histgram.png', dpi=300)
	plt.close()
