import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

factors =['%Infection']

shRNA = {
"LacZ":1,
"Ovol2":2,
"Rnf215":3,
"Cldn34c4":4,
"Rd3":5,
"Cebpg":6,
"Rbm26":7,
"P2ry2":8
}

df_analysis = pd.read_csv("./Figure.3_raw_data_for_shRNA_screening/Table.S4_stats.csv")
df_analysis["order"] = df_analysis["Target"].replace(shRNA)
df_analysis = df_analysis.sort_values(by="order", ascending=True)

for values in factors:
	df_values = df_analysis.loc[:,["Target",values]].dropna()
	LacZ = df_values[df_values["Target"] == "LacZ"][values].tolist()
	for target in ["Ovol2","Rnf215","Cldn34c4","Rd3","Cebpg","Rbm26","P2ry2"]:
		Target = df_values[df_values["Target"] == target][values].tolist()
		try:
			u = stats.mannwhitneyu(LacZ, Target)
			t = stats.ttest_ind(LacZ, Target)
			print(f'({values}-{target}):mannwhitney={u},ttest={t}')
			with open("./Figure.S3E_Infectious_Efficiency/Figure.S3E_InfectioousEfficiency_stats.txt", "a", encoding='utf-8', newline='\n') as f:
				f.write(f'({values}-{target}):mannwhitney={u},ttest={t} \n')
		except TypeError:
			print("not appropriate")
	sns.set(style="ticks", font_scale=2.5, font='arial',
		rc = {'figure.figsize':(15,8), 'axes.linewidth':1.5})
	plt.rcParams['figure.subplot.bottom'] = 0.2
	plt.rcParams['figure.subplot.right'] = 0.55
	plt.rcParams['figure.subplot.left'] = 0.1
	sns.violinplot(
		data=df_values,
		x="Target", y=values,
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
		data=df_values,
		x="Target", y=values,
		linewidth = 3,
		edgecolor = "black",
		size = 20,
		palette="mako"
		)
	plt.xlabel("shRNA", fontsize=32, rotation=0)
	plt.ylabel(f'{values}', fontsize=32, rotation=90)
	plt.xticks(fontsize=32, rotation=30, ha="right")
	plt.yticks(fontsize=32, rotation=0, ha="right")
	plt.ylim([
	df_values[values].min() - df_values[values].std(),
	df_values[values].max() + df_values[values].std()
	])
	sns.despine()
	plt.tight_layout()
	plt.savefig(f'./Figure.S3E_Infectious_Efficiency/Figure.S3E_InfectioousEfficiency_{values}.eps', dpi=300)
	plt.savefig(f'./Figure.S3E_Infectious_Efficiency/Figure.S3E_InfectioousEfficiency_{values}.png', dpi=300)
	plt.close()


