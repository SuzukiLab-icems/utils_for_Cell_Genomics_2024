import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot(data):
	plt.rcParams['figure.subplot.top'] = 0.95
	plt.rcParams['figure.subplot.bottom'] = 0.15
	plt.rcParams['figure.subplot.right'] = 0.98
	plt.rcParams['figure.subplot.left'] = 0.82
	sns.set(
		style="ticks",
		font_scale=2.5,
		font='arial',
		rc = {'figure.figsize':(10,8), 'axes.linewidth':1.5}
		)
	sns.heatmap(
		data=data,
		vmin=0,
		center=1.3,
		vmax=2,
		cmap="GnBu"
		)
	plt.xlabel("Species", fontsize=34, rotation=0)
	plt.ylabel("GO terms", fontsize=34)
	plt.xticks(fontsize=34, rotation=20, ha="right")
	plt.yticks(fontsize=34, rotation=0, ha="right")
	plt.tight_layout()
	plt.savefig('IO/Figure.5D.eps',format='eps',dpi=300)
	plt.savefig('IO/Figure.5D.png',format='png',dpi=300)
	plt.close()

def main():
	df = pd.read_excel('IO/Figure.5D_GO_summary_table.xlsx', sheet_name=None)
	df_human = df["H.S.DIRECT"].loc[:,["GO_id","GO_terms","Benjamini"]]
	df_human["species"] = "Human"
	df_mouse = df["M.S.DIRECT"].loc[:,["GO_id","GO_terms","Benjamini"]]
	df_mouse["species"] = "Mouse"
	df_analysis = pd.concat([df_mouse,df_human])
	df_analysis = df_analysis.pivot_table(values="Benjamini",columns="species",index="GO_terms").loc[:,["Mouse","Human"]].fillna(1)
	df_analysis = -1*np.log10(df_analysis)
	df_analysis = df_analysis.sort_values(by="Mouse", ascending=False)
	df_analysis = df_analysis[df_analysis["Mouse"] > 0].head(10)
	plot(data=df_analysis)

if __name__ == "__main__":
	main()

