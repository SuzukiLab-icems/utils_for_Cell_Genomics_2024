import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def data_processing(file):
	df = pd.read_excel(file, sheet_name=None)
	#first
	first = df["GO_1st"]
	first["-Log10(p-value)"] = -1 * np.log10(first["PValue"])
	first = first.loc[:,["GO_term","Fold Enrichment","-Log10(p-value)"]]
	first["round"] = "1st"
	#second
	second = df["GO_2nd"]
	second["-Log10(p-value)"] = -1 * np.log10(second["PValue"])
	second = second.loc[:,["GO_term","Fold Enrichment","-Log10(p-value)"]]
	second["round"] = "2nd"
	#merged for pivot table
	df_merge = pd.concat([first,second])
	p_heatmap = df_merge.pivot_table(values="-Log10(p-value)", index="GO_term", columns="round").fillna(-1 * np.log10(1))
	p_heatmap = p_heatmap.sort_values("2nd", ascending = False).reset_index()
	es_heatmap = df_merge.pivot_table(values="Fold Enrichment", index="GO_term", columns="round").fillna(0)
	es_heatmap["Enrichment Score"] = np.log2(0.1 + es_heatmap["2nd"]) - np.log2(0.1 + es_heatmap["1st"])
	es_label = es_heatmap.loc[:,["Enrichment Score"]].reset_index()
	analysis_table = p_heatmap.merge(es_label, how="left", on="GO_term")
	return analysis_table

def plot(table, font_size, width, height):
	data = table.head(40)
	plt.rcParams['figure.subplot.top'] = 0.95
	plt.rcParams['figure.subplot.bottom'] = 0.15
	plt.rcParams['figure.subplot.right'] = 0.98
	plt.rcParams['figure.subplot.left'] = 0.82
	sns.set(
		style="ticks",
		font_scale=2.5,
		font='arial',
		rc = {'figure.figsize':(width, height), 'axes.linewidth':1.5}
		)
	sns.heatmap(
		data=data,
		vmin=table["2nd"].min(),
		center=table["2nd"].mean(),
		vmax=table["2nd"].max(),
		cmap="GnBu"
		)
	plt.xlabel("round of screen", fontsize=font_size, rotation=270)
	plt.ylabel("GO terms", fontsize=font_size, rotation=180)
	plt.xticks(fontsize=font_size, rotation=270)
	plt.yticks(fontsize=font_size, rotation=180, ha="right")
	plt.tight_layout()

if __name__ == "__main__":
	file = "Figure.S2I_summary.xlsx"
	analysis_table = data_processing(file)
	heatmap_table = analysis_table.loc[:,["GO_term","1st","2nd"]].set_index("GO_term")
	plot(table = heatmap_table, font_size=70, width=39, height=60)
	plt.savefig('Figure.S2I_GOanalysis.eps',dpi=300)
	plt.savefig('Figure.S2I_GOanalysis.png',dpi=300)
	plt.close()
	
