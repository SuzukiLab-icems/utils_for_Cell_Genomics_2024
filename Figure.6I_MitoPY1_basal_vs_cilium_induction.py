import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

file = 'Figure.6I_table.xlsx'

"""Figure.6I(MitoPY1 Ctrl vs. Cilium Induction)"""
df_6I = pd.read_excel(file, sheet_name=None)["Figure.6I"]
df_6I["normalzied MitoPY1"] = 10000 * df_6I["MitoPY1"] / df_6I["Hoechst"]

font_size = 32
x_size = 14
y_size = 7
sns.set(style="ticks", font_scale=2.5, font='arial',
		rc = {'figure.figsize':(x_size,y_size), 'axes.linewidth':1.5})
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.55
plt.rcParams['figure.subplot.left'] = 0.1

for serum in ["-","+"]:
	df_serum = df_6I[df_6I["serum"] == serum]
	WT = df_serum[df_serum["RD3"] == "WT"]["normalzied MitoPY1"].tolist()
	KO = df_serum[df_serum["RD3"] == "KO"]["normalzied MitoPY1"].tolist()
	u = stats.mannwhitneyu(WT, KO)
	t = stats.ttest_ind(WT, KO)
	with open('Figure.6I_Figure.6I_vs_Cilium_Induction.txt', "a", encoding='utf-8', newline='\n') as f:
		f.write(f'Serum{serum}:WT vs. KO:mannwhitney={u};ttest={t} \n')

for rd3 in ["WT","KO"]:
	df_rd3 = df_6I[df_6I["RD3"] == rd3]
	ctrl = df_rd3[df_rd3["serum"] == "-"]["normalzied MitoPY1"].tolist()
	cilium = df_rd3[df_rd3["serum"] == "+"]["normalzied MitoPY1"].tolist()
	u = stats.mannwhitneyu(ctrl, cilium)
	t = stats.ttest_ind(ctrl, cilium)
	with open('Figure.6I_vs_Cilium_Induction.txt', "a", encoding='utf-8', newline='\n') as f:
		f.write(f'RD3{rd3}:Ctrl vs. Cilium Induction:mannwhitney={u};ttest={t} \n')

sns.barplot(
	data = df_6I,
	x="condition",
	y="normalzied MitoPY1",
	order = ["WT(+)","WT(-)","KO(+)","KO(-)"],
	palette = ["bisque","coral","bisque","coral"],
	ci=None,
	edgecolor='darkgrey',
	linewidth=1,
	)
sns.stripplot(
	data = df_6I,
	x="condition",
	y="normalzied MitoPY1",
	order = ["WT(+)","WT(-)","KO(+)","KO(-)"],
	color = "ivory",
	s=30,
	jitter = 0.10,
	edgecolor='darkgrey', linewidth=1,
	)

plt.xticks(fontsize=font_size, rotation=0, ha="right")
plt.xlabel("Cilium Induction", fontsize=font_size, rotation=0)
plt.ylabel("Normalized MitoPY1", fontsize=font_size, rotation=90)
plt.yticks(fontsize=font_size, rotation=0, ha="right")
plt.ylim([
	0,
	df_6I["normalzied MitoPY1"].max() + df_6I["normalzied MitoPY1"].std()
])
plt.tight_layout()
plt.savefig("Figure.6I.eps", dpi=300)
plt.savefig("Figure.6I.png", dpi=300)
plt.close()
