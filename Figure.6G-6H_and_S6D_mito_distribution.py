import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

file = 'summary_for_analysis.csv'
font_size = 46
x_size = 12
y_size = 12.5

df = pd.read_csv(file)

df_diff = df.pivot_table(values = "Intensity", index = ["Position","rep","Serum"], columns = ["RD3"])
df_diff["diff"] = df_diff["KO"] - df_diff["WT"]
df_diff = df_diff.reset_index()

sns.set(style="ticks", font_scale=2.5, font='arial',
		rc = {'figure.figsize':(x_size,y_size), 'axes.linewidth':1.5})
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.55
plt.rcParams['figure.subplot.left'] = 0.1

"""Density"""
for serum in np.unique(df["Serum"].tolist()):
	df_serum = df[df["Serum"] == serum]
	df_WT = df_serum[df_serum["RD3"] == "WT"]
	df_KO = df_serum[df_serum["RD3"] == "KO"]
	#p-value
	for pos in np.unique(df["Position"].tolist()):
		WT = df_WT[df_WT["Position"] == pos]["Intensity"].tolist()
		KO = df_KO[df_KO["Position"] == pos]["Intensity"].tolist()
		u = stats.mannwhitneyu(WT, KO)
		t = stats.ttest_ind(WT, KO)
		with open(f'Stats_Figure.6G_and_S6D_Serum{serum}.txt', "a", encoding='utf-8', newline='\n') as f:
			f.write(f'Serum({serum});Position_{pos};WT vs. KO;mannwhitney={u};ttest={t} \n')
	df_serum["pos"] = df_serum["Position"] -1
	sns.lineplot(
		data = df_serum,
		x="pos",
		y="Intensity",
		hue = "RD3",
		hue_order = ["WT", "KO"],
		palette = "mako",
		lw=5,
		ci = None,
		linestyle='dashdot'
		)
	sns.stripplot(
		data = df_serum,
		x="Position",
		y="Intensity",
		hue = "RD3",
		hue_order = ["WT", "KO"],
		palette = "mako",
		s = 25,
		jitter = 0.1,
		edgecolor='lightgrey', linewidth=0.5,
		)
	plt.xlabel("Distance from gTubulin", fontsize=font_size, rotation=0)
	plt.ylabel("Positional Density", fontsize=font_size, rotation=90)
	plt.xticks(fontsize=font_size, rotation=0, ha="right")
	plt.yticks(fontsize=font_size, rotation=0, ha="right")
	plt.ylim([
		df_serum["Intensity"].min() - df_serum["Intensity"].std(),
		df_serum["Intensity"].max() + df_serum["Intensity"].std()
	])
	plt.tight_layout()
	plt.savefig(f'Figure.6G_and_S6D_Serum{serum}.eps', dpi=300)
	plt.savefig(f'Figure.6G_and_S6D_Serum{serum}.png', dpi=300)
	plt.close()

"""Difference"""
for pos in np.unique(df_diff["Position"]):
	df_stat = df_diff[df_diff["Position"] == pos]
	Ctrl = df_stat[df_stat["Serum"] == "-"]["diff"].tolist()
	Starvation = df_stat[df_stat["Serum"] == "+"]["diff"].tolist()
	u = stats.mannwhitneyu(Ctrl, Starvation)
	t = stats.ttest_ind(Ctrl, Starvation)
	with open(f'Stats_Figure.6H_distribution_difference.txt', "a", encoding='utf-8', newline='\n') as f:
		f.write(f'Position({pos});Ctrl-Starvation;mannwhitney={u};ttest={t} \n')
sns.barplot(
	data = df_diff,
	x="Position",
	y="diff",
	hue = "Serum",
	hue_order = ["+", "-"],
	palette = ["bisque","coral"],
	ci = None,
	edgecolor='lightgrey', linewidth=0.5,
	)
sns.stripplot(
	data = df_diff,
	x="Position",
	y="diff",
	hue = "Serum",
	hue_order = ["+", "-"],
	palette = ["bisque","coral"],
	s = 25,
	jitter = 0.1,
	edgecolor='lightgrey', linewidth=0.5,
	)
plt.xlabel("Distance from gTubulin", fontsize=font_size, rotation=0)
plt.ylabel("Difference(KO-WT)", fontsize=font_size, rotation=90)
plt.xticks(fontsize=font_size, rotation=0, ha="right")
plt.yticks(fontsize=font_size, rotation=0, ha="right")
plt.ylim([
	df_diff["diff"].min() - df_diff["diff"].std(),
	df_diff["diff"].max() + df_diff["diff"].std()
])
plt.tight_layout()
plt.savefig('Figure.6H_difference.eps', dpi=300)
plt.savefig('Figure.6H_difference.png', dpi=300)
plt.close()
