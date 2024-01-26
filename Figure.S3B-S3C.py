import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

sns.set(
	style="ticks",
	font_scale=2.5,
	font='arial',
	rc = {'figure.figsize':(8,5), 'axes.linewidth':1.5}
	)
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25

df = pd.read_excel("Figure.S3B-S3C_230912_Summary.xlsx",sheet_name=None)["for_analysis"]

for gene in np.unique(df["gene"].tolist()):
	data = df[df["gene"] == gene]
	sns.stripplot(
		data=data,
		x="weeks",y="-2∆Ct",
		s=25,
		color="mediumorchid",
		linewidth=0.5,
		edgecolor="black")
	plt.xlabel(f'week',fontsize=34)
	plt.ylabel(f'2-∆Ct(vs. Actb)',fontsize=34)
	plt.xticks(fontsize=34,rotation=0)
	plt.yticks(fontsize=34,rotation=0)
	plt.ylim([0.0000001, 0.001])
	plt.yscale('log')
	plt.tight_layout()
	sns.despine()
	plt.savefig(f'Figure.S3B-S3C_{gene}_expression.eps',dpi=300)
	plt.savefig(f'Figure.S3B-S3C_{gene}_expression.png',dpi=300)
	plt.close()


