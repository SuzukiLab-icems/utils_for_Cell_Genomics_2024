import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter('ignore')

sns.set(
	style="ticks",
	font_scale=2.5,
	font='arial',
	rc = {'figure.figsize':(7.5,4.5), 'axes.linewidth':1.5}
	)
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25
palette = ["bisque","sandybrown", "coral", "indianred"]

df = pd.read_excel("Figure.S10A-S10G_settings_Summary.xlsx",sheet_name=None)["summary"]
plzf = [df.columns[n] for n in np.arange(len(df.columns)) if df.columns.str.contains("Type")[n] == True]

for p in plzf:
	sns.stripplot(
		data=df,
		x="dose(mg/kg)",y=p,
		s=20, edgecolor="black", linewidth=0.5,
		palette = palette
	)
	plt.xlabel("Busulfan(mg/kg)",fontsize=34)
	plt.ylabel(f'%Population',fontsize=34)
	plt.xticks(fontsize=34,rotation=0)
	plt.yticks(fontsize=34,rotation=0)
	plt.ylim([0,1.08*df[p].max()])
	plt.tight_layout()
	sns.despine()
	plt.savefig(f'Figure.S10F-S10G/Figure.S10F-S10G_Settings.poplation.{p}.eps', dpi=300)
	plt.savefig(f'Figure.S10F-S10G/Figure.S10F-S10G_Settings.poplation.{p}.png', dpi=300)
	plt.close()
