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
	rc = {'figure.figsize':(12,6), 'axes.linewidth':1.5}
	)
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25
palette = ["bisque","sandybrown", "coral", "indianred"]

df = pd.read_excel("Figure.S10A-S10G_settings_Summary.xlsx",sheet_name=None)["summary"]

sns.stripplot(
	data=df,
	x="dose(mg/kg)",y="Gcna1+",
	s=25, edgecolor="black", linewidth=0.5,
	palette = palette
)

plt.xlabel("Busulfan(mg/kg)",fontsize=34)
plt.ylabel(f'%Population of Gcna1+',fontsize=34)
plt.xticks(fontsize=34,rotation=0)
plt.yticks(fontsize=34,rotation=0)
plt.ylim([0,40])
plt.tight_layout()
sns.despine()
plt.savefig("Figure.S10C/Figure.S10C_Settings.poplation.gcna1+.eps", dpi=300)
plt.savefig("Figure.S10C/Figure.S10C_Settings.poplation.gcna1+.png", dpi=300)
plt.close()
