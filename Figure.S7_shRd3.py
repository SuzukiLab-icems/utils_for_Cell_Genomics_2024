import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import warnings
warnings.simplefilter('ignore')

sns.set(
	style="ticks",
	font_scale=2.5,
	font='arial',
	rc = {'figure.figsize':(15,10), 'axes.linewidth':1.5}
	)
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25
palette = ["bisque","indianred"]

df = pd.read_excel("Figure.S7-S8_summary_table.xlsx",sheet_name=None)
for dict in df.keys():
	table = df[dict]
	table = table[(table['NAC'] == '-') & (table['injection_quality(0:×,1:∆,2:○)'] != 0)]
	start = [n for n,k in enumerate(table.columns) if k == 'rep'][0] + 1
	with open(f'Figure.S7/Figure.S7_stats_{dict}.txt', 'w') as f:
		for s in table.columns[start:]:
			if table[s].dtype != 'float64': continue
			table = table[table['NAC'] == '-']
			sns.swarmplot(
				data=table,
				x='Week',y=s,
				hue='shRNA',
				s=35, edgecolor="black", linewidth=0.5,
				palette = palette
				)
			plt.xlabel('Week',fontsize=34)
			plt.ylabel(f'{s}',fontsize=34)
			plt.xticks(fontsize=34,rotation=0)
			plt.yticks(fontsize=34,rotation=0)
			plt.ylim([0,1.05*table[s].max()])
			plt.legend(
				bbox_to_anchor=(0.98, 1.05),
				loc='upper left',
				frameon=False,
				fontsize=34,
				fancybox=False,
				edgecolor="black")
			plt.tight_layout()
			sns.despine()
			plt.savefig(f'Figure.S7/Figure.S7_shRd3_{dict}_{s}.eps', dpi=300)
			plt.savefig(f'Figure.S7/Figure.S7_shRd3_{dict}_{s}.png', dpi=300)
			plt.close()
			for day in np.unique(table['Week']):
				df_pnd = table[table['Week'] == day]
				shLacZ = df_pnd[df_pnd['shRNA'] == 'shLacZ'].dropna()[s].tolist()
				shRd3 = df_pnd[df_pnd['shRNA'] == 'shRd3'].dropna()[s].tolist()
				statistic,p_value = stats.ttest_ind(shLacZ, shRd3)
				_ = f.write(f'{dict}|{s}|{day}week|t-value:{statistic}|p-value:{p_value}\n')
		f.close()
