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
	rc = {'figure.figsize':(12,10), 'axes.linewidth':1.5}
	)
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25
palette = ["thistle","orchid","thistle","orchid"]
order = ['shLacZ:NAC(-)','shLacZ:NAC(+)','shRd3:NAC(-)','shRd3:NAC(+)']

df = pd.read_excel("Figure.S7-S8_summary_table.xlsx",sheet_name=None)
for dict in df.keys():
	table = df[dict]
	table = table[(table['injection_quality(0:×,1:∆,2:○)'] == 2) & (table['Week'] == 8)]
	table['class'] = table['shRNA'] + str(':NAC(') + table['NAC'] + str(')')
	start = [n for n,k in enumerate(table.columns) if k == 'rep'][0] + 1
	with open(f'./Figure.S8/Figure.S8_stats_NAC_{dict}.txt', 'w') as f:
		for s in table.columns[start:]:
			if table[s].dtype != 'float64': continue
			sns.swarmplot(
				data=table,
				x='class',y=s,
				order=order,
				s=35, color='white',edgecolor="black", linewidth=0.5
				)
			sns.violinplot(
				data=table,
				x='class',
				y=s,
				order=order,
				palette=palette,
				inner = None
				)
			plt.xlabel('class',fontsize=34)
			plt.ylabel(f'{s}',fontsize=34)
			plt.xticks(fontsize=34,rotation=30)
			plt.yticks(fontsize=34,rotation=0)
			if s == 'mouse_weight(mg)': plt.ylim([0,1.1*table[s].max()])
			plt.legend(
				bbox_to_anchor=(0.98, 1.05),
				loc='upper left',
				frameon=False,
				fontsize=34,
				fancybox=False,
				edgecolor="black")
			plt.tight_layout()
			sns.despine()
			plt.savefig(f'./Figure.S8/Figure.S8_shRd3_NAC_{dict}_{s}.eps', dpi=300)
			plt.savefig(f'./Figure.S8/Figure.S8_shRd3_NAC_{dict}_{s}.png', dpi=300)
			plt.close()
			'''stats'''
			for A in order:
				value_A = table[table['class'] == A][s].tolist()
				for B in order:
					value_B = table[table['class'] == B][s].tolist()
					statistic,p_value = stats.ttest_ind(value_A, value_B)
					_ = f.write(f'{dict}|{s}|{A} vs. {B}|t-value:{statistic}|p-value:{p_value}\n')
		f.close()
	'''Index(shRd3/shLacZ) Comparison Between NAC(-), NAC(+)'''
	if dict == 'weight_number':
		data = df[dict]
		data = data[(data['injection_quality(0:×,1:∆,2:○)'] == 2) & (data['Week'] == 8)].dropna()
		df_data = data.pivot_table(values='weight_index', columns=['shRNA','NAC'], index='rep')
		df_data['NAC(+)'], df_data['NAC(-)'] = \
			df_data[('shRd3', '+')] / df_data[('shLacZ', '+')], df_data[('shRd3', '-')] / df_data[('shLacZ', '-')]
		pos, neg = \
			df_data.loc[:,['NAC(+)']].dropna().rename(columns={"NAC(+)":"Testis Weight(shRd3/shLacZ)"}),\
			df_data.loc[:,['NAC(-)']].dropna().rename(columns={"NAC(-)":"Testis Weight(shRd3/shLacZ)"})
		pos["NAC"], neg["NAC"] = "NAC(+)", "NAC(-)"
		'''ratio'''
		df_plot = pd.DataFrame(pd.concat([pos,neg],axis=0))
		sns.violinplot(
			data=df_plot,
			x='NAC',
			y="Testis Weight(shRd3/shLacZ)",
			order = ['NAC(-)','NAC(+)'],
			palette = ["thistle","orchid"],
			inner = None
			)
		sns.swarmplot(
			data=df_plot,
			x="NAC",y="Testis Weight(shRd3/shLacZ)",
			order = ['NAC(-)','NAC(+)'],
			color='white',s=35,linewidth=0.5,edgecolor='black'
			)
		plt.xlabel('NAC',fontsize=34)
		plt.ylabel(f'Testis Weight(shRd3/shLacZ)',fontsize=34)
		plt.xticks(fontsize=34,rotation=30)
		plt.yticks(fontsize=34,rotation=0)
		plt.ylim([-0.5,1.5])
		plt.legend(
			bbox_to_anchor=(0.98, 1.05),
			loc='upper left',
			frameon=False,
			fontsize=34,
			fancybox=False,
			edgecolor="black")
		plt.tight_layout()
		sns.despine()
		plt.savefig(f'./Figure.S8/Figure.S8_Index(shRd3_vs_shLacZ)_Comparison_Between_NAC(-),NAC(+).eps', dpi=300)
		plt.savefig(f'./Figure.S8/Figure.S8_Index(shRd3_vs_shLacZ)_Comparison_Between_NAC(-),NAC(+).png', dpi=300)
		plt.close()
		with open(f'./Figure.S8/Figure.S8_stats_Index(shRd3_vs_shLacZ)_Comparison_Between_NAC(-),NAC(+).txt', 'w') as f:
			statistic,p_value = stats.ttest_ind(
				neg['Testis Weight(shRd3/shLacZ)'].tolist(),
				pos['Testis Weight(shRd3/shLacZ)'].tolist()
			)
			_ = f.write(f'Index(shRd3_vs_shLacZ)_Comparison_Between_NAC(-),NAC(+)|{dict}|t-value:{statistic}|p-value:{p_value}\n')
		f.close()

