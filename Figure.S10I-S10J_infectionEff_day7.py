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
        rc = {'figure.figsize':(7,8.5), 'axes.linewidth':1.5}
        )
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25

df = pd.read_excel('Figure.S10H-S10J_Summary_InfectionEff_day7.xlsx',sheet_name=None)['day7']

pnd = 29
days = 7


'''population'''
pop = df[df['PND'] == pnd]
for m in ['%CDH1','%KIT']:
	sns.stripplot(
		data=pop,
		order=['-','+'],
		y=m,x="busulfan",
		s=30, edgecolor="black", linewidth=0.5,
		palette=["skyblue","dodgerblue"]
	)
	plt.xlabel('Busulfan',fontsize=34)
	plt.ylabel(f'{m}',fontsize=34)
	plt.xticks(fontsize=34)
	plt.yticks(fontsize=34)
	plt.tight_layout()
	sns.despine()
	plt.savefig(f'Figure.S10I-S10J/Figure.S10H-S10J_busulfan_infection_population_{m}_{days}.eps')
	plt.savefig(f'Figure.S10I-S10J/Figure.S10H-S10J_busulfan_infection_population_{m}_{days}.png')
	plt.close()
	
'''infection'''
inf = df[df['class'] != 'CTRL']
for m in ['CDH1:%GFP(+)','KIT:%GFP(+)']:
	sns.stripplot(
		data=inf,
		order=['-','+'],
		y=m,x="busulfan",
		s=30, edgecolor="black", linewidth=0.5,
		palette=["skyblue","dodgerblue"]
	)
	plt.xlabel('Busulfan',fontsize=34)
	plt.ylabel(f'{m}',fontsize=34)
	plt.xticks(fontsize=34)
	plt.yticks(fontsize=34)
	plt.tight_layout()
	sns.despine()
	plt.savefig(f'Figure.S10I-S10J/Figure.S10H-S10J_busulfan_infection_population_{m}_{days}.eps')
	plt.savefig(f'Figure.S10I-S10J/Figure.S10H-S10J_busulfan_infection_population_{m}_{days}.png')
	plt.close()
