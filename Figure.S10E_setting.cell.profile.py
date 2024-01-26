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
        rc = {'figure.figsize':(15,12), 'axes.linewidth':1.5}
        )
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.right'] = 0.95
plt.rcParams['figure.subplot.left'] = 0.25

cell_types = [
"Spermatogonia, Secondary Spermatocytes",
"Primary Spermatocytes",
"Round Spermatids",
"Elongating Spermatids",
"Somatic Cells",
"Activated Somatic Cells"
]

df = pd.read_excel("Figure.S10A-S10G_settings_Summary.xlsx",sheet_name=None)["summary"]
df = df.set_index("dose(mg/kg)").loc[:,cell_types].groupby(by="dose(mg/kg)").mean().transpose()

sns.heatmap(
	data=df,
	annot=True, fmt="1.1f",
	vmin=0,center=15,
	cmap="YlGnBu",
	linecolor="white",
	linewidth=0.25

)
plt.xlabel('dose(mg/kg)',fontsize=34)
plt.ylabel('cell types',fontsize=34)
plt.xticks(fontsize=34)
plt.yticks(fontsize=34)
plt.tight_layout()
plt.savefig(f'Figure.S10E/Figure.S10E_population_matrix.png')
plt.close()
