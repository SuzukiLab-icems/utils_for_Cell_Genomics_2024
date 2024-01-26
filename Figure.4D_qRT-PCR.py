import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns
from matplotlib.font_manager import FontProperties
from datetime import datetime
import matplotlib.ticker as ticker

df = pd.read_csv("Figure.4D_qRT-PCR.csv")
df_Rd3 = df[df['gene'] == 'Rd3']

sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(12,10), 'axes.linewidth':1.5})
sns.color_palette("viridis_r", as_cmap=True)

sns.regplot(
    data = df_Rd3,
    x="week",y="(-2∆Ct)",
    order=3,
    x_jitter=.3,
    color = 'lavender',
	scatter_kws={'edgecolor':'darkslateblue','s':500},
	line_kws={'color':'mediumorchid','linewidth':5,'ls':'-.'},
    ci=None)

plt.tick_params(labelsize = 38, width = 1.5)
plt.xticks([1,2,3,4,5,6,7,8],fontsize=38, rotation=0, ha ='right')
plt.xlabel("Week", fontsize=38)
plt.ylabel("-∆Ct**2(Rd3/Actb)", fontsize=38)
plt.ylim([-0.005,0.065])
plt.tight_layout()
sns.despine()

plt.tight_layout()
plt.savefig("Figure.4D.eps",format='eps',dpi=300)
plt.savefig("Figure.4D.png",format='png',dpi=300)
