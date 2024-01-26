import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns

#Data Loading
df = pd.read_csv("Figure.S3D_KDeff_data_for_parse.csv")

'''Ovol2'''
df_Ovol2 = df[df['target'] == 'Ovol2']
df_Ovol2 = df_Ovol2.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Ovol2['normLacZ'] = df_Ovol2['sh.LacZ'] / df_Ovol2['sh.LacZ']
df_Ovol2['normOvol2'] = df_Ovol2['sh.Ovol2'] / df_Ovol2['sh.LacZ']

df_normLacZ = df_Ovol2.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normOvol2 = df_Ovol2.loc[:,['normOvol2']].rename(columns={'normOvol2':'(-2**Ct)'})
df_normOvol2['shrna'] = 'sh.Ovol2'

df_Ovol2 = pd.concat([df_normLacZ,df_normOvol2])

'''Rnf215'''
df_Rnf215 = df[df['target'] == 'Rnf215']
df_Rnf215 = df_Rnf215.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Rnf215['normLacZ'] = df_Rnf215['sh.LacZ'] / df_Rnf215['sh.LacZ']
df_Rnf215['normRnf215'] = df_Rnf215['sh.Rnf215'] / df_Rnf215['sh.LacZ']

df_normLacZ = df_Rnf215.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normRnf215 = df_Rnf215.loc[:,['normRnf215']].rename(columns={'normRnf215':'(-2**Ct)'})
df_normRnf215['shrna'] = 'sh.Rnf215'

df_Rnf215 = pd.concat([df_normLacZ,df_normRnf215])

'''Cldn34c4'''
df_Cldn34c4 = df[df['target'] == 'Cldn34c4']
df_Cldn34c4 = df_Cldn34c4.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Cldn34c4['normLacZ'] = df_Cldn34c4['sh.LacZ'] / df_Cldn34c4['sh.LacZ']
df_Cldn34c4['normCldn34c4'] = df_Cldn34c4['sh.Cldn34c4'] / df_Cldn34c4['sh.LacZ']

df_normLacZ = df_Cldn34c4.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normCldn34c4 = df_Cldn34c4.loc[:,['normCldn34c4']].rename(columns={'normCldn34c4':'(-2**Ct)'})
df_normCldn34c4['shrna'] = 'sh.Cldn34c4'

df_Cldn34c4 = pd.concat([df_normLacZ,df_normCldn34c4])

'''Rd3'''
df_Rd3 = df[df['target'] == 'Rd3']
df_Rd3 = df_Rd3.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Rd3['normLacZ'] = df_Rd3['sh.LacZ'] / df_Rd3['sh.LacZ']
df_Rd3['normRd3'] = df_Rd3['sh.Rd3'] / df_Rd3['sh.LacZ']

df_normLacZ = df_Rd3.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normRd3 = df_Rd3.loc[:,['normRd3']].rename(columns={'normRd3':'(-2**Ct)'})
df_normRd3['shrna'] = 'sh.Rd3'

df_Rd3 = pd.concat([df_normLacZ,df_normRd3])

'''Cebpg'''
df_Cebpg = df[df['target'] == 'Cebpg']
df_Cebpg = df_Cebpg.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Cebpg['normLacZ'] = df_Cebpg['sh.LacZ'] / df_Cebpg['sh.LacZ']
df_Cebpg['normCebpg'] = df_Cebpg['sh.Cebpg'] / df_Cebpg['sh.LacZ']

df_normLacZ = df_Cebpg.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normCebpg = df_Cebpg.loc[:,['normCebpg']].rename(columns={'normCebpg':'(-2**Ct)'})
df_normCebpg['shrna'] = 'sh.Cebpg'

df_Cebpg = pd.concat([df_normLacZ,df_normCebpg])

'''Rbm26'''
df_Rbm26 = df[df['target'] == 'Rbm26']
df_Rbm26 = df_Rbm26.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_Rbm26['normLacZ'] = df_Rbm26['sh.LacZ'] / df_Rbm26['sh.LacZ']
df_Rbm26['normRbm26'] = df_Rbm26['sh.Rbm26'] / df_Rbm26['sh.LacZ']

df_normLacZ = df_Rbm26.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normRbm26 = df_Rbm26.loc[:,['normRbm26']].rename(columns={'normRbm26':'(-2**Ct)'})
df_normRbm26['shrna'] = 'sh.Rbm26'

df_Rbm26 = pd.concat([df_normLacZ,df_normRbm26])

'''P2ry2'''
df_P2ry2 = df[df['target'] == 'P2ry2']
df_P2ry2 = df_P2ry2.pivot_table(values='(-2**Ct)', columns=['shrna'], index='stats_replicates')
df_P2ry2['normLacZ'] = df_P2ry2['sh.LacZ'] / df_P2ry2['sh.LacZ']
df_P2ry2['normP2ry2'] = df_P2ry2['sh.P2ry2'] / df_P2ry2['sh.LacZ']

df_normLacZ = df_P2ry2.loc[:,['normLacZ']].rename(columns={'normLacZ':'(-2**Ct)'})
df_normLacZ['shrna'] = 'sh.LacZ'
df_normP2ry2 = df_P2ry2.loc[:,['normP2ry2']].rename(columns={'normP2ry2':'(-2**Ct)'})
df_normP2ry2['shrna'] = 'sh.P2ry2'

df_P2ry2 = pd.concat([df_normLacZ,df_normP2ry2])

#Definition
sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(40,8), 'axes.linewidth':1.5})
sns.color_palette("magma", as_cmap=True)

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(1, 7)
labels = ['sh.LacZ', 'sh.Target']

#Ovol2
sns.violinplot(data = df_Ovol2, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax1)
sns.swarmplot(data = df_Ovol2, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax1)

#ax1.set_xlim([-1000,24000])
ax1.set_ylim([-0.05,1.05])
ax1.set_ylabel("Relative Fold Change", fontsize=40)
ax1.set_xlabel("sh.Ovol2", fontsize=40)
ax1.tick_params(labelsize = 40, width = 1.5 )
ax1.set_xticklabels(labels, rotation=30, ha ='right')

#Rnf215
sns.violinplot(data = df_Rnf215, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax2)

sns.swarmplot(data = df_Rnf215, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax2)

#ax2.set_xlim([-1000,24000])
ax2.set_ylim([-0.05,1.05])
ax2.set_ylabel("Relative Fold Change", fontsize=40)
ax2.set_xlabel("sh.Rnf215", fontsize=40)
ax2.tick_params(labelsize = 40, width = 1.5 )
ax2.set_xticklabels(labels, rotation=30, ha ='right')

#Cldn34c4
sns.violinplot(data = df_Cldn34c4, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax3)

sns.swarmplot(data = df_Cldn34c4, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax3)

#ax3.set_xlim([-1000,24000])
ax3.set_ylim([-0.05,1.05])
ax3.set_ylabel("Relative Fold Change", fontsize=40)
ax3.set_xlabel("sh.Cldn34c4", fontsize=40)
ax3.tick_params(labelsize = 40, width = 1.5 )
ax3.set_xticklabels(labels, rotation=30, ha ='right')

#Rd3
sns.violinplot(data = df_Rd3, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax4)

sns.swarmplot(data = df_Rd3, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax4)

#ax4.set_xlim([-1000,24000])
ax4.set_ylim([-0.05,1.05])
ax4.set_ylabel("Relative Fold Change", fontsize=40)
ax4.set_xlabel("sh.Rd3", fontsize=40)
ax4.tick_params(labelsize = 40, width = 1.5 )
ax4.set_xticklabels(labels, rotation=30, ha ='right')

#Cebpg
sns.violinplot(data = df_Cebpg, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax5)

sns.swarmplot(data = df_Cebpg, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax5)

#ax5.set_xlim([-1000,24000])
ax5.set_ylim([-0.05,1.05])
ax5.set_ylabel("Relative Fold Change", fontsize=40)
ax5.set_xlabel("sh.Cebpg", fontsize=40)
ax5.tick_params(labelsize = 40, width = 1.5 )
ax5.set_xticklabels(labels, rotation=30, ha ='right')

#Rbm26
sns.violinplot(data = df_Rbm26, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax6)

sns.swarmplot(data = df_Rbm26, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax6)

#ax6.set_xlim([-1000,24000])
ax6.set_ylim([-0.05,1.05])
ax6.set_ylabel("Relative Fold Change", fontsize=40)
ax6.set_xlabel("sh.Rbm26", fontsize=40)
ax6.tick_params(labelsize = 40, width = 1.5 )
ax6.set_xticklabels(labels, rotation=30, ha ='right')

#P2ry2
sns.violinplot(data = df_P2ry2, x="shrna", y="(-2**Ct)", inner=None, scale="width", bw=0.25, linewidth = 3, edgecolor = "black", size = 20, color="whitesmoke",
    saturation=1, ax=ax7)

sns.swarmplot(data = df_P2ry2, x="shrna", y="(-2**Ct)", linewidth = 3,
    edgecolor = "black", size = 20, palette="magma", ax=ax7)

#ax7.set_xlim([-1000,24000])
ax7.set_ylim([-0.05,1.05])
ax7.set_ylabel("Relative Fold Change", fontsize=40)
ax7.set_xlabel("sh.P2ry2", fontsize=40)
ax7.tick_params(labelsize = 40, width = 1.5 )
ax7.set_xticklabels(labels, rotation=30, ha ='right')

plt.tight_layout()
sns.despine()

plt.savefig("Figure.S3D.eps",format='eps',dpi=400)
plt.savefig("Figure.S3D.png",format='png',dpi=400)
