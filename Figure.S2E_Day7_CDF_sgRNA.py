#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns

#import files
df3 = pd.read_table("./F2-3_1_D7N3/mageck_result/F2-3_1.count.txt", index_col=0)
df4 = pd.read_table("./F2-4_1_D7N9/mageck_result/F2-4_1.count.txt", index_col=0)
df5 = pd.read_table("./INPUT/INPUT.count.txt", index_col=0)

df3 = df3.loc[:, ["sample1"]]
df4 = df4.loc[:, ["sample1"]]
df5 = df5.loc[:, ["sample1"]]

df3_del_nocount = df3[df3.sample1 !=0 ]
df4_del_nocount = df4[df4.sample1 !=0 ]
df5_del_nocount = df5[df5.sample1 !=0 ]

df3_fix = df3_del_nocount.dropna()
df4_fix = df4_del_nocount.dropna()
df5_fix = df5_del_nocount.dropna()

df3_fix.to_csv('Figure.S2E_Day7_N3_count.csv')
df4_fix.to_csv('Figure.S2E_Day7_N9_count.csv')
df5_fix.to_csv('Figure.S2E_INPUT_count.csv')

#calculate log2CPM
sum_df3 = df3_fix.sum()
sum_df4 = df4_fix.sum()
sum_df5 = df5_fix.sum()

print(sum_df3, sum_df4, sum_df5)

df3_CPM = 10**6 * df3_fix / sum_df3
df3_log2_CPM = np.log2(df3_CPM)
df3_log2_CPM.to_csv('Figure.S2E_Day7_N3_normalized_result.csv')

df4_CPM = 10**6 * df4_fix / sum_df4
df4_log2_CPM = np.log2(df4_CPM)
df4_log2_CPM.to_csv('Figure.S2E_Day7_N9_normalized_result.csv')

df5_CPM = 10**6 * df5_fix / sum_df5
df5_log2_CPM = np.log2(df5_CPM)
df5_log2_CPM.to_csv('Figure.S2E_INPUT_normalized_result.csv')

df_x1 = pd.read_csv("Figure.S2E_Day7_N3_normalized_result.csv", index_col=0)
df_x2 = pd.read_csv("Figure.S2E_Day7_N9_normalized_result.csv", index_col=0)
df_x5 = pd.read_csv("Figure.S2E_INPUT_normalized_result.csv", index_col=0)

print(df_x1['sample1'].count(), df_x2['sample1'].count(), df_x5['sample1'].count())

df_x1['condition'] = 'Day7_N3'
df_x2['condition'] = 'Day7_N9'
df_x5['condition'] = 'INPUT'

#df_merge = pd.concat([df_x1, df_x2, df_x5])
df_merge = pd.concat([df_x2, df_x5])
df_merge = df_merge.rename(columns={'sample1': 'log2CPM'})

df_merge.to_csv("Figure.S2E_merge.csv")

sns.set(style="ticks", font_scale=2.5, font='helvetica',
    rc = {'figure.figsize':(8,10), 'axes.linewidth':1.5})
sns.color_palette("viridis_r", as_cmap=True)

sns.ecdfplot(data=df_merge,
    palette="viridis_r",
    x="log2CPM",
    hue="condition",
    lw=6,
    stat="count",
    legend=False)

plt.xticks(fontsize=34)
plt.yticks(fontsize=34)

plt.ylim([0,65000])
plt.xlim([-5,10])
plt.xlabel("log2CPM", fontsize=34)
plt.ylabel("count", fontsize=34)
plt.title('Day7', fontsize=34)
plt.tight_layout()
sns.despine()

plt.savefig("Figure.S2E_Day7_sgRNA.png",format='png',dpi=300)
plt.savefig("Figure.S2E_Day7_sgRNA.eps",format='eps',dpi=300)

