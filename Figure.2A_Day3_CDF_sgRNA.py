#!/usr/bin/python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn as sns

#import files
df1 = pd.read_table("./F2-1_1_D3N3/mageck_result/F2-1_1.count.txt", index_col=0)
df2 = pd.read_table("./F2-2_1_D3N9/mageck_result/F2-2_1.count.txt", index_col=0)
df5 = pd.read_table("./INPUT/INPUT.count.txt", index_col=0)

df1 = df1.loc[:, ["sample1"]]
df2 = df2.loc[:, ["sample1"]]
df5 = df5.loc[:, ["sample1"]]

df1_del_nocount = df1[df1.sample1 !=0 ]
df2_del_nocount = df2[df2.sample1 !=0 ]
df5_del_nocount = df5[df5.sample1 !=0 ]

df1_fix = df1_del_nocount.dropna()
df2_fix = df2_del_nocount.dropna()
df5_fix = df5_del_nocount.dropna()

df1_fix.to_csv('Figure.2A_Day3_N3_count.csv')
df2_fix.to_csv('Figure.2A_Day3_N9_count.csv')
df5_fix.to_csv('Figure.2A_INPUT_count.csv')

#calculate log2CPM
sum_df1 = df1_fix.sum()
sum_df2 = df2_fix.sum()
sum_df5 = df5_fix.sum()

print(sum_df1, sum_df2, sum_df5)

df1_CPM = 10**6 * df1_fix / sum_df1
df1_log2_CPM = np.log2(df1_CPM)
df1_log2_CPM.to_csv('Figure.2A_Day3_N3_normalized_result.csv')

df2_CPM = 10**6 * df2_fix / sum_df2
df2_log2_CPM = np.log2(df2_CPM)
df2_log2_CPM.to_csv('Figure.2A_Day3_N9_normalized_result.csv')

df5_CPM = 10**6 * df5_fix / sum_df5
df5_log2_CPM = np.log2(df5_CPM)
df5_log2_CPM.to_csv('Figure.2A_INPUT_normalized_result.csv')

#cumulative distribution fraction
df_x1 = pd.read_csv("Figure.2A_Day3_N3_normalized_result.csv")
df_x2 = pd.read_csv("Figure.2A_Day3_N9_normalized_result.csv")
df_x5 = pd.read_csv("Figure.2A_INPUT_normalized_result.csv")

print(df_x1['sample1'].count(), df_x2['sample1'].count(), df_x5['sample1'].count())

df_x1['condition'] = 'Day3_N3'
df_x2['condition'] = 'Day3_N9'
df_x5['condition'] = 'INPUT'

df_merge = pd.concat([df_x1, df_x2, df_x5])
df_merge = df_merge.rename(columns={'sample1': 'log2CPM'})

df_merge.to_csv("Figure.2A_merge.csv")
df_merge = pd.read_csv("Figure.2A_merge.csv")

sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(8,8), 'axes.linewidth':1.5})
sns.color_palette("viridis_r", as_cmap=True)

sns.ecdfplot(data=df_merge,
    palette="viridis_r",
    x="log2CPM",
    hue="condition",
    lw=3,
    stat="count",
    legend=False)

plt.tick_params(labelsize = 30, width = 1.5 )
plt.ylim([0,65000])
plt.xlim([-5,10])
plt.xlabel("log2CPM", fontsize=30)
plt.ylabel("count", fontsize=30)
plt.title('Day3', fontsize=30)
plt.tight_layout()
sns.despine()

plt.legend(labels=["INPUT", "9 Testis","3 Testis"],
    loc='upper left',
    frameon=False,
    fontsize=30,
    fancybox=False,
    edgecolor="black")

plt.savefig("Figure.2A_Day3_sgRNA.eps",format='eps',dpi=300)
plt.savefig("Figure.2A_Day3_sgRNA.png",format='png',dpi=300)

#plt.show()

exit=0
