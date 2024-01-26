import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df_A1 = pd.read_table("./PoolA/PoolA_1st.txt", sep="\t")
df_A2 = pd.read_table("./PoolA/PoolA_2nd.txt", sep="\t")
df_A3 = pd.read_table("./PoolA/PoolA_3rd.txt", sep="\t")
df_B1 = pd.read_table("./PoolB/PoolB_1st.txt", sep="\t")
df_B2 = pd.read_table("./PoolB/PoolB_2nd.txt", sep="\t")
df_B3 = pd.read_table("./PoolB/PoolB_3rd.txt", sep="\t")

df_A1 = df_A1.rename(columns = { df_A1.columns[2]: '1st' })
df_A2 = df_A2.rename(columns = { df_A2.columns[2]: '2nd' })
df_A3 = df_A3.rename(columns = { df_A3.columns[2]: '3rd' })
df_B1 = df_B1.rename(columns = { df_B1.columns[2]: '1st' })
df_B2 = df_B2.rename(columns = { df_B2.columns[2]: '2nd' })
df_B3 = df_B3.rename(columns = { df_B3.columns[2]: '3rd' })

df_A1 = df_A1.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df_A2 = df_A2.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df_A3 = df_A3.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df_B1 = df_B1.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df_B2 = df_B2.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()
df_B3 = df_B3.groupby('Gene').apply(lambda x: x.sum()).drop('Gene',axis=1).reset_index()

df_A1["Log2(CPM+1)_A_1st"] = np.log2(1 + 1000000 * (df_A1["1st"] / df_A1["1st"].sum()))
df_A2["Log2(CPM+1)_A_2nd"] = np.log2(1 + 1000000 * (df_A2["2nd"] / df_A2["2nd"].sum()))
df_A3["Log2(CPM+1)_A_3rd"] = np.log2(1 + 1000000 * (df_A3["3rd"] / df_A3["3rd"].sum()))
df_B1["Log2(CPM+1)_B_1st"] = np.log2(1 + 1000000 * (df_B1["1st"] / df_B1["1st"].sum()))
df_B2["Log2(CPM+1)_B_2nd"] = np.log2(1 + 1000000 * (df_B2["2nd"] / df_B2["2nd"].sum()))
df_B3["Log2(CPM+1)_B_3rd"] = np.log2(1 + 1000000 * (df_B3["3rd"] / df_B3["3rd"].sum()))

df_A_merge = df_A1.merge(df_A2.loc[:, ["Gene", "Log2(CPM+1)_A_2nd"]], on="Gene")
df_A_merge = df_A_merge.merge(df_A3.loc[:, ["Gene", "Log2(CPM+1)_A_3rd"]], on="Gene")
df_B_merge = df_B1.merge(df_B2.loc[:, ["Gene", "Log2(CPM+1)_B_2nd"]], on="Gene")
df_B_merge = df_B_merge.merge(df_B3.loc[:, ["Gene", "Log2(CPM+1)_B_3rd"]], on="Gene")
df_merge = df_A_merge.merge(df_B_merge.loc[:, ["Gene", "Log2(CPM+1)_B_1st", "Log2(CPM+1)_B_2nd", "Log2(CPM+1)_B_3rd"]], on="Gene")
df_merge = df_merge.loc[:, ["Gene",
    "Log2(CPM+1)_A_1st", "Log2(CPM+1)_A_2nd", "Log2(CPM+1)_A_3rd",
    "Log2(CPM+1)_B_1st", "Log2(CPM+1)_B_2nd", "Log2(CPM+1)_B_3rd"
    ]
]
df_merge["A_Log2FC(2nd.vs.1st)"] = df_merge["Log2(CPM+1)_A_2nd"] - df_merge["Log2(CPM+1)_A_1st"]
df_merge["B_Log2FC(2nd.vs.1st)"] = df_merge["Log2(CPM+1)_B_2nd"] - df_merge["Log2(CPM+1)_B_1st"]
df_merge["A_Log2FC(3rd.vs.2nd)"] = df_merge["Log2(CPM+1)_A_3rd"] - df_merge["Log2(CPM+1)_A_2nd"]
df_merge["B_Log2FC(3rd.vs.2nd)"] = df_merge["Log2(CPM+1)_B_3rd"] - df_merge["Log2(CPM+1)_B_2nd"]
df_sig_A = df_merge[df_merge["A_Log2FC(2nd.vs.1st)"] >= 1.5]
df_sig_B = df_merge[df_merge["B_Log2FC(2nd.vs.1st)"] >= 1.5]
df_sig = pd.concat([df_sig_A, df_sig_B])
df_sig.to_csv("Figure.2G_sig.csv")
df_sig = pd.read_csv("Figure.2G_sig.csv")

df_sig_A_2nd = df_sig.loc[:, ["Gene", "Log2(CPM+1)_A_2nd", "A_Log2FC(3rd.vs.2nd)"]].rename(columns = { "Log2(CPM+1)_A_2nd":"Log2(CPM+1)", "A_Log2FC(3rd.vs.2nd)":"Log2FC"})
df_sig_A_3rd = df_sig.loc[:, ["Gene", "Log2(CPM+1)_A_3rd", "A_Log2FC(3rd.vs.2nd)"]].rename(columns = { "Log2(CPM+1)_A_3rd":"Log2(CPM+1)", "A_Log2FC(3rd.vs.2nd)":"Log2FC"})
df_sig_B_2nd = df_sig.loc[:, ["Gene", "Log2(CPM+1)_B_2nd", "B_Log2FC(3rd.vs.2nd)"]].rename(columns = { "Log2(CPM+1)_B_2nd":"Log2(CPM+1)", "B_Log2FC(3rd.vs.2nd)":"Log2FC"})
df_sig_B_3rd = df_sig.loc[:, ["Gene", "Log2(CPM+1)_B_3rd", "B_Log2FC(3rd.vs.2nd)"]].rename(columns = { "Log2(CPM+1)_B_3rd":"Log2(CPM+1)", "B_Log2FC(3rd.vs.2nd)":"Log2FC"})

df_sig_A_2nd["Pool"] = "A"
df_sig_A_3rd["Pool"] = "A"
df_sig_B_2nd["Pool"] = "B"
df_sig_B_3rd["Pool"] = "B"

df_sig_2nd = pd.concat([df_sig_A_2nd, df_sig_B_2nd])
df_sig_3rd = pd.concat([df_sig_A_3rd, df_sig_B_3rd])

df_sig_2nd = df_sig_2nd.sort_values(by="Log2FC", ascending=False)
df_sig_3rd = df_sig_3rd.sort_values(by="Log2FC", ascending=False)

df_sig_2nd["Screen"] = "1st"
df_sig_3rd["Screen"] = "2nd"

df_sig = pd.concat([df_sig_2nd, df_sig_3rd])
df_sig = df_sig[df_sig["Log2FC"] != 0]
df_sig = df_sig[df_sig["Gene"] != "NonTargetingControlGuideForMouse_0854"] ##library noise
df_sig = df_sig.replace("NonTargetingControlGuideForMouse", "sgNTC", regex=True)
df_sig = df_sig.replace("Gm13305/Il11ra1/Il11ra2/", "Gm13305/Il11ra1/Il11ra2", regex=True)
df_sig = df_sig.replace("Hba-a1/Hba-a2/", "Hba-a1/Hba-a2", regex=True)
df_sig.to_csv("Figure.2G_Significant_Gene_Stats.csv")
df_sig = pd.read_csv("Figure.2G_Significant_Gene_Stats.csv")

df_sig_high = df_sig[(df_sig["Screen"] == "2nd") & (df_sig["Log2FC"] > 0)]
df_sig_high["Exp"] = "High"
df_sig_low = df_sig[(df_sig["Screen"] == "2nd") & (df_sig["Log2FC"] <= 0)]
df_sig_low["Exp"] = "Low"
df_sig_1st = df_sig[df_sig["Screen"] == "1st"]
df_sig_FC = pd.concat([df_sig_1st, df_sig_high, df_sig_low])
df_sig_FC.to_csv("Figure.2G_sig_FC.csv")

df_sig = pd.read_csv("Figure.2G_Significant_Gene_Stats.csv")
df_sig_FC = pd.read_csv("Figure.2G_sig_FC.csv")

sns.set(style="ticks", font_scale=2.5, font='arial',
    rc = {'figure.figsize':(40,10), 'axes.linewidth':1.5})
sns.color_palette("viridis", as_cmap=True)
color_dict = dict({"High":'orchid',
                   "Low":'thistle'
                   }
)
kwargs = {"edgecolor":"black", "linewidth":1.5}

fig,ax = plt.subplots()

ax = sns.barplot(x="Gene", y="Log2(CPM+1)",
    hue="Screen", data=df_sig,
    palette=["royalblue", "midnightblue"],
    linewidth=3,
    edgecolor="black",
    zorder=0)

ax2 = ax.twinx()
ax2.axhline(c='mediumvioletred', lw=3, ls='--', alpha = 0.6, zorder=1)
sns.scatterplot(data=df_sig_FC, x="Gene", y="Log2FC",
    palette=color_dict,
    hue="Exp",
    size="Exp",
    sizes=(400,800),
    legend=False,
    **kwargs,
    zorder=2,
    ax=ax2)

df_1st = df_sig[df_sig["Screen"] == "1st"]
ax.set_ylabel('Log2(CPM+1)', fontsize=30)
ax2.set_ylabel('Foldchange(2nd.vs.1st)', fontsize=30)
ax.set_xticklabels(labels=df_1st["Gene"], rotation=315, ha ='left')
ax.legend(
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    frameon=False,
    fontsize=30,
    fancybox=False,
    edgecolor="black")
plt.tight_layout()
plt.savefig("Figure.2G_Summary.eps",format='eps',dpi=300)
plt.savefig("Figure.2G_Summary.png",format='png',dpi=300)
