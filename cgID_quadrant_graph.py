import pandas as pd
from pathlib import Path
import warnings
import streamlit as st
from scipy import stats
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.express as px
warnings.filterwarnings("ignore")

pickle_path = Path(r'C:\Users\abhin\OneDrive\Documents\GitHub\Research\processed_data.pk1')

table = pd.read_pickle(pickle_path)

data = []
data2 = []
pospos = [0] * 100
posposcg = [str('empty')] * 100
posposrcc= [0] * 100
posposnorm = [0] * 100
posneg = [0] * 100
posnegcg = [str('empty')] * 100
posnegrcc = [0] * 100
posnegnorm = [0] * 100
negneg = [0] * 100
negnegcg = [str('empty')] * 100
negnegrcc= [0] * 100
negnegnorm = [0] * 100
negpos = [0] * 100
negposcg = [str('empty')] * 100
negposrcc= [0] * 100
negposnorm = [0] * 100

# def get_gene_info(cg_sites, df):
#     gene_symbols = []

#     for cg_site in cg_sites:
#         if cg_site in df.index:
#             row = df.loc[cg_site]
#             gene_symbols.append('' if pd.isna(row['Gene_Symbol']) else row['Gene_Symbol'])
#         else:
#             gene_symbols.append(None)

#     return gene_symbols

i = 0
for (columnName, columnData) in table.iteritems():
    if (i >= 17):
        young_rcc = table.loc[((table["rcc"] == "rcc") & (table["age_at_initial_pathologic_diagnosis"] <= 50)), columnName].dropna()
        old_rcc = table.loc[((table["rcc"] == "rcc") & (table["age_at_initial_pathologic_diagnosis"] > 50)), columnName].dropna()
        young_normal = table.loc[((table["rcc"] == "normal") & (table["age_at_initial_pathologic_diagnosis"] <= 50)), columnName].dropna()
        old_normal = table.loc[((table["rcc"] == "normal") & (table["age_at_initial_pathologic_diagnosis"] > 50)), columnName].dropna()
        if not young_rcc.empty and not old_rcc.empty:
            _, p_val = stats.mannwhitneyu(young_rcc, old_rcc)
        if (_, p_val <= 0.01):
            rcc_old = old_rcc.mean()
            rcc_young = young_rcc.mean()
            normal_old = old_normal.mean()
            normal_young = young_normal.mean()
            rcc_diff = rcc_old - rcc_young
            normal_diff = normal_old - normal_young
            data.append((normal_diff, rcc_diff, columnName))
            if (rcc_diff > 0 and normal_diff <= 0):
                curr = rcc_diff - normal_diff
                for i in range(0, 100):
                    if (curr > posneg[i]):
                        posneg.insert(i, curr)
                        posneg.pop()
                        posnegcg.insert(i, columnName)
                        posnegcg.pop()
                        posnegrcc.insert(i, rcc_diff)
                        posnegnorm.insert(i, normal_diff)
                        posnegrcc.pop()
                        posnegnorm.pop()
                        break
            elif (rcc_diff > 0 and normal_diff > 0):
                curr = rcc_diff + normal_diff
                for i in range(0, 100):
                    if (curr > pospos[i]):
                        pospos.insert(i, curr)
                        pospos.pop()
                        posposcg.insert(i, columnName)
                        posposcg.pop()
                        posposrcc.insert(i, rcc_diff)
                        posposnorm.insert(i, normal_diff)
                        posposrcc.pop()
                        posposnorm.pop()
                        break
            elif (rcc_diff <= 0 and normal_diff <=0):
                curr = rcc_diff + normal_diff
                for i in range(0, 100):
                    if (curr < negneg[i]):
                        negneg.insert(i, curr)
                        negneg.pop()
                        negnegcg.insert(i, columnName)
                        negnegcg.pop()
                        negnegrcc.insert(i, rcc_diff)
                        negnegnorm.insert(i, normal_diff)
                        negnegrcc.pop()
                        negnegnorm.pop()
                        break
            else:
                curr = normal_diff - rcc_diff
                for i in range(0, 100):
                    if (curr > negpos[i]):
                        negpos.insert(i, curr)
                        negpos.pop()
                        negposcg.insert(i, columnName)
                        negposcg.pop()
                        negposrcc.insert(i, rcc_diff)
                        negposnorm.insert(i, normal_diff)
                        negposrcc.pop()
                        negposnorm.pop()
                        break
    i = i + 1

for j in range(0, 10):
    data2.append((posposnorm[j], posposrcc[j], posposcg[j]))
    data2.append((posnegnorm[j], posnegrcc[j], posnegcg[j]))
    data2.append((negnegnorm[j], negnegrcc[j], negnegcg[j]))
    data2.append((negposnorm[j], negposrcc[j], negposcg[j]))

df = pd.DataFrame(data, columns=["normal_diff", "rcc_diff", "cgID"])
df_sampled = df.sample(frac=0.05, random_state=42)
df2 = pd.DataFrame(data2, columns=["normal_diff", "rcc_diff", "cgID"])

# gene_symbols, _, _ = get_gene_info(df_sampled["cgID"], df2)


fig = px.scatter(
    df_sampled,
    x="normal_diff",
    y="rcc_diff",
    title="Normal and RCC Differences",
    labels={"normal_diff": "Normal (Old - Young)", "rcc_diff": "RCC (Old - Young)"},
    range_x=[-0.15, 0.15],
    range_y=[-0.25, 0.25],
    marginal_x="histogram",
    marginal_y="histogram",
    opacity=0.4,
    color_discrete_sequence=["orange"],
    hover_name="cgID",
    hover_data={"normal_diff": True, "rcc_diff": True},
    height = 1000,
    width = 1000
)

fig.update_traces(
    marker=dict(color="orange", opacity=1.0), 
    selector=dict(type='histogram'),
)

fig2 = px.scatter(
    df2,
    x="normal_diff",
    y="rcc_diff",
    title="Normal and RCC Differences",
    labels={"normal_diff": "Normal (Old - Young)", "rcc_diff": "RCC (Old - Young)"},
    marginal_x="histogram",
    marginal_y="histogram",
    range_x=[-0.15, 0.15],
    range_y=[-0.25, 0.25],
    opacity=0.4,
    color_discrete_sequence=["orange"],
    hover_name="cgID",
    hover_data={"normal_diff": True, "rcc_diff": True},
    height = 1000,
    width = 1000
)

fig2.update_traces(
    marker=dict(color="black", opacity=1.0),
    selector=dict(type='histogram'),
)

fig.add_traces(fig2.data)

for x_val in [-0.2, -0.1, 0, 0.1, 0.2]:
    fig.add_shape(
        type="line",
        x0=x_val,
        x1=x_val,
        y0=-0.25,
        y1=0.25,
        line=dict(color="lightgray", width=1),
    )

fig.add_shape(
    type="line",
    x0=-0.15,
    x1=0.15,
    y0=-0.25,
    y1=-0.25,
    line=dict(color="black", width=2),
)

fig.add_shape(
    type="line",
    x0=-0.15,
    x1=0.15,
    y0=0.25,
    y1=0.25,
    line=dict(color="black", width=2),
)

fig.add_shape(
    type="line",
    x0=-0.15,
    x1=-0.15,
    y0=-0.25,
    y1=0.25,
    line=dict(color="black", width=2),
)

fig.add_shape(
    type="line",
    x0=0.15,
    x1=0.15,
    y0=-0.25,
    y1=0.25,
    line=dict(color="black", width=2),
)

fig.add_shape(
    type="line",
    x0=0,
    x1=0,
    y0=-1,
    y1=1,
    line=dict(color="black", width=2),
)

fig.add_shape(
    type="line",
    x0=-1,
    x1=1,
    y0=0,
    y1=0,
    line=dict(color="black", width=2),
)

st.plotly_chart(fig)