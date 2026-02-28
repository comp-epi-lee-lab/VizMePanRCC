import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

st.set_page_config(layout="wide")

@st.cache_data
def load_data():
    pickle_path = Path(__file__).parent/'data'/'pickle_file.pk1'
    return pd.read_pickle(pickle_path)

table = load_data()
table = table[table['subtype'] == 'kirc']
df = table.dropna()

df_normals = table[table['rcc'] == 'rcc'].copy()
cpg_columns = [col for col in df_normals.columns if col.startswith('cg')]

print(f"Number of CpG sites: {len(cpg_columns)}")
print(f"Number of normal samples: {len(df_normals)}")

df_normals['mean_methylation'] = df_normals[cpg_columns].mean(axis=1, skipna=True)

df_normals = df_normals.dropna(subset=['age_at_initial_pathologic_diagnosis', 'mean_methylation'])

def assign_age_group(age):
    if 20 <= age < 40:
        return '20-39'
    elif 40 <= age < 50:
        return '40-49'
    elif 50 <= age < 60:
        return '50-59'
    elif 60 <= age < 70:
        return '60-69'
    elif 70 <= age < 200:
        return '70+'
    else:
        return None

df_normals['age_group'] = df_normals['age_at_initial_pathologic_diagnosis'].apply(assign_age_group)
df_normals = df_normals.dropna(subset=['age_group'])

age_groups_order = ['20-39', '40-49', '50-59', '60-69', '70+']
age_groups_present = [ag for ag in age_groups_order if ag in df_normals['age_group'].values]

fig = go.Figure()

for i, age_group in enumerate(age_groups_present):
    data = df_normals[df_normals['age_group'] == age_group]['mean_methylation']
    count = len(data)
    
    fig.add_trace(go.Box(
        y=data,
        x=[i] * count,
        name=age_group,
        boxpoints="all",
        jitter=0.2,
        marker=dict(color="#8977ad"),
        showlegend=False
    ))

tick_labels = [f"{ag}<br>(N = {len(df_normals[df_normals['age_group'] == ag])})" 
               for ag in age_groups_present]

fig.update_layout(
    title='Mean Methylation Ratio Across All CpGs (kirc) by Age Group (Tumor Samples)',
    title_font=dict(size=30),
    yaxis_title='<b>Mean Methylation Ratio</b>',
    yaxis_title_font=dict(size=22),
    xaxis_title='<b>Age Group</b>',
    xaxis_title_font=dict(size=22),
    xaxis=dict(
        tickmode='array',
        tickvals=list(range(len(age_groups_present))),
        ticktext=tick_labels,
        tickfont=dict(size=20)
    ),
    yaxis=dict(
        tickmode='linear',
        dtick=0.01,
        showgrid=False,
        tickfont=dict(size=22)
    ),
    plot_bgcolor='rgba(0, 0, 0, 0)',
    font=dict(
        family="Arial",
        size=24
    ),
    margin=dict(r=20, b=170),
    height=700,
    width=1200
)

fig.update_layout(yaxis_range=[0.42, 0.53])
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

import scipy.stats as stats
valid = df_normals.dropna(subset=['age_at_initial_pathologic_diagnosis', 'mean_methylation'])
valid = valid[valid['age_at_initial_pathologic_diagnosis'] <= 90]
r, p_value = stats.pearsonr(valid['age_at_initial_pathologic_diagnosis'], valid['mean_methylation'])

print(f"Pearson correlation coefficient: {r:.4f}")
print(f"P-value: {p_value:.6e}")

st.plotly_chart(fig)