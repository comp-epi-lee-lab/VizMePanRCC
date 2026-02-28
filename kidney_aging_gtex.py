import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

st.set_page_config(layout="wide")

@st.cache_data
def load_data():
    meth_path = Path(__file__).parent / 'data' / 'GTEx_Kidney.meth.csv'
    meth_df = pd.read_csv(meth_path, sep='\t')
    meth_df = meth_df.set_index('cgID')
    
    anno_path = Path(__file__).parent / 'data' / 'GTEx_Kidney.anno.csv'
    anno_df = pd.read_csv(anno_path, sep='\t', index_col=0)
    
    return meth_df, anno_df

meth_df, anno_df = load_data()

st.write(f"Number of CpG sites: {len(meth_df)}")
st.write(f"Number of samples: {len(meth_df.columns)}")

mean_methylation = meth_df.mean(axis=0)
results_df = pd.DataFrame({
    'sample_id': mean_methylation.index,
    'mean_methylation': mean_methylation.values
})

age_dict = {}
for col in anno_df.columns:
    age_range = anno_df.loc['age', col]
    if pd.notna(age_range):
        if '-' in str(age_range):
            low, high = map(int, age_range.split('-'))
            age_dict[col] = (low + high) / 2
        else:
            age_dict[col] = None
    else:
        age_dict[col] = None

results_df['age'] = results_df['sample_id'].map(age_dict)
results_df = results_df.dropna(subset=['age', 'mean_methylation'])

st.write(f"Number of samples with valid age and methylation: {len(results_df)}")

normal_20_40 = results_df[(results_df['age'] >= 20) & (results_df['age'] < 40)]['mean_methylation']
normal_40_50 = results_df[(results_df['age'] >= 40) & (results_df['age'] < 50)]['mean_methylation']
normal_50_60 = results_df[(results_df['age'] >= 50) & (results_df['age'] < 60)]['mean_methylation']
normal_60_70 = results_df[(results_df['age'] >= 60) & (results_df['age'] < 70)]['mean_methylation']

count_20_40 = len(normal_20_40)
count_40_50 = len(normal_40_50)
count_50_60 = len(normal_50_60)
count_60_70 = len(normal_60_70)

fig = go.Figure()

fig.add_trace(go.Box(
    y=normal_20_40,
    x=[0] * count_20_40,
    boxpoints="all",
    jitter=0.2,
    marker=dict(color="darkturquoise"),
    showlegend=False
))

fig.add_trace(go.Box(
    y=normal_40_50,
    x=[1] * count_40_50,
    boxpoints="all",
    jitter=0.2,
    marker=dict(color="darkturquoise"),
    showlegend=False
))

fig.add_trace(go.Box(
    y=normal_50_60,
    x=[2] * count_50_60,
    boxpoints="all",
    jitter=0.2,
    marker=dict(color="darkturquoise"),
    showlegend=False
))

fig.add_trace(go.Box(
    y=normal_60_70,
    x=[3] * count_60_70,
    boxpoints="all",
    jitter=0.2,
    marker=dict(color="darkturquoise"),
    showlegend=False
))

fig.update_layout(
    title='Mean Methylation Ratio Across All CpGs by Age Group (GTEx Normal Kidney)',
    title_font=dict(size=30),
    yaxis_title='<b>Mean Methylation Ratio</b>',
    yaxis_title_font=dict(size=22),
    xaxis_title='<b>Age Group</b>',
    xaxis_title_font=dict(size=22),
    xaxis=dict(
        tickmode='array',
        tickvals=[0, 1, 2, 3, 4],
        ticktext=[
            f"20-39<br>(N = {count_20_40})",
            f"40-49<br>(N = {count_40_50})",
            f"50-59<br>(N = {count_50_60})",
            f"60-69<br>(N = {count_60_70})",
        ],
        tickfont=dict(size=20)
    ),
    yaxis=dict(
        tickmode='linear',
        dtick=0.05,
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

fig.update_layout(yaxis_range=[0.57, 0.63])
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_xaxes(categoryorder='category ascending')

st.write("\n**Summary Statistics by Age Group:**")
if count_20_40 > 0:
    st.write(f"20-39: N={count_20_40}, Mean={normal_20_40.mean():.4f}, Median={normal_20_40.median():.4f}")
if count_40_50 > 0:
    st.write(f"40-49: N={count_40_50}, Mean={normal_40_50.mean():.4f}, Median={normal_40_50.median():.4f}")
if count_50_60 > 0:
    st.write(f"50-59: N={count_50_60}, Mean={normal_50_60.mean():.4f}, Median={normal_50_60.median():.4f}")
if count_60_70 > 0:
    st.write(f"60-69: N={count_60_70}, Mean={normal_60_70.mean():.4f}, Median={normal_60_70.median():.4f}")

st.write("\n**Statistical Tests (Mann-Whitney U test):**")

if count_20_40 > 0 and count_40_50 > 0:
    _, p_val = stats.mannwhitneyu(normal_20_40, normal_40_50)
    st.write(f"20-39 vs 40-49: p-value = {p_val:.5e}" + ('*' if p_val < 0.05 else ''))

if count_40_50 > 0 and count_50_60 > 0:
    _, p_val = stats.mannwhitneyu(normal_40_50, normal_50_60)
    st.write(f"40-49 vs 50-59: p-value = {p_val:.5e}" + ('*' if p_val < 0.05 else ''))

if count_50_60 > 0 and count_60_70 > 0:
    _, p_val = stats.mannwhitneyu(normal_50_60, normal_60_70)
    st.write(f"50-59 vs 60-69: p-value = {p_val:.5e}" + ('*' if p_val < 0.05 else ''))

st.write("\n**Correlation with Age:**")
results_corr = results_df[results_df['age'] < 70]
correlation, p_corr = stats.spearmanr(results_corr['age'], results_corr['mean_methylation'])
st.write(f"Spearman correlation: r = {correlation:.4f}, p-value = {p_corr:.5e}" + ('*' if p_corr < 0.05 else ''))

st.plotly_chart(fig)