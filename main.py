import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

csv_path = Path(r'data/kipan_clin_meth.tsv')
pickle_path = Path(r'data/pickle_file.pk1')

st.set_page_config(page_title=None, page_icon=None, layout="wide", initial_sidebar_state="auto", menu_items=None)

if(pickle_path.is_file() == False or csv_path.stat().st_mtime > pickle_path.stat().st_mtime):

    table = pd.read_csv(csv_path, sep='\t')
    table.to_pickle(pickle_path)

    large_20_differences = [-float('inf')] * 20
    cg_large_20_differences = [str('empty')] * 20

    i = 0
    for (columnName, columnData) in table.items():
        if (i >= 17):
            difference = abs(table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), columnName].mean() - table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), columnName].mean())
            for j in range(0, 20):
                if (difference > large_20_differences[j]):
                    large_20_differences.insert(j, difference)
                    large_20_differences.pop()
                    cg_large_20_differences.insert(j, columnName)
                    cg_large_20_differences.pop()
                    break
        i = i + 1

table = pd.read_pickle(pickle_path)
df = table.dropna()

col1, col2 = st.columns([1, 3])
search_pressed = False

with col1:
    cg_value = st.text_input('Search CG Value', value='cg01774894', placeholder='cgXXXXXXXX')
    age = st.checkbox("Age", value=True, key='age')
    gender = st.checkbox("Gender", value=True, key='gender')
    lts = st.checkbox("Long Term Survivorship", value=True, key='lts')
    stage = st.checkbox("Stage", value=True, key='stage')
    search = st.button('Search')
    if search:
        search_pressed = True

with col2:
    if search_pressed == True:
        if cg_value not in table.columns:
            st.warning('No data found for the provided CG value.')
        else:
            df = table.sort_values(cg_value, ascending=False)
            df = table.dropna(subset=[cg_value])

            if age:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                under50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'rcc')).sum()
                under50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'normal')).sum()
                over50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'rcc')).sum()
                over50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'normal')][cg_value],
                                    x=[0]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'normal')][cg_value],
                                    x=[1]*over50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="#3fa6a8"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[2]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[3]*over50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="#8977ad"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value) + ' (young [<= 50] vs old [>50])',
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>age and condition</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "young normal<br>(N = " + str(under50_normal_count) + ")",
                            "old normal<br>(N = " + str(over50_normal_count) + ")",
                            "early-onset RCC<br>(N = " + str(under50_rcc_count) + ")",
                            "late-onset RCC<br>(N = " + str(over50_rcc_count) + ")",
                        ],
                        tickfont=dict(size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20, b=320),
                    height=750,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.1])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                old_rcc = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                young_rcc = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                old_normal = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value].dropna()
                young_normal = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value].dropna()

                difference_rcc_age = old_rcc.mean() - \
                                young_rcc.mean()
                _, p_value_rcc_age = stats.mannwhitneyu(
                    young_rcc, old_rcc
                )
                difference_normal_age = old_normal.mean() - \
                                young_normal.mean()
                _, p_value_normal_age = stats.mannwhitneyu(
                    young_normal, old_normal
                )
                difference_under50 = young_normal.mean() - \
                                            young_rcc.mean()
                _, p_value_under50 = stats.mannwhitneyu(
                    young_normal, young_rcc
                )

                difference_over50 = old_normal.mean() - \
                                            old_rcc.mean()
                _, p_value_over50 = stats.mannwhitneyu(
                    old_normal, old_rcc
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'mean differences & p-values:',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'RCC old - RCC young: {difference_rcc_age:.5f} | p-value: {p_value_rcc_age:.5e}' + ('*' if p_value_rcc_age < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'normal old - normal young: {difference_normal_age:.5f} | p-value: {p_value_normal_age:.5e}' + ('*' if p_value_normal_age < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )

                significance1 = 'ns'
                if p_value_normal_age < 0.0005:
                    significance1 = '***'
                elif p_value_normal_age < 0.005:
                    significance1 = '**'
                elif p_value_normal_age < 0.05:
                    significance1 = '*'

                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=1,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=0,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=1,
                    y0=1.03,
                    x1=1,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=0.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance1
                )

                significance2 = 'ns'
                if p_value_rcc_age < 0.0005:
                    significance2 = '***'
                elif p_value_rcc_age < 0.005:
                    significance2 = '**'
                elif p_value_rcc_age < 0.05:
                    significance2 = '*'

                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=3,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=2,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=3,
                    y0=1.03,
                    x1=3,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=2.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance2
                )

                st.plotly_chart(fig)

            if age:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                under50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'rcc')).sum()
                under50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'normal')).sum()
                over50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'rcc')).sum()
                over50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'normal')][cg_value],
                                    x=[0]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'normal')][cg_value],
                                    x=[2]*over50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="#3fa6a8"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 50) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[1]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 50) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[3]*over50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="#8977ad"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value) + ' (young [<= 50] vs old [>50])',
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>age and condition</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 2, 1, 3],
                        ticktext=[
                            "young normal<br>(N = " + str(under50_normal_count) + ")",
                            "old normal<br>(N = " + str(over50_normal_count) + ")",
                            "early-onset RCC<br>(N = " + str(under50_rcc_count) + ")",
                            "late-onset RCC<br>(N = " + str(over50_rcc_count) + ")",
                        ],
                        tickfont=dict(size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20, b=320),
                    height=750,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.1])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                old_rcc = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                young_rcc = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                old_normal = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value].dropna()
                young_normal = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value].dropna()

                difference_rcc_age = old_rcc.mean() - \
                                young_rcc.mean()
                _, p_value_rcc_age = stats.mannwhitneyu(
                    young_rcc, old_rcc
                )
                difference_normal_age = old_normal.mean() - \
                                young_normal.mean()
                _, p_value_normal_age = stats.mannwhitneyu(
                    young_normal, old_normal
                )
                difference_under50 = young_normal.mean() - \
                                            young_rcc.mean()
                _, p_value_under50 = stats.mannwhitneyu(
                    young_normal, young_rcc
                )

                difference_over50 = old_normal.mean() - \
                                            old_rcc.mean()
                _, p_value_over50 = stats.mannwhitneyu(
                    old_normal, old_rcc
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'mean differences & p-values:',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'normal young - RCC young: {difference_under50:.5f} | p-value: {p_value_under50:.5e}' + ('*' if p_value_under50 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'normal old - RCC old: {difference_over50:.5f} | p-value: {p_value_over50:.5e}' + ('*' if p_value_over50 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )

                significance1 = 'ns'
                if p_value_under50 < 0.0005:
                    significance1 = '***'
                elif p_value_under50 < 0.005:
                    significance1 = '**'
                elif p_value_under50 < 0.05:
                    significance1 = '*'

                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=1,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=0,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=1,
                    y0=1.03,
                    x1=1,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=0.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance1
                )

                significance2 = 'ns'
                if p_value_over50 < 0.0005:
                    significance2 = '***'
                elif p_value_over50 < 0.005:
                    significance2 = '**'
                elif p_value_over50 < 0.05:
                    significance2 = '*'

                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=3,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=2,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=3,
                    y0=1.03,
                    x1=3,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=2.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance2
                )

                st.plotly_chart(fig)


            if age:
                rcc_40 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_40 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_60 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_60 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_70 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_70 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_80 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_80 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_90 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_90 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "normal"), cg_value].dropna()

                under40_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40), cg_value].count()
                under50_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50), cg_value].count()
                under60_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60), cg_value].count()
                under70_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70), cg_value].count()
                under80_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80), cg_value].count()
                under90_count = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90), cg_value].count()

                under40_rcc_mean = rcc_40.mean()
                under40_normal_mean = normal_40.mean()
                under50_rcc_mean = rcc_50.mean()
                under50_normal_mean = normal_50.mean()
                under60_rcc_mean = rcc_60.mean()
                under60_normal_mean = normal_60.mean()
                under70_rcc_mean = rcc_70.mean()
                under70_normal_mean = normal_70.mean()
                under80_rcc_mean = rcc_80.mean()
                under80_normal_mean = normal_80.mean()
                under90_rcc_mean = rcc_90.mean()
                under90_normal_mean = normal_90.mean()
                mean_values_rcc = [ under40_rcc_mean, under50_rcc_mean, under60_rcc_mean, under70_rcc_mean, under80_rcc_mean, under90_rcc_mean]
                mean_values_normal = [under40_normal_mean, under50_normal_mean, under60_normal_mean, under70_normal_mean, under80_normal_mean, under90_normal_mean]
                labels = ["20-39<br>(N = " + str(under40_count) + ")", "40-49<br>(N = " + str(under50_count) + ")", "50-59<br>(N = " + str(under60_count) + ")", "60-69<br>(N = " + str(under70_count) + ")", "70-79<br>(N = " + str(under80_count) + ")", "80+<br>(N = " + str(under90_count) + ")"]
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=labels, y=mean_values_normal, mode='lines+markers', name='Normal', line=dict(color='darkturquoise'), marker=dict(size=12)))
                fig.add_trace(go.Scatter(x=labels, y=mean_values_rcc, mode='lines+markers', name='RCC', line=dict(color='mediumpurple'), marker=dict(size=12)))
                fig.update_layout(
                    title='age line plot of ' + str(cg_value) + ' (cohorts by 10)',
                    title_font=dict(size=30),
                    xaxis_title='<b>age ranges</b>',
                    yaxis_title='<b>mean methylation ratio</b>',
                    xaxis_title_font=dict(size=22),
                    yaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickfont=dict(size=22),
                    ),
                    font=dict(
                        family="Arial",
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22),
                    ),
                    legend=dict(orientation="h", font=dict(size=22), y=-0.27),
                    margin=dict(r=80, b=150),
                    height=600,
                    width=1050
                )
                fig.update_layout(yaxis_range=[0,1.0])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                st.plotly_chart(fig)
                

            if age:
                rcc_40 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_40 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_60 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_60 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_70 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_70 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_80 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_80 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "normal"), cg_value].dropna()
                rcc_90 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "rcc"), cg_value].dropna()
                normal_90 = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "normal"), cg_value].dropna()
               
                under40_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')).sum()
                under40_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')).sum()
                under50_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) &( df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')).sum()
                under50_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')).sum()
                under60_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')).sum()
                under60_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')).sum()
                under70_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')).sum()
                under70_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')).sum()
                under80_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'rcc')).sum()
                under80_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'normal')).sum()
                under90_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 80)) & (df['rcc'] == 'rcc')).sum()
                under90_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 80)) & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')][cg_value], x=[0]*under40_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')][cg_value], x=[1]*under40_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')][cg_value], x=[2]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')][cg_value], x=[3]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')][cg_value], x=[4]*under60_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')][cg_value], x=[5]*under60_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')][cg_value], x=[6]*under70_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')][cg_value], x=[7]*under70_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'normal')][cg_value], x=[8]*under80_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'rcc')][cg_value], x=[9]*under80_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 80)) & (df['rcc'] == 'normal')][cg_value], x=[10]*under90_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 80)) & (df['rcc'] == 'rcc')][cg_value], x=[11]*under90_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value) + ' (cohorts by 10)',
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>age ranges</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
                        ticktext=[
                            "20-39 normal (N = " + str(under40_normal_count) + ")",
                            "20-39 RCC (N = " + str(under40_rcc_count) + ")",
                            "40-49 normal (N = " + str(under50_normal_count) + ")",
                            "40-49 RCC (N = " + str(under50_rcc_count) + ")",
                            "50-59 normal (N = " + str(under60_normal_count) + ")",
                            "50-59 RCC (N = " + str(under60_rcc_count) + ")",
                            "60-69 normal (N = " + str(under70_normal_count) + ")",
                            "60-69 RCC (N = " + str(under70_rcc_count) + ")",
                            "70-79 normal (N = " + str(under80_normal_count) + ")",
                            "70-79 RCC (N = " + str(under80_rcc_count) + ")",
                            "80-89 normal (N = " + str(under90_normal_count) + ")",
                            "80-89 RCC (N = " + str(under90_rcc_count) + ")",
                        ],
                        tickfont=dict(family="Arial", size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=120, b=500),
                    height = 950,
                    width = 1100
                )
                fig.update_layout(yaxis_range=[0,1.0])
                fig.update_xaxes(tickangle=270, showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                _, p_value_under40 = stats.mannwhitneyu(
                rcc_40, normal_40
                )
                _, p_value_under50 = stats.mannwhitneyu(
                rcc_50, normal_50
                )
                _, p_value_under60 = stats.mannwhitneyu(
                rcc_60, normal_60
                )
                _, p_value_under70 = stats.mannwhitneyu(
                rcc_70, normal_70
                )
                _, p_value_under80 = stats.mannwhitneyu(
                rcc_80, normal_80
                )
                _, p_value_under90 = stats.mannwhitneyu(
                rcc_90, normal_90
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.74,
                    xref='paper',
                    yref='paper',
                    text=f'p-values (RCC vs normal):',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.84,
                    xref='paper',
                    yref='paper',
                    text=f'20-39: {p_value_under40:.5e}' + ('*' if p_value_under40 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.94,
                    xref='paper',
                    yref='paper',
                    text=f'40-49: {p_value_under50:.5e}' + ('*' if p_value_under50 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.04,
                    xref='paper',
                    yref='paper',
                    text=f'50-59: {p_value_under60:.5e}' + ('*' if p_value_under60 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.14,
                    xref='paper',
                    yref='paper',
                    text=f'60-69: {p_value_under70:.5e}' + ('*' if p_value_under70 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.24,
                    xref='paper',
                    yref='paper',
                    text=f'70-79: {p_value_under80:.5e}' + ('*' if p_value_under80 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.34,
                    xref='paper',
                    yref='paper',
                    text=f'80-89: {p_value_under90:.5e}' + ('*' if p_value_under90 < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)

            if gender:
                male_rcc_count = ((df['gender'] == 'male') & (df['rcc'] == 'rcc')).sum()
                male_normal_count = ((df['gender'] == 'male') & (df['rcc'] == 'normal')).sum()
                female_rcc_count = ((df['gender'] == 'female') & (df['rcc'] == 'rcc')).sum()
                female_normal_count = ((df['gender'] == 'female') & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
                normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
                rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
                normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]
                fig.add_trace(go.Box(y=normal_male[cg_value], x=[0] * len(normal_male),
                                    boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=normal_female[cg_value], x=[1] * len(normal_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=rcc_male[cg_value], x=[2] * len(rcc_male),
                     boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=rcc_female[cg_value], x=[3] * len(rcc_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='gender plot of ' + str(cg_value),
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>gender and condition</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "normal male<br>(N = " + str(male_normal_count) + ")",
                            "normal female<br>(N = " + str(female_normal_count) + ")",
                            "RCC male<br>(N = " + str(male_rcc_count) + ")",
                            "RCC female<br>(N = " + str(female_rcc_count) + ")",
                        ],
                        tickfont=dict(size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ), 
                    margin=dict(r=20, b=330, l=20),       
                    height=750,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.1])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

                rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
                normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
                rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
                normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]

                rcc_male_values = rcc_male[cg_value].dropna()
                rcc_female_values = rcc_female[cg_value].dropna()
                normal_male_values = normal_male[cg_value].dropna()
                normal_female_values = normal_female[cg_value].dropna()

                difference_rcc_male_vs_female = rcc_female_values.mean() - \
                                rcc_male_values.mean()
                _, p_value_rcc_male_vs_female = stats.mannwhitneyu(
                    rcc_male_values, rcc_female_values
                )
                difference_normal_male_vs_female = normal_female_values.mean() - \
                                normal_male_values.mean()
                _, p_value_normal_male_vs_female = stats.mannwhitneyu(
                    normal_male_values, normal_female_values
                )
                difference_male_vs_normal = normal_male_values.mean() - \
                                            rcc_male_values.mean()
                _, p_value_male_vs_normal = stats.mannwhitneyu(
                    normal_male_values,
                    rcc_male_values
                )

                difference_female_vs_normal = normal_female_values.mean() - \
                                            rcc_female_values.mean()
                _, p_value_female_vs_normal = stats.mannwhitneyu(
                    normal_female_values, rcc_female_values
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'mean differences & p-values:',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'RCC female - RCC male: {difference_rcc_male_vs_female:.5f} | p-value: {p_value_rcc_male_vs_female:.5e}' + ('*' if p_value_rcc_male_vs_female < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'normal female - normal male: {difference_normal_male_vs_female:.5f} | p-value: {p_value_normal_male_vs_female:.5e}' + ('*' if p_value_normal_male_vs_female < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )

                significance1 = 'ns'
                if p_value_normal_male_vs_female < 0.0005:
                    significance1 = '***'
                elif p_value_normal_male_vs_female < 0.005:
                    significance1 = '**'
                elif p_value_normal_male_vs_female < 0.05:
                    significance1 = '*'

                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=1,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=0,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=1,
                    y0=1.03,
                    x1=1,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=0.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance1
                )

                significance2 = 'ns'
                if p_value_rcc_male_vs_female < 0.0005:
                    significance2 = '***'
                elif p_value_rcc_male_vs_female < 0.005:
                    significance2 = '**'
                elif p_value_rcc_male_vs_female < 0.05:
                    significance2 = '*'

                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=3,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=2,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=3,
                    y0=1.03,
                    x1=3,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=2.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance2
                )

                st.plotly_chart(fig)

            if gender:
                male_rcc_count = ((df['gender'] == 'male') & (df['rcc'] == 'rcc')).sum()
                male_normal_count = ((df['gender'] == 'male') & (df['rcc'] == 'normal')).sum()
                female_rcc_count = ((df['gender'] == 'female') & (df['rcc'] == 'rcc')).sum()
                female_normal_count = ((df['gender'] == 'female') & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
                normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
                rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
                normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]
                fig.add_trace(go.Box(y=normal_male[cg_value], x=[0] * len(normal_male),
                                    boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=normal_female[cg_value], x=[2] * len(normal_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
                fig.add_trace(go.Box(y=rcc_male[cg_value], x=[1] * len(rcc_male),
                     boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=rcc_female[cg_value], x=[3] * len(rcc_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='gender plot of ' + str(cg_value),
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>gender and condition</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 2, 1, 3],
                        ticktext=[
                            "normal male<br>(N = " + str(male_normal_count) + ")",
                            "normal female<br>(N = " + str(female_normal_count) + ")",
                            "RCC male<br>(N = " + str(male_rcc_count) + ")",
                            "RCC female<br>(N = " + str(female_rcc_count) + ")",
                        ],
                        tickfont=dict(size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ), 
                    margin=dict(r=20, b=330, l=20),       
                    height=750,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.1])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

                rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
                normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
                rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
                normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]

                rcc_male_values = rcc_male[cg_value].dropna()
                rcc_female_values = rcc_female[cg_value].dropna()
                normal_male_values = normal_male[cg_value].dropna()
                normal_female_values = normal_female[cg_value].dropna()

                difference_rcc_male_vs_female = rcc_female_values.mean() - \
                                rcc_male_values.mean()
                _, p_value_rcc_male_vs_female = stats.mannwhitneyu(
                    rcc_male_values, rcc_female_values
                )
                difference_normal_male_vs_female = normal_female_values.mean() - \
                                normal_male_values.mean()
                _, p_value_normal_male_vs_female = stats.mannwhitneyu(
                    normal_male_values, normal_female_values
                )
                difference_male_vs_normal = normal_male_values.mean() - \
                                            rcc_male_values.mean()
                _, p_value_male_vs_normal = stats.mannwhitneyu(
                    normal_male_values,
                    rcc_male_values
                )

                difference_female_vs_normal = normal_female_values.mean() - \
                                            rcc_female_values.mean()
                _, p_value_female_vs_normal = stats.mannwhitneyu(
                    normal_female_values, rcc_female_values
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'mean differences & p-values:',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'normal male - RCC male: {difference_male_vs_normal:.5f} | p-value: {p_value_male_vs_normal:.5e}' + ('*' if p_value_male_vs_normal < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'normal female - RCC female: {difference_female_vs_normal:.5f} | p-value: {p_value_female_vs_normal:.5e}' + ('*' if p_value_female_vs_normal < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )

                significance1 = 'ns'
                if p_value_male_vs_normal < 0.0005:
                    significance1 = '***'
                elif p_value_male_vs_normal < 0.005:
                    significance1 = '**'
                elif p_value_male_vs_normal < 0.05:
                    significance1 = '*'

                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=1,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=0,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=1,
                    y0=1.03,
                    x1=1,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=0.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance1
                )

                significance2 = 'ns'
                if p_value_female_vs_normal < 0.0005:
                    significance2 = '***'
                elif p_value_female_vs_normal < 0.005:
                    significance2 = '**'
                elif p_value_female_vs_normal < 0.05:
                    significance2 = '*'

                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=3,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=2,
                    y0=1.03,
                    x1=2,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=3,
                    y0=1.03,
                    x1=3,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=2.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance2
                )

                st.plotly_chart(fig)

            if lts:
                undercount = (df['days_to_death'] <= 1825 ).sum()
                overcount = (df['days_to_death'] > 1825).sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[df['days_to_death'] <= 1825][cg_value], x=[0] * undercount, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False)
                )
                fig.add_trace(
                    go.Box(y=df[df['days_to_death'] > 1825][cg_value], x=[1] * overcount, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False)
                )
                
                fig.update_layout(
                    title='long term survivorship plot of ' + str(cg_value),
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>years until death</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1],
                        ticktext=[
                            "under 5 years<br>(N = " + str(undercount) + ")",
                            "over 5 years<br>(N = " + str(overcount) + ")",
                        ],
                        tickfont=dict(size=20)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20, b=170),
                    height=600,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.1])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )

                less_days_to_death = table.loc[table["days_to_death"] <= 1825, cg_value].dropna()
                more_days_to_death = table.loc[table["days_to_death"] > 1825, cg_value].dropna()
                difference = less_days_to_death.mean() - more_days_to_death.mean()
                _, p_value = stats.mannwhitneyu(less_days_to_death, more_days_to_death)

                fig.add_annotation(
                    x=0.5,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'Mean difference (under 5 years - over 5 years): {difference:.5f} | p-value: {p_value:.5e}' + ('*' if p_value < 0.05 else ''),
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )

                significance = 'ns'
                if p_value < 0.0005:
                    significance = '***'
                elif p_value < 0.005:
                    significance = '**'
                elif p_value < 0.05:
                    significance = '*'

                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=1,
                    y1=1.03,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=0,
                    y0=1.03,
                    x1=0,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_shape(
                    type="line",
                    x0=1,
                    y0=1.03,
                    x1=1,
                    y1=1.01,
                    line=dict(color="black", width=2)
                )
                fig.add_annotation(
                    x=0.5,
                    y=1.08,
                    font=dict(size=24),
                    showarrow=False,
                    text=significance
                )

                st.plotly_chart(fig)

            if stage:
                s1count = (df['stage'] == 'stage i').sum()
                s2count = (df['stage'] == 'stage ii').sum()
                s3count = (df['stage'] == 'stage iii').sum()
                s4count = (df['stage'] == 'stage iv').sum()
                normalcount = (df['rcc'] == 'normal').sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[cg_value], x=df['stage'], name="stage", boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False)
                )

                fig.add_trace(
                    go.Box(y=df[df['rcc'] == 'normal'][cg_value], name="normal", x=df[df['rcc'] == 'normal']['rcc'], boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False)
                )
                fig.update_layout(
                    title='stage plot of ' + str(cg_value),
                    title_font=dict(size=30),
                    yaxis_title='<b>methylation ratio</b>',
                    yaxis_title_font=dict(size=22),
                    xaxis_title='<b>stage</b>',
                    xaxis_title_font=dict(size=22),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3, 4],
                        ticktext=[
                            "normal<br>(N = " + str(normalcount) + ")",
                            "stage i<br>(N = " + str(s1count) + ")",
                            "stage ii<br>(N = " + str(s2count) + ")",
                            "stage iii<br>(N = " + str(s3count) + ")",
                            "stage iv<br>(N = " + str(s4count) + ")"
                        ],
                        tickfont=dict(size=19)
                    ),
                    font=dict(
                        family="Arial",
                        size=24
                    ), 
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.2,
                        showgrid = False,
                        tickfont=dict(size=22)
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    margin=dict(r=20),
                    height=600,
                    width=1000
                )
                fig.update_layout(yaxis_range=[0,1.0])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)

            

            