import numpy as np
import scipy.stats as stats
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import warnings
import plotly.express as px
warnings.filterwarnings("ignore")

req_cols = ["vital_status", "age_at_initial_pathologic_diagnosis", "days_to_death", "gender", "rcc", "stage", "cg05423956", "cg03544320", "cg01541443", "cg00000029", "cg00001583", "cg00004773", "cg00015699", "cg00020794"]
cg_values = ["cg05423956", "cg03544320", "cg01541443", "cg00000029", "cg00001583", "cg00004773", "cg00015699", "cg00020794"]
table = pd.read_csv(
    'C:/Users/abhin/Downloads/kipan_clin_meth_20221210.tsv (1)/kipan_clin_meth_20221210.tsv', sep='\t', usecols=req_cols)
col1, col2 = st.columns([1, 3])
search_pressed = False

with col1:
    cg_value = st.text_input('Search CG Value', placeholder='cgXXXXXXXX')
    stage = st.checkbox("Stage")
    gender = st.checkbox("Gender")
    age = st.checkbox("Age")
    lts = st.checkbox("Long Term Survivorship")
    search = st.button('Search')
    if search:
        search_pressed = True

with col2:
    if search_pressed == True:
        if cg_value not in cg_values: 
            st.warning('No data found for the provided CG value.')
        else:

            if stage:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
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
                    go.Box(y=df[df['rcc'] == 'normal'][cg_value], name="normal", x=df[df['rcc'] == 'normal']['rcc'], boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False)
                )
                fig.update_layout(
                    title='stage plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3, 4],
                        ticktext=[
                            "normal (N = " + str(normalcount) + ")",
                            "stage i (N = " + str(s1count) + ")",
                            "stage ii (N = " + str(s2count) + ")",
                            "stage iii (N = " + str(s3count) + ")",
                            "stage iv (N = " + str(s4count) + ")"
                        ]
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1,
                        showgrid = False,
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20),
                )
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)

            if gender:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                male_rcc_count = ((df['gender'] == 'male') & (df['rcc'] == 'rcc')).sum()
                male_normal_count = ((df['gender'] == 'male') & (df['rcc'] == 'normal')).sum()
                female_rcc_count = ((df['gender'] == 'female') & (df['rcc'] == 'rcc')).sum()
                female_normal_count = ((df['gender'] == 'female') & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
                normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
                rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
                normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]
                fig.add_trace(go.Box(y=rcc_male[cg_value], x=[0] * len(rcc_male),
                     boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=normal_male[cg_value], x=[1] * len(normal_male),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=rcc_female[cg_value], x=[2] * len(rcc_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=normal_female[cg_value], x=[3] * len(normal_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='gender plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "rcc male (N = " + str(male_rcc_count) + ")",
                            "normal male (N = " + str(male_normal_count) + ")",
                            "rcc female (N = " + str(female_rcc_count) + ")",
                            "normal female (N = " + str(female_normal_count) + ")",
                        ],
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1,
                        showgrid = False
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ), 
                    margin=dict(r=20, b=120),       
                )
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                difference_rcc_male_vs_female = table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value].mean() - \
                                table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_rcc_male_vs_female = stats.mannwhitneyu(
                    table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value],
                    table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value]
                )
                difference_normal_male_vs_female = table.loc[(table["gender"] == "male") & (table["rcc"] == "normal"), cg_value].mean() - \
                                table.loc[(table["gender"] == "female") & (table["rcc"] == "normal"), cg_value].mean()
                _, p_value_normal_male_vs_female = stats.mannwhitneyu(
                    table.loc[(table["gender"] == "male") & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["gender"] == "female") & (table["rcc"] == "normal"), cg_value]
                )
                difference_male_vs_normal = table.loc[(table["gender"] == "male") & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_male_vs_normal = stats.mannwhitneyu(
                    table.loc[(table["gender"] == "male") & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value]
                )

                difference_female_vs_normal = table.loc[(table["gender"] == "female") & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_female_vs_normal = stats.mannwhitneyu(
                    table.loc[(table["gender"] == "female") & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value]
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'rcc mean difference: {difference_rcc_male_vs_female:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'rcc p-value: {p_value_rcc_male_vs_female:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.25,
                    xref='paper',
                    yref='paper',
                    text=f'normal mean difference: {difference_normal_male_vs_female:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.25,
                    xref='paper',
                    yref='paper',
                    text=f'normal p-value: {p_value_normal_male_vs_female:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'male mean difference: {difference_male_vs_normal:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'male p-value: {p_value_male_vs_normal:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'female mean difference: {difference_female_vs_normal:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'female p-value: {p_value_female_vs_normal:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)

            if age:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                under55_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] <= 55) & (df['rcc'] == 'rcc')).sum()
                under55_normal_count = ((df['age_at_initial_pathologic_diagnosis'] <= 55) & (df['rcc'] == 'normal')).sum()
                over55_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'rcc')).sum()
                over55_normal_count = ((df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 55) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[0]*under55_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 55) & (df['rcc'] == 'normal')][cg_value],
                                    x=[1]*under55_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[2]*over55_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'normal')][cg_value],
                                    x=[3]*over55_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "<=55 rcc (N = " + str(under55_rcc_count) + ")",
                            "<=55 normal (N = " + str(under55_normal_count) + ")",
                            ">55 rcc (N = " + str(over55_rcc_count) + ")",
                            ">55 normal (N = " + str(over55_normal_count) + ")",
                        ]
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1,
                        showgrid = False
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20, b=120),
                )
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                difference_rcc_age = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "rcc"), cg_value].mean() - \
                                table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_rcc_age = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "rcc"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "rcc"), cg_value]
                )
                difference_normal_age = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "normal"), cg_value].mean() - \
                                table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "normal"), cg_value].mean()
                _, p_value_normal_age = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "normal"), cg_value]
                )
                difference_under55 = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_under55 = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "rcc"), cg_value]
                )

                difference_over55 = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_over55 = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "rcc"), cg_value])
                fig.add_annotation(
                    x=0.28,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'rcc mean difference: {difference_rcc_age:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'rcc p-value: {p_value_rcc_age:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.25,
                    xref='paper',
                    yref='paper',
                    text=f'normal mean difference: {difference_normal_age:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.25,
                    xref='paper',
                    yref='paper',
                    text=f'normal p-value: {p_value_normal_age:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'under 55 mean difference: {difference_under55:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'under 55 p-value: {p_value_under55:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'over 55 mean difference: {difference_over55:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'over 55 p-value: {p_value_over55:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)

            if lts:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                undercount = (df['days_to_death'] <= 1825 ).sum()
                overcount = (df['days_to_death'] > 1825).sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[df['days_to_death'] <= 1825][cg_value], x=[0] * undercount, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False)
                )
                fig.add_trace(
                    go.Box(y=df[df['days_to_death'] > 1825][cg_value], x=[1] * overcount, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False)
                )
                fig.update_layout(
                    title='long term survivorship plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1],
                        ticktext=[
                            "under 5 years (N = " + str(undercount) + ")",
                            "over 5 years (N = " + str(overcount) + ")",
                        ]
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1,
                        showgrid = False
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20),
                )
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                difference = table.loc[table["days_to_death"] <= 1825, cg_value].mean() - table.loc[table["days_to_death"] > 1825, cg_value].mean()
                _, p_value = stats.mannwhitneyu(table.loc[table["days_to_death"] <= 1825, cg_value], table.loc[table["days_to_death"] > 1825, cg_value])
                fig.add_annotation(
                    x=0.28,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'Mean difference: {difference:.4f}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.2,
                    xref='paper',
                    yref='paper',
                    text=f'p-value: {p_value:.4e}',
                    showarrow=False,
                    font=dict(size=12),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)