import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

csv_path = Path(r'C:\Users\abhin\Downloads\kipan_clin_meth_20221210.tsv (1)\kipan_clin_meth_20221210.tsv')
pickle_path = Path(r'C:\Users\abhin\OneDrive\Documents\GitHub\Research\processed_data.pk1')

large_20_differences = [0.15642470275754355, 0.1553588536675063, 0.15482110244865024, 0.1544712550474568, 0.1533709537114315, 0.15219904231398718, 0.15074757851385756, 0.14662371028511223, 0.14650689524987492, 0.14342406505768346, 0.14025356892466606, 0.13877998112299295, 0.13861957326956653, 0.13688944108480505, 0.13586140448081035, 0.13453964654821216, 0.13434813435909798, 0.13305454319644583, 0.13005207422393905, 0.12843642705731367]
cg_large_20_differences = ["cg02598441", "cg06296331", "cg04039555", "cg01774894", "cg00249511", "cg03045635", "cg05982757", "cg08940787", "cg00939495", "cg04597433", "cg02632185", "cg05721365", "cg06113789", "cg08827307", "cg00576279", "cg02637318", "cg04387835", "cg08360726", "cg01382860", "cg06546677"]

if(pickle_path.is_file() == False or csv_path.stat().st_mtime > pickle_path.stat().st_mtime):

    table = pd.read_csv(csv_path, sep='\t')
    table.to_pickle(pickle_path)

    large_20_differences = [-float('inf')] * 20
    cg_large_20_differences = [str('empty')] * 20

    i = 0
    for (columnName, columnData) in table.iteritems():
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
        if cg_value not in table.columns:
            st.warning('No data found for the provided CG value.')
        else:
            df = table.sort_values(cg_value, ascending=False)
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
                        dtick=0.1,
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
                        dtick=0.1,
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
                fig.update_layout(yaxis_range=[0,1.0])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                difference_rcc_male_vs_female = table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value].mean() - \
                                table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_rcc_male_vs_female = stats.mannwhitneyu(
                    table.loc[(table["gender"] == "male") & (table["rcc"] == "rcc"), cg_value],
                    table.loc[(table["gender"] == "female") & (table["rcc"] == "rcc"), cg_value]
                )
                difference_normal_male_vs_female = table.loc[(table["gender"] == "female") & (table["rcc"] == "normal"), cg_value].mean() - \
                                table.loc[(table["gender"] == "male") & (table["rcc"] == "normal"), cg_value].mean()
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
                    text=f'RCC female - RCC male: {difference_rcc_male_vs_female:.5f} | p-value: {p_value_rcc_male_vs_female:.5e}',
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
                    text=f'normal female - normal male: {difference_normal_male_vs_female:.5f} | p-value: {p_value_normal_male_vs_female:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.7,
                    xref='paper',
                    yref='paper',
                    text=f'normal male - RCC male: {difference_male_vs_normal:.5f} | p-value: {p_value_male_vs_normal:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.8,
                    xref='paper',
                    yref='paper',
                    text=f'normal female - RCC female: {difference_female_vs_normal:.5f} | p-value: {p_value_female_vs_normal:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
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
                            dtick=0.1,
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
                    fig.update_layout(yaxis_range=[0,1.0])
                    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                    fig.update_xaxes(
                        categoryorder='category ascending'
                    )
                    difference = table.loc[table["days_to_death"] <= 1825, cg_value].mean() - table.loc[table["days_to_death"] > 1825, cg_value].mean()
                    _, p_value = stats.mannwhitneyu(table.loc[table["days_to_death"] <= 1825, cg_value], table.loc[table["days_to_death"] > 1825, cg_value])
                    fig.add_annotation(
                        x=0.5,
                        y=-0.4,
                        xref='paper',
                        yref='paper',
                        text=f'Mean difference (under 5 years - over 5 years): {difference:.5f} | p-value: {p_value:.5e}',
                        showarrow=False,
                        font=dict(size=22),
                        align='center',
                        xanchor='center',
                        yanchor='top'
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
                        dtick=0.1,
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
                fig.update_layout(yaxis_range=[0,1.0])
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                difference_rcc_age = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value].mean() - \
                                table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_rcc_age = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value]
                )
                difference_normal_age = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value].mean() - \
                                table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value].mean()
                _, p_value_normal_age = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value]
                )
                difference_under50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_under50 = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value]
                )

                difference_over50 = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value].mean() - \
                                            table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value].mean()
                _, p_value_over50 = stats.mannwhitneyu(
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value],
                    table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value])
                rcc_young_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "rcc"), cg_value].mean()
                rcc_old_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "rcc"), cg_value].mean()
                normal_young_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 50) & (table["rcc"] == "normal"), cg_value].mean()
                normal_old_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] > 50) & (table["rcc"] == "normal"), cg_value].mean()
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
                    text=f'RCC old - RCC young: {difference_rcc_age:.5f} | p-value: {p_value_rcc_age:.5e}',
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
                    text=f'normal old - normal young: {difference_normal_age:.5f} | p-value: {p_value_normal_age:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.7,
                    xref='paper',
                    yref='paper',
                    text=f'normal young - RCC young: {difference_under50:.5f} | p-value: {p_value_under50:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.8,
                    xref='paper',
                    yref='paper',
                    text=f'normal old - RCC old: {difference_over50:.5f} | p-value: {p_value_over50:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)

            if age:
                under40_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40))).sum()
                under50_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) &( df['age_at_initial_pathologic_diagnosis'] < 50))).sum()
                under60_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60))).sum()
                under70_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70))).sum()
                under80_count = (((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80))).sum()
                under90_count = (((df['age_at_initial_pathologic_diagnosis'] >= 80))).sum()

                under40_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "rcc"), cg_value].mean()
                under40_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40) & (table["rcc"] == "normal"), cg_value].mean()
                under50_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "rcc"), cg_value].mean()
                under50_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50) & (table["rcc"] == "normal"), cg_value].mean()
                under60_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "rcc"), cg_value].mean()
                under60_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60) & (table["rcc"] == "normal"), cg_value].mean()
                under70_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "rcc"), cg_value].mean()
                under70_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70) & (table["rcc"] == "normal"), cg_value].mean()
                under80_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "rcc"), cg_value].mean()
                under80_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80) & (table["rcc"] == "normal"), cg_value].mean()
                under90_rcc_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "rcc"), cg_value].mean()
                under90_normal_mean = table.loc[(table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90) & (table["rcc"] == "normal"), cg_value].mean()
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
                        dtick=0.1,
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
                            "80+ normal (N = " + str(under90_normal_count) + ")",
                            "80+ RCC (N = " + str(under90_rcc_count) + ")",
                        ],
                        tickfont=dict(size=16)
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1,
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
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                _, p_value_under40 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 20) & (table["age_at_initial_pathologic_diagnosis"] < 40)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under50 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 40) & (table["age_at_initial_pathologic_diagnosis"] < 50)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under60 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 50) & (table["age_at_initial_pathologic_diagnosis"] < 60)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under70 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 60) & (table["age_at_initial_pathologic_diagnosis"] < 70)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under80 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 70) & (table["age_at_initial_pathologic_diagnosis"] < 80)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under90 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 80)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 80)) & (table["rcc"] == "normal"), cg_value]
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.5,
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
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'20-39: {p_value_under40:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.7,
                    xref='paper',
                    yref='paper',
                    text=f'40-49: {p_value_under50:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.8,
                    xref='paper',
                    yref='paper',
                    text=f'50-59: {p_value_under60:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-0.9,
                    xref='paper',
                    yref='paper',
                    text=f'60-69: {p_value_under70:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.00,
                    xref='paper',
                    yref='paper',
                    text=f'70-79: {p_value_under80:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.5,
                    y=-1.1,
                    xref='paper',
                    yref='paper',
                    text=f'80+: {p_value_under90:.5e}',
                    showarrow=False,
                    font=dict(size=22),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)