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
            difference = abs(table.loc[(table["age_at_initial_pathologic_diagnosis"] <= 55) & (table["rcc"] == "rcc"), columnName].mean() - table.loc[(table["age_at_initial_pathologic_diagnosis"] > 55) & (table["rcc"] == "rcc"), columnName].mean())
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
                fig.add_trace(go.Box(y=rcc_female[cg_value], x=[1] * len(rcc_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=normal_male[cg_value], x=[2] * len(normal_male),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=normal_female[cg_value], x=[3] * len(normal_female),
                                    boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='gender plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "RCC male (N = " + str(male_rcc_count) + ")",
                            "RCC female (N = " + str(female_rcc_count) + ")",
                            "normal male (N = " + str(male_normal_count) + ")",
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
                    margin=dict(r=20, b=145, l=20),       
                )
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
                    x=0.201,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'mean difference of RCC female and RCC male: {difference_rcc_male_vs_female:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.748,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'p-value of RCC female and RCC male: {p_value_rcc_male_vs_female:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.218,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'mean difference of normal female and normal male: {difference_normal_male_vs_female:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7629,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'p-value of normal female and normal male: {p_value_normal_male_vs_female:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.205,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'mean difference of RCC male and normal male: {difference_male_vs_normal:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.748,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'p-value of RCC male and normal male: {p_value_male_vs_normal:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.22,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'mean difference of RCC female and normal female: {difference_female_vs_normal:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7617,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'p-value of RCC female and normal female: {p_value_female_vs_normal:.4e}',
                    showarrow=False,
                    font=dict(size=13),
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
                        font=dict(size=13),
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
                        font=dict(size=13),
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
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'rcc')][cg_value],
                                    x=[1]*over55_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] <= 55) & (df['rcc'] == 'normal')][cg_value],
                                    x=[2]*under55_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] > 55) & (df['rcc'] == 'normal')][cg_value],
                                    x=[3]*over55_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value) + ' (young vs old)',
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1, 2, 3],
                        ticktext=[
                            "<=55 rcc (N = " + str(under55_rcc_count) + ")",
                            ">55 rcc (N = " + str(over55_rcc_count) + ")",
                            "<=55 normal (N = " + str(under55_normal_count) + ")",
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
                    margin=dict(r=20, b=150),
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
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'rcc mean difference: {difference_rcc_age:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.3,
                    xref='paper',
                    yref='paper',
                    text=f'rcc p-value: {p_value_rcc_age:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'normal mean difference: {difference_normal_age:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'normal p-value: {p_value_normal_age:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'under 55 mean difference: {difference_under55:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'under 55 p-value: {p_value_under55:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.28,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'over 55 mean difference: {difference_over55:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.7,
                    y=-0.6,
                    xref='paper',
                    yref='paper',
                    text=f'over 55 p-value: {p_value_over55:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)

            if age:
                data = go.Bar(x = cg_large_20_differences, y = large_20_differences)
                layout = go.Layout(
                    yaxis=dict(
                        range=[0, 1],
                        dtick=0.05,
                    ),
                    height = 800,
                    width = 950
                )
                fig = go.Figure(data = [data], layout=layout)
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_layout(
                    title='20 largest mean differences in methylation values of rcc patients (<=55 vs >55)'
                )
                st.plotly_chart(fig)

            if age:
                under20_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 0) & (df['age_at_initial_pathologic_diagnosis'] < 20)) & (df['rcc'] == 'rcc')).sum()
                under20_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 0) & (df['age_at_initial_pathologic_diagnosis'] < 20)) & (df['rcc'] == 'normal')).sum()
                under30_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 30)) & (df['rcc'] == 'rcc')).sum()
                under30_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 30)) & (df['rcc'] == 'normal')).sum()
                under40_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 30) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')).sum()
                under40_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 30) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')).sum()
                under50_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) &( df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')).sum()
                under50_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')).sum()
                under60_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')).sum()
                under60_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')).sum()
                under70_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')).sum()
                under70_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')).sum()
                under80_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'rcc')).sum()
                under80_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'normal')).sum()
                under90_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 80) & (df['age_at_initial_pathologic_diagnosis'] < 90)) & (df['rcc'] == 'rcc')).sum()
                under90_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 80) & (df['age_at_initial_pathologic_diagnosis'] < 90)) & (df['rcc'] == 'normal')).sum()
                under100_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 90) & (df['age_at_initial_pathologic_diagnosis'] < 100)) & (df['rcc'] == 'rcc')).sum()
                under100_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 90) & (df['age_at_initial_pathologic_diagnosis'] < 100)) & (df['rcc'] == 'normal')).sum()
                fig = go.Figure()
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 0) & (df['age_at_initial_pathologic_diagnosis'] < 20)) & (df['rcc'] == 'rcc')][cg_value], x=[0]*under20_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 0) & (df['age_at_initial_pathologic_diagnosis'] < 20)) & (df['rcc'] == 'normal')][cg_value], x=[1]*under20_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 30)) & (df['rcc'] == 'rcc')][cg_value], x=[2]*under30_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 30)) & (df['rcc'] == 'normal')][cg_value], x=[3]*under30_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 30) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')][cg_value], x=[4]*under40_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 30) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')][cg_value], x=[5]*under40_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')][cg_value], x=[6]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')][cg_value], x=[7]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')][cg_value], x=[8]*under60_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')][cg_value], x=[9]*under60_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')][cg_value], x=[10]*under70_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')][cg_value], x=[11]*under70_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'rcc')][cg_value], x=[12]*under80_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['age_at_initial_pathologic_diagnosis'] < 80)) & (df['rcc'] == 'normal')][cg_value], x=[13]*under80_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 80) & (df['age_at_initial_pathologic_diagnosis'] < 90)) & (df['rcc'] == 'rcc')][cg_value], x=[14]*under90_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 80) & (df['age_at_initial_pathologic_diagnosis'] < 90)) & (df['rcc'] == 'normal')][cg_value], x=[15]*under90_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 90) & (df['age_at_initial_pathologic_diagnosis'] < 100)) & (df['rcc'] == 'rcc')][cg_value], x=[16]*under100_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="cyan"), showlegend=False))
                fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 90) & (df['age_at_initial_pathologic_diagnosis'] < 100)) & (df['rcc'] == 'normal')][cg_value], x=[17]*under100_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="purple"), showlegend=False))
                fig.update_layout(
                    title='age plot of ' + str(cg_value) + ' (cohorts by 10)',
                        xaxis=dict(
                            tickmode='array',
                            tickvals=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                            ticktext=[
                                "<20 RCC (N = " + str(under20_rcc_count) + ")",
                                "<20 normal (N = " + str(under20_normal_count) + ")",
                                "<30 RCC (N = " + str(under30_rcc_count) + ")",
                                "<30 normal (N = " + str(under30_normal_count) + ")",
                                "<40 RCC (N = " + str(under40_rcc_count) + ")",
                                "<40 normal (N = " + str(under40_normal_count) + ")",
                                "<50 RCC (N = " + str(under50_rcc_count) + ")",
                                "<50 normal (N = " + str(under50_normal_count) + ")",
                                "<60 RCC (N = " + str(under60_rcc_count) + ")",
                                "<60 normal (N = " + str(under60_normal_count) + ")",
                                "<70 RCC (N = " + str(under70_rcc_count) + ")",
                                "<70 normal (N = " + str(under70_normal_count) + ")",
                                "<80 RCC (N = " + str(under80_rcc_count) + ")",
                                "<80 normal (N = " + str(under80_normal_count) + ")",
                                "<90 RCC (N = " + str(under90_rcc_count) + ")",
                                "<90 normal (N = " + str(under90_normal_count) + ")",
                                "<100 RCC (N = " + str(under100_rcc_count) + ")",
                                "<100 normal (N = " + str(under100_normal_count) + ")",
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
                        margin=dict(r=80, b=270),
                        height = 800,
                        width = 1000
                        )
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                _, p_value_under40 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 30) & (table["age_at_initial_pathologic_diagnosis"] < 40)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 30) & (table["age_at_initial_pathologic_diagnosis"] < 40)) & (table["rcc"] == "normal"), cg_value]
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
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 80) & (table["age_at_initial_pathologic_diagnosis"] < 90)) & (table["rcc"] == "normal"), cg_value]
                )
                _, p_value_under100 = stats.mannwhitneyu(
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 90) & (table["age_at_initial_pathologic_diagnosis"] < 100)) & (table["rcc"] == "rcc"), cg_value],
                table.loc[((table["age_at_initial_pathologic_diagnosis"] >= 90) & (table["age_at_initial_pathologic_diagnosis"] < 100)) & (table["rcc"] == "normal"), cg_value]
                )
                fig.add_annotation(
                    x=0.14,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 40 (RCC vs normal): {p_value_under40:.4f}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.16,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 50 (RCC vs normal): {p_value_under50:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.16,
                    y=-0.45,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 60 (RCC vs normal): {p_value_under60:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.16,
                    y=-0.5,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 70 (RCC vs normal): {p_value_under70:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.65,
                    y=-0.35,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 80 (RCC vs normal): {p_value_under80:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.65,
                    y=-0.4,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 90 (RCC vs normal): {p_value_under90:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                fig.add_annotation(
                    x=0.66,
                    y=-0.45,
                    xref='paper',
                    yref='paper',
                    text=f'p-value under 100 (RCC vs normal): {p_value_under100:.4e}',
                    showarrow=False,
                    font=dict(size=13),
                    align='center',
                    xanchor='center',
                    yanchor='top'
                )
                st.plotly_chart(fig)