import numpy as np
import scipy.stats as stats
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import warnings
import plotly.express as px
warnings.filterwarnings("ignore")

req_cols = ["vital_status", "age_at_initial_pathologic_diagnosis", "days_to_death", "gender", "rcc", "stage", "cg05423956", "cg03544320", "cg01541443"]
cg_values = ["cg05423956", "cg03544320", "cg01541443"]
table = pd.read_csv(
    #'C:/Users/abhin/Downloads/kipan_clin_meth_20221210.tsv (1)/kipan_clin_meth_20221210.tsv', sep='\t', usecols=req_cols)
    '/Users/hayanlee/Project/Research/panRCC/TCGA.KIPAN/kipan_clin_meth.tsv', sep='\t', usecols=req_cols)
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
                    go.Box(y=df[cg_value], x=df['stage'], name="stage", boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"))
                )

                fig.add_trace(
                    go.Box(y=df[df['rcc'] == 'normal'][cg_value], name="normal", x=df[df['rcc'] == 'normal']['rcc'], boxpoints="all", jitter=0.2, marker=dict(color="cyan"))
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
                        dtick=0.1
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                )
                for i in range(1, 10):
                    fig.add_hline(y=i/10, line_color='black', line_width=1, opacity=0.2)
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)

            if gender:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                malecount = (df['gender'] == 'male').sum()
                femalecount = (df['gender'] == 'female').sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[cg_value], x=df['gender'], name="gender", boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"))
                )
                fig.update_layout(
                    title='gender plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1],
                        ticktext=[
                            "male (N = " + str(malecount) + ")",
                            "female (N = " + str(femalecount) + ")",
                        ]
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20),
                )
                for i in range(1, 10):
                    fig.add_hline(y=i/10, line_color='black', line_width=1, opacity=0.2)
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)

            if lts:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                undercount = (df['days_to_death'] <= 1825 ).sum()
                overcount = (df['days_to_death'] > 1825).sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[cg_value], x=df['days_to_death']<=1825, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"))
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
                        dtick=0.1
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20),
                )
                for i in range(1, 10):
                    fig.add_hline(y=i/10, line_color='black', line_width=1, opacity=0.2)
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)

            if age:
                df = table.sort_values(cg_value, ascending=False)
                df.head()
                under55count = (df['age_at_initial_pathologic_diagnosis'] <= 55 ).sum()
                over55count = (df['age_at_initial_pathologic_diagnosis'] > 55).sum()
                fig = go.Figure()
                fig.add_trace(
                    go.Box(y=df[cg_value], x=df['age_at_initial_pathologic_diagnosis']<=55, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"))
                )
                fig.update_layout(
                    title='age plot of ' + str(cg_value),
                    xaxis=dict(
                        tickmode='array',
                        tickvals=[0, 1],
                        ticktext=[
                            "age under 55 (N = " + str(under55count) + ")",
                            "age over 55 (N = " + str(over55count) + ")",
                        ]
                    ),
                    yaxis=dict(
                        tickmode='linear',
                        dtick=0.1
                    ),
                    plot_bgcolor='rgba(0, 0, 0, 0)',
                    font=dict(
                        family="Arial",
                        size=24
                    ),
                    margin=dict(r=20),
                )
                for i in range(1, 10):
                    fig.add_hline(y=i/10, line_color='black', line_width=1, opacity=0.2)
                fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
                fig.update_xaxes(
                    categoryorder='category ascending'
                )
                st.plotly_chart(fig)
