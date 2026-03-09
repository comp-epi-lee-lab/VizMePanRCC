import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

st.set_page_config(layout="wide")
st.title("VizMePanRCC")

@st.cache_data
def load_data():
    pickle_path = Path(__file__).parent / 'data' / 'pickle_file.pk1'
    return pd.read_pickle(pickle_path)


@st.cache_data
def load_gtex():
    meth = pd.read_csv(
        Path(__file__).parent / 'data' / 'GTEx_Kidney.meth.csv', sep='\t'
    ).set_index('cgID')
    anno = pd.read_csv(
        Path(__file__).parent / 'data' / 'GTEx_Kidney.anno.csv', sep='\t', index_col=0
    )
    return meth, anno


@st.cache_data
def load_epic():
    epic_path = Path(__file__).parent / 'data' / 'EPIC-8v2-0_A1.sim.v6.sorted.h.tsv'
    df = pd.read_csv(epic_path, sep='\t', usecols=['Name', 'UCSC_RefGene_Name',
                     'Relation_to_UCSC_CpG_Island', 'GencodeV41_Group'], index_col='Name')
    df['UCSC_RefGene_Name'] = df['UCSC_RefGene_Name'].fillna('').apply(
        lambda x: '; '.join(dict.fromkeys(g for g in x.split(';') if g)) if x else ''
    )
    return df


@st.cache_data
def load_cortex_expr():
    path = Path(__file__).parent / 'data' / 'gene_tpm_v11_kidney_cortex.gct.simple.tsv'
    return pd.read_csv(path, sep='\t', index_col='Description')


@st.cache_data
def load_medulla_expr():
    path = Path(__file__).parent / 'data' / 'gene_tpm_v11_kidney_medulla.gct.simple.tsv'
    return pd.read_csv(path, sep='\t', index_col='Description')


table = load_data()
meth_df, anno_df = load_gtex()
epic_anno = load_epic()
cortex_expr = load_cortex_expr()
medulla_expr = load_medulla_expr()

@st.cache_data
def compute_all_tests(subtype):
    if subtype == 'All':
        df = table[table['rcc'] == 'rcc'].copy()
    else:
        df = table[(table['rcc'] == 'rcc') & (table['subtype'] == subtype)].copy()
    df = df.dropna(subset=['age_at_initial_pathologic_diagnosis'])

    young_rcc = df[df['age_at_initial_pathologic_diagnosis'] < 50]
    old_rcc = df[df['age_at_initial_pathologic_diagnosis'] >= 50]

    cpg_columns = table.columns[17:]
    results = []
    for col in cpg_columns:
        y = young_rcc[col].dropna()
        o = old_rcc[col].dropna()
        if len(y) < 3 or len(o) < 3:
            continue
        _, p_val = stats.mannwhitneyu(y, o, alternative='two-sided')
        mean_diff = y.mean() - o.mean()
        results.append({'cpg': col, 'p_value': p_val, 'mean_diff': mean_diff,
                        'mean_young': y.mean(), 'mean_old': o.mean()})

    if not results:
        return pd.DataFrame(), df

    res_df = pd.DataFrame(results)
    res_df['fdr'] = stats.false_discovery_control(res_df['p_value'].values)
    return res_df, df


def filter_results(res_df, top_n):
    return res_df.sort_values('fdr').head(top_n).reset_index(drop=True)


def build_gtex_sample_df(cpg):
    if cpg not in meth_df.index:
        return pd.DataFrame()
    meth_vals = meth_df.loc[cpg]
    rows = []
    for sample_id, meth_val in meth_vals.items():
        if pd.isna(meth_val):
            continue
        age_bin = anno_df.loc['age', sample_id] if 'age' in anno_df.index else None
        smoker = anno_df.loc['smoker_status', sample_id] if 'smoker_status' in anno_df.index else None
        if pd.isna(smoker) or str(smoker).strip() == '':
            smoker = 'NA'
        rows.append({'sample_id': sample_id, 'methylation': meth_val,
                     'age_bin': age_bin, 'smoker_status': str(smoker)})
    return pd.DataFrame(rows).dropna(subset=['age_bin', 'methylation'])


def make_expr_figure(gene_name):
    """Return a Plotly figure of cortex vs medulla TPM for each gene in gene_name.
    gene_name may be '; '-separated for multi-gene CpGs."""
    genes = [g.strip() for g in gene_name.split(';') if g.strip()]
    fig = go.Figure()
    colors = {'Cortex': 'steelblue', 'Medulla': 'coral'}
    found_any = False
    for gene in genes:
        for tissue, df_expr in [('Cortex', cortex_expr), ('Medulla', medulla_expr)]:
            if gene in df_expr.index:
                vals = df_expr.loc[gene].values.astype(float)
                vals = vals[~np.isnan(vals)]
                fig.add_trace(go.Box(
                    y=vals,
                    name=f'{gene} — {tissue}',
                    boxpoints='all', jitter=0.3,
                    marker_color=colors[tissue],
                    line_color=colors[tissue],
                    showlegend=True,
                ))
                found_any = True
    fig.update_layout(
        yaxis_title='<b>TPM</b>',
        yaxis=dict(showgrid=False, tickfont=dict(size=16)),
        xaxis=dict(tickfont=dict(size=16)),
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(family='Arial', size=18),
        title_font=dict(size=24),
        height=550,
        width=900,
        margin=dict(r=20, b=80),
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    return fig, found_any


SMOKER_COLORS = {
    'Never': 'darkturquoise',
    'Former': 'goldenrod',
    'Current': 'tomato',
    'NA': 'lightgrey',
}
AGE_BIN_ORDER = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89']

tab1, tab2, tab3 = st.tabs(["CpG Explorer", "eoRCC Markers", "GTEx Normal Kidney"])

with tab1:
    st.header("CpG Explorer")

    ctrl_col, _ = st.columns([1, 3])
    with ctrl_col:
        selected_subtype_t1 = st.selectbox(
            'Select RCC Subtype', options=['All', 'kirc', 'kirp', 'kich'],
            index=1, key='t1_subtype'
        )
        cg_value = st.text_input(
            'Search CG Value', placeholder='cgXXXXXXXX',
            value='cg02275016', key='t1_cg'
        )
        age = st.checkbox("Age", key='t1_age', value=True)
        lts = st.checkbox("Long Term Survivorship", key='t1_lts', value=False)
        stage = st.checkbox("Stage", key='t1_stage', value=False)
        gender = st.checkbox("Gender", key='t1_gender', value=False)

    if cg_value not in table.columns:
        st.warning('No data found for the provided CG value.')
    else:
        if selected_subtype_t1 == "All":
            df = table.copy()
        elif selected_subtype_t1 == "kich":
            df_kich = table[table["subtype"] == "kich"]
            df_normals = table[(table["subtype"].isin(["kirc", "kirp"])) & (table["rcc"] == "normal")]
            df = pd.concat([df_kich, df_normals], ignore_index=True)
        else:
            df = table[table["subtype"] == selected_subtype_t1]
        df = df.sort_values(cg_value, ascending=False)
        df = df.dropna(subset=[cg_value])

        if age:
            under50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] < 50) & (df['rcc'] == 'rcc')).sum()
            under50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] < 50) & (df['rcc'] == 'normal')).sum()
            over50_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['rcc'] == 'rcc')).sum()
            over50_normal_count = ((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['rcc'] == 'normal')).sum()
            fig = go.Figure()
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] < 50) & (df['rcc'] == 'normal')][cg_value],
                                x=[0]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['rcc'] == 'normal')][cg_value],
                                x=[1]*over50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="#3fa6a8"), showlegend=False))
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] < 50) & (df['rcc'] == 'rcc')][cg_value],
                                x=[2]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['rcc'] == 'rcc')][cg_value],
                                x=[3]*over50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="#8977ad"), showlegend=False))
            fig.update_layout(
                title='age plot of ' + str(cg_value) + ' (young [< 50] vs old [>= 50])',
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
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                plot_bgcolor='rgba(0, 0, 0, 0)',
                font=dict(family="Arial", size=24),
                margin=dict(r=20, b=320),
                height=750,
                width=1000
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_xaxes(categoryorder='category ascending')
            old_rcc = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 50) & (df["rcc"] == "rcc"), cg_value].dropna()
            young_rcc = df.loc[(df["age_at_initial_pathologic_diagnosis"] < 50) & (df["rcc"] == "rcc"), cg_value].dropna()
            old_normal = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 50) & (df["rcc"] == "normal"), cg_value].dropna()
            young_normal = df.loc[(df["age_at_initial_pathologic_diagnosis"] < 50) & (df["rcc"] == "normal"), cg_value].dropna()
            difference_rcc_age = old_rcc.mean() - young_rcc.mean()
            _, p_value_rcc_age = stats.mannwhitneyu(young_rcc, old_rcc)
            difference_normal_age = old_normal.mean() - young_normal.mean()
            _, p_value_normal_age = stats.mannwhitneyu(young_normal, old_normal)
            difference_under50 = young_normal.mean() - young_rcc.mean()
            _, p_value_under50 = stats.mannwhitneyu(young_normal, young_rcc)
            difference_over50 = old_normal.mean() - old_rcc.mean()
            _, p_value_over50 = stats.mannwhitneyu(old_normal, old_rcc)
            fig.add_annotation(x=0.5, y=-0.4, xref='paper', yref='paper',
                text='mean differences & p-values:', showarrow=False, font=dict(size=22),
                align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.5, xref='paper', yref='paper',
                text=f'RCC old - RCC young: {difference_rcc_age:.5f} | p-value: {p_value_rcc_age:.5e}' + ('*' if p_value_rcc_age < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.6, xref='paper', yref='paper',
                text=f'normal old - normal young: {difference_normal_age:.5f} | p-value: {p_value_normal_age:.5e}' + ('*' if p_value_normal_age < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.7, xref='paper', yref='paper',
                text=f'normal young - RCC young: {difference_under50:.5f} | p-value: {p_value_under50:.5e}' + ('*' if p_value_under50 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.8, xref='paper', yref='paper',
                text=f'normal old - RCC old: {difference_over50:.5f} | p-value: {p_value_over50:.5e}' + ('*' if p_value_over50 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            st.plotly_chart(fig)

        if age:
            rcc_40 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 20) & (df["age_at_initial_pathologic_diagnosis"] < 40) & (df["rcc"] == "rcc"), cg_value].dropna()
            normal_40 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 20) & (df["age_at_initial_pathologic_diagnosis"] < 40) & (df["rcc"] == "normal"), cg_value].dropna()
            rcc_50 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 40) & (df["age_at_initial_pathologic_diagnosis"] < 50) & (df["rcc"] == "rcc"), cg_value].dropna()
            normal_50 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 40) & (df["age_at_initial_pathologic_diagnosis"] < 50) & (df["rcc"] == "normal"), cg_value].dropna()
            rcc_60 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 50) & (df["age_at_initial_pathologic_diagnosis"] < 60) & (df["rcc"] == "rcc"), cg_value].dropna()
            normal_60 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 50) & (df["age_at_initial_pathologic_diagnosis"] < 60) & (df["rcc"] == "normal"), cg_value].dropna()
            rcc_70 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 60) & (df["age_at_initial_pathologic_diagnosis"] < 70) & (df["rcc"] == "rcc"), cg_value].dropna()
            normal_70 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 60) & (df["age_at_initial_pathologic_diagnosis"] < 70) & (df["rcc"] == "normal"), cg_value].dropna()
            rcc_o70 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 70) & (df["rcc"] == "rcc"), cg_value].dropna()
            normal_o70 = df.loc[(df["age_at_initial_pathologic_diagnosis"] >= 70) & (df["rcc"] == "normal"), cg_value].dropna()

            under40_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')).sum()
            under40_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')).sum()
            under50_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')).sum()
            under50_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')).sum()
            under60_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')).sum()
            under60_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')).sum()
            under70_rcc_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')).sum()
            under70_normal_count = (((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')).sum()
            over70_rcc_count = ((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['rcc'] == 'rcc')).sum()
            over70_normal_count = ((df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['rcc'] == 'normal')).sum()

            mean_values_rcc = [rcc_40.mean(), rcc_50.mean(), rcc_60.mean(), rcc_70.mean(), rcc_o70.mean()]
            mean_values_normal = [normal_40.mean(), normal_50.mean(), normal_60.mean(), normal_70.mean(), normal_o70.mean()]
            labels = [
                "20-39<br>(N = " + str(under40_normal_count) + "|" + str(under40_rcc_count) + ")",
                "40-49<br>(N = " + str(under50_normal_count) + "|" + str(under50_rcc_count) + ")",
                "50-59<br>(N = " + str(under60_normal_count) + "|" + str(under60_rcc_count) + ")",
                "60-69<br>(N = " + str(under70_normal_count) + "|" + str(under70_rcc_count) + ")",
                "70+<br>(N = " + str(over70_normal_count) + "|" + str(over70_rcc_count) + ")"
            ]
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=labels, y=mean_values_normal, mode='lines+markers', name='Normal',
                                     line=dict(color='darkturquoise'), marker=dict(size=12)))
            fig.add_trace(go.Scatter(x=labels, y=mean_values_rcc, mode='lines+markers', name='RCC',
                                     line=dict(color='mediumpurple'), marker=dict(size=12)))
            fig.update_layout(
                title='age line plot of ' + str(cg_value),
                title_font=dict(size=30),
                xaxis_title='<b>age ranges</b>',
                yaxis_title='<b>mean methylation ratio</b>',
                xaxis_title_font=dict(size=22),
                yaxis_title_font=dict(size=22),
                xaxis=dict(tickfont=dict(size=22)),
                font=dict(family="Arial"),
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                legend=dict(orientation="h", font=dict(size=22), y=-0.4),
                margin=dict(r=80, b=150),
                height=600,
                width=600
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig)

            fig = go.Figure()
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'normal')][cg_value], x=[0]*under40_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 20) & (df['age_at_initial_pathologic_diagnosis'] < 40)) & (df['rcc'] == 'rcc')][cg_value], x=[1]*under40_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'normal')][cg_value], x=[2]*under50_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 40) & (df['age_at_initial_pathologic_diagnosis'] < 50)) & (df['rcc'] == 'rcc')][cg_value], x=[3]*under50_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'normal')][cg_value], x=[4]*under60_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 50) & (df['age_at_initial_pathologic_diagnosis'] < 60)) & (df['rcc'] == 'rcc')][cg_value], x=[5]*under60_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'normal')][cg_value], x=[6]*under70_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[((df['age_at_initial_pathologic_diagnosis'] >= 60) & (df['age_at_initial_pathologic_diagnosis'] < 70)) & (df['rcc'] == 'rcc')][cg_value], x=[7]*under70_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['rcc'] == 'normal')][cg_value], x=[8]*over70_normal_count, boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[(df['age_at_initial_pathologic_diagnosis'] >= 70) & (df['rcc'] == 'rcc')][cg_value], x=[9]*over70_rcc_count, boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
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
                        "70+ normal (N = " + str(over70_normal_count) + ")",
                        "70+ RCC (N = " + str(over70_rcc_count) + ")",
                    ],
                    tickfont=dict(size=16)
                ),
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                plot_bgcolor='rgba(0, 0, 0, 0)',
                font=dict(family="Arial", size=24),
                margin=dict(r=120, b=500),
                height=950,
                width=1100
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_xaxes(categoryorder='category ascending')
            _, p_value_under40 = stats.mannwhitneyu(rcc_40, normal_40)
            _, p_value_under50 = stats.mannwhitneyu(rcc_50, normal_50)
            _, p_value_under60 = stats.mannwhitneyu(rcc_60, normal_60)
            _, p_value_under70 = stats.mannwhitneyu(rcc_70, normal_70)
            _, p_value_over70 = stats.mannwhitneyu(rcc_o70, normal_o70)
            fig.add_annotation(x=0.5, y=-0.5, xref='paper', yref='paper',
                text='p-values (RCC vs normal):', showarrow=False, font=dict(size=22),
                align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.6, xref='paper', yref='paper',
                text=f'20-39: {p_value_under40:.5e}' + ('*' if p_value_under40 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.7, xref='paper', yref='paper',
                text=f'40-49: {p_value_under50:.5e}' + ('*' if p_value_under50 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.8, xref='paper', yref='paper',
                text=f'50-59: {p_value_under60:.5e}' + ('*' if p_value_under60 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.9, xref='paper', yref='paper',
                text=f'60-69: {p_value_under70:.5e}' + ('*' if p_value_under70 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-1.00, xref='paper', yref='paper',
                text=f'70+: {p_value_over70:.5e}' + ('*' if p_value_over70 < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            st.plotly_chart(fig)

        if lts:
            undercount = (df['days_to_death'] <= 1825).sum()
            overcount = (df['days_to_death'] > 1825).sum()
            fig = go.Figure()
            fig.add_trace(go.Box(y=df[df['days_to_death'] <= 1825][cg_value], x=[0]*undercount,
                                  boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=df[df['days_to_death'] > 1825][cg_value], x=[1]*overcount,
                                  boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.update_layout(
                title='long term survivorship plot of ' + str(cg_value),
                title_font=dict(size=30),
                yaxis_title='<b>methylation ratio</b>',
                yaxis_title_font=dict(size=22),
                xaxis_title='<b>years until death</b>',
                xaxis_title_font=dict(size=22),
                xaxis=dict(
                    tickmode='array', tickvals=[0, 1],
                    ticktext=["under 5 years<br>(N = " + str(undercount) + ")",
                              "over 5 years<br>(N = " + str(overcount) + ")"],
                    tickfont=dict(size=20)
                ),
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                plot_bgcolor='rgba(0, 0, 0, 0)',
                font=dict(family="Arial", size=24),
                margin=dict(r=20, b=170),
                height=600,
                width=1000
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_xaxes(categoryorder='category ascending')
            less_days_to_death = df.loc[df["days_to_death"] <= 1825, cg_value].dropna()
            more_days_to_death = df.loc[df["days_to_death"] > 1825, cg_value].dropna()
            difference = less_days_to_death.mean() - more_days_to_death.mean()
            _, p_value = stats.mannwhitneyu(less_days_to_death, more_days_to_death)
            fig.add_annotation(x=0.5, y=-0.4, xref='paper', yref='paper',
                text=f'Mean difference (under 5 years - over 5 years): {difference:.5f} | p-value: {p_value:.5e}' + ('*' if p_value < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            st.plotly_chart(fig)

        if stage:
            s1count = (df['stage'] == 'stage i').sum()
            s2count = (df['stage'] == 'stage ii').sum()
            s3count = (df['stage'] == 'stage iii').sum()
            s4count = (df['stage'] == 'stage iv').sum()
            normalcount = (df['rcc'] == 'normal').sum()
            fig = go.Figure()
            fig.add_trace(go.Box(y=df[cg_value], x=df['stage'], name="stage",
                                  boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=df[df['rcc'] == 'normal'][cg_value], name="normal",
                                  x=df[df['rcc'] == 'normal']['rcc'],
                                  boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.update_layout(
                title='stage plot of ' + str(cg_value),
                title_font=dict(size=30),
                yaxis_title='<b>methylation ratio</b>',
                yaxis_title_font=dict(size=22),
                xaxis_title='<b>stage</b>',
                xaxis_title_font=dict(size=22),
                xaxis=dict(
                    tickmode='array', tickvals=[0, 1, 2, 3, 4],
                    ticktext=[
                        "normal<br>(N = " + str(normalcount) + ")",
                        "stage i<br>(N = " + str(s1count) + ")",
                        "stage ii<br>(N = " + str(s2count) + ")",
                        "stage iii<br>(N = " + str(s3count) + ")",
                        "stage iv<br>(N = " + str(s4count) + ")"
                    ],
                    tickfont=dict(size=19)
                ),
                font=dict(family="Arial", size=24),
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                plot_bgcolor='rgba(0, 0, 0, 0)',
                margin=dict(r=20),
                height=600,
                width=1000
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_xaxes(categoryorder='category ascending')
            st.plotly_chart(fig)

        if gender:
            male_rcc_count = ((df['gender'] == 'male') & (df['rcc'] == 'rcc')).sum()
            male_normal_count = ((df['gender'] == 'male') & (df['rcc'] == 'normal')).sum()
            female_rcc_count = ((df['gender'] == 'female') & (df['rcc'] == 'rcc')).sum()
            female_normal_count = ((df['gender'] == 'female') & (df['rcc'] == 'normal')).sum()
            rcc_male = df[(df['gender'] == 'male') & (df['rcc'] == 'rcc')]
            normal_male = df[(df['gender'] == 'male') & (df['rcc'] == 'normal')]
            rcc_female = df[(df['gender'] == 'female') & (df['rcc'] == 'rcc')]
            normal_female = df[(df['gender'] == 'female') & (df['rcc'] == 'normal')]
            fig = go.Figure()
            fig.add_trace(go.Box(y=normal_male[cg_value], x=[0]*len(normal_male),
                                  boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=normal_female[cg_value], x=[1]*len(normal_female),
                                  boxpoints="all", jitter=0.2, marker=dict(color="darkturquoise"), showlegend=False))
            fig.add_trace(go.Box(y=rcc_male[cg_value], x=[2]*len(rcc_male),
                                  boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.add_trace(go.Box(y=rcc_female[cg_value], x=[3]*len(rcc_female),
                                  boxpoints="all", jitter=0.2, marker=dict(color="mediumpurple"), showlegend=False))
            fig.update_layout(
                title='gender plot of ' + str(cg_value),
                title_font=dict(size=30),
                yaxis_title='<b>methylation ratio</b>',
                yaxis_title_font=dict(size=22),
                xaxis_title='<b>gender and condition</b>',
                xaxis_title_font=dict(size=22),
                xaxis=dict(
                    tickmode='array', tickvals=[0, 1, 2, 3],
                    ticktext=[
                        "normal male<br>(N = " + str(male_normal_count) + ")",
                        "normal female<br>(N = " + str(female_normal_count) + ")",
                        "RCC male<br>(N = " + str(male_rcc_count) + ")",
                        "RCC female<br>(N = " + str(female_rcc_count) + ")",
                    ],
                    tickfont=dict(size=20)
                ),
                yaxis=dict(tickmode='linear', dtick=0.1, showgrid=False, tickfont=dict(size=22)),
                plot_bgcolor='rgba(0, 0, 0, 0)',
                font=dict(family="Arial", size=24),
                margin=dict(r=20, b=330, l=20),
                height=750,
                width=1000
            )
            fig.update_layout(yaxis_range=[0, 1.0])
            fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            rcc_male_values = rcc_male[cg_value].dropna()
            rcc_female_values = rcc_female[cg_value].dropna()
            normal_male_values = normal_male[cg_value].dropna()
            normal_female_values = normal_female[cg_value].dropna()
            difference_rcc_male_vs_female = rcc_female_values.mean() - rcc_male_values.mean()
            _, p_value_rcc_male_vs_female = stats.mannwhitneyu(rcc_male_values, rcc_female_values)
            difference_normal_male_vs_female = normal_female_values.mean() - normal_male_values.mean()
            _, p_value_normal_male_vs_female = stats.mannwhitneyu(normal_male_values, normal_female_values)
            difference_male_vs_normal = normal_male_values.mean() - rcc_male_values.mean()
            _, p_value_male_vs_normal = stats.mannwhitneyu(normal_male_values, rcc_male_values)
            difference_female_vs_normal = normal_female_values.mean() - rcc_female_values.mean()
            _, p_value_female_vs_normal = stats.mannwhitneyu(normal_female_values, rcc_female_values)
            fig.add_annotation(x=0.5, y=-0.4, xref='paper', yref='paper',
                text='mean differences & p-values:', showarrow=False, font=dict(size=22),
                align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.5, xref='paper', yref='paper',
                text=f'RCC female - RCC male: {difference_rcc_male_vs_female:.5f} | p-value: {p_value_rcc_male_vs_female:.5e}' + ('*' if p_value_rcc_male_vs_female < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.6, xref='paper', yref='paper',
                text=f'normal female - normal male: {difference_normal_male_vs_female:.5f} | p-value: {p_value_normal_male_vs_female:.5e}' + ('*' if p_value_normal_male_vs_female < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.7, xref='paper', yref='paper',
                text=f'normal male - RCC male: {difference_male_vs_normal:.5f} | p-value: {p_value_male_vs_normal:.5e}' + ('*' if p_value_male_vs_normal < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            fig.add_annotation(x=0.5, y=-0.8, xref='paper', yref='paper',
                text=f'normal female - RCC female: {difference_female_vs_normal:.5f} | p-value: {p_value_female_vs_normal:.5e}' + ('*' if p_value_female_vs_normal < 0.05 else ''),
                showarrow=False, font=dict(size=22), align='center', xanchor='center', yanchor='top')
            st.plotly_chart(fig, use_container_width=False)


# ── Tab 2: eoRCC Markers (eorcc_markers.py) ───────────────────────────────────
with tab2:
    st.header("eoRCC Markers: early-onset (<50) vs late-onset (>=50) RCC")

    ctrl_col2, _ = st.columns([1, 3])
    with ctrl_col2:
        selected_subtype_t2 = st.selectbox(
            'Select RCC Subtype', options=['All', 'kirc', 'kirp', 'kich'],
            index=1, key='t2_subtype'
        )
        top_n_t2 = st.number_input(
            'Top N markers to display', min_value=1, max_value=200,
            value=20, step=1, key='t2_top_n'
        )

    all_results_t2, rcc_df_t2 = compute_all_tests(selected_subtype_t2)
    sig_df_t2 = filter_results(all_results_t2, top_n_t2)

    if sig_df_t2.empty:
        st.warning("No significant markers found with the current thresholds. Try lowering the FDR or methylation difference threshold.")
    else:
        n_young_t2 = (rcc_df_t2['age_at_initial_pathologic_diagnosis'] < 50).sum()
        n_old_t2 = (rcc_df_t2['age_at_initial_pathologic_diagnosis'] >= 50).sum()
        st.markdown(f"**{len(sig_df_t2)} significant eoRCC markers found** | eoRCC (<50): N={n_young_t2} | loRCC (>=50): N={n_old_t2}")

        annotated_t2 = sig_df_t2.join(epic_anno, on='cpg')
        display_df_t2 = annotated_t2.rename(columns={
            'cpg': 'CpG ID',
            'mean_young': 'Mean eoRCC (<50)',
            'mean_old': 'Mean loRCC (>=50)',
            'mean_diff': 'Mean Diff (young - old)',
            'fdr': 'FDR p-value',
            'p_value': 'Raw p-value',
            'UCSC_RefGene_Name': 'Gene',
            'Relation_to_UCSC_CpG_Island': 'CpG Island Relation',
            'GencodeV41_Group': 'GencodeV41 Group',
        })[['CpG ID', 'Gene', 'CpG Island Relation', 'GencodeV41 Group',
            'Mean eoRCC (<50)', 'Mean loRCC (>=50)', 'Mean Diff (young - old)',
            'Raw p-value', 'FDR p-value']]
        st.dataframe(display_df_t2.style.format({
            'Mean eoRCC (<50)': '{:.4f}',
            'Mean loRCC (>=50)': '{:.4f}',
            'Mean Diff (young - old)': '{:.4f}',
            'Raw p-value': '{:.3e}',
            'FDR p-value': '{:.3e}'
        }), use_container_width=True)

        st.divider()
        selected_cpg_t2 = st.selectbox(
            'Select a CpG marker to visualize',
            options=sig_df_t2['cpg'].tolist(), index=0, key='t2_cpg'
        )

        gene_t2 = annotated_t2.loc[annotated_t2['cpg'] == selected_cpg_t2, 'UCSC_RefGene_Name'].values[0]
        cpg_label_t2 = f'{selected_cpg_t2} ({gene_t2})' if gene_t2 else selected_cpg_t2

        row_t2 = sig_df_t2[sig_df_t2['cpg'] == selected_cpg_t2].iloc[0]
        young_vals_t2 = rcc_df_t2.loc[rcc_df_t2['age_at_initial_pathologic_diagnosis'] < 50, selected_cpg_t2].dropna()
        old_vals_t2 = rcc_df_t2.loc[rcc_df_t2['age_at_initial_pathologic_diagnosis'] >= 50, selected_cpg_t2].dropna()
        _, raw_p_t2 = stats.mannwhitneyu(young_vals_t2, old_vals_t2, alternative='two-sided')
        fdr_p_t2 = row_t2['fdr']
        mean_diff_t2 = row_t2['mean_diff']

        annotation_text_t2 = (
            f'eoRCC (<50) mean: {young_vals_t2.mean():.4f} | loRCC (>=50) mean: {old_vals_t2.mean():.4f}<br>'
            f'Mean diff (young - old): {mean_diff_t2:.5f}<br>'
            f'Mann-Whitney p-value: {raw_p_t2:.3e}' + (' *' if raw_p_t2 < 0.05 else '') +
            f' | FDR p-value: {fdr_p_t2:.3e}' + (' *' if fdr_p_t2 < 0.05 else '')
        )

        common_layout_t2 = dict(
            yaxis=dict(range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                       tickfont=dict(size=18), title='<b>methylation ratio</b>',
                       title_font=dict(size=20)),
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(family='Arial', size=20),
            title_font=dict(size=26),
            margin=dict(r=20, b=260),
            height=700,
            width=600,
        )

        col_left_t2, col_right_t2 = st.columns(2)

        with col_left_t2:
            fig_box_t2 = go.Figure()
            fig_box_t2.add_trace(go.Box(
                y=young_vals_t2, x=[0]*len(young_vals_t2),
                boxpoints='all', jitter=0.2,
                marker=dict(color='mediumpurple'), line=dict(color='mediumpurple'), showlegend=False
            ))
            fig_box_t2.add_trace(go.Box(
                y=old_vals_t2, x=[1]*len(old_vals_t2),
                boxpoints='all', jitter=0.2,
                marker=dict(color='#8977ad'), line=dict(color='#8977ad'), showlegend=False
            ))
            fig_box_t2.update_layout(
                title=f'RCC: {cpg_label_t2}',
                xaxis=dict(
                    tickmode='array', tickvals=[0, 1],
                    ticktext=[
                        f'eoRCC (<50)<br>(N={len(young_vals_t2)})',
                        f'loRCC (>=50)<br>(N={len(old_vals_t2)})'
                    ],
                    tickfont=dict(size=18)
                ),
                **common_layout_t2
            )
            fig_box_t2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_box_t2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_box_t2.add_annotation(
                x=0.5, y=-0.35, xref='paper', yref='paper',
                text=annotation_text_t2, showarrow=False, font=dict(size=16),
                align='center', xanchor='center', yanchor='top'
            )
            st.plotly_chart(fig_box_t2, use_container_width=False)

        with col_right_t2:
            if selected_subtype_t2 == 'All':
                normal_df_t2 = table.loc[
                    (table['rcc'] != 'rcc') & table['age_at_initial_pathologic_diagnosis'].notna()
                ]
            else:
                normal_df_t2 = table.loc[
                    (table['rcc'] != 'rcc') & (table['subtype'] == selected_subtype_t2) &
                    table['age_at_initial_pathologic_diagnosis'].notna()
                ]
            normal_young_t2 = normal_df_t2.loc[normal_df_t2['age_at_initial_pathologic_diagnosis'] < 50, selected_cpg_t2].dropna()
            normal_old_t2 = normal_df_t2.loc[normal_df_t2['age_at_initial_pathologic_diagnosis'] >= 50, selected_cpg_t2].dropna()
            normal_annotation_t2 = (
                f'Normal (<50) mean: {normal_young_t2.mean():.4f} | Normal (>=50) mean: {normal_old_t2.mean():.4f}<br>'
                f'N young={len(normal_young_t2)} | N old={len(normal_old_t2)}'
            )
            fig_norm_t2 = go.Figure()
            fig_norm_t2.add_trace(go.Box(
                y=normal_young_t2, x=[0]*len(normal_young_t2),
                boxpoints='all', jitter=0.2,
                marker=dict(color='darkturquoise'), line=dict(color='darkturquoise'), showlegend=False
            ))
            fig_norm_t2.add_trace(go.Box(
                y=normal_old_t2, x=[1]*len(normal_old_t2),
                boxpoints='all', jitter=0.2,
                marker=dict(color='darkturquoise'), line=dict(color='darkturquoise'), showlegend=False
            ))
            fig_norm_t2.update_layout(
                title=f'Normal: {cpg_label_t2}',
                xaxis=dict(
                    tickmode='array', tickvals=[0, 1],
                    ticktext=[
                        f'Normal (<50)<br>(N={len(normal_young_t2)})',
                        f'Normal (>=50)<br>(N={len(normal_old_t2)})'
                    ],
                    tickfont=dict(size=18)
                ),
                yaxis=dict(range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                           tickfont=dict(size=18), title='<b>methylation ratio</b>',
                           title_font=dict(size=20)),
                plot_bgcolor='rgba(0,0,0,0)',
                font=dict(family='Arial', size=20),
                title_font=dict(size=26),
                margin=dict(r=20, b=260),
                height=700,
                width=600,
            )
            fig_norm_t2.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_norm_t2.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_norm_t2.add_annotation(
                x=0.5, y=-0.35, xref='paper', yref='paper',
                text=normal_annotation_t2, showarrow=False, font=dict(size=16),
                align='center', xanchor='center', yanchor='top'
            )
            st.plotly_chart(fig_norm_t2, use_container_width=False)

        if gene_t2:
            st.divider()
            fig_expr_t2, found_t2 = make_expr_figure(gene_t2)
            if found_t2:
                fig_expr_t2.update_layout(title=f'GTEx Kidney Gene Expression — {gene_t2}')
                st.plotly_chart(fig_expr_t2, use_container_width=False, key='expr_t2')
            else:
                st.info(f"No expression data found for **{gene_t2}**.")


# ── Tab 3: GTEx Normal Kidney (eorcc_markers_gtex.py) ────────────────────────
with tab3:
    st.header("GTEx Normal Kidney by Age")

    ctrl_col3, _ = st.columns([1, 3])
    with ctrl_col3:
        selected_subtype_t3 = st.selectbox(
            'Select RCC Subtype', options=['All', 'kirc', 'kirp', 'kich'],
            index=1, key='t3_subtype'
        )
        top_n_t3 = st.number_input(
            'Top N markers to display', min_value=1, max_value=200,
            value=20, step=1, key='t3_top_n'
        )

    all_results_t3, rcc_df_t3 = compute_all_tests(selected_subtype_t3)
    sig_df_t3 = filter_results(all_results_t3, top_n_t3)

    if sig_df_t3.empty:
        st.warning("No significant markers found. Try a different subtype or increase top N.")
    else:
        n_young_t3 = (rcc_df_t3['age_at_initial_pathologic_diagnosis'] < 50).sum()
        n_old_t3 = (rcc_df_t3['age_at_initial_pathologic_diagnosis'] >= 50).sum()
        st.markdown(
            f"**{len(sig_df_t3)} significant eoRCC markers** | "
            f"eoRCC (<50): N={n_young_t3} | loRCC (>=50): N={n_old_t3}"
        )

        display_df_t3 = sig_df_t3.join(epic_anno, on='cpg').rename(columns={
            'cpg': 'CpG ID',
            'mean_young': 'Mean eoRCC (<50)',
            'mean_old': 'Mean loRCC (>=50)',
            'mean_diff': 'Mean Diff (young - old)',
            'fdr': 'FDR p-value',
            'p_value': 'Raw p-value',
            'UCSC_RefGene_Name': 'Gene',
            'Relation_to_UCSC_CpG_Island': 'CpG Island Relation',
            'GencodeV41_Group': 'GencodeV41 Group',
        })[['CpG ID', 'Gene', 'CpG Island Relation', 'GencodeV41 Group',
            'Mean eoRCC (<50)', 'Mean loRCC (>=50)', 'Mean Diff (young - old)',
            'Raw p-value', 'FDR p-value']]
        st.dataframe(
            display_df_t3.style.format({
                'Mean eoRCC (<50)': '{:.4f}', 'Mean loRCC (>=50)': '{:.4f}',
                'Mean Diff (young - old)': '{:.4f}',
                'Raw p-value': '{:.3e}', 'FDR p-value': '{:.3e}'
            }),
            use_container_width=True
        )

        st.divider()
        selected_cpg_t3 = st.selectbox(
            'Select a CpG marker to visualize',
            options=sig_df_t3['cpg'].tolist(), index=0, key='t3_cpg'
        )

        gene_t3 = display_df_t3.loc[display_df_t3['CpG ID'] == selected_cpg_t3, 'Gene'].values[0]
        cpg_label_t3 = f'{selected_cpg_t3} ({gene_t3})' if gene_t3 else selected_cpg_t3

        sample_df_t3 = build_gtex_sample_df(selected_cpg_t3)

        if sample_df_t3.empty:
            st.warning(f"{selected_cpg_t3} not found in GTEx methylation data.")
        else:
            present_bins_t3 = [b for b in AGE_BIN_ORDER if b in sample_df_t3['age_bin'].values]
            present_statuses_t3 = [s for s in ['Never', 'Former', 'Current', 'NA']
                                    if s in sample_df_t3['smoker_status'].values]

            fig_gtex = go.Figure()
            for status in present_statuses_t3:
                status_df = sample_df_t3[sample_df_t3['smoker_status'] == status]
                x_vals, y_vals, text_vals = [], [], []
                for age_bin in present_bins_t3:
                    group = status_df[status_df['age_bin'] == age_bin]
                    x_vals.extend([age_bin] * len(group))
                    y_vals.extend(group['methylation'].tolist())
                    text_vals.extend(group['sample_id'].tolist())
                fig_gtex.add_trace(go.Box(
                    x=x_vals, y=y_vals,
                    name=status,
                    marker_color=SMOKER_COLORS[status],
                    line_color=SMOKER_COLORS[status],
                    fillcolor=SMOKER_COLORS[status],
                    opacity=0.7,
                    boxpoints='all', jitter=0.4, pointpos=0,
                    marker=dict(size=6, line=dict(width=0.5, color='black')),
                    text=text_vals,
                    hovertemplate=(
                        'Sample: %{text}<br>'
                        'Methylation: %{y:.4f}<br>'
                        'Smoker: ' + status +
                        '<extra></extra>'
                    ),
                ))

            tick_labels_t3 = {
                age_bin: f'{age_bin}<br>(N={len(sample_df_t3[sample_df_t3["age_bin"] == age_bin])})'
                for age_bin in present_bins_t3
            }
            fig_gtex.update_layout(
                boxmode='group',
                title=f'GTEx Normal Kidney — {cpg_label_t3}',
                title_font=dict(size=26),
                yaxis=dict(
                    range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                    tickfont=dict(size=18),
                    title='<b>methylation ratio</b>', title_font=dict(size=20)
                ),
                xaxis=dict(
                    categoryorder='array', categoryarray=present_bins_t3,
                    tickmode='array',
                    tickvals=present_bins_t3,
                    ticktext=[tick_labels_t3[b] for b in present_bins_t3],
                    tickfont=dict(size=18), title='<b>Age group</b>', title_font=dict(size=20)
                ),
                legend=dict(title='Smoker status', font=dict(size=16)),
                plot_bgcolor='rgba(0,0,0,0)',
                font=dict(family='Arial', size=20),
                margin=dict(r=20, b=120),
                height=700,
                width=1100,
            )
            fig_gtex.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_gtex.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig_gtex, use_container_width=False)

            if gene_t3:
                st.divider()
                fig_expr_t3, found_t3 = make_expr_figure(gene_t3)
                if found_t3:
                    fig_expr_t3.update_layout(title=f'GTEx Kidney Gene Expression — {gene_t3}')
                    st.plotly_chart(fig_expr_t3, use_container_width=False, key='expr_t3')
                else:
                    st.info(f"No expression data found for **{gene_t3}**.")
