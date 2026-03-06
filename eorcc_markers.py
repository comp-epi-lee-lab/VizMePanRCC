import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

st.set_page_config(layout="wide")
st.title("eoRCC Markers: early-onset (<50) vs late-onset (>=50) RCC")

@st.cache_data
def load_data():
    pickle_path = Path(__file__).parent / 'data' / 'pickle_file.pk1'
    return pd.read_pickle(pickle_path)


@st.cache_data
def load_epic():
    epic_path = Path(__file__).parent / 'data' / 'EPIC-8v2-0_A1.sim.v6.sorted.h.tsv'
    df = pd.read_csv(epic_path, sep='\t', usecols=['Name', 'UCSC_RefGene_Name',
                     'Relation_to_UCSC_CpG_Island', 'GencodeV41_Group'], index_col='Name')
    # Deduplicate semicolon-separated gene names
    df['UCSC_RefGene_Name'] = df['UCSC_RefGene_Name'].fillna('').apply(
        lambda x: '; '.join(dict.fromkeys(g for g in x.split(';') if g)) if x else ''
    )
    return df


table = load_data()
epic_anno = load_epic()

@st.cache_data
def compute_all_tests(subtype):
    """Heavy step: run MW tests on every CpG. Cached only by subtype."""
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
    fdr_vals = stats.false_discovery_control(res_df['p_value'].values)
    res_df['fdr'] = fdr_vals
    return res_df, df


def filter_results(res_df, top_n):
    """Cheap step: rank by FDR and return top N."""
    sig = res_df.sort_values('fdr')
    return sig.head(top_n).reset_index(drop=True)

with st.sidebar:
    subtype_options = ['All', 'kirc', 'kirp', 'kich']
    selected_subtype = st.selectbox('Select RCC Subtype', options=subtype_options, index=1)
    top_n = st.number_input('Top N markers to display', min_value=1, max_value=200,
                            value=20, step=1)
    run_pressed = st.button('Run Analysis')

run_pressed = True

all_results, rcc_df = compute_all_tests(selected_subtype)
sig_df = filter_results(all_results, top_n)

if sig_df.empty:
    st.warning("No significant markers found with the current thresholds. Try lowering the FDR or methylation difference threshold.")
else:
    n_young = (rcc_df['age_at_initial_pathologic_diagnosis'] < 50).sum()
    n_old = (rcc_df['age_at_initial_pathologic_diagnosis'] >= 50).sum()

    st.markdown(f"**{len(sig_df)} significant eoRCC markers found** | eoRCC (<50): N={n_young} | loRCC (>=50): N={n_old}")

    annotated = sig_df.join(epic_anno, on='cpg')
    display_df = annotated.rename(columns={
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
    display_df = display_df.style.format({
        'Mean eoRCC (<50)': '{:.4f}',
        'Mean loRCC (>=50)': '{:.4f}',
        'Mean Diff (young - old)': '{:.4f}',
        'Raw p-value': '{:.3e}',
        'FDR p-value': '{:.3e}'
    })
    st.dataframe(display_df, use_container_width=True)

    st.divider()
    selected_cpg = st.selectbox(
        'Select a CpG marker to visualize',
        options=sig_df['cpg'].tolist(),
        index=0
    )

    row = sig_df[sig_df['cpg'] == selected_cpg].iloc[0]
    young_vals = rcc_df.loc[rcc_df['age_at_initial_pathologic_diagnosis'] < 50, selected_cpg].dropna()
    old_vals = rcc_df.loc[rcc_df['age_at_initial_pathologic_diagnosis'] >= 50, selected_cpg].dropna()
    _, raw_p = stats.mannwhitneyu(young_vals, old_vals, alternative='two-sided')
    fdr_p = row['fdr']
    mean_diff = row['mean_diff']

    annotation_text = (
        f'eoRCC (<50) mean: {young_vals.mean():.4f} | loRCC (>=50) mean: {old_vals.mean():.4f}<br>'
        f'Mean diff (young - old): {mean_diff:.5f}<br>'
        f'Mann-Whitney p-value: {raw_p:.3e}' + (' *' if raw_p < 0.05 else '') +
        f' | FDR p-value: {fdr_p:.3e}' + (' *' if fdr_p < 0.05 else '')
    )

    common_layout = dict(
        yaxis=dict(range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                   tickfont=dict(size=18), title='<b>methylation ratio</b>',
                   title_font=dict(size=20)),
        xaxis=dict(
            tickmode='array', tickvals=[0, 1],
            ticktext=[
                f'eoRCC (<50)<br>(N={len(young_vals)})',
                f'loRCC (>=50)<br>(N={len(old_vals)})'
            ],
            tickfont=dict(size=18)
        ),
        plot_bgcolor='rgba(0,0,0,0)',
        font=dict(family='Arial', size=20),
        title_font=dict(size=26),
        margin=dict(r=20, b=260),
        height=700,
        width=600,
    )

    col_left, col_right = st.columns(2)

    with col_left:
        fig_box = go.Figure()
        fig_box.add_trace(go.Box(
            y=young_vals, x=[0] * len(young_vals),
            boxpoints='all', jitter=0.2,
            marker=dict(color='mediumpurple'),
            line=dict(color='mediumpurple'),
            showlegend=False
        ))
        fig_box.add_trace(go.Box(
            y=old_vals, x=[1] * len(old_vals),
            boxpoints='all', jitter=0.2,
            marker=dict(color='#8977ad'),
            line=dict(color='#8977ad'),
            showlegend=False
        ))
        fig_box.update_layout(
            title=f'RCC Samples: {selected_cpg}',
            **common_layout
        )
        fig_box.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_box.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_box.add_annotation(
            x=0.5, y=-0.35, xref='paper', yref='paper',
            text=annotation_text,
            showarrow=False, font=dict(size=16),
            align='center', xanchor='center', yanchor='top'
        )
        st.plotly_chart(fig_box, use_container_width=False)

    with col_right:
        if selected_subtype == 'All':
            normal_df = table.loc[
                (table['rcc'] != 'rcc') & table['age_at_initial_pathologic_diagnosis'].notna()
            ]
        else:
            normal_df = table.loc[
                (table['rcc'] != 'rcc') & (table['subtype'] == selected_subtype) & table['age_at_initial_pathologic_diagnosis'].notna()
            ]
        normal_young = normal_df.loc[normal_df['age_at_initial_pathologic_diagnosis'] < 50, selected_cpg].dropna()
        normal_old = normal_df.loc[normal_df['age_at_initial_pathologic_diagnosis'] >= 50, selected_cpg].dropna()
        normal_annotation_text = (
            f'Normal (<50) mean: {normal_young.mean():.4f} | Normal (>=50) mean: {normal_old.mean():.4f}<br>'
            f'N young={len(normal_young)} | N old={len(normal_old)}'
        )
        normal_layout = dict(
            yaxis=dict(range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                       tickfont=dict(size=18), title='<b>methylation ratio</b>',
                       title_font=dict(size=20)),
            xaxis=dict(
                tickmode='array', tickvals=[0, 1],
                ticktext=[
                    f'Normal (<50)<br>(N={len(normal_young)})',
                    f'Normal (>=50)<br>(N={len(normal_old)})'
                ],
                tickfont=dict(size=18)
            ),
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(family='Arial', size=20),
            title_font=dict(size=26),
            margin=dict(r=20, b=260),
            height=700,
            width=600,
        )
        fig_norm = go.Figure()
        fig_norm.add_trace(go.Box(
            y=normal_young, x=[0] * len(normal_young),
            boxpoints='all', jitter=0.2,
            marker=dict(color='darkturquoise'),
            line=dict(color='darkturquoise'),
            showlegend=False
        ))
        fig_norm.add_trace(go.Box(
            y=normal_old, x=[1] * len(normal_old),
            boxpoints='all', jitter=0.2,
            marker=dict(color='darkturquoise'),
            line=dict(color='darkturquoise'),
            showlegend=False
        ))
        fig_norm.update_layout(
            title=f'Normal samples: {selected_cpg}',
            **normal_layout
        )
        fig_norm.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_norm.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_norm.add_annotation(
            x=0.5, y=-0.35, xref='paper', yref='paper',
            text=normal_annotation_text,
            showarrow=False, font=dict(size=16),
            align='center', xanchor='center', yanchor='top'
        )
        st.plotly_chart(fig_norm, use_container_width=False)
