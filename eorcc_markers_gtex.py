import pandas as pd
import streamlit as st
import plotly.graph_objs as go
import scipy.stats as stats
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

st.set_page_config(layout="wide")
st.title("GTEx Normal Kidney by Age")

SMOKER_COLORS = {
    'Never':   'darkturquoise',
    'Former':  'goldenrod',
    'Current': 'tomato',
    'NA':      'lightgrey',
}
AGE_BIN_ORDER = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89']


@st.cache_data
def load_pickle():
    pickle_path = Path(__file__).parent / 'data' / 'pickle_file.pk1'
    return pd.read_pickle(pickle_path)


@st.cache_data
def load_gtex():
    meth = pd.read_csv(
        Path(__file__).parent / 'data' / 'GTEx_Kidney.meth.csv',
        sep='\t'
    ).set_index('cgID')
    anno = pd.read_csv(
        Path(__file__).parent / 'data' / 'GTEx_Kidney.anno.csv',
        sep='\t', index_col=0
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


table = load_pickle()
meth_df, anno_df = load_gtex()
epic_anno = load_epic()


@st.cache_data
def compute_all_tests(subtype):
    if subtype == 'All':
        df = table[table['rcc'] == 'rcc'].copy()
    else:
        df = table[(table['rcc'] == 'rcc') & (table['subtype'] == subtype)].copy()
    df = df.dropna(subset=['age_at_initial_pathologic_diagnosis'])

    young_rcc = df[df['age_at_initial_pathologic_diagnosis'] < 50]
    old_rcc   = df[df['age_at_initial_pathologic_diagnosis'] >= 50]

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
    """Return a DataFrame with methylation, age_bin, smoker_status for each GTEx sample."""
    if cpg not in meth_df.index:
        return pd.DataFrame()
    meth_vals = meth_df.loc[cpg]
    rows = []
    for sample_id, meth_val in meth_vals.items():
        if pd.isna(meth_val):
            continue
        age_bin = anno_df.loc['age', sample_id] if 'age' in anno_df.index else None
        smoker  = anno_df.loc['smoker_status', sample_id] if 'smoker_status' in anno_df.index else None
        if pd.isna(smoker) or str(smoker).strip() == '':
            smoker = 'NA'
        rows.append({'sample_id': sample_id, 'methylation': meth_val,
                     'age_bin': age_bin, 'smoker_status': str(smoker)})
    return pd.DataFrame(rows).dropna(subset=['age_bin', 'methylation'])


# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    subtype_options = ['All', 'kirc', 'kirp', 'kich']
    selected_subtype = st.selectbox('Select RCC Subtype', options=subtype_options, index=1)
    top_n = st.number_input('Top N markers to display', min_value=1, max_value=200,
                            value=20, step=1)

all_results, rcc_df = compute_all_tests(selected_subtype)
sig_df = filter_results(all_results, top_n)

if sig_df.empty:
    st.warning("No significant markers found. Try a different subtype or increase top N.")
    st.stop()

n_young = (rcc_df['age_at_initial_pathologic_diagnosis'] < 50).sum()
n_old   = (rcc_df['age_at_initial_pathologic_diagnosis'] >= 50).sum()
st.markdown(
    f"**{len(sig_df)} significant eoRCC markers** | "
    f"eoRCC (<50): N={n_young} | loRCC (>=50): N={n_old}"
)

display_df = sig_df.join(epic_anno, on='cpg').rename(columns={
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
    display_df.style.format({
        'Mean eoRCC (<50)': '{:.4f}', 'Mean loRCC (>=50)': '{:.4f}',
        'Mean Diff (young - old)': '{:.4f}',
        'Raw p-value': '{:.3e}', 'FDR p-value': '{:.3e}'
    }),
    use_container_width=True
)

st.divider()
selected_cpg = st.selectbox('Select a CpG marker to visualize',
                             options=sig_df['cpg'].tolist(), index=0)

# ── GTEx data for selected CpG ───────────────────────────────────────────────
sample_df = build_gtex_sample_df(selected_cpg)

if sample_df.empty:
    st.warning(f"{selected_cpg} not found in GTEx methylation data.")
    st.stop()

# Keep only age bins present in the data, in chronological order
present_bins = [b for b in AGE_BIN_ORDER if b in sample_df['age_bin'].values]

# ── Build figure ─────────────────────────────────────────────────────────────
present_statuses = [s for s in ['Never', 'Former', 'Current', 'NA']
                    if s in sample_df['smoker_status'].values]

fig = go.Figure()

for status in present_statuses:
    status_df = sample_df[sample_df['smoker_status'] == status]
    x_vals, y_vals, text_vals = [], [], []
    for age_bin in present_bins:
        group = status_df[status_df['age_bin'] == age_bin]
        x_vals.extend([age_bin] * len(group))
        y_vals.extend(group['methylation'].tolist())
        text_vals.extend(group['sample_id'].tolist())
    fig.add_trace(go.Box(
        x=x_vals,
        y=y_vals,
        name=status,
        marker_color=SMOKER_COLORS[status],
        line_color=SMOKER_COLORS[status],
        fillcolor=SMOKER_COLORS[status],
        opacity=0.7,
        boxpoints='all',
        jitter=0.4,
        pointpos=0,
        marker=dict(size=6, line=dict(width=0.5, color='black')),
        text=text_vals,
        hovertemplate=(
            'Sample: %{text}<br>'
            'Methylation: %{y:.4f}<br>'
            'Smoker: ' + status +
            '<extra></extra>'
        ),
    ))

# Build x-axis tick labels with total N per age bin
tick_labels = {
    age_bin: f'{age_bin}<br>(N={len(sample_df[sample_df["age_bin"] == age_bin])})'
    for age_bin in present_bins
}

fig.update_layout(
    boxmode='group',
    title=f'GTEx Normal Kidney — {selected_cpg}',
    title_font=dict(size=26),
    yaxis=dict(
        range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
        tickfont=dict(size=18),
        title='<b>methylation ratio</b>', title_font=dict(size=20)
    ),
    xaxis=dict(
        categoryorder='array', categoryarray=present_bins,
        tickmode='array',
        tickvals=present_bins,
        ticktext=[tick_labels[b] for b in present_bins],
        tickfont=dict(size=18), title='<b>Age group</b>', title_font=dict(size=20)
    ),
    legend=dict(title='Smoker status', font=dict(size=16)),
    plot_bgcolor='rgba(0,0,0,0)',
    font=dict(family='Arial', size=20),
    margin=dict(r=20, b=120),
    height=700,
    width=1100,
)
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

st.plotly_chart(fig, use_container_width=False)
