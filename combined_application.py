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
    return (
        res_df[res_df['p_value'] < 0.05]
        .sort_values('mean_diff', ascending=True)
        .head(top_n)
        .reset_index(drop=True)
    )


@st.cache_data
def compute_tcga_age_correlation(subtype):
    if subtype == 'All':
        df = table[table['rcc'] == 'rcc'].copy()
    else:
        df = table[(table['rcc'] == 'rcc') & (table['subtype'] == subtype)].copy()
    df = df.dropna(subset=['age_at_initial_pathologic_diagnosis'])
    ages = df['age_at_initial_pathologic_diagnosis']

    cpg_columns = table.columns[17:]
    results = []
    for col in cpg_columns:
        vals = df[col]
        mask = vals.notna()
        if mask.sum() < 3:
            continue
        r, p_val = stats.pearsonr(vals[mask], ages[mask])
        if np.isnan(p_val):
            continue
        results.append({
            'cpg': col, 'r': r, 'p_value': p_val,
            'n_samples': int(mask.sum()),
            'mean_methylation': float(vals[mask].mean()),
        })
    if not results:
        return pd.DataFrame(), df
    return pd.DataFrame(results), df


def filter_tcga_age_results(res_df, top_n):
    return (
        res_df[res_df['p_value'] < 0.05]
        .sort_values('r', ascending=True)
        .head(top_n)
        .reset_index(drop=True)
    )


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
    """Return a Plotly figure of cortex TPM for each gene in gene_name.
    gene_name may be '; '-separated for multi-gene CpGs."""
    genes = [g.strip() for g in gene_name.split(';') if g.strip()]
    fig = go.Figure()
    found_any = False
    for gene in genes:
        if gene in cortex_expr.index:
            vals = cortex_expr.loc[gene].values.astype(float)
            vals = vals[~np.isnan(vals)]
            fig.add_trace(go.Box(
                y=vals,
                name=gene,
                boxpoints='all', jitter=0.3,
                marker_color='steelblue',
                line_color='steelblue',
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


EXPR_TARGET_GENES = sorted([
    g for g in cortex_expr.index
    if (g.startswith('CD163') or g.startswith('IL') or g.startswith('CCL') or
        g.startswith('CXCL') or g.startswith('CXCR') or g.startswith('CX3') or
        g == 'CHI3L1')
])


def _build_expr_age_subset():
    anno_donor_map = {col: '-'.join(col.split('-')[:2]) for col in anno_df.columns}
    cortex_donor_map = {col: '-'.join(col.split('-')[:2]) for col in cortex_expr.columns}
    overlap_donors = set(anno_donor_map.values()) & set(cortex_donor_map.values())
    donor_to_age = {
        donor: anno_df.loc['age', sample]
        for sample, donor in anno_donor_map.items()
        if donor in overlap_donors
    }
    cortex_samples = {
        col: cortex_donor_map[col]
        for col in cortex_expr.columns
        if cortex_donor_map[col] in overlap_donors
    }
    sub = cortex_expr[list(cortex_samples.keys())].T.copy()
    sub['age_bin'] = [donor_to_age[cortex_samples[s]] for s in sub.index]
    sub['age_group'] = sub['age_bin'].apply(
        lambda ab: 'young' if int(ab.split('-')[0]) < 50 else 'old'
    )
    return sub


expr_age_df = _build_expr_age_subset()


EXPR_RATIO_EPS = 1e-6


def compute_expr_age_tests():
    young = expr_age_df[expr_age_df['age_group'] == 'young']
    old = expr_age_df[expr_age_df['age_group'] == 'old']
    age_numeric_all = expr_age_df['age_bin'].apply(lambda ab: int(ab.split('-')[0])).astype(float)
    results = []
    for gene in EXPR_TARGET_GENES:
        if gene not in expr_age_df.columns:
            continue
        y_vals = young[gene].dropna().values.astype(float)
        o_vals = old[gene].dropna().values.astype(float)
        if len(y_vals) < 2 or len(o_vals) < 2:
            continue
        if max(y_vals.max(), o_vals.max()) < 0.01:
            continue
        try:
            _, p_val = stats.ttest_ind(y_vals, o_vals, equal_var=False)
            p_val = float(p_val)
        except Exception:
            p_val = float('nan')
        ratio = float(o_vals.mean()) / (float(y_vals.mean()) + EXPR_RATIO_EPS)

        pearson_r, pearson_p, spearman_r, spearman_p = np.nan, np.nan, np.nan, np.nan
        gene_vals = expr_age_df[gene].dropna()
        if len(gene_vals) >= 3:
            ages_for_gene = age_numeric_all.loc[gene_vals.index].values
            try:
                pearson_r, pearson_p = stats.pearsonr(ages_for_gene, gene_vals.values.astype(float))
                pearson_r, pearson_p = float(pearson_r), float(pearson_p)
                spearman_r, spearman_p = stats.spearmanr(ages_for_gene, gene_vals.values.astype(float))
                spearman_r, spearman_p = float(spearman_r), float(spearman_p)
            except Exception:
                pearson_r, pearson_p, spearman_r, spearman_p = np.nan, np.nan, np.nan, np.nan

        results.append({
            'gene': gene,
            'mean_young': float(y_vals.mean()),
            'mean_old': float(o_vals.mean()),
            'ratio': ratio,
            'significant_ratio': ratio > 2,
            'p_value': p_val,
            'age_pearson_r': pearson_r,
            'age_pearson_p': pearson_p,
            'age_spearman_r': spearman_r,
            'age_spearman_p': spearman_p,
        })
    res_df = pd.DataFrame(results)
    res_df = res_df[res_df['p_value'] < 0.20]
    return res_df.sort_values('age_pearson_r', ascending=False).reset_index(drop=True)


@st.cache_data
def compute_gene_age_correlation(genes):
    age_numeric = expr_age_df['age_bin'].apply(lambda ab: int(ab.split('-')[0])).values.astype(float)
    results = []
    for gene in genes:
        if gene not in expr_age_df.columns:
            continue
        vals = expr_age_df[gene].dropna().values.astype(float)
        if len(vals) < 3 or vals.max() < 0.01:
            continue
        try:
            r, p = stats.pearsonr(vals, age_numeric[:len(vals)])
        except Exception:
            continue
        if np.isnan(p):
            continue
        results.append({'gene': gene, 'gene_r': r, 'gene_p_value': p})
    return pd.DataFrame(results).set_index('gene') if results else pd.DataFrame(columns=['gene_r', 'gene_p_value'])


def attach_gene_correlation(ranked_df, gene_col='Gene'):
    all_genes = set()
    for g in ranked_df[gene_col]:
        all_genes.update(x.strip() for x in str(g).split(';') if x.strip())
    gene_corr = compute_gene_age_correlation(tuple(sorted(all_genes)))

    def pick_best(gene_str):
        candidates = [x.strip() for x in str(gene_str).split(';') if x.strip()]
        matches = gene_corr.loc[gene_corr.index.intersection(candidates)]
        if matches.empty:
            return pd.Series({'Linked Gene': None, 'Gene r': np.nan, 'Gene p-value': np.nan})
        best = matches.sort_values('gene_p_value').iloc[0]
        return pd.Series({
            'Linked Gene': matches.sort_values('gene_p_value').index[0],
            'Gene r': best['gene_r'],
            'Gene p-value': best['gene_p_value'],
        })

    extra = ranked_df[gene_col].apply(pick_best)
    return pd.concat([ranked_df, extra], axis=1)


SMOKER_COLORS = {
    'Never': 'darkturquoise',
    'Former': 'goldenrod',
    'Current': 'tomato',
    'NA': 'lightgrey',
}
AGE_BIN_ORDER = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79', '80-89']

tab1, tab2, tab3, tab4, tab5 = st.tabs(["TCGA CpG Explorer", "TCGA eoRCC CpG Markers", "GTEx CpG Age Association", "GTEx CpG Smoke Association", "GTEx Gene Expression w/ Aging"])

with tab1:
    st.header("CpG Explorer")
    st.caption("Dataset: TCGA PanRCC methylation (TCGA-KIRC, TCGA-KIRP, TCGA-KICH)")

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
                title='age plot of ' + str(cg_value) + ' (young [< 50 years old] vs old [>= 50 years old])',
                title_font=dict(size=30),
                yaxis_title='<b>methylation ratio</b>',
                yaxis_title_font=dict(size=22),
                xaxis_title='<b>age and condition</b>',
                xaxis_title_font=dict(size=22),
                xaxis=dict(
                    tickmode='array',
                    tickvals=[0, 1, 2, 3],
                    ticktext=[
                        "young normal<br>(<50 years old)<br>(N = " + str(under50_normal_count) + ")",
                        "old normal<br>(>=50 years old)<br>(N = " + str(over50_normal_count) + ")",
                        "early-onset RCC<br>(<50 years old)<br>(N = " + str(under50_rcc_count) + ")",
                        "late-onset RCC<br>(>=50 years old)<br>(N = " + str(over50_rcc_count) + ")",
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
                "20-39 years old<br>(N = " + str(under40_normal_count) + "|" + str(under40_rcc_count) + ")",
                "40-49 years old<br>(N = " + str(under50_normal_count) + "|" + str(under50_rcc_count) + ")",
                "50-59 years old<br>(N = " + str(under60_normal_count) + "|" + str(under60_rcc_count) + ")",
                "60-69 years old<br>(N = " + str(under70_normal_count) + "|" + str(under70_rcc_count) + ")",
                "70+ years old<br>(N = " + str(over70_normal_count) + "|" + str(over70_rcc_count) + ")"
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
                        "20-39 years old normal (N = " + str(under40_normal_count) + ")",
                        "20-39 years old RCC (N = " + str(under40_rcc_count) + ")",
                        "40-49 years old normal (N = " + str(under50_normal_count) + ")",
                        "40-49 years old RCC (N = " + str(under50_rcc_count) + ")",
                        "50-59 years old normal (N = " + str(under60_normal_count) + ")",
                        "50-59 years old RCC (N = " + str(under60_rcc_count) + ")",
                        "60-69 years old normal (N = " + str(under70_normal_count) + ")",
                        "60-69 years old RCC (N = " + str(under70_rcc_count) + ")",
                        "70+ years old normal (N = " + str(over70_normal_count) + ")",
                        "70+ years old RCC (N = " + str(over70_rcc_count) + ")",
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
    st.caption("Dataset: TCGA PanRCC methylation (TCGA-KIRC, TCGA-KIRP, TCGA-KICH)")

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
            'Raw p-value']]
        st.dataframe(display_df_t2.style.format({
            'Mean eoRCC (<50)': '{:.4f}',
            'Mean loRCC (>=50)': '{:.4f}',
            'Mean Diff (young - old)': '{:.4f}',
            'Raw p-value': '{:.3e}',
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
            f'eoRCC (<50 years old) mean: {young_vals_t2.mean():.4f} | loRCC (>=50 years old) mean: {old_vals_t2.mean():.4f}<br>'
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
                        f'eoRCC (<50 years old)<br>(N={len(young_vals_t2)})',
                        f'loRCC (>=50 years old)<br>(N={len(old_vals_t2)})'
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
                f'Normal (<50 years old) mean: {normal_young_t2.mean():.4f} | Normal (>=50 years old) mean: {normal_old_t2.mean():.4f}<br>'
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
                        f'Normal (<50 years old)<br>(N={len(normal_young_t2)})',
                        f'Normal (>=50 years old)<br>(N={len(normal_old_t2)})'
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
            st.info("This shows overall gene expression (TPM) across all samples. For gene TPM broken down by age group, go to the **GTEx Normal Kidney by Age** tab.")
            fig_expr_t2, found_t2 = make_expr_figure(gene_t2)
            if found_t2:
                fig_expr_t2.update_layout(title=f'GTEx Kidney Gene Expression — {gene_t2}')
                st.plotly_chart(fig_expr_t2, use_container_width=False, key='expr_t2')
            else:
                st.info(f"No expression data found for **{gene_t2}**.")


# ── Tab 3: GTEx Normal Kidney (eorcc_markers_gtex.py) ────────────────────────
with tab3:
    st.header("GTEx Normal Kidney by Age")
    st.caption("Methylation: GTEx v11 Kidney Cortex (β values) | Expression: GTEx v11 Kidney Cortex & Medulla (TPM)")

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
        custom_cg_t3 = st.text_input(
            'Add custom cgIDs (comma-separated)',
            placeholder='cg12345678, cg87654321',
            key='t3_custom_cg'
        )

    all_results_t3, rcc_df_t3 = compute_tcga_age_correlation(selected_subtype_t3)
    sig_df_t3 = filter_tcga_age_results(all_results_t3, top_n_t3)

    if sig_df_t3.empty:
        st.warning("No significant markers found. Try a different subtype or increase top N.")
    else:
        st.markdown(
            f"**{len(sig_df_t3)} significant CpGs** (Pearson r with age, p<0.05) | "
            f"N={len(rcc_df_t3)} RCC patients, age range "
            f"{rcc_df_t3['age_at_initial_pathologic_diagnosis'].min():.0f}-"
            f"{rcc_df_t3['age_at_initial_pathologic_diagnosis'].max():.0f}"
        )

        display_df_t3 = sig_df_t3.join(epic_anno, on='cpg').rename(columns={
            'cpg': 'CpG ID',
            'r': 'Pearson r (meth vs age)',
            'p_value': 'Raw p-value',
            'mean_methylation': 'Mean Methylation',
            'n_samples': 'N samples',
            'UCSC_RefGene_Name': 'Gene',
            'Relation_to_UCSC_CpG_Island': 'CpG Island Relation',
            'GencodeV41_Group': 'GencodeV41 Group',
        })
        display_df_t3 = attach_gene_correlation(display_df_t3, gene_col='Gene')
        display_df_t3 = display_df_t3[['CpG ID', 'Gene', 'CpG Island Relation', 'GencodeV41 Group',
            'Pearson r (meth vs age)', 'Raw p-value', 'Mean Methylation', 'N samples',
            'Linked Gene', 'Gene r', 'Gene p-value']]
        st.caption(
            "Linked Gene / Gene r / Gene p-value use only 5 GTEx donors with matched "
            "RNA-seq + methylation age data — interpret with caution."
        )
        st.dataframe(
            display_df_t3.style.format({
                'Pearson r (meth vs age)': '{:.3f}',
                'Raw p-value': '{:.3e}',
                'Mean Methylation': '{:.4f}',
                'Gene r': '{:.3f}',
                'Gene p-value': '{:.3e}',
            }),
            use_container_width=True
        )

        st.divider()
        custom_cgs_t3 = [c.strip() for c in custom_cg_t3.split(',') if c.strip()] if custom_cg_t3 else []
        base_options_t3 = sig_df_t3['cpg'].tolist()
        cpg_options_t3 = base_options_t3 + [c for c in custom_cgs_t3 if c not in base_options_t3]
        selected_cpg_t3 = st.selectbox(
            'Select a CpG marker to visualize',
            options=cpg_options_t3, index=0, key='t3_cpg'
        )

        match_t3 = display_df_t3.loc[display_df_t3['CpG ID'] == selected_cpg_t3, 'Gene']
        if not match_t3.empty:
            gene_t3 = match_t3.values[0]
        elif selected_cpg_t3 in epic_anno.index:
            gene_t3 = epic_anno.loc[selected_cpg_t3, 'UCSC_RefGene_Name']
        else:
            gene_t3 = ''
        cpg_label_t3 = f'{selected_cpg_t3} ({gene_t3})' if gene_t3 else selected_cpg_t3

        row_t3 = sig_df_t3[sig_df_t3['cpg'] == selected_cpg_t3]
        if not row_t3.empty and selected_cpg_t3 in rcc_df_t3.columns:
            r_t3 = row_t3['r'].values[0]
            p_t3 = row_t3['p_value'].values[0]
            ages_t3 = rcc_df_t3['age_at_initial_pathologic_diagnosis']
            meth_t3 = rcc_df_t3[selected_cpg_t3]
            valid_t3 = meth_t3.notna()

            fig_scatter_t3 = go.Figure()
            fig_scatter_t3.add_trace(go.Scatter(
                x=ages_t3[valid_t3], y=meth_t3[valid_t3], mode='markers',
                marker=dict(color='mediumpurple', size=8, line=dict(width=0.5, color='black')),
                showlegend=False,
            ))
            if valid_t3.sum() >= 2:
                slope, intercept = np.polyfit(ages_t3[valid_t3], meth_t3[valid_t3], 1)
                x_line = np.array([ages_t3[valid_t3].min(), ages_t3[valid_t3].max()])
                fig_scatter_t3.add_trace(go.Scatter(
                    x=x_line, y=slope * x_line + intercept, mode='lines',
                    line=dict(color='firebrick', dash='dash'), showlegend=False,
                ))
            fig_scatter_t3.update_layout(
                title=f'TCGA RCC: Methylation vs Age — {cpg_label_t3}<br><sup>Pearson r={r_t3:.3f}, p={p_t3:.3e}</sup>',
                xaxis_title='<b>Age at diagnosis</b>', yaxis_title='<b>methylation ratio</b>',
                yaxis=dict(range=[0, 1.0], showgrid=False),
                plot_bgcolor='rgba(0,0,0,0)', font=dict(family='Arial', size=18),
                height=600, width=900, margin=dict(r=20, b=80),
            )
            fig_scatter_t3.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_scatter_t3.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig_scatter_t3, use_container_width=False, key='scatter_t3')

        st.divider()
        sample_df_t3 = build_gtex_sample_df(selected_cpg_t3)

        if sample_df_t3.empty:
            st.warning(f"{selected_cpg_t3} not found in GTEx methylation data.")
        else:
            present_bins_t3 = [b for b in AGE_BIN_ORDER if b in sample_df_t3['age_bin'].values]

            fig_gtex = go.Figure()
            x_vals, y_vals, text_vals = [], [], []
            for age_bin in present_bins_t3:
                group = sample_df_t3[sample_df_t3['age_bin'] == age_bin]
                x_vals.extend([age_bin] * len(group))
                y_vals.extend(group['methylation'].tolist())
                text_vals.extend(group['sample_id'].tolist())
            fig_gtex.add_trace(go.Box(
                x=x_vals, y=y_vals,
                marker_color='steelblue',
                line_color='steelblue',
                boxpoints='all', jitter=0.4, pointpos=0,
                marker=dict(size=6, line=dict(width=0.5, color='black')),
                text=text_vals,
                hovertemplate='Sample: %{text}<br>Methylation: %{y:.4f}<extra></extra>',
                showlegend=False,
            ))

            tick_labels_t3 = {
                age_bin: f'{age_bin} years old<br>(N={len(sample_df_t3[sample_df_t3["age_bin"] == age_bin])})'
                for age_bin in present_bins_t3
            }
            fig_gtex.update_layout(
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


# ── Tab 4: GTEx Smoker Status ─────────────────────────────────────────────────
with tab4:
    st.header("GTEx Normal Kidney by Smoker Status")
    st.caption("Methylation: GTEx v11 Kidney Cortex (β values) | Expression: GTEx v11 Kidney Cortex & Medulla (TPM)")

    ctrl_col4, _ = st.columns([1, 3])
    with ctrl_col4:
        selected_subtype_t4 = st.selectbox(
            'Select RCC Subtype', options=['All', 'kirc', 'kirp', 'kich'],
            index=1, key='t4_subtype'
        )
        top_n_t4 = st.number_input(
            'Top N markers to display', min_value=1, max_value=200,
            value=20, step=1, key='t4_top_n'
        )
        custom_cg_t4 = st.text_input(
            'Add custom cgIDs (comma-separated)',
            placeholder='cg12345678, cg87654321',
            key='t4_custom_cg'
        )

    all_results_t4, rcc_df_t4 = compute_tcga_age_correlation(selected_subtype_t4)
    sig_df_t4 = filter_tcga_age_results(all_results_t4, top_n_t4)

    if sig_df_t4.empty:
        st.warning("No significant markers found. Try a different subtype or increase top N.")
    else:
        st.markdown(
            f"**{len(sig_df_t4)} significant CpGs** (Pearson r with age, p<0.05) | "
            f"N={len(rcc_df_t4)} RCC patients, age range "
            f"{rcc_df_t4['age_at_initial_pathologic_diagnosis'].min():.0f}-"
            f"{rcc_df_t4['age_at_initial_pathologic_diagnosis'].max():.0f}"
        )

        display_df_t4 = sig_df_t4.join(epic_anno, on='cpg').rename(columns={
            'cpg': 'CpG ID',
            'r': 'Pearson r (meth vs age)',
            'p_value': 'Raw p-value',
            'mean_methylation': 'Mean Methylation',
            'n_samples': 'N samples',
            'UCSC_RefGene_Name': 'Gene',
            'Relation_to_UCSC_CpG_Island': 'CpG Island Relation',
            'GencodeV41_Group': 'GencodeV41 Group',
        })
        display_df_t4 = attach_gene_correlation(display_df_t4, gene_col='Gene')
        display_df_t4 = display_df_t4[['CpG ID', 'Gene', 'CpG Island Relation', 'GencodeV41 Group',
            'Pearson r (meth vs age)', 'Raw p-value', 'Mean Methylation', 'N samples',
            'Linked Gene', 'Gene r', 'Gene p-value']]
        st.caption(
            "Linked Gene / Gene r / Gene p-value use only 5 GTEx donors with matched "
            "RNA-seq + methylation age data — interpret with caution."
        )
        st.dataframe(
            display_df_t4.style.format({
                'Pearson r (meth vs age)': '{:.3f}',
                'Raw p-value': '{:.3e}',
                'Mean Methylation': '{:.4f}',
                'Gene r': '{:.3f}',
                'Gene p-value': '{:.3e}',
            }),
            use_container_width=True
        )

        st.divider()
        custom_cgs_t4 = [c.strip() for c in custom_cg_t4.split(',') if c.strip()] if custom_cg_t4 else []
        base_options_t4 = sig_df_t4['cpg'].tolist()
        cpg_options_t4 = base_options_t4 + [c for c in custom_cgs_t4 if c not in base_options_t4]
        selected_cpg_t4 = st.selectbox(
            'Select a CpG marker to visualize',
            options=cpg_options_t4, index=0, key='t4_cpg'
        )

        match_t4 = display_df_t4.loc[display_df_t4['CpG ID'] == selected_cpg_t4, 'Gene']
        if not match_t4.empty:
            gene_t4 = match_t4.values[0]
        elif selected_cpg_t4 in epic_anno.index:
            gene_t4 = epic_anno.loc[selected_cpg_t4, 'UCSC_RefGene_Name']
        else:
            gene_t4 = ''
        cpg_label_t4 = f'{selected_cpg_t4} ({gene_t4})' if gene_t4 else selected_cpg_t4

        row_t4 = sig_df_t4[sig_df_t4['cpg'] == selected_cpg_t4]
        if not row_t4.empty and selected_cpg_t4 in rcc_df_t4.columns:
            r_t4 = row_t4['r'].values[0]
            p_t4 = row_t4['p_value'].values[0]
            ages_t4 = rcc_df_t4['age_at_initial_pathologic_diagnosis']
            meth_t4 = rcc_df_t4[selected_cpg_t4]
            valid_t4 = meth_t4.notna()

            fig_scatter_t4 = go.Figure()
            fig_scatter_t4.add_trace(go.Scatter(
                x=ages_t4[valid_t4], y=meth_t4[valid_t4], mode='markers',
                marker=dict(color='mediumpurple', size=8, line=dict(width=0.5, color='black')),
                showlegend=False,
            ))
            if valid_t4.sum() >= 2:
                slope, intercept = np.polyfit(ages_t4[valid_t4], meth_t4[valid_t4], 1)
                x_line = np.array([ages_t4[valid_t4].min(), ages_t4[valid_t4].max()])
                fig_scatter_t4.add_trace(go.Scatter(
                    x=x_line, y=slope * x_line + intercept, mode='lines',
                    line=dict(color='firebrick', dash='dash'), showlegend=False,
                ))
            fig_scatter_t4.update_layout(
                title=f'TCGA RCC: Methylation vs Age — {cpg_label_t4}<br><sup>Pearson r={r_t4:.3f}, p={p_t4:.3e}</sup>',
                xaxis_title='<b>Age at diagnosis</b>', yaxis_title='<b>methylation ratio</b>',
                yaxis=dict(range=[0, 1.0], showgrid=False),
                plot_bgcolor='rgba(0,0,0,0)', font=dict(family='Arial', size=18),
                height=600, width=900, margin=dict(r=20, b=80),
            )
            fig_scatter_t4.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_scatter_t4.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig_scatter_t4, use_container_width=False, key='scatter_t4')

        st.divider()
        sample_df_t4 = build_gtex_sample_df(selected_cpg_t4)

        if sample_df_t4.empty:
            st.warning(f"{selected_cpg_t4} not found in GTEx methylation data.")
        else:
            present_bins_t4 = [b for b in AGE_BIN_ORDER if b in sample_df_t4['age_bin'].values]
            present_statuses_t4 = [s for s in ['Never', 'Former', 'Current', 'NA']
                                    if s in sample_df_t4['smoker_status'].values]

            fig_gtex_t4 = go.Figure()
            for status in present_statuses_t4:
                status_df = sample_df_t4[sample_df_t4['smoker_status'] == status]
                x_vals, y_vals, text_vals = [], [], []
                for age_bin in present_bins_t4:
                    group = status_df[status_df['age_bin'] == age_bin]
                    x_vals.extend([age_bin] * len(group))
                    y_vals.extend(group['methylation'].tolist())
                    text_vals.extend(group['sample_id'].tolist())
                fig_gtex_t4.add_trace(go.Box(
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

            tick_labels_t4 = {
                age_bin: f'{age_bin} years old<br>(N={len(sample_df_t4[sample_df_t4["age_bin"] == age_bin])})'
                for age_bin in present_bins_t4
            }
            fig_gtex_t4.update_layout(
                boxmode='group',
                title=f'GTEx Normal Kidney by Smoker Status — {cpg_label_t4}',
                title_font=dict(size=26),
                yaxis=dict(
                    range=[0, 1.0], tickmode='linear', dtick=0.1, showgrid=False,
                    tickfont=dict(size=18),
                    title='<b>methylation ratio</b>', title_font=dict(size=20)
                ),
                xaxis=dict(
                    categoryorder='array', categoryarray=present_bins_t4,
                    tickmode='array',
                    tickvals=present_bins_t4,
                    ticktext=[tick_labels_t4[b] for b in present_bins_t4],
                    tickfont=dict(size=18), title='<b>Age group</b>', title_font=dict(size=20)
                ),
                legend=dict(title='Smoker status', font=dict(size=16)),
                plot_bgcolor='rgba(0,0,0,0)',
                font=dict(family='Arial', size=20),
                margin=dict(r=20, b=120),
                height=700,
                width=1100,
            )
            fig_gtex_t4.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_gtex_t4.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig_gtex_t4, use_container_width=False)

            if gene_t4:
                st.divider()
                fig_expr_t4, found_t4 = make_expr_figure(gene_t4)
                if found_t4:
                    fig_expr_t4.update_layout(title=f'GTEx Kidney Gene Expression — {gene_t4}')
                    st.plotly_chart(fig_expr_t4, use_container_width=False, key='expr_t4')
                else:
                    st.info(f"No expression data found for **{gene_t4}**.")


# ── Tab 5: Gene Expression by Age ─────────────────────────────────────────────
with tab5:
    st.header("Gene Expression by Age (Cortex Only)")
    st.caption("GTEx v11 Kidney Cortex TPM | Age linked via overlapping donors with GTEx methylation anno")

    n_young_t5 = (expr_age_df['age_group'] == 'young').sum()
    n_old_t5 = (expr_age_df['age_group'] == 'old').sum()
    st.info(
        f"Analysis uses {len(expr_age_df)} donors with matched age data from the methylation anno file. "
        f"Young (<50): N={n_young_t5} | Old (≥50): N={n_old_t5}. "
        f"Screened to p<0.20 (relaxed from the app's usual p<0.05) given the small sample size."
    )

    expr_results_t5 = compute_expr_age_tests()

    if expr_results_t5.empty:
        st.warning("No gene expression data found for target gene families.")
    else:
        st.dataframe(
            expr_results_t5.rename(columns={
                'gene': 'Gene',
                'mean_young': 'Mean Young (<50)',
                'mean_old': 'Mean Old (≥50)',
                'ratio': 'Old/Young Ratio',
                'significant_ratio': 'Significant (ratio>2)',
                'p_value': 'p-value',
                'age_pearson_r': 'Pearson r (vs age)',
                'age_pearson_p': 'Pearson p-value',
                'age_spearman_r': 'Spearman r (vs age)',
                'age_spearman_p': 'Spearman p-value',
            }).style.format({
                'Mean Young (<50)': '{:.3f}',
                'Mean Old (≥50)': '{:.3f}',
                'Old/Young Ratio': '{:.3f}',
                'p-value': '{:.3e}',
                'Pearson r (vs age)': '{:.3f}',
                'Pearson p-value': '{:.3e}',
                'Spearman r (vs age)': '{:.3f}',
                'Spearman p-value': '{:.3e}',
            }),
            use_container_width=True
        )

        st.divider()
        selected_gene_t5 = st.selectbox(
            'Select gene to visualize',
            options=expr_results_t5['gene'].tolist(),
            index=0, key='t5_gene'
        )

        young_t5 = expr_age_df[expr_age_df['age_group'] == 'young']
        old_t5 = expr_age_df[expr_age_df['age_group'] == 'old']
        young_vals_t5 = young_t5[selected_gene_t5].dropna().values.astype(float)
        old_vals_t5 = old_t5[selected_gene_t5].dropna().values.astype(float)
        young_ages_t5 = young_t5['age_bin'].values
        old_ages_t5 = old_t5['age_bin'].values
        row_t5 = expr_results_t5[expr_results_t5['gene'] == selected_gene_t5].iloc[0]

        fig_t5 = go.Figure()
        fig_t5.add_trace(go.Box(
            y=young_vals_t5, x=[0] * len(young_vals_t5),
            boxpoints='all', jitter=0.2,
            marker=dict(color='mediumpurple'), line=dict(color='mediumpurple'),
            text=list(young_ages_t5),
            hovertemplate='%{text}<br>TPM: %{y:.3f}<extra></extra>',
            showlegend=False,
        ))
        fig_t5.add_trace(go.Box(
            y=old_vals_t5, x=[1] * len(old_vals_t5),
            boxpoints='all', jitter=0.2,
            marker=dict(color='tomato'), line=dict(color='tomato'),
            text=list(old_ages_t5),
            hovertemplate='%{text}<br>TPM: %{y:.3f}<extra></extra>',
            showlegend=False,
        ))
        fig_t5.update_layout(
            title=f'Gene Expression by Age — {selected_gene_t5}',
            title_font=dict(size=26),
            yaxis_title='<b>TPM</b>',
            xaxis=dict(
                tickmode='array', tickvals=[0, 1],
                ticktext=[
                    f'Young (<50 years old)<br>(N={len(young_vals_t5)})',
                    f'Old (≥50 years old)<br>(N={len(old_vals_t5)})',
                ],
                tickfont=dict(size=18),
            ),
            yaxis=dict(showgrid=False, tickfont=dict(size=16)),
            plot_bgcolor='rgba(0,0,0,0)',
            font=dict(family='Arial', size=18),
            margin=dict(r=20, b=200),
            height=600,
            width=600,
        )
        fig_t5.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_t5.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
        fig_t5.add_annotation(
            x=0.5, y=-0.3, xref='paper', yref='paper',
            text=(
                f'Mean young: {row_t5["mean_young"]:.3f} | Mean old: {row_t5["mean_old"]:.3f}<br>'
                f'Old/Young Ratio: {row_t5["ratio"]:.3f}'
                + (' *significant (ratio>2)*' if row_t5['significant_ratio'] else '') + '<br>'
                f't-test p-value: {row_t5["p_value"]:.3e}'
            ),
            showarrow=False, font=dict(size=16),
            align='center', xanchor='center', yanchor='top'
        )
        st.plotly_chart(fig_t5, use_container_width=False)

        age_numeric_t5 = expr_age_df['age_bin'].apply(lambda ab: int(ab.split('-')[0])).astype(float)
        gene_vals_t5 = expr_age_df[selected_gene_t5].dropna()
        if len(gene_vals_t5) >= 2:
            ages_t5 = age_numeric_t5.loc[gene_vals_t5.index]

            fig_scatter_t5 = go.Figure()
            fig_scatter_t5.add_trace(go.Scatter(
                x=ages_t5, y=gene_vals_t5, mode='markers',
                marker=dict(color='mediumpurple', size=10, line=dict(width=0.5, color='black')),
                showlegend=False,
            ))
            slope, intercept = np.polyfit(ages_t5, gene_vals_t5, 1)
            x_line = np.array([ages_t5.min(), ages_t5.max()])
            fig_scatter_t5.add_trace(go.Scatter(
                x=x_line, y=slope * x_line + intercept, mode='lines',
                line=dict(color='firebrick', dash='dash'), showlegend=False,
            ))
            fig_scatter_t5.update_layout(
                title=(
                    f'GTEx Cortex: Expression vs Age — {selected_gene_t5}<br>'
                    f'<sup>Pearson r={row_t5["age_pearson_r"]:.3f}, p={row_t5["age_pearson_p"]:.3e}</sup>'
                ),
                xaxis_title='<b>Age (bin lower bound)</b>', yaxis_title='<b>TPM</b>',
                yaxis=dict(showgrid=False),
                plot_bgcolor='rgba(0,0,0,0)', font=dict(family='Arial', size=18),
                height=600, width=900, margin=dict(r=20, b=80),
            )
            fig_scatter_t5.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            fig_scatter_t5.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
            st.plotly_chart(fig_scatter_t5, use_container_width=False, key='scatter_t5')
