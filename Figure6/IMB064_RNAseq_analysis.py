#!/usr/bin/env python3
"""
IMB064 brain RNA-seq — analysis code (Figure 1; Supplementary Tables S1-S3)
==========================================================================
Single-file, computation-only pipeline that reproduces every value reported
to readers in Figure 1 and Supplementary Tables S1-S3. No figure-drawing code.

Upstream step (input, not included): differential expression was performed
with DESeq2 on the raw counts. The per-region DESeq2 result tables (gene
symbols as index; columns include log2FoldChange, pvalue, padj) are the INPUT
to this script. A gene is a DEG when padj < 0.05 and |log2FoldChange| > 0.5.

Sub-commands
------------
  deg       DEG identification + counts per region        -> Fig. 1A
  reactome  Reactome over-representation (gProfiler g:SCS) -> Fig. 1B / Table S2
  overlap   Overlap with 5 ASD gene sets (Fisher's exact) -> Fig. 1C / Table S3
  modules   Gandal WGCNA module enrichment (Fisher + BH)  -> Fig. 1D / Table S1

Run `python IMB064_RNAseq_analysis.py <sub-command> --help` for arguments.

Data sources — publicly available supplementary files
------------------------------------------------------
Human->mouse gene-name conversion  (used by `overlap` and `modules`)
  * Ensembl BioMart one-to-one orthologs. Ensembl release 113, accessed
    June 2026, via BioMart (https://www.ensembl.org/biomart/martview).
    Dataset: Human genes (hsapiens_gene_ensembl). Attributes: "Gene name"
    (external_gene_name) and "Mouse gene name"
    (mmusculus_homolog_associated_gene_name), restricted to one-to-one
    orthologs (homology type ortholog_one2one). Saved as a two-column CSV
    with header `human_gene,mouse_gene` (17,052 ortholog pairs). Passed via
    --biomart.

ASD gene sets  (used by `overlap`)
  1. Satterstrom et al. 2020, Cell 180(3):568-584. DOI 10.1016/j.cell.2019.12.036.
     Supplement mmc2.xlsx, sheet "102_ASD" — 102 high-confidence ASD risk
     genes (human symbols -> mouse).
  2. Gandal et al. 2022, Nature 611:532-539. DOI 10.1038/s41586-022-05377-7.
     Supplementary Data 3, sheet "DEGene_Statistics" — human ASD cortex DEGs
     filtered to WholeCortex_ASD_FDR < 0.05 (human symbols -> mouse).
  3. Mooney et al. 2025, Molecular Neurobiology 62:10614-10634.
     DOI 10.1007/s12035-025-04900-x. Supplement, sheet "Table 2_DEGs_BTBR" —
     BTBR prefrontal-cortex DEGs (mouse symbols; no conversion).
  4-5. Fazel Darbandi et al. 2024, "Five autism-associated transcriptional
     regulators target shared loci proximal to brain-expressed genes,"
     Cell Reports (2024), article 114329. DOI 10.1016/j.celrep.2024.114329
     (PMC11235582).
       * CRISPRi targets  -> supplement mmc8.xlsx, sheets "Table S7F" (Tbr1)
         and "Table S7G" (Arid1b); union of genes with padj < 0.05 in either
         CRISPRi RNA-seq (mouse gene_name; no conversion).
       * Enhancer-proximal -> supplement mmc4.xlsx, sheet "Human_5TRa_Wide";
         protein-coding genes (gene_list field, Gene_type == protein_coding)
         near DISTAL peaks (prox_distal == 'distal'); human symbols -> mouse.

Gandal WGCNA module assignments  (used by `modules`)
  * Gandal et al. 2022 (as above), Supplementary Data 5, sheet "Gene_Level"
    (column WGCNA_module). The unassigned "grey" module ('M0_grey') is
    excluded from testing and from the BH multiple-testing correction.

Reactome  (used by `reactome`)
  * Queried live via the gProfiler web service (g:GOSt, g:SCS correction,
    sources=['REAC']); no local file.

Environment: Python 3.x with pandas, numpy, scipy, statsmodels, openpyxl,
gprofiler-official. The `reactome` sub-command needs network access to the
gProfiler web service; the others run offline.
"""

import argparse
import json
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def load_deg_table(deg_file, padj_cutoff=0.05, lfc_cutoff=0.5):
    """Return (full DataFrame, significant-DEG DataFrame) from a DESeq2 table."""
    df = pd.read_csv(deg_file, sep='\t', index_col=0)
    sig = df[(df['padj'] < padj_cutoff) & (df['log2FoldChange'].abs() > lfc_cutoff)].copy()
    return df, sig


def load_biomart(biomart_file):
    """Ensembl BioMart one-to-one human->mouse ortholog dict."""
    df = pd.read_csv(biomart_file)
    mapping = dict(zip(df['human_gene'], df['mouse_gene']))
    print(f"Loaded {len(mapping)} one-to-one orthologs from BioMart")
    return mapping


def to_mouse(gene, mapping):
    return mapping.get(gene) if isinstance(gene, str) and gene else None


# ---------------------------------------------------------------------------
# 1. DEG identification and counts (Fig. 1A)
# ---------------------------------------------------------------------------
def cmd_deg(args):
    os.makedirs(args.output_dir, exist_ok=True)
    regions = {'cortex': args.cortex, 'hippocampus': args.hippo, 'amygdala': args.amygdala}
    summary = []
    for name, path in regions.items():
        full, sig = load_deg_table(path, args.padj_cutoff, args.lfc_cutoff)
        sig = sig.copy()
        sig['direction'] = ['up' if x > 0 else 'down' for x in sig['log2FoldChange']]
        n_up = int((sig['direction'] == 'up').sum())
        n_down = int((sig['direction'] == 'down').sum())
        summary.append({'region': name, 'background_genes': len(full),
                        'DEGs_total': len(sig), 'up': n_up, 'down': n_down})
        sig.to_csv(os.path.join(args.output_dir, f'{name}_DEGs.csv'))
        print(f"  {name}: {n_up} up, {n_down} down ({len(sig)} total DEGs)")
    pd.DataFrame(summary).to_csv(os.path.join(args.output_dir, 'deg_counts_summary.csv'), index=False)
    print(f"\nSaved DEG counts + per-region DEG tables to {args.output_dir}/")


# ---------------------------------------------------------------------------
# 2. Reactome over-representation (Fig. 1B / Table S2)
# ---------------------------------------------------------------------------
def cmd_reactome(args):
    from gprofiler import GProfiler  # lazy import (needs network)
    os.makedirs(args.output_dir, exist_ok=True)
    _, sig = load_deg_table(args.deg_file, args.padj_cutoff, args.lfc_cutoff)
    up = sig[sig['log2FoldChange'] > 0]
    down = sig[sig['log2FoldChange'] < 0]
    print(f"DEGs: {len(sig)} total, {len(up)} up, {len(down)} down")

    gp = GProfiler(return_dataframe=True)

    def run(gene_list, label):
        res = gp.profile(organism=args.organism, query=list(gene_list),
                         sources=['REAC'], significance_threshold_method='g_SCS',
                         user_threshold=0.05, no_evidences=False)
        res['direction'] = label
        print(f"  {label}: {len(res)} significant Reactome pathways")
        return res

    print("\nRunning Reactome over-representation analysis (gProfiler g:SCS)...")
    run(sig.index.tolist(), 'All').to_csv(os.path.join(args.output_dir, 'gprofiler_reactome_all.csv'), index=False)
    run(up.index.tolist(), 'Up').to_csv(os.path.join(args.output_dir, 'gprofiler_reactome_up.csv'), index=False)
    run(down.index.tolist(), 'Down').to_csv(os.path.join(args.output_dir, 'gprofiler_reactome_down.csv'), index=False)
    print(f"\nSaved Reactome results (Table S2) to {args.output_dir}/gprofiler_reactome_*.csv")


# ---------------------------------------------------------------------------
# 3. ASD gene-set overlap (Fig. 1C / Table S3)
# ---------------------------------------------------------------------------
def fisher_overlap(degs, geneset, background, name):
    """One-sided Fisher's exact test (enrichment) + Woolf 95% CI."""
    gs = geneset & background
    a = len(degs & gs)
    b = len(degs - gs)
    c = len(gs - degs)
    d = len(background - degs - gs)
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
    a2, b2, c2, d2 = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log(a2 * d2 / (b2 * c2))
    se = np.sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
    result = {'geneset': name, 'overlap': a, 'geneset_in_bg': len(gs),
              'odds_ratio': round(odds_ratio, 2),
              'ci_low': round(np.exp(log_or - 1.96 * se), 2),
              'ci_high': round(np.exp(log_or + 1.96 * se), 2),
              'p_value': p_value}
    print(f"  {name}: {a}/{len(gs)}, OR={odds_ratio:.2f}, p={p_value:.2e}")
    return result


def cmd_overlap(args):
    os.makedirs(args.output_dir, exist_ok=True)
    full, sig = load_deg_table(args.deg_file, args.padj_cutoff, args.lfc_cutoff)
    background, degs = set(full.index), set(sig.index)
    print(f"Background: {len(background)} genes, DEGs: {len(degs)}")
    mapping = load_biomart(args.biomart)
    results = []

    if args.satterstrom:
        sat = pd.read_excel(args.satterstrom, sheet_name='102_ASD')
        genes = set(filter(None, [to_mouse(str(g), mapping) for g in sat.iloc[:, 0].dropna() if isinstance(g, str)]))
        results.append(fisher_overlap(degs, genes, background, 'Satterstrom 2020 — 102 ASD risk genes'))

    if args.gandal:
        gandal = pd.read_excel(args.gandal, sheet_name='DEGene_Statistics')
        gandal = gandal[gandal['WholeCortex_ASD_FDR'] < 0.05]
        genes = set(filter(None, [to_mouse(str(g), mapping) for g in gandal['external_gene_name'].dropna() if isinstance(g, str)]))
        results.append(fisher_overlap(degs, genes, background, 'Gandal 2022 — Human ASD cortex DEGs'))

    if args.mooney:
        mooney = pd.read_excel(args.mooney, sheet_name='Table 2_DEGs_BTBR')
        results.append(fisher_overlap(degs, set(mooney['GeneSymbol'].dropna()), background,
                                      'Mooney 2025 — BTBR cortex DEGs'))

    if args.fd_mmc4 and args.fd_mmc8:
        # CRISPRi union: padj < 0.05 in Tbr1 (S7F) or Arid1b (S7G); mouse gene_name
        s7f = pd.read_excel(args.fd_mmc8, sheet_name='Table S7F')
        s7g = pd.read_excel(args.fd_mmc8, sheet_name='Table S7G')
        crispri = (set(s7f[s7f['padj'] < 0.05]['gene_name'].dropna())
                   | set(s7g[s7g['padj'] < 0.05]['gene_name'].dropna()))
        results.append(fisher_overlap(degs, crispri, background,
                                      'Fazel Darbandi 2024 — CRISPRi targets (Tbr1/Arid1b)'))
        # Enhancer-proximal: protein-coding genes near DISTAL Human_5TRa_Wide peaks -> mouse
        tra = pd.read_excel(args.fd_mmc4, sheet_name='Human_5TRa_Wide')
        distal = tra[tra['prox_distal'] == 'distal']
        enh_human = set()
        for gl in distal['gene_list'].dropna():
            for part in str(gl).split(','):
                fields = part.split(':')  # Sym:Ensembl:Strand:Gene_type:Distance
                if len(fields) >= 5 and fields[3] == 'protein_coding' and fields[0].strip():
                    enh_human.add(fields[0].strip())
        enh = set(filter(None, [to_mouse(g, mapping) for g in enh_human]))
        results.append(fisher_overlap(degs, enh, background,
                                      'Fazel Darbandi 2024 — Enhancer-proximal genes'))

    pd.DataFrame(results).to_csv(os.path.join(args.output_dir, 'asd_geneset_overlap_results.csv'), index=False)
    with open(os.path.join(args.output_dir, 'asd_geneset_overlap_results.json'), 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved overlap results (Table S3) to {args.output_dir}/asd_geneset_overlap_results.csv")


# ---------------------------------------------------------------------------
# 4. WGCNA module enrichment (Fig. 1D / Table S1)
# ---------------------------------------------------------------------------
def load_modules(module_file, mapping):
    """Gandal Supp Data 5 module assignments; exclude the 'M0_grey' unassigned module."""
    df = pd.read_excel(module_file, sheet_name='Gene_Level')
    df['mouse_gene'] = df['external_gene_name'].apply(lambda g: mapping.get(g) if isinstance(g, str) else None)
    modules = {}
    for mod in df['WGCNA_module'].dropna().unique():
        if str(mod) == 'M0_grey' or str(mod).startswith('M0_') or mod == 'M0':
            continue  # unassigned grey module — not a real co-expression module
        modules[mod] = set(df[df['WGCNA_module'] == mod]['mouse_gene'].dropna())
    print(f"Loaded {len(modules)} WGCNA modules (excluding M0/unassigned)")
    return modules


def cmd_modules(args):
    os.makedirs(args.output_dir, exist_ok=True)
    full, sig = load_deg_table(args.deg_file, args.padj_cutoff, args.lfc_cutoff)
    background, degs = set(full.index), set(sig.index)
    print(f"Background: {len(background)} genes, DEGs: {len(degs)}")
    mapping = load_biomart(args.biomart)
    modules = load_modules(args.gandal_modules, mapping)

    rows = []
    for mod_name, mod_genes in sorted(modules.items()):
        mod_in_bg = mod_genes & background
        a = len(degs & mod_in_bg)
        b = len(degs - mod_in_bg)
        c = len(mod_in_bg - degs)
        d = len(background - degs - mod_in_bg)
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append({'module': mod_name, 'module_size': len(mod_in_bg),
                     'overlap': a, 'odds_ratio': round(odds_ratio, 2), 'p_value': p_value})
    res = pd.DataFrame(rows)
    reject, fdr, _, _ = multipletests(res['p_value'], method='fdr_bh')
    res['fdr'] = fdr
    res['significant'] = reject
    print(f"Tested {len(res)} modules, {int(res['significant'].sum())} significant (FDR < 0.05)")

    out = os.path.join(args.output_dir, 'module_enrichment_results.csv')
    res.to_csv(out, index=False)
    for _, r in res[res['significant']].sort_values('odds_ratio', ascending=False).iterrows():
        print(f"  {r['module']:20s} OR={r['odds_ratio']:5.2f}  {r['overlap']}/{r['module_size']}  "
              f"p={r['p_value']:.2e}  FDR={r['fdr']:.2e}")
    print(f"\nSaved module enrichment results (Table S1) to {out}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description='IMB064 brain RNA-seq analysis (Figure 1; Tables S1-S3).',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)

    d = sub.add_parser('deg', help='DEG identification + counts (Fig. 1A)')
    d.add_argument('--cortex', required=True)
    d.add_argument('--hippo', required=True)
    d.add_argument('--amygdala', required=True)
    d.add_argument('--output_dir', required=True)
    d.add_argument('--padj_cutoff', type=float, default=0.05)
    d.add_argument('--lfc_cutoff', type=float, default=0.5)
    d.set_defaults(func=cmd_deg)

    r = sub.add_parser('reactome', help='Reactome over-representation (Fig. 1B / Table S2)')
    r.add_argument('--deg_file', required=True, help='cortex DESeq2 table')
    r.add_argument('--output_dir', required=True)
    r.add_argument('--organism', default='mmusculus')
    r.add_argument('--padj_cutoff', type=float, default=0.05)
    r.add_argument('--lfc_cutoff', type=float, default=0.5)
    r.set_defaults(func=cmd_reactome)

    o = sub.add_parser('overlap', help='ASD gene-set overlap (Fig. 1C / Table S3)')
    o.add_argument('--deg_file', required=True, help='cortex DESeq2 table')
    o.add_argument('--biomart', required=True, help='BioMart one-to-one CSV (human_gene,mouse_gene)')
    o.add_argument('--satterstrom', help='Satterstrom 2020 mmc2.xlsx')
    o.add_argument('--gandal', help='Gandal 2022 Supp Data 3 (DEGene_Statistics)')
    o.add_argument('--mooney', help='Mooney 2025 Supp Table 3')
    o.add_argument('--fd_mmc4', help='Fazel Darbandi 2024 mmc4.xlsx (5TRa enhancer loci)')
    o.add_argument('--fd_mmc8', help='Fazel Darbandi 2024 mmc8.xlsx (CRISPRi DEG, Table S7F/S7G)')
    o.add_argument('--output_dir', required=True)
    o.add_argument('--padj_cutoff', type=float, default=0.05)
    o.add_argument('--lfc_cutoff', type=float, default=0.5)
    o.set_defaults(func=cmd_overlap)

    m = sub.add_parser('modules', help='WGCNA module enrichment (Fig. 1D / Table S1)')
    m.add_argument('--deg_file', required=True, help='cortex DESeq2 table')
    m.add_argument('--biomart', required=True, help='BioMart one-to-one CSV (human_gene,mouse_gene)')
    m.add_argument('--gandal_modules', required=True, help='Gandal Supp Data 5 (Gene_Level sheet)')
    m.add_argument('--output_dir', required=True)
    m.add_argument('--padj_cutoff', type=float, default=0.05)
    m.add_argument('--lfc_cutoff', type=float, default=0.5)
    m.set_defaults(func=cmd_modules)
    return p


def main():
    args = build_parser().parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
