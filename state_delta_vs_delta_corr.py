#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import pearsonr, spearmanr

# -------------------------- Utils --------------------------
def log(msg: str):
    print(f"[INFO] {msg}", flush=True)

def get_bulk_log_mean_counts(adata):
    """计算 Bulk 级别的 log1p(mean counts): log1p(mean(expm1(X)))"""
    if adata.n_obs == 0: return None
    X = adata.X
    if issparse(X):
        Xc = X.copy()
        Xc.data = np.expm1(Xc.data)
        mean_counts = np.asarray(Xc.mean(axis=0)).ravel()
    else:
        mean_counts = np.expm1(np.asarray(X)).mean(axis=0).ravel()
    return np.log1p(mean_counts)

def corr_safe(a, b, method="pearson"):
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    if a.size < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]): return np.nan
    return spearmanr(a, b)[0] if method == "spearman" else pearsonr(a, b)[0]

def plot_heatmap(mat, row_labels, col_labels, title, out_png, out_pdf, vmin=-1, vmax=1):
    mat = np.asarray(mat, dtype=float)
    nrows, ncols = mat.shape
    fig_w, fig_h = max(6, 0.55 * ncols + 3), max(4, 0.45 * nrows + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="coolwarm", aspect="auto")
    ax.set_xticks(np.arange(ncols)); ax.set_yticks(np.arange(nrows))
    ax.set_xticklabels(col_labels, rotation=45, ha="right"); ax.set_yticklabels(row_labels)
    ax.set_title(title)
    for i in range(nrows):
        for j in range(ncols):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "NA", ha="center", va="center", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Correlation")
    plt.tight_layout()
    fig.savefig(out_png, dpi=220); fig.savefig(out_pdf); plt.close(fig)
    log(f"Saved: {out_png}")

# -------------------------- Main --------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred-pert", required=True, help="预测的扰动后文件 (pred.h5ad)")
    parser.add_argument("--pred-unpert", required=True, help="预测的空扰动文件 (pred_unpert.h5ad)")
    parser.add_argument("--base-h5ad", required=True, help="原始 939 基因基础文件")
    parser.add_argument("--tag", required=True)
    parser.add_argument("--real-samples", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log(f"Running Delta-vs-Delta analysis for {args.tag}...")
    ad_p = sc.read_h5ad(args.pred_pert)
    ad_u = sc.read_h5ad(args.pred_unpert)
    ad_b = sc.read_h5ad(args.base_h5ad)

    # 1. 筛选 Stroma 细胞 (全列扫描)
    def get_stroma_mask(ad):
        mask = np.zeros(ad.n_obs, dtype=bool)
        for col in ad.obs.columns:
            if ad.obs[col].astype(str).str.contains("Stroma", case=False).any():
                mask |= ad.obs[col].astype(str).str.contains("Stroma", case=False)
        return mask

    ad_p_s = ad_p[get_stroma_mask(ad_p)].copy()
    ad_u_s = ad_u[get_stroma_mask(ad_u)].copy()
    ad_b_s = ad_b[get_stroma_mask(ad_b)].copy()

    # 2. 基因对齐 (18k -> 939)
    base_genes = ad_b_s.var_names.tolist()
    def align(ad_source):
        pred_map = {g: i for i, g in enumerate(ad_source.var_names)}
        aligned_X = np.zeros((ad_source.n_obs, len(base_genes)))
        raw_X = ad_source.X.toarray() if issparse(ad_source.X) else ad_source.X
        for j, g in enumerate(base_genes):
            if g in pred_map: aligned_X[:, j] = raw_X[:, pred_map[g]]
        return sc.AnnData(X=aligned_X, obs=ad_source.obs, var=ad_b_s.var)

    ad_p_aligned = align(ad_p_s)
    ad_u_aligned = align(ad_u_s)

    # 3. 计算样本列表
    ctrl_samples = ["RCC3_1", "RCC3_2", "RCC3_4", "RCC3_6", "RCC3_9", "RCC3_10", "RCC5_4", "RCC5_6", "RCC5_16", "RCC5_18"]
    real_samples = [s.strip() for s in args.real_samples.split(",")]
    available = set(ad_b_s.obs['batch_var'].unique())
    rows, cols = [s for s in ctrl_samples if s in available], [s for s in real_samples if s in available]

    # 4. 计算 Delta
    # Real Delta = Real_Perturb - Real_Control_Baseline
    data_baseline = get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'].isin(rows)])
    real_deltas = {s: get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'] == s]) - data_baseline for s in cols}
    
    # Sim Delta = Pred_Perturb - Pred_Unperturbed
    # 这里我们按 batch 计算模拟的扰动效应
    sim_deltas = {}
    for s in rows:
        vec_p = get_bulk_log_mean_counts(ad_p_aligned[ad_p_aligned.obs['batch_var'] == s])
        vec_u = get_bulk_log_mean_counts(ad_u_aligned[ad_u_aligned.obs['batch_var'] == s])
        sim_deltas[s] = vec_p - vec_u

    # 5. 绘图
    for m in ["pearson", "spearman"]:
        mat = np.array([[corr_safe(sim_deltas[r], real_deltas[c], m) for c in cols] for r in rows])
        plot_heatmap(mat, rows, cols, f"STATE {args.tag} | Delta-vs-Delta | {m.capitalize()}", 
                     f"{args.outdir}/delta_corr_{args.tag}_{m}.png", f"{args.outdir}/delta_corr_{args.tag}_{m}.pdf")

if __name__ == "__main__":
    main()
