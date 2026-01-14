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
    """计算 Bulk 级别的 log1p(mean counts)"""
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
    """安全计算相关性"""
    if a is None or b is None: return np.nan
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if a.size < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]): return np.nan
    return spearmanr(a, b)[0] if method == "spearman" else pearsonr(a, b)[0]

def plot_heatmap(mat, row_labels, col_labels, title, out_png, out_pdf, vmin=-1, vmax=1):
    if mat.size == 0 or np.all(np.isnan(mat)):
        log(f"跳过绘图 {title}: 无有效数据点。")
        return
    
    mat = np.asarray(mat, dtype=float)
    nrows, ncols = mat.shape
    fig_w, fig_h = max(6, 0.55 * ncols + 3), max(4, 0.45 * nrows + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap="coolwarm", aspect="auto")
    
    ax.set_xticks(np.arange(ncols))
    ax.set_yticks(np.arange(nrows))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticklabels(row_labels)
    ax.set_title(title)

    for i in range(nrows):
        for j in range(ncols):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "NA", ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Correlation")
    plt.tight_layout()
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)
    log(f"Saved: {out_png}")

# -------------------------- Main --------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pred-pert", required=True, help="预测的扰动后文件 (pred.h5ad)")
    parser.add_argument("--pred-unpert", required=True, help="预测的空扰动文件 (pred_unpert.h5ad)")
    parser.add_argument("--base-h5ad", required=True, help="原始 939 基因基础文件")
    parser.add_argument("--tag", required=True, help="输出标签前缀")
    parser.add_argument("--real-samples", required=True, help="逗号分隔的真实样本ID")
    parser.add_argument("--ctrl-samples", required=True, help="逗号分隔的对照样本ID (用于仿真基线)")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log(f"Loading data for Split Analysis: {args.tag}...")
    ad_p = sc.read_h5ad(args.pred_pert)
    ad_u = sc.read_h5ad(args.pred_unpert)
    ad_b = sc.read_h5ad(args.base_h5ad)

    # 定义要分析的细胞类型
    analyze_types = ['Stroma', 'Tumor']

    # 解析样本列表
    ctrl_samples = [s.strip() for s in args.ctrl_samples.split(",")]
    real_samples = [s.strip() for s in args.real_samples.split(",")]

    log(f"Control Samples ({len(ctrl_samples)}): {ctrl_samples}")
    log(f"Real Samples ({len(real_samples)}): {real_samples}")

    # ------------------ 循环分析 Stroma 和 Tumor ------------------
    for cell_type in analyze_types:
        log(f"\n>>> Processing Sub-population: {cell_type} <<<")
        
        # 1. 筛选特定细胞类型
        def get_mask(ad, target_type):
            mask = np.zeros(ad.n_obs, dtype=bool)
            for col in ad.obs.columns:
                if ad.obs[col].dtype.name == 'category' or ad.obs[col].dtype == object:
                    try:
                        # 严格匹配（包含 target_type 字符串）
                        series_str = ad.obs[col].astype(str).str.strip()
                        col_mask = series_str.str.contains(target_type, case=False, regex=True)
                        if col_mask.any():
                            mask |= col_mask
                    except Exception:
                        continue
            return mask

        mask_p = get_mask(ad_p, cell_type)
        mask_u = get_mask(ad_u, cell_type)
        mask_b = get_mask(ad_b, cell_type)

        log(f"  Cells kept ({cell_type}): Pred={mask_p.sum()}, Base={mask_b.sum()}")

        if mask_p.sum() == 0 or mask_b.sum() == 0:
            log(f"  [WARN] No {cell_type} cells found. Skipping.")
            continue

        # 切片
        ad_p_s = ad_p[mask_p].copy()
        ad_u_s = ad_u[mask_u].copy()
        ad_b_s = ad_b[mask_b].copy()

        # 2. 基因对齐 (以 base 为准)
        base_genes = ad_b_s.var_names.tolist()
        def align(ad_source):
            pred_map = {g: i for i, g in enumerate(ad_source.var_names)}
            aligned_X = np.zeros((ad_source.n_obs, len(base_genes)))
            raw_X = ad_source.X.toarray() if issparse(ad_source.X) else ad_source.X
            for j, g in enumerate(base_genes):
                if g in pred_map: 
                    aligned_X[:, j] = raw_X[:, pred_map[g]]
            return sc.AnnData(X=aligned_X, obs=ad_source.obs, var=ad_b_s.var)

        ad_p_aligned = align(ad_p_s)
        ad_u_aligned = align(ad_u_s)

        # 3. 匹配样本
        # 确保 batch_var
        if 'batch_var' not in ad_b_s.obs.columns:
            if 'dataset_name' in ad_b_s.obs.columns:
                ad_b_s.obs['batch_var'] = ad_b_s.obs['dataset_name']
                ad_p_aligned.obs['batch_var'] = ad_p_aligned.obs['dataset_name']
                ad_u_aligned.obs['batch_var'] = ad_u_aligned.obs['dataset_name']
            else:
                log("  [ERROR] 'batch_var' not found.")
                continue

        available = set(ad_b_s.obs['batch_var'].unique())
        rows = [s for s in ctrl_samples if s in available]
        cols = [s for s in real_samples if s in available]

        if not rows or not cols:
            log("  [ERROR] No matching samples found in data.")
            continue

        # 4. 计算 Delta (均只针对当前 cell_type)
        # Real Delta = Real_Perturb(Type) - Real_Control_Baseline(Type)
        # 注意：这里的基线仅由 Control 样本中的该类型细胞构成
        data_baseline = get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'].isin(rows)])
        
        real_deltas = {}
        for s in cols:
            val = get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'] == s])
            if val is not None and data_baseline is not None:
                real_deltas[s] = val - data_baseline
            else:
                real_deltas[s] = None

        # Sim Delta = Pred_Perturb(Type) - Pred_Unpert(Type)
        sim_deltas = {}
        for s in rows:
            vec_p = get_bulk_log_mean_counts(ad_p_aligned[ad_p_aligned.obs['batch_var'] == s])
            vec_u = get_bulk_log_mean_counts(ad_u_aligned[ad_u_aligned.obs['batch_var'] == s])
            if vec_p is not None and vec_u is not None:
                sim_deltas[s] = vec_p - vec_u
            else:
                sim_deltas[s] = None

        # 5. 绘图
        for m in ["pearson", "spearman"]:
            # 构建矩阵
            mat = np.full((len(rows), len(cols)), np.nan)
            for i, r in enumerate(rows):
                for j, c in enumerate(cols):
                    if r in sim_deltas and c in real_deltas:
                        mat[i, j] = corr_safe(sim_deltas[r], real_deltas[c], m)
            
            # 文件名带上 cell_type
            fname_base = f"delta_corr_{args.tag}_{cell_type}_{m}"
            out_png = os.path.join(args.outdir, fname_base + ".png")
            out_pdf = os.path.join(args.outdir, fname_base + ".pdf")
            
            title = f"{args.tag} [{cell_type}] | Delta-vs-Delta | {m.capitalize()}"
            plot_heatmap(mat, rows, cols, title, out_png, out_pdf)

    log("\nAll Split Analyses Completed.")

if __name__ == "__main__":
    main()
