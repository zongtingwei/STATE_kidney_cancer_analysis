#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import pearsonr, spearmanr

# -------------------------- 样本列表硬编码 --------------------------
SARCOMOID_CONTROL = [
    "RCC3_7", "RCC3_8", "RCC3_11", "RCC3_12", "RCC3_13", "RCC3_14", "RCC3_15", "RCC3_16",
    "RCC5_7", "RCC5_8", "RCC5_9", "RCC5_10", "RCC5_11", "RCC5_12"
]

SARCOMOID_PD1CTLA4 = ["RCC4_3", "RCC4_4"]
SARCOMOID_PD1VEGFR = ["RCC4_5", "RCC4_6"]

# -------------------------- 工具函数 --------------------------
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
    """安全计算相关性"""
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if a.size < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]): return np.nan
    return spearmanr(a, b)[0] if method == "spearman" else pearsonr(a, b)[0]

def plot_heatmap(mat, row_labels, col_labels, title, out_png, out_pdf):
    """绘制热图并保存为 PNG 和 PDF"""
    if mat.size == 0 or np.all(np.isnan(mat)):
        log(f"跳过绘图 {title}: 无有效数据点。")
        return
    fig_w, fig_h = max(6, 0.55 * len(col_labels) + 3), max(4, 0.45 * len(row_labels) + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, vmin=-1, vmax=1, cmap="coolwarm", aspect="auto")
    
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticklabels(row_labels)
    ax.set_title(title)

    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "NA", ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Correlation")
    plt.tight_layout()
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)
    log(f"已保存热图: {out_png}")

# -------------------------- 核心逻辑 --------------------------
def main():
    parser = argparse.ArgumentParser(description="STATE Sarcomoid Delta-vs-Delta Analysis")
    parser.add_argument("--pred-pert", required=True, help="扰动后的预测文件 (pred.h5ad)")
    parser.add_argument("--pred-unpert", required=True, help="空扰动的预测文件 (pred_unpert.h5ad)")
    parser.add_argument("--base-h5ad", required=True, help="原始 939 基因基础文件")
    parser.add_argument("--task", choices=["ctla4", "vegfr"], required=True, help="分析任务")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log(f"开始 Sarcomoid 样本分析任务: {args.task.upper()}")
    
    # 加载数据
    ad_p = sc.read_h5ad(args.pred_pert)
    ad_u = sc.read_h5ad(args.pred_unpert)
    ad_b = sc.read_h5ad(args.base_h5ad)

    # 1. Stroma 细胞筛选
    mask = ad_b.obs.stack().str.contains("Stroma", case=False).unstack().any(axis=1)
    ad_b_s = ad_b[mask].copy()
    ad_p_s = ad_p[ad_p.obs_names.isin(ad_b_s.obs_names)].copy()
    ad_u_s = ad_u[ad_u.obs_names.isin(ad_b_s.obs_names)].copy()
    log(f"找到 Stroma 细胞数量: {ad_b_s.n_obs}")

    # 2. 基因空间还原 (18k -> 939)
    base_genes = ad_b_s.var_names.tolist()
    def align(source):
        pred_map = {g: i for i, g in enumerate(source.var_names)}
        aligned_X = np.zeros((source.n_obs, len(base_genes)))
        raw_X = source.X.toarray() if issparse(source.X) else source.X
        for j, g in enumerate(base_genes):
            if g in pred_map: aligned_X[:, j] = raw_X[:, pred_map[g]]
        return sc.AnnData(X=aligned_X, obs=source.obs, var=ad_b_s.var)

    ad_p_aligned = align(ad_p_s)
    ad_u_aligned = align(ad_u_s)

    # 3. 确定行与列
    available = set(ad_b_s.obs['batch_var'].unique())
    rows = [s for s in SARCOMOID_CONTROL if s in available]
    
    if args.task == "ctla4":
        cols = [s for s in SARCOMOID_PD1CTLA4 if s in available]
        tag = "Sarcomoid_PD1_CTLA4"
    else:
        cols = [s for s in SARCOMOID_PD1VEGFR if s in available]
        tag = "Sarcomoid_PD1_VEGFR"

    log(f"行 (模拟样本): {rows}")
    log(f"列 (真实样本): {cols}")

    if not rows or not cols:
        log("[ERROR] 样本匹配失败，请检查 batch_var 内容。")
        return

    # 4. 计算差异向量 (Delta)
    # 真实 Delta = 真实扰动样本 - Sarcomoid对照组全局均值
    data_baseline = get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'].isin(rows)])
    real_deltas = {s: get_bulk_log_mean_counts(ad_b_s[ad_b_s.obs['batch_var'] == s]) - data_baseline for s in cols}
    
    # 模拟 Delta = 预测扰动 - 预测空扰动 (针对每一个对照样本)
    sim_deltas = {s: get_bulk_log_mean_counts(ad_p_aligned[ad_p_aligned.obs['batch_var'] == s]) - 
                     get_bulk_log_mean_counts(ad_u_aligned[ad_u_aligned.obs['batch_var'] == s]) for s in rows}

    # 5. 绘图
    for m in ["pearson", "spearman"]:
        mat = np.array([[corr_safe(sim_deltas[r], real_deltas[c], m) for c in cols] for r in rows])
        out_p = os.path.join(args.outdir, f"sarcomoid_{tag}_{m}.png")
        out_f = os.path.join(args.outdir, f"sarcomoid_{tag}_{m}.pdf")
        plot_heatmap(mat, rows, cols, f"Sarcomoid {tag} | Delta-vs-Delta | {m.capitalize()}", out_p, out_f)

if __name__ == "__main__":
    main()
