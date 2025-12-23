#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import anndata as ad
import sys

def load_adata_safely(file_path):
    """安全读取 h5ad，跳过可能导致报错的 uns"""
    try:
        return sc.read_h5ad(file_path)
    except Exception:
        print(f"[INFO] 标准读取失败，尝试安全提取核心数据: {file_path}")
        adata_backed = ad.read_h5ad(file_path, backed='r')
        return ad.AnnData(
            X=adata_backed.X[:], 
            obs=adata_backed.obs.copy(), 
            var=adata_backed.var.copy()
        )

def main():
    ap = argparse.ArgumentParser(description="将 h5ad 基因维度对齐到模型要求的 18080 维")
    ap.add_argument("--in-h5ad", required=True, help="输入的 939 基因 h5ad")
    ap.add_argument("--out-h5ad", required=True, help="输出的 18080 基因 h5ad")
    ap.add_argument("--gene-list", required=True, help="模型基因列表 gene_names.csv")
    args = ap.parse_args()

    # 1. 加载目标基因列表 (18080 维)
    target_genes = pd.read_csv(args.gene_list, header=None)[0].astype(str).tolist()
    n_target = len(target_genes)
    print(f"[INFO] 目标基因数: {n_target}")

    # 2. 加载原始数据
    adata = load_adata_safely(args.in_h5ad)
    print(f"[INFO] 原始数据维度: {adata.shape}")

    # 3. 准备索引映射
    # 找出哪些目标基因在原始数据中是存在的
    existing_genes = set(adata.var_names)
    intersect_genes = [g for g in target_genes if g in existing_genes]
    print(f"[INFO] 共有 {len(intersect_genes)} 个基因在 CosMx 面板中命中。")

    # 4. 创建一个全零的稀疏矩阵 (n_obs, 18080)
    new_X = sp.lil_matrix((adata.n_obs, n_target), dtype=adata.X.dtype)

    # 5. 填充数据
    # 为了效率，我们按列进行填充
    for i, gene in enumerate(target_genes):
        if gene in existing_genes:
            # 获取该基因在原始 adata 中的索引
            old_idx = adata.var_names.get_loc(gene)
            # 这里的 adata.X[:, old_idx] 需要处理成 1D 数组
            col_data = adata.X[:, old_idx]
            if sp.issparse(col_data):
                col_data = col_data.toarray()
            new_X[:, i] = col_data.reshape(-1, 1)

    # 转为 CSR 格式以提高性能
    new_X = new_X.tocsr()

    # 6. 构建新的 AnnData
    new_var = pd.DataFrame(index=target_genes)
    new_var['gene_name'] = target_genes
    
    new_adata = ad.AnnData(
        X=new_X,
        obs=adata.obs.copy(),
        var=new_var,
        uns=adata.uns.copy() if hasattr(adata, 'uns') else {}
    )

    # 7. 保存结果
    new_adata.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"[OK] 已成功对齐并保存至: {args.out_h5ad}")

if __name__ == "__main__":
    main()
