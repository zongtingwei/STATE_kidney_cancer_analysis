#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import scanpy as sc
import sys
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Create PD1+CTLA4 and PD1+VEGFR designs (Stroma + Tumor) based on strict celltype matching.")
    ap.add_argument("--base-h5ad", required=True, help="Path to base h5ad file")
    ap.add_argument("--out-pd1ctla4", required=True, help="Output path for PD1+CTLA4 design")
    ap.add_argument("--out-pd1vegfr", required=True, help="Output path for PD1+VEGFR design")
    ap.add_argument("--control-label", default="non-targeting")
    # 注意：请确认你的 h5ad 里细胞类型列名是 'celltype' 还是 'cell_type'
    # 根据你提供的 VirtualTissue 代码片段，这里默认设为 'celltype'
    ap.add_argument("--celltype-col", default="celltype", help="Column for cell type matching")
    ap.add_argument("--annotation-col", default="annotation", help="Column to filter Stroma/Tumor")
    args = ap.parse_args()

    print(f"Loading {args.base_h5ad}...")
    adata = sc.read_h5ad(args.base_h5ad)
    
    # ---------------------------------------------------------
    # 1. 筛选 Stroma + Tumor
    # ---------------------------------------------------------
    if args.annotation_col in adata.obs.columns:
        print(f"Filtering column '{args.annotation_col}' for ['Stroma', 'Tumor']...")
        
        # 去除首尾空格
        adata.obs[args.annotation_col] = adata.obs[args.annotation_col].astype(str).str.strip()
        
        # 筛选
        mask = adata.obs[args.annotation_col].isin(['Stroma', 'Tumor'])
        n_before = adata.n_obs
        adata = adata[mask].copy()
        n_after = adata.n_obs
        
        print(f"Cells kept: {n_after} (dropped {n_before - n_after})")
        
        if n_after == 0:
            print("[ERROR] No cells left after filtering!")
            sys.exit(1)
    else:
        print(f"[WARN] Column '{args.annotation_col}' not found! Using ALL cells.")

    # ---------------------------------------------------------
    # 2. 定义严格的细胞匹配逻辑 (这部分是你日志里缺失的)
    # ---------------------------------------------------------
    if args.celltype_col not in adata.obs.columns:
        print(f"[ERROR] Celltype column '{args.celltype_col}' not found in adata.obs!")
        print(f"Available columns: {adata.obs.columns.tolist()}")
        sys.exit(1)

    print(f"Using celltype column: '{args.celltype_col}'")
    
    # 获取 series 以便向量化操作
    celltypes = adata.obs[args.celltype_col].astype(str)

    # === 核心修改：严格匹配逻辑 ===
    # 1. CD8 T cell (Exact match)
    mask_cd8 = (celltypes == 'CD8 T cell')

    # 2. Regulatory T cell (Exact match)
    mask_treg = (celltypes == 'Regulatory T cell')

    # 3. Endothelium (Substring match, 'endothelium' in ct)
    # case=False 忽略大小写, regex=False 纯字符串匹配
    mask_endo = celltypes.str.contains('endothelium', case=False, regex=False)

    # === 打印统计 (这就是你刚才没看到的部分) ===
    print(f"  Found {mask_cd8.sum()} 'CD8 T cell'")
    print(f"  Found {mask_treg.sum()} 'Regulatory T cell'")
    print(f"  Found {mask_endo.sum()} cells containing 'endothelium'")

    # ---------------------------------------------------------
    # 3. Generate PD1 + CTLA4 Design
    # ---------------------------------------------------------
    print("\nGenerating PD1 + CTLA4 design...")
    ad1 = adata.copy()
    ad1.obs["target_gene"] = args.control_label
    
    # 批量赋值 (比 for 循环快得多)
    if mask_cd8.any():
        ad1.obs.loc[mask_cd8, "target_gene"] = "PDCD1"
    if mask_treg.any():
        ad1.obs.loc[mask_treg, "target_gene"] = "CTLA4"

    # 打印统计以确认
    print("Design 1 Summary:")
    print(ad1.obs["target_gene"].value_counts())
    
    ad1.write_h5ad(args.out_pd1ctla4, compression="gzip")
    print(f"[OK] wrote design: {args.out_pd1ctla4}")

    # ---------------------------------------------------------
    # 4. Generate PD1 + VEGFR (KDR) Design
    # ---------------------------------------------------------
    print("\nGenerating PD1 + VEGFR design...")
    ad2 = adata.copy()
    ad2.obs["target_gene"] = args.control_label
    
    # 批量赋值
    if mask_cd8.any():
        ad2.obs.loc[mask_cd8, "target_gene"] = "PDCD1"
    if mask_endo.any():
        ad2.obs.loc[mask_endo, "target_gene"] = "KDR"
    
    # 打印统计以确认
    print("Design 2 Summary:")
    print(ad2.obs["target_gene"].value_counts())

    ad2.write_h5ad(args.out_pd1vegfr, compression="gzip")
    print(f"[OK] wrote design: {args.out_pd1vegfr}")

if __name__ == "__main__":
    main()
