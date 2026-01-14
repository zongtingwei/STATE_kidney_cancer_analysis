#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import scanpy as sc
import sys

def match_any(s: str, patterns):
    s_low = str(s).lower()
    return any(re.search(p, s_low) for p in patterns)

def main():
    ap = argparse.ArgumentParser(description="Create PD1+CTLA4 and PD1+VEGFR designs (Stroma + Tumor).")
    ap.add_argument("--base-h5ad", required=True, help="Path to base h5ad file")
    ap.add_argument("--out-pd1ctla4", required=True, help="Output path for PD1+CTLA4 design")
    ap.add_argument("--out-pd1vegfr", required=True, help="Output path for PD1+VEGFR design")
    ap.add_argument("--control-label", default="non-targeting")
    ap.add_argument("--celltype-col", default="cell_type", help="Column for cell type matching")
    # 新增 annotation-col 参数，默认 annotation
    ap.add_argument("--annotation-col", default="annotation", help="Column to filter Stroma/Tumor")
    args = ap.parse_args()

    print(f"Loading {args.base_h5ad}...")
    adata = sc.read_h5ad(args.base_h5ad)
    
    # ---------------------------------------------------------
    # 关键修改：筛选 Stroma + Tumor
    # ---------------------------------------------------------
    if args.annotation_col in adata.obs.columns:
        print(f"Filtering column '{args.annotation_col}' for ['Stroma', 'Tumor']...")
        
        # 1. 去除首尾空格 (解决 'Tumor ' 问题)
        adata.obs[args.annotation_col] = adata.obs[args.annotation_col].astype(str).str.strip()
        
        # 2. 筛选
        mask = adata.obs[args.annotation_col].isin(['Stroma', 'Tumor'])
        n_before = adata.n_obs
        adata = adata[mask].copy()
        n_after = adata.n_obs
        
        print(f"Cells kept: {n_after} (dropped {n_before - n_after})")
        
        if n_after == 0:
            print("[ERROR] No cells left after filtering! Check your annotation column values.")
            sys.exit(1)
    else:
        print(f"[WARN] Column '{args.annotation_col}' not found! Using ALL cells.")

    # ---------------------------------------------------------
    # 提取细胞类型 (注意：必须在筛选后提取，保证长度一致)
    # ---------------------------------------------------------
    ct = adata.obs[args.celltype_col].astype(str).values

    # 定义正则匹配模式
    # - PDCD1 perturb on CD8 T cell
    # - CTLA4 perturb on Treg
    # - KDR perturb on endothelium
    cd8_pat  = [r"\bcd8\b", r"cd8 t", r"cd8\+? t"]
    treg_pat = [r"\btreg\b", r"regulatory t", r"treg cell"]
    endo_pat = [r"endotheli", r"endothelial"]

    # -------- 1. Generate PD1 + CTLA4 Design --------
    print("Generating PD1 + CTLA4 design...")
    ad1 = adata.copy()
    ad1.obs["target_gene"] = args.control_label
    
    # 获取列索引以便快速修改
    # (如果数据量巨大，建议改用 boolean masking，但维持原逻辑对于 <200k 细胞通常没问题)
    target_idx = ad1.obs.columns.get_loc("target_gene")

    for i, cti in enumerate(ct):
        if match_any(cti, cd8_pat):
            ad1.obs.iloc[i, target_idx] = "PDCD1"
        elif match_any(cti, treg_pat):
            ad1.obs.iloc[i, target_idx] = "CTLA4"

    ad1.write_h5ad(args.out_pd1ctla4, compression="gzip")
    print(f"[OK] wrote design: {args.out_pd1ctla4}")

    # -------- 2. Generate PD1 + VEGFR (KDR) Design --------
    print("Generating PD1 + VEGFR design...")
    ad2 = adata.copy()
    ad2.obs["target_gene"] = args.control_label
    
    target_idx2 = ad2.obs.columns.get_loc("target_gene")

    for i, cti in enumerate(ct):
        if match_any(cti, cd8_pat):
            ad2.obs.iloc[i, target_idx2] = "PDCD1"
        elif match_any(cti, endo_pat):
            ad2.obs.iloc[i, target_idx2] = "KDR"

    ad2.write_h5ad(args.out_pd1vegfr, compression="gzip")
    print(f"[OK] wrote design: {args.out_pd1vegfr}")

if __name__ == "__main__":
    main()
