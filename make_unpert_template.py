#!/usr/bin/env python3
import scanpy as sc
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="现有的推理模板 (如 val_template_pd1ctla4_18k.h5ad)")
    parser.add_argument("--output", required=True, help="生成的空扰动模板")
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    
    # 将所有扰动标签重置为对照组标签
    print("Setting all target_gene to 'non-targeting'...")
    adata.obs['target_gene'] = "non-targeting"
    if 'guide_id' in adata.obs.columns:
        adata.obs['guide_id'] = "non-targeting"
    
    adata.write_h5ad(args.output, compression="gzip")
    print(f"Success! Saved to {args.output}")

if __name__ == "__main__":
    main()
