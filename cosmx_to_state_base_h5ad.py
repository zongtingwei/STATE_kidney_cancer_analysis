import argparse
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import anndata as ad
import sys

def pick_first_existing(obs, candidates):
    """从 obs 列名中挑选第一个存在的列"""
    for c in candidates:
        if c in obs.columns:
            return c
    return None

def load_adata_safely(file_path):
    """
    安全读取 h5ad 文件。
    如果遇到 encoding_type='null' 等版本不兼容错误，
    则尝试只读取核心数据 (X, obs, var)，跳过 uns。
    """
    try:
        print(f"[INFO] 正在尝试标准读取: {file_path}")
        return sc.read_h5ad(file_path)
    except Exception as e:
        print(f"[WARNING] 标准读取失败，报错信息: {e}")
        print(f"[INFO] 尝试绕过 uns 元数据进行读取...")
        
        try:
            # 这种方式会尝试以“后端映射”方式打开，通常能绕过一些解析错误
            adata_backed = ad.read_h5ad(file_path, backed='r')
            # 提取核心组件创建新对象
            new_adata = ad.AnnData(
                X=adata_backed.X[:], 
                obs=adata_backed.obs.copy(), 
                var=adata_backed.var.copy()
            )
            print("[SUCCESS] 成功提取核心数据 (X, obs, var)，已丢弃可能导致错误的 uns 字典。")
            return new_adata
        except Exception as e2:
            print(f"[ERROR] 即使尝试安全读取也失败了: {e2}")
            sys.exit(1)

def main():
    ap = argparse.ArgumentParser(
        description="Convert CosMx h5ad to STATE/cell-load friendly h5ad WITHOUT gene alignment."
    )
    ap.add_argument("--in-h5ad", required=True, help="输入路径")
    ap.add_argument("--out-h5ad", required=True, help="输出路径")

    ap.add_argument("--celltype-col", default=None, help="细胞类型所在的 obs 列名")
    ap.add_argument("--sample-col", default=None, help="样本/批次所在的 obs 列名")

    ap.add_argument("--control-label", default="non-targeting")
    ap.add_argument("--force-csr", action="store_true", help="强制将矩阵转为 CSR 格式")
    ap.add_argument("--add-gene-name-col", action="store_true", help="在 var 中添加 gene_name 列")
    args = ap.parse_args()

    # 使用安全读取函数
    adata = load_adata_safely(args.in_h5ad)

    # ----------- 选择细胞类型和批次列 -----------
    celltype_col = args.celltype_col or pick_first_existing(
        adata.obs,
        ["celltype", "cell_type", "CellType", "cell_type_final", "annotation", "celltype_major"]
    )
    if celltype_col is None:
        raise ValueError(
            "无法自动找到细胞类型列。请使用 --celltype-col 指定，例如 --celltype-col celltype"
        )

    sample_col = args.sample_col or pick_first_existing(
        adata.obs,
        ["dataset_name", "batch_var", "sample", "sample_id", "orig.ident", "Run_name", "fov"]
    )

    # ----------- 标准化 STATE 要求的字段 -----------
    print(f"[INFO] 正在标准化列: cell_type={celltype_col}, batch_var={sample_col}")
    
    adata.obs["cell_type"] = adata.obs[celltype_col].astype(str)

    if sample_col is None:
        adata.obs["batch_var"] = "batch0"
    else:
        adata.obs["batch_var"] = adata.obs[sample_col].astype(str)

    # 兼容性字段
    if "batch" not in adata.obs.columns:
        adata.obs["batch"] = adata.obs["batch_var"].astype(str)
    
    # 扰动标签 (Perturbation label)
    if "target_gene" not in adata.obs.columns:
        adata.obs["target_gene"] = args.control_label
    else:
        adata.obs["target_gene"] = adata.obs["target_gene"].astype(str).fillna(args.control_label)

    # guide_id 兼容性
    if "guide_id" not in adata.obs.columns:
        adata.obs["guide_id"] = adata.obs["target_gene"].astype(str)

    # ----------- 矩阵格式转换 -----------
    if args.force_csr:
        print("[INFO] 正在转换矩阵为 CSR 格式...")
        if sp.issparse(adata.X):
            adata.X = adata.X.tocsr()
        else:
            adata.X = sp.csr_matrix(adata.X)

    # ----------- 基因名安全处理 -----------
    adata.var_names_make_unique()
    if args.add_gene_name_col and ("gene_name" not in adata.var.columns):
        adata.var["gene_name"] = adata.var_names.astype(str)

    # ----------- 写入文件 -----------
    print(f"[INFO] 正在保存到: {args.out_h5ad}")
    # 移除旧的 uns 以防保存时再次触发错误
    adata.uns = {} 
    adata.write_h5ad(args.out_h5ad, compression="gzip")
    
    print(f"[OK] 转换完成！")
    print(f"     保留基因数: {adata.n_vars}")
    print(f"     保留细胞数: {adata.n_obs}")

if __name__ == "__main__":
    main()
