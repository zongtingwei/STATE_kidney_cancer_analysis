# STATE_kidney_cancer_analysis
STATE_kidney_cancer_analysis
## 📖 Overview
Nichecompass_Analysis

## 🚀 How to run
```bash
python cosmx_to_state_base_h5ad.py \
  --in-h5ad  /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/adata.h5ad \
  --out-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --celltype-col celltype \
  --sample-col dataset_name \
  --force-csr \
  --add-gene-name-col
```
### Build the pd1ctla4 & pd1vegfr datasets
```bash
cd /data5/zongtingwei/vcc/state
conda activate state
```

```bash
python make_kidney_therapy_designs.py \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --out-pd1ctla4 /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1ctla4.h5ad \
  --out-pd1vegfr /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1vegfr.h5ad
```

```bash
uv run state tx preprocess_infer \
  --adata /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1ctla4.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1ctla4.h5ad \
  --control-condition "non-targeting" \
  --pert-col "target_gene" \
  --seed 42
```

```bash
uv run state tx preprocess_infer \
  --adata /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1vegfr.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1vegfr.h5ad \
  --control-condition "non-targeting" \
  --pert-col "target_gene" \
  --seed 42
```

