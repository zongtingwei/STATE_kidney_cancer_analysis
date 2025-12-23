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
#### PD1+CTLA4 939 to 18080genes
```bash
python align_to_model_genes.py \
  --in-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1ctla4.h5ad \
  --out-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1ctla4_18k.h5ad \
  --gene-list /data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/gene_names.csv
```
#### PD1+VEGFR 939 to 18080genes
```bash
python align_to_model_genes.py \
  --in-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1vegfr.h5ad \
  --out-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1vegfr_18k.h5ad \
  --gene-list /data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/gene_names.csv
```

#### Preprocess
```bash
uv run state tx preprocess_infer \
  --adata /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1ctla4_18k.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1ctla4_18k.h5ad \
  --control-condition "non-targeting" \
  --pert-col "target_gene" \
  --seed 42
```

```bash
uv run state tx preprocess_infer \
  --adata /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_design_pd1vegfr_18k.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1vegfr_18k.h5ad \
  --control-condition "non-targeting" \
  --pert-col "target_gene" \
  --seed 42
```

#### Inference
```bash
uv run state tx infer \
  --output "/data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad" \
  --model-dir "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints" \
  --checkpoint "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints/best.ckpt" \
  --adata "/public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1ctla4_18k.h5ad" \
  --pert-col "target_gene"
```

```bash
uv run state tx infer \
  --output "/data5/zongtingwei/vcc/state/pd1vegfr/pred_vegfr.h5ad" \
  --model-dir "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints" \
  --checkpoint "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints/best.ckpt" \
  --adata "/public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1vegfr_18k.h5ad" \
  --pert-col "target_gene"
```

#### draw heatmaps




