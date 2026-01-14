# STATE_kidney_cancer_analysis
STATE_kidney_cancer_analysis
## 📖 Overview
STATE_kidney_cancer_analysis

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


#### make PD1+CTLA4 pred_unperturbed_datasets
```bash
cd /public2/YouYuning/zongtingwei/cifmv2
conda activate cifmv2
```
```bash
python make_unpert_template.py \
  --input /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1ctla4_18k.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_CONTROL_ctla4_18k.h5ad
```

```bash
uv run state tx infer \
  --output "/data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad" \
  --model-dir "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints" \
  --checkpoint "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints/best.ckpt" \
  --adata "/public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_CONTROL_ctla4_18k.h5ad" \
  --pert-col "target_gene"
```
#### make PD1+VEGFR pred_unperturbed_datasets
```bash
python make_unpert_template.py \
  --input /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_pd1vegfr_18k.h5ad \
  --output /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_CONTROL_vegfr_18k.h5ad
```

```bash
uv run state tx infer \
  --output "/data5/zongtingwei/vcc/state/pd1vegfr/pred_unpert.h5ad" \
  --model-dir "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints" \
  --checkpoint "/data5/zongtingwei/vcc/state/pd1ctla4_checkpoint/checkpoints/best.ckpt" \
  --adata "/public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/val_template_CONTROL_vegfr_18k.h5ad" \
  --pert-col "target_gene"
```



#### draw PD1+CTLA4 heatmaps
```bash
python state_delta_vs_delta_corr.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag PD1_CTLA4_Delta \
  --real-samples "RCC4_1,RCC4_2" \
  --outdir ./state_results_delta_ctla4
```

#### draw PD1+VEGFR heatmaps
```bash
python state_delta_vs_delta_corr.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1vegfr/pred_vegfr.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1vegfr/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag PD1_VEGFR_Delta \
  --real-samples "RCC4_7,RCC4_8,RCC4_17,RCC4_18,RCC5_1,RCC5_2" \
  --outdir ./state_results_delta_vegfr
```

#### draw PD1+CTLA4 Sarcomoid heatmaps
```bash
python state_sarcomoid_analysis.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --task ctla4 \
  --outdir ./sarcomoid_heatmaps_ctla4
```


#### draw PD1+VEGFR Sarcomoid heatmaps
```bash
python state_sarcomoid_analysis.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1vegfr/pred_vegfr.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1vegfr/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --task vegfr \
  --outdir ./sarcomoid_heatmaps_vegfr
```

#### draw PD1 + CTLA4 non-sarcomoid split heatmaps
```bash
python state_delta_vs_delta_corr_split_debug.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag NonSarc_PD1_CTLA4 \
  --ctrl-samples "RCC3_1,RCC3_2,RCC3_3,RCC3_4,RCC3_5,RCC3_6,RCC3_9,RCC3_10,RCC5_4,RCC5_5,RCC5_6,RCC5_15,RCC5_16,RCC5_18" \
  --real-samples "RCC4_1,RCC4_2" \
  --outdir ./state_results_split_nonsarc_ctla4
```

#### draw PD1 + VEGFR non-sarcomoid split heatmaps
```bash
python state_delta_vs_delta_corr_split_debug.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag NonSarc_PD1_CTLA4 \
  --ctrl-samples "RCC3_1,RCC3_2,RCC3_3,RCC3_4,RCC3_5,RCC3_6,RCC3_9,RCC3_10,RCC5_4,RCC5_5,RCC5_6,RCC5_15,RCC5_16,RCC5_18" \
  --real-samples "RCC4_1,RCC4_2" \
  --outdir ./state_results_split_nonsarc_ctla4
```

#### draw PD1 + CTLA4 sarcomoid split heatmaps
```bash
python state_delta_vs_delta_corr_split_debug.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1ctla4/pred.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1ctla4/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag Sarc_PD1_CTLA4 \
  --ctrl-samples "RCC3_7,RCC3_8,RCC3_11,RCC3_12,RCC3_13,RCC3_14,RCC3_15,RCC3_16,RCC5_7,RCC5_8,RCC5_9,RCC5_10,RCC5_11,RCC5_12" \
  --real-samples "RCC4_3,RCC4_4" \
  --outdir ./state_results_split_sarc_ctla4
```

#### draw PD1 + VEGFR sarcomoid split heatmaps
```bash
python state_delta_vs_delta_corr_split_debug.py \
  --pred-pert /data5/zongtingwei/vcc/state/pd1vegfr/pred_vegfr.h5ad \
  --pred-unpert /data5/zongtingwei/vcc/state/pd1vegfr/pred_unpert.h5ad \
  --base-h5ad /public2/YouYuning/zongtingwei/datasets/cifmv2_datasets/kidneyccrcc_cosmx_soupir/kidney_cosmx_state_base.h5ad \
  --tag Sarc_PD1_VEGFR \
  --ctrl-samples "RCC3_7,RCC3_8,RCC3_11,RCC3_12,RCC3_13,RCC3_14,RCC3_15,RCC3_16,RCC5_7,RCC5_8,RCC5_9,RCC5_10,RCC5_11,RCC5_12" \
  --real-samples "RCC4_5,RCC4_6" \
  --outdir ./state_results_split_sarc_vegfr
```

