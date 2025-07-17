# Nuclear Segmentation and Morphometric Profiling

This module features deep learning–based pipelines for segmenting nuclei and extracting morphometric features from cardiac tissue sections in embryonic wild-type mice. The goal is to identify structural nuclear abnormalities associated with Lem2 mediated Dilated cardiomyopathy.

---

## ⚙️ Setup & Requirements

Install all dependencies via the provided `requirements.txt` file:

```bash
cd NuclearSegmentation/
pip install -r requirements.txt
```

---

# 🚀 Reproducing the Analysis

## 1. Benchmarking Segmentation Models
Benchmarking was performed to identify the optimal segmentation model for cardiac nuclear imaging:

- Exploring_2D_models.ipynb:
Explores the 2D slices/frames of a 3D image and assesses the performance of the Pre-trained Stardist2D model. It also briefly compares Stardist2D with Cellpose.

- Exploring_stardist3D.ipynb
Test StarDist 3D on a small stack to explore its feasibility. This model was not used in the final analysis.

## 2. Training the StarDist Model
Custom training was performed on a manually annotated dataset to optimise segmentation performance for embryonic mouse hearts.

Notebook:
TrainingAndMorphometry/SD_training_Top30Micron_500epochs.ipynb

Training Dataset:
31 2D frames (2048x2048) taken from a single 3D image. Annotated using the Stardist2D plugin in Fiji followed by manual refinement of predicted labels.

Model Output:
Trained model saved to:
```bash
models/stardist_trained_Top30Micron_500epoch_v0/
├── config.json
├── weights_best.h5
├── thresholds.json
└── weights_last.h5
```
📌 This model was used for all subsequent morphometric profiling in the study.

## 3. Nuclear shape descriptors were extracted from segmented masks using the trained StarDist model.

Notebook:
TrainingAndMorphometry/SD_2dMorphometry_Top30Micron.ipynb

Workflow:

1. Load high-resolution cardiac tissue images.

2. Apply segmentation using weights_best.h5 from trained StarDist model.

3. Extract morphometric features of the 2D nuclear cross-sections such as:
- Area
- Solidity (abnormally shaped nuclei with invaginations have solidity < 1)
- Euler Number (Ruptered Nuclear sections would have Euler Number < 1 while it's >1 for merged nuclear sections)

Export results for downstream analysis and flag the abnormally shaped nuclear sections based on the relevant threshold.

## 4. Experimental Denoising with Noise2Void (N2V)
Preliminary experiments were conducted using the Noise2Void framework to reduce imaging noise and improve segmentation quality. However, these denoising steps were not incorporated into the final pipeline due to limited validation and inconclusive benefits.

Training Notebook:
Experiments/n2v_20250610.ipynb

Trained Model Output:
(Excluded from GitHub due to size)

Noise Analysis Notebook:
Experiments/ExploringNoiseProfile.ipynb

⚠️ Note:
The N2V model is still under development and was not used in the final morphometric or segmentation results reported in the study. Future work may involve deeper benchmarking and integration.
   



