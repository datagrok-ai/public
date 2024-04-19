# EDA changelog

## 1.1.20 (2024-04-16)

Improve dimensionality reduction dialog.

## 1.1.19 (2024-04-15)

* PLS components computation
* Update UI for `Multivariate Analysis` - main feature & demo app

## 1.1.18 (2024-04-15)

Fixed webGPU numeric distance with 0 range

## 1.1.17 (2024-04-14)

### Features
* Add webGPU UMAP implementation.
* Add webGPU Sparse matrix calculation.
* Add webGPU option for MCL.

## 1.1.16 (2024-04-05)

Add webGPU KNN calculation option to UMAP.

## 1.1.15 (2024-02-27)

Add support for post processing functions in dimensionality reduction. The method `Dimensionality Reduction` (see Top Menu > ML > Dimensionality Reduction) now allows to apply post processing functions to the resulting embeddings.

## 1.1.14 (2024-02-20)

Fix NW distance function

## 1.1.13 (2024-02-12)

Fix PCA & PLS applying to table from local file

## 1.1.12 (2024-02-12)

Add MCL clustering. The method is available in Top Menu > ML > Cluster > MCL.

## 1.1.11 (2023-12-28)

Add missing values imputation using the KNN method.

## 1.1.10 (2023-12-28)

Add option to pass random seed to dimensionality reduction methods. This allows to reproduce results of dimensionality reduction.

## 1.1.9 (2023-12-22)

Improvements to multi column dimensionality reduction.

## 1.1.8 (2023-12-18)

Removed separate methods of dimensionality reduction and substituted with a single method `Dimensionality Reduction` (see Top Menu > ML > Dimensionality Reduction) that supports t-SNE and UMAP. The method allows to use multiple columns (like number, string, Molecule, Macromolecule, etc.) as features with different distance functions. The method also allows to cluster resulting embeddings using DBSCAN algorithm and color resulting scatterplot according to clusters.

## 1.1.7 (2023-12-14)

Add DBSCAN clustering to Top Menu > ML > Cluster > DBSCAN.

## 1.1.6 (2023-10-31)

Methods descriptions are updated.

## 1.1.5 (2023-10-23)

Anova top-menu item is improved.

## 1.1.4 (2023-10-11)

Add one-way ANOVA (see Top Menu > ML > Analysis of Variances (ANOVA)...)

## 1.1.3 (2023-08-24)

* The method PCA is replaced to Top Menu > ML > Dimensionality Reduction
* The method UMPA is added (see Top Menu > ML > Dimensionality Reduction)
* The method t-SNE is added (see Top Menu > ML > Dimensionality Reduction)
* The method SPE is added (see Top Menu > ML > Dimensionality Reduction)

## 1.1.2 (2023-07-27)

* Move Multivariate Analysis using partial least squares (PLS) regression to Top Menu | ML
* In PLS and PCA, just numerical columns can be selected as features
