# EDA changelog

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

Implemented ANOVA.

### Features

* One-way ANOVA (see Top Menu > ML > Analysis of Variances (ANOVA)...)

## 1.1.3 (2023-08-24)

This release is centered around enhancing user-friendliness and addressing concerns.

### Features

* The method PCA is replaced to Top Menu > ML > Dimensionality Reduction
* The method UMPA is added (see Top Menu > ML > Dimensionality Reduction)
* The method t-SNE is added (see Top Menu > ML > Dimensionality Reduction)
* The method SPE is added (see Top Menu > ML > Dimensionality Reduction)

## 1.1.2 (2023-07-27)

This release focuses on improving usability.

### Features

* Move Multivariate Analysis using partial least squares (PLS) regression to Top Menu | ML
* In PLS and PCA, just numerical columns can be selected as features
