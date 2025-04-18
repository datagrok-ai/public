---
title: "Data science"
sidebar_position: 4
---

Datagrok was built for data scientists, by data scientists. Our goal is to let scientist focus on science, not
infrastructure.

Out of the box, Datagrok provides all tools necessary for
[data ingestion](../../../access/files/files.md),
[transformation](../../../transform/transform.md),
[visualization](../../../visualize/viewers/viewers.md),
[analysis](use-cases/eda.md),
[modeling](../../../learn/learn.md), as well as [deploying models](../../../learn/learn.md#deployment)
and scientific analyses. Scripts and models can be written in any language, such as R or Python.

## Data munging

80% of data science time is spent retrieving and cleaning the data. By natively integrating with the data retrieval
mechanisms and having built-in collaboration features, we drastically reduce that time.

## Reproducibility

* Version control for data
* Version control for code
* Capturing metadata about OS, dependencies, CPU loads, etc

## Data provenance

[Data provenance](../../../govern/audit/data-provenance.md) is the ability to fully understand everything that the result depends
upon.

## Data pipelines

Data pipelines is a core component of the Datagrok platform designed to let end users define
jobs that would get data from disparate data sources, clean or merge the data if needed, run transformations, build
interactive dashboards based on the retrieved data, and publish these dashboards.

## Statistical hypothesis testing

Available hypothesis tests:

* [Welch's t-test](https://en.wikipedia.org/wiki/Welch%27s_t-test): #\{x.TTest}
* [Kolmogorov–Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test): #\{x.KSTest}

Return p-values.

Tests are available on Context Panel in "Commands" section for two selected numerical columns without missing values. Or
from Functions browser "Help | Actions", see "Math or Statistics"
sections.

## Normalization

Available normalizations:

* [min-max](https://en.wikipedia.org/wiki/Feature_scaling)
* [z-scores](https://en.wikipedia.org/wiki/Standard_score)

Any numerical columns can be normalized via Context Panel in "Commands" section. Or from Functions browser "Help |
Actions", see "Math" section.

## Interactive methods

### [Clustering](../../../explore/cluster-data.md)

Performs clustering using k-means algorithm. An interactive visualization lets you see clustering results in real-time,
which makes the process a lot more intuitive.

### [Missing Value Imputation](../../../explore/missing-values-imputation.md)

Allows to do fast and simple missing values imputation using k-nearest neighbours algorithm.

### [Random Data Generation](../../../transform/random-data.md)

Generate columns with random data with different distributions (Normal, Log-Normal, Binomial, Poisson, Uniform), using
the specified parameters.

### Multivariate analysis

[Multivariate Analysis](../../../explore/multivariate-analysis.md) plugin implements partial least squares (PLS)
algorithm. It is an easy-to-interpret, commonly used approach for multidimensional data analysis. It shows the following
on the interactive viewers: scores, explained variance, correlation loadings, predicted vs. reference, and regression
coefficients.

### Predictive modeling

Train models, apply them, compare performance characteristics, deploy, share. Currently, there two ways to train models:

* Using build-in plugin for modelling based on R Caret via OpenCPU. It allows to train models:
  * SVM (linear or radial)
  * Random Forests
  * GBM
* Using [H2O](https://h2o.ai). A model can be built using H2O UI and than exported into the platform in POJO format and
  used from the "Model browser". Supports the following models:
  * Deep Learning (Neural Networks)
  * Distributed Random Forest (DRF)
  * Generalized Linear Model (GLM)
  * Gradient Boosting Machine (GBM)
  * Naïve Bayes Classifier
  * K-Means Clustering
  * Principal Component Analysis (PCA)

Trained models can be shared with other users. In addition to making them discoverable and reusable, the platform might
also [suggest applying models to the datasets](../../../govern/catalog/self-learning-platform.md) (
potentially opened by other users) when it deduces that the input dataset is of the same structure as the dataset the
model was trained on.

See also:

* [Life sciences](use-cases/life-sciences.md)
