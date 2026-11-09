# EDA

[![Docs: Predictive Modeling](https://img.shields.io/badge/docs-Predictive%20Modeling-blue.svg)](https://datagrok.ai/help/learn/)
[![Docs: Dimensionality Reduction](https://img.shields.io/badge/docs-Dimensionality%20Reduction-green.svg)](https://datagrok.ai/help/explore/dim-reduction)
[![Docs: Cluster Data](https://img.shields.io/badge/docs-Cluster%20Data-yellow.svg)](https://datagrok.ai/help/explore/cluster-data)
[![Docs: Multivariate Analysis](https://img.shields.io/badge/docs-Multivariate%20Analysis-red.svg)](https://datagrok.ai/help/explore/multivariate-analysis)
[![Docs: Group Comparison](https://img.shields.io/badge/docs-Group%20Comparison-purple.svg)](https://datagrok.ai/help/explore/group-comparison)
[![Tutorial: Multivariate Analysis](https://img.shields.io/badge/tutorial-Multivariate%20Analysis-olive.svg)](https://public.datagrok.ai/apps/tutorials/Tutorials/MachineLearning/MultivariateAnalysis)

EDA is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai) platform. It provides the following exploratory data analysis and machine learning tools:

* Dimensionality reduction
  * Principal component analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis))
  * Uniform manifold approximation and projection ([UMAP](https://arxiv.org/abs/1802.03426))
  * t-distributed stochastic neighbor embedding ([t-SNE](https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding))
  * Stochastic proximity embedding ([SPE](https://pubmed.ncbi.nlm.nih.gov/12820129/))
* Multivariate analysis
  * Partial least squares regression ([PLS](https://en.wikipedia.org/wiki/Partial_least_squares_regression))
  * Quadratic PLS regression
* Supervised machine learning
  * [Linear regression](https://en.wikipedia.org/wiki/Linear_regression)
  * [Softmax](https://en.wikipedia.org/wiki/Multinomial_logistic_regression) classifier
  * [PLS](https://en.wikipedia.org/wiki/Partial_least_squares_regression) regression
  * [Gradient Boosting](https://en.wikipedia.org/wiki/Gradient_boosting)
  * [Support vector machine](https://en.wikipedia.org/wiki/Support-vector_machine) (SVM) classification and regression
  * [Random forest](https://en.wikipedia.org/wiki/Random_forest)
  * [Extra trees](https://en.wikipedia.org/wiki/Random_forest#ExtraTrees) (extremely randomized trees)
  * [Hist gradient boosting](https://en.wikipedia.org/wiki/Gradient_boosting) (histogram-based)
  * [CatBoost](https://en.wikipedia.org/wiki/CatBoost)
  * [Decision tree](https://en.wikipedia.org/wiki/Decision_tree_learning)
  * [k-nearest neighbors](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) classification and regression
  * [Naive Bayes](https://en.wikipedia.org/wiki/Naive_Bayes_classifier) classifier
  * [Linear discriminant analysis](https://en.wikipedia.org/wiki/Linear_discriminant_analysis) (LDA) classifier
  * [Neural network](https://en.wikipedia.org/wiki/Multilayer_perceptron) (multilayer perceptron)
  * [RealMLP](https://github.com/dholzmueller/pytabkit) (tuned-default tabular MLP, pytabkit)
  * [TabM](https://github.com/yandex-research/tabm) (parameter-efficient MLP ensemble)

  The machine learning methods, PCA and PLS run on Crux (scikit-learn, XGBoost and CatBoost ports and the RealMLP / TabM neural tabular models compiled to WebAssembly). Classifiers return the class probabilities (`proba_<class>` columns) after the prediction; SVM only with Class probabilities checked (Platt scaling: an internal 5-fold cross-validation, about 5 times the training time), which Forge sets from its Predict probability. Every classifier's training function declares `meta.mlprobabilities` (`true`, or the name of the input that turns the probabilities on); Forge shows AUC-ROC and the ROC curve from them. `Eda:classificationReport` gives Forge AUC-ROC, average precision, the per-class values and the ROC and precision-recall curves (Forge computes balanced accuracy and MCC itself).

  Models saved before the move to Crux are still applied, with the same predictions as before (prediction only, no class probabilities): PLS regression models with stored feature names, Linear regression models from EDA 1.8.0 and Softmax, XGBoost and SVM models from EDA 1.9.0. Older models must be retrained.

  Hyperparameters: XGBoost has Iterations (boosting rounds), Rate (eta), Max depth, L2 / L1 (lambda / alpha on the leaf weights), Row sample and Column sample (the share of rows / features sampled for each tree, default 1) and Model seed. A Model seed input exists only where Crux uses one: Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost, Decision Tree, Neural Network, RealMLP, TabM, XGBoost (only with a sample below 1) and Linear Regression (its solvers are deterministic with these settings). PLS Regression, Softmax, SVM, k-nearest neighbors, Naive Bayes and LDA have no seed (SVM with Class probabilities uses the fixed seed 0, so its probabilities are reproducible).

  Inputs that depend on others: `meta.mlenabled` on a training function lists, per input, the conditions under which it is used (JSON `{"<input>": <condition> | [<condition>, ...]}`; a condition's clauses must all hold, an array means any of them; a clause is `"<other input>": [values]`, `"<other input>": {"gt": n}` / `{"lt": n}` or `"task": ["classification" | "regression"]`). Forge disables an input whose rule fails and keeps its value. SVM: Gamma for the RBF, polynomial and sigmoid kernels, Degree for polynomial, Offset for polynomial and sigmoid, Epsilon for regression; Linear Regression: Iterations only with L1 > 0; XGBoost: Model seed only with Row sample or Column sample below 1.

  Search spaces: every training function declares `meta.mlparallel` (fits may run in parallel: EDA runs Crux in a pool of up to four workers, half the logical CPUs) and, every method but Naive Bayes, `meta.mlspace`: the ranges of its settings worth searching, in the `ml.tune` space format (`int` / `float` with `low`, `high` and `log`, `categorical` with `choices`) plus the end that gives the simpler model (`simpler`). Forge's Tune samples these ranges (narrowed by the user's Search options); its first trial is the form's current values when they lie inside the ranges, and the seed is never searched. For example, XGBoost searches Iterations 10-100, Rate 0.03-0.5, Max depth 1-10, L2 0.1-10 and L1 0.01-10 (log scale but for the depth), Row sample and Column sample 0.5-1.
* Group comparison
  * Two-sample [t-test](https://en.wikipedia.org/wiki/Student%27s_t-test) (Welch and Student)
  * One-way analysis of variances ([ANOVA](https://en.wikipedia.org/wiki/One-way_analysis_of_variance))
  * Control comparisons against a single control ([Dunnett's test](https://en.wikipedia.org/wiki/Dunnett%27s_test) and Holm-corrected Welch's t-tests)
* Missing data imputation
  * k-nearest neighbors method ([KNN](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm))
* Multi-objective optimization
  * [Pareto front](https://en.wikipedia.org/wiki/Pareto_front) viewer
  * [Probabilistic MPO](https://pmc.ncbi.nlm.nih.gov/articles/PMC4716604/)
