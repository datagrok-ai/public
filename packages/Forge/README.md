# Forge

[![Docs: Predictive Modeling](https://img.shields.io/badge/docs-Predictive%20Modeling-blue.svg)](https://datagrok.ai/help/learn/)

Forge is a [package](https://datagrok.ai/help/develop/#packages) for the [Datagrok](https://datagrok.ai) platform. It
brings predictive modeling to the platform: train a model on a table, apply it to new data, and manage your models in
one catalog. On a server with Forge, **ML | Models**, **Browse > Platform > Predictive models** and a table's **Models**
pane open Forge; without it the platform offers no predictive modeling tool. Forge needs Datagrok
**<minimum platform version: TBD>** or later.

* Training (**ML | Models | Train Model...**)
  * Regression for a numerical target, classification for a text or yes/no target
  * Every installed method: EDA's Linear Regression, Softmax, PLS Regression, XGBoost, SVM, Random Forest, Extra Trees,
    Hist Gradient Boosting, CatBoost, Decision Tree, k-Nearest Neighbors, Naive Bayes, LDA, Neural Network, RealMLP and
    TabM; Chem's Chemprop for molecules, where it is installed; scripts that follow the same contract
  * Fast methods retrain as you change the inputs; the results show performance, plots and insights
  * **Auto** finds the best model: it tries the suitable methods (Default: the quick ones; Custom: the ones you check)
    and trains the best one. **Search** is **Default** (the recommended settings) or **Custom** (the metric, the time or
    trial limit and the methods to try)
  * **Auto-tune** of one method searches better settings on your data
  * Preparation: one-hot encoding, skipping unique categories, skipping or imputing missing values
  * Validation: k-fold cross-validation or a hold-out split
  * Class probabilities and a threshold for classifiers that give them
  * **Trained models** cards to compare the trainings of the session
* Saving: a model keeps only a summary of its training data, a link to its source or an uploaded copy. Training data
  never leaves the browser unless you choose the copy or train with a method that runs on the server
* Applying (**ML | Models | Apply Model...**, a model's **Apply...**, the table's **Models** pane): the model's
  features are matched to columns of the same or a close name and kind, and any of them can be chosen by hand; rows
  with missing values are skipped or imputed; the prediction is added as a column
* **Models** view: search, sort, cards, details, performance, activity (training runs and applications) and comparison
  of models; edit, share and delete them
* Sharing: like any other entity, or through a Space by dragging the model onto it
* Scripts and Datagrok AI: `Forge:listMethods`, `listModels`, `getModel`, `trainModel`, `evaluateModel`, `applyModel`,
  `updateModel`, `shareModel` and `deleteModel`; `grok.ml.applyModel` applies a Forge model

Models, training runs and applications are stored in the `forge` schema of the platform's domain databases. A method is
a set of functions in any package that share `meta.mlname` and declare their role in `meta.mlrole` (`train`, `apply`,
`isApplicable`, optionally `isInteractive` and `visualize`); see [ARCHITECTURE.md](ARCHITECTURE.md), "Engine contract".
