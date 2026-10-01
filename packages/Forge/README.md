# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version is the foundation. It shows the modeling engines available on your server and the model catalog,
which stays empty for now. Training and applying models arrive in the next versions.

## Open Forge

* **Browse > Apps > Forge**
* **ML | Forge | Models**

The view has two sections:

* **Engines**: every engine Forge found, with its package, kind, roles and the hyperparameters it accepts
* **Models**: the model catalog with the name, engine, task, target, data storage mode, number of training
  rows and creation date of each model

The built-in **ML | Models** tools keep working next to Forge. The two do not share models.

## Engines

Forge does not train models itself. Engines live in other packages and in scripts:

* **EDA**: Linear Regression, Softmax, PLS Regression, XGBoost, SVM
* **Chem**: Chemprop, on servers where it is published
* **Scripts** that follow the same contract

An engine appears only when its package is installed.

To add an engine, register functions that share the same `meta.mlname` (the engine name) and set
`meta.mlrole` on each of them:

| Role | Required | Purpose |
|---|---|---|
| `train` | yes | Trains a model on a table and a target column, returns the model |
| `apply` | yes | Applies a trained model to a table, returns the predictions |
| `isApplicable` | yes | Tells whether the engine can learn from the given data |
| `isInteractive` | no | Tells whether training is fast enough to rerun on every change |
| `visualize` | no | Shows engine-specific views of a trained model |

The inputs of the `train` function, except the table and the target column, become the engine hyperparameters.

## Data storage

Every model records what it keeps of its training data:

* **None**: only a summary of the data, enough to check that a new table fits the model
* **Reference**: a link to the source of the data (a file, a query or a script); the data itself is not uploaded
* **Copy**: a copy of the training table, uploaded to the server

Training data never leaves your browser unless you choose **Copy**.
