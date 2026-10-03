# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version trains and saves models with XGBoost and lists them in the catalog. Applying a saved model to new data
arrives in the next version.

## Train a model

Open a table, then **ML | Forge | Train...**. The **Predictive model** view opens with the inputs on the left:

* **Table**: the table to learn from, the current one by default
* **Target**: the column to predict, the last column by default. A numerical target makes a regression model;
  a text or boolean target makes a classifier
* **Features**: the columns to learn from. By default, all numerical columns except the target and except integer
  columns whose values are all different, such as row numbers and ids
* **Method**: the machine learning method, XGBoost in this version
* the method's settings (hyperparameters), with their default values

Hover an input to see what it is for. The inputs are checked as you change them. A problem is shown on the input it
concerns: the input turns red and its tooltip says what to fix, for example that the target is also checked as a
feature. **Train** is available only when the selection can be trained, and not while a training runs.

Click **Train**. A progress bar shows the training in the task bar; you can cancel it there. The **Results** pane
then shows a table of metrics with a **Train** and a **Validation** column:

* **Train**: the quality of the model on the rows it was trained on
* **Validation**: the quality on rows the model has not seen. Forge splits the table into five parts at random,
  trains five times on four parts and predicts the fifth, and measures all these predictions together
  (5-fold cross-validation). Each training draws a new random split; its seed is shown and saved with the model

Hover a metric name to see what it means:

| Task | Metrics |
|---|---|
| Regression | **MSE**, **RMSE** (mean squared error and its root; lower is better), **MAE** (mean absolute error; lower is better), **R squared** (the share of the target's variation the model explains; 1 is perfect) |
| Classification | **Accuracy** (the share of correct predictions), **F1** (balances precision and sensitivity; for more than two classes, the average over the classes) |
| Two classes | also **Sensitivity**, **Specificity**, **Precision**, **Negative Predicted Value**, computed for the positive class named under the table (the first class in alphabetical order) |

Training runs in the browser's main thread: on large tables the page pauses for a few seconds.

## Save a model

After training, click **Save** at the top left of the view, enter a **Name** and an optional **Description**, and
click **OK**. **Save** is available once a model is trained; after saving it greys out until the next training, so
one training is saved once. Forge saves:

* the trained model
* the method, the target, the features in training order, the hyperparameters, the metrics, and how the data was
  split
* a summary of the training data: the number of rows, statistics of each column and a checksum

The training data itself never leaves your browser. Every model catalog open in this browser tab shows the new model
right away.

## Model catalog

Open the catalog from **Browse > Apps > Forge** or **ML | Forge | Models**. The view has two sections:

* **Methods**: every method Forge found, with its package, method type (function or script), roles and the
  hyperparameters it accepts
* **Models**: the model catalog with the name, method, task, target, data storage mode, number of training
  rows and creation date of each model. The refresh icon reloads it; the trash icon (**Delete model**) removes the
  chosen model after a confirmation and is available only while a model is chosen. The catalog also reloads itself
  after a model is saved or deleted

The built-in **ML | Models** tools keep working next to Forge. The two do not share models.

## Methods

Forge does not train models itself. Methods live in other packages and in scripts:

* **EDA**: Linear Regression, Softmax, PLS Regression, XGBoost, SVM
* **Chem**: Chemprop, on servers where it is published
* **Scripts** that follow the same contract

A method appears only when its package is installed. This version trains with XGBoost only; the other methods are
listed in the catalog view.

To add a method, register functions that share the same `meta.mlname` (the method name) and set
`meta.mlrole` on each of them:

| Role | Required | Purpose |
|---|---|---|
| `train` | yes | Trains a model on a table and a target column, returns the model |
| `apply` | yes | Applies a trained model to a table, returns the predictions |
| `isApplicable` | yes | Tells whether the method can learn from the given data |
| `isInteractive` | no | Tells whether training is fast enough to rerun on every change |
| `visualize` | no | Shows method-specific views of a trained model |

The inputs of the `train` function, except the table and the target column, become the method's hyperparameters;
their default values and descriptions (shown as tooltips) come from the function's annotation.

## Data storage

Every model records what it keeps of its training data:

* **None**: only a summary of the data, enough to check that a new table fits the model
* **Reference**: a link to the source of the data (a file, a query or a script); the data itself is not uploaded
* **Copy**: a copy of the training table, uploaded to the server

This version saves every model with **None**. **Reference** and **Copy** arrive later. Training data never leaves
your browser unless you choose **Copy**.
