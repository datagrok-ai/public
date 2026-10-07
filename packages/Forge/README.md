# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version trains models with XGBoost, saves them, lists them in the catalog with their details, performance and
activity, compares them, and applies them to new data, from the menu, the model's context menu or from scripts.

## Train a model

Open a table, then **ML | Forge | Train...**. The **Predictive model** view opens with the inputs on the left, in two
groups you can fold with their chevrons: **Data** (Table, Target, Features, Missing values) and **Method** (Method and
its settings). Both start open, and a group with a problem opens itself.

* **Table**: the table to learn from, the current one by default
* **Target**: the column to predict, the last column by default. A numerical target makes a regression model;
  a text or boolean target makes a classifier. Rows without a target value are left out of the training; the
  tooltip of **Target** says how many
* **Features**: the columns to learn from. By default, all numerical columns except the target and except integer
  columns whose values are all different, such as row numbers and ids. Dates and very large whole numbers (bigint
  columns) cannot be features
* **Missing values**: a choice shown only when a checked feature has empty cells; its tooltip lists the columns and
  how many cells are empty. **Skip rows** (the default) leaves those rows out. **Impute** fills the empty cells from
  the most similar rows (k nearest neighbors, from the EDA package) on a copy, with the settings **Neighbors** and
  **Distance** shown under it; your table keeps its empty cells
* **Method**: the machine learning method, XGBoost in this version
* the method's settings (hyperparameters), with their default values

Hover an input to see what it is for. The inputs are checked as you change them. A problem is shown on the input it
concerns: the input turns red and its tooltip says what to fix, for example that the target is also checked as a
feature. **Train** is available only when the selection can be trained, and not while a training runs.

Click **Train**. A progress bar shows the training in the task bar; you can cancel it there. The **Results** pane
then shows a grid of metrics with a **Metric**, a **Train** and a **Validation** column:

* **Train**: the quality of the model on the rows it was trained on
* **Validation**: the quality on rows the model has not seen. Forge splits the table into five parts at random,
  trains five times on four parts and predicts the fifth, and measures all these predictions together
  (5-fold cross-validation). Each training draws a new random split; its seed is shown and saved with the model

Under the grid, a bullet list gives `Rows: <used> used, <skipped> skipped (missing values)` when rows were left out
because of missing values, `Validation: 5-fold cross-validation on <rows> rows`, `Seed: <seed>` (the seed can be
selected, and the copy icon right after it, **Copy the seed**, copies it) and, for two classes, `Positive class:
<class>`.

Hover a column header or a metric name to see what it means; hover a value to see it in full precision:

| Task | Metrics |
|---|---|
| Regression | **MSE**, **RMSE** (mean squared error and its root; lower is better), **MAE** (mean absolute error; lower is better), **R squared** (the share of the target's variation the model explains; 1 is perfect) |
| Classification | **Accuracy** (the share of correct predictions), **F1** (balances precision and sensitivity; for more than two classes, the average over the classes) |
| Two classes | also **Sensitivity**, **Specificity**, **Precision**, **Negative Predicted Value**, computed for the positive class named under the table (the first class in alphabetical order) |

Training runs in the browser's main thread: on large tables the page pauses for a few seconds.

## Save a model

After training, click **Save** at the top left of the view, enter a **Name** (prefilled; **OK** is unavailable while
it is empty), an optional **Description** and optional **Tags** (type a tag in full and press Enter to make it a chip;
a comma inside a tag splits it into two tags; the box offers no browser autofill list, since a value picked there
would make no chip), and click **OK**. **Save** is available once a model is trained; after saving it greys out
until the next training, so one training is saved once. Forge saves:

* the trained model
* the method, the target, the features in training order, the hyperparameters, the metrics, and how the data was
  split
* a summary of the training data: the number of rows, statistics of each column and a checksum

The training data itself never leaves your browser. Every model catalog open in this browser tab shows the new model
right away.

## Apply a model

Open a table, then **ML | Forge | Apply...** (or use the play icon in the catalog). The **Apply predictive model**
dialog has:

* **Table**: the table to add the prediction to
* **Model**: every model you can see, up to 10000; the models whose features all have a close column in the table
  come first, newest first. Models with the same name show their creation time (with seconds, or a number such as
  `#2`, when that is the same too)
* **Columns**: one line, `N of M matched`, that opens with its chevron into one row per feature of the model, in
  training order. Each row is prefilled with the table column of the same or a close name and the same kind of data
  (numbers for a numerical feature; the same type and semantic type otherwise). A feature without such a column stays
  empty. Only columns of the right kind are offered. Hover a row to see what the feature needs and whether the chosen
  column has empty cells. The rows start folded when every feature has a column; otherwise they start open and the
  line is red
* **Missing values**: a choice shown only when a chosen column has empty cells. **Skip rows** (the default) gives
  those rows no prediction; **Impute** fills the empty cells on a copy (**Neighbors**, **Distance**, needs the EDA
  package) and predicts every row except those it cannot fill (for example, a row whose features are all empty)
* **More options** (folded): **Batch size**, the rows predicted in one step (10000 by default); lower it for heavy
  methods

Opened for a chosen model (the catalog's play icon, a model's **Apply...**), the dialog starts on the catalog's
**Applicable to** table while it is open, otherwise on the current table if the model fits it, otherwise on the first
open table it fits, otherwise on the current table (or the first open one).

A problem (an empty row, a column of the wrong kind, a column used twice, a column removed from the table meanwhile,
a model Forge cannot apply) marks the input it concerns, and **OK** stays unavailable with the first problem as its
tooltip. **OK** closes the dialog and adds the column `<target> (predicted)` (`<target> (predicted 2)` and so on when
the name is taken) to the table, tagged with the model's id (`forge.model`). A progress bar in the task bar shows the
batches; cancel it there to stop before the next batch. A balloon names the new column and the number of skipped
rows. Every application, completed, failed or cancelled, is recorded with the model.

## From scripts

```javascript
await grok.functions.call('Forge:applyModel',
  {model: 'Iris species', table: grok.shell.t, columnNamesMap: {}, showProgress: true});
```

`model` is a model id or a name that only one model you can see has. Features are matched to columns of the same name
(ignoring case); `columnNamesMap` maps the others (`{'Sepal.Length': 'sepal length'}`). A feature without a column, or
a column of the wrong kind, fails the call with a message and no dialog. Rows with missing values get no prediction.

## Model catalog

Open the catalog from **Browse > Apps > Forge** or **ML | Forge | Models**. The view has two sections:

* **Methods**: a grid of every method Forge found, with its package, method type (function or script), roles and the
  hyperparameters it accepts; hover a header, or a **Method type**, **Roles** or **Hyperparameters** cell, for what it
  means
* **Models**: the model catalog with the name, method, task, target, data storage mode, number of training
  rows, tags and creation date of each model. The refresh icon reloads it; the play icon (**Apply model**) opens the
  **Apply predictive model** dialog for the chosen model and shows the table after applying; the trash icon
  (**Delete model**) removes the chosen model after a confirmation. Both are available only while a model is chosen.
  The Compare icon (**Compare in a new view**) is available while two or more models are selected. The icons sit
  next to the **Models** title, blue while available and grey while not. The catalog also reloads itself after a model
  is saved, edited or deleted, keeping the chosen and the selected models; the context panel never keeps a deleted
  model (it shows the chosen model, a comparison of the selected ones that are left, or nothing)

**Applicable to**, on the line under the **Models** title, is a table input over the open tables (its folder icon
opens a file): choose one to see only the models whose features all have a close column in it; empty shows every
model, and closing the chosen table empties it. Nothing is uploaded for that.

Click a model to see it in the context panel. Its title has the model icon, the favorites star, the name and the
model's commands, then come the panes:

* **Details**: the description, **Author**, **Created**, **Updated**, **Table** (the training table and its rows),
  **Last run** (the newest application, with `(failed)` or `(cancelled)` when it did not complete; `Never` before the
  first one), **Applications**, **Features**, **Target**, **Method**, **Task**, **Applicable to** (the open tables
  that fit), and **Tags**, which you can edit right there
* **Performance**: the metrics saved with the model (the **Metric / Train / Validation** grid of the training) and
  how it was validated
* **Activity**: a grid of the model's applications, newest first (up to 100): **When**, **Who**, **Table**, **Rows**,
  **Prediction column**, **Status**, **Source** (`ui` or `api`) and **Duration (ms)**; hover a header, a status or a
  source for what it means
* **Sharing**: the groups the model is shared with and **Share...** (see "Share a model")
* **History**: every change of the model's record (who, when, what)

Wherever the platform shows a model (the Domains view, Browse, the title in the context panel), its menu has Forge's
**Apply...** (the Apply dialog) and **Download** (the model file, `<model name>.bin`) next to the platform's own
commands. The platform's **Edit...** changes every field of the record and cannot open a saved model yet (a platform
issue), and its **Delete** removes the record only, leaving the model file behind. Right-clicking a row of the catalog
offers **Apply...**, **Edit model...** (Name, Description and Tags), **Download** and **Delete model** (with its file
and its applications), and **Compare** for a selection of two or more models.

With two or more models selected in the catalog, the context panel shows them under the title `Compare N models`
side by side as forms across the panel's width, one per model: the name, method, task, target, training rows,
creation date and the validation metrics. The forms come from the **Forms** viewer of the PowerGrid package; without
PowerGrid, the panel shows the comparison as a table and says `Install the PowerGrid package to see the comparison
as forms.` With fewer selected models the panel shows the chosen model again, or, when no model is chosen, the one
selected model; with neither it is emptied.

**Compare** opens the table view **Compare models**: one row per model with its name, description, method, task,
target, training rows and creation date, and the train and validation value of every metric the models have, with the
Forms viewer next to the grid.

In **Browse > Platform > Domains > forge > model**, the card view shows each model as the built-in tool did: the open
tables it fits, what it predicts, by which features, the method and the creation date.

A prediction column's context panel has a **Predicted by** pane with the card of the model that made it.

## Share a model

Every saved model belongs to its author, who can see, apply, edit, delete and share it; nobody else sees it until it
is shared. The author shares a model with **Share...** in its **Sharing** pane, and the pane lists the new groups as
soon as the share is done. A user a model is shared with sees it in the catalog and applies it, and the application
is recorded with the model. Every user trains, saves and deletes models of their own the same way.

Training runs are not shared with their model: each user sees only their own. On a Datagrok 1.28.0 server an
administrator, too, sees only the models shared with them, not every user's (a platform issue).

The built-in **ML | Models** tools keep working next to Forge. The two do not share models.

## Methods

Forge does not train models itself. Methods live in other packages and in scripts:

* **EDA**: Linear Regression, Softmax, PLS Regression, XGBoost, SVM
* **Chem**: Chemprop, on servers where it is published
* **Scripts** that follow the same contract

A method appears only when its package is installed. This version trains with XGBoost only; the other methods are
listed in the catalog view. **Impute** uses the EDA package's `knnImpute` function; without it, **Missing values**
offers **Skip rows** only. Imputing a large table takes a while.

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
