# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version trains models with every installed method, retrains them as you change the inputs, saves them with or
without their training data, lists them in the catalog with their details, performance and
activity, compares them, and applies them to new data, from the menu, the model's context menu or from scripts.

## Train a model

Open a table, then **ML | Forge | Train...**. The **Predictive model** view opens with the inputs on the left, in two
groups you can fold with their chevrons: **Data** (Table, Target, Features, Missing values) and **Method** (Method and
its settings). Both start open, and a group with a problem opens itself. **Results** is on the right; drag the
splitter between them to resize.

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
* **Method**: the machine learning method. The list holds the methods that can learn from the chosen target and
  features (a text target, for example, leaves out the regression methods). Until you choose one, Forge suggests a
  method as the built-in tool did: XGBoost for a classifier, PLS Regression for a regression with five or more
  numerical features, Linear Regression for any other regression (Chemprop for a molecule feature, where it is
  installed). The method you choose stays while it can learn from the data; when a change rules it out, a balloon
  says `<Method> cannot be used with this selection; <suggested> is chosen.` When no method can learn from the
  data, **Method** is empty and red: `No method can learn from this selection. Check the features and the target.`* the method's settings (hyperparameters), with their default values. Values you change are kept per method while
  the view is open, also when you switch the method or the table

Hover an input to see what it is for. The inputs are checked as you change them. A problem is shown on the input it
concerns: the input turns red and its tooltip says what to fix, for example that the target is also checked as a
feature, or that a setting is out of its range (`Value must be less than 100`). Nothing trains while there is a
problem.

The model trains by itself: 200 ms after each change that leaves a valid selection, a new training starts, and a
change during a training stops it (it is not recorded). A method that says it is too slow for this data (for example
SVM on 20000 rows) does not retrain on every change; **Results** then reads `<Method> takes a while on <N> rows, so it
does not retrain on every change. Click Train.` and a **Train** button appears under the inputs, at their right edge.
It trains the current selection (its tooltip: `Train the model on the current selection.`). **Train** is never shown
for a method that retrains by itself. For a slow method it is disabled, with the reason as its tooltip, while the
selection has a problem (the target, the features, or a setting out of its range), is checked, a training runs, or
after a training until you change an input. When the problem is a setting out of its range (a hyperparameter or a
missing-values setting), **Results** reads `Fix the settings.` A progress bar `Training <Method> model` shows each training in
the task bar; you can cancel it there (`Training was cancelled.`). While a training runs, a small loader is
shown next to the **Results** header (the previous results stay while the model retrains by itself). The **Results** pane then shows a grid of metrics
with a **Metric**, a **Train** and a **Validation** column:

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
| Regression | **MSE**, **RMSE** (mean squared error and its root; lower is better), **MAE** (mean absolute error; lower is better), **R2** (the share of the target's variation the model explains; 1 is perfect) |
| Classification | **Accuracy** (the share of correct predictions), **F1** (balances precision and sensitivity; for more than two classes, the average over the classes) |
| Two classes | also **Sensitivity**, **Specificity**, **Precision**, **Negative Predicted Value**, computed for the positive class named under the table (the first class in alphabetical order) |

Training runs in the browser's main thread: on large tables the page pauses for a few seconds.

## Save a model

After training, click **Save** at the top left of the view, enter a **Name** (prefilled; **OK** is unavailable while
it is empty), an optional **Description** and optional **Tags** (type a tag in full and press Enter to make it a chip;
a comma inside a tag splits it into two tags; the box offers no browser autofill list, since a value picked there
would make no chip), choose the **Data storage** (see "Data storage"), and click **OK**. **Save** is available once a
model is trained and no training runs; after saving it greys out until the next training, so one training is saved
once. Forge saves:

* the trained model
* the method, the target, the features in training order, the hyperparameters, the metrics, and how the data was
  split
* a summary of the training data: the number of rows, statistics of each column and a checksum
* with **Reference**, a link to where the table came from; with **Copy**, an uploaded copy of the training columns

A method whose model file itself holds training rows (SVM keeps some of them as support vectors) is recorded with
**Contains training rows** `true` in the model's details; the dialog shows no warning.

A balloon `Model "<name>" saved.` follows, with the links **Apply...** (the Apply dialog for the new model on the
training table) and **Show in the catalog** (opens the catalog, or brings an open one to the front). Every model
catalog open in this browser tab shows the new model right away.

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
  **Data storage** (`None`, `Reference` or `Copy`), **Data source** (with a reference: the file path, the query name
  or `a script`), **Data copy** (with a copy: the uploaded table, or `missing` when it is gone), **Last run** (the newest application, with `(failed)` or `(cancelled)` when it did not complete; `Never` before the
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

A method appears only when its package is installed, and every installed method trains in the **Predictive model**
view. Forge itself never passes an empty cell to a method: **Missing values** skips or fills them first. Methods
that run in the browser (all of EDA's) keep the training data there; a method that runs on the server (a Python or R
script, a function marked `meta.mlserver: true`, or Chemprop) receives it. Chemprop and script
methods follow the same contract but were not available to test on this version's server. **Impute** uses the EDA package's `knnImpute` function; without it, **Missing values**
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

Every model records what it keeps of its training data. You choose it in **Data storage** of the **Save model**
dialog; the line under the choice says what it means:

* **Reference**: a link to the source of the data (a file, a query or a script the platform recorded when the table
  was opened); the data itself is not uploaded. Line: `A link to <file path, query name or a script> is
  saved; the data stays where it is.` Offered, and preselected, only when the table's origin is known: a table opened
  from **Browse > Files**, a query or a script has one; a table built by a script in memory or opened from a local
  file has none
* **None**: only a summary of the data, enough to check that a new table fits the model. Line: `Only a summary of
  the data is saved.` The default when there is no reference
* **Copy**: the feature and target columns, every row, uploaded as the table `<model name> (training data)`. Line:
  `The training columns (<N> rows) will be uploaded to the server.` Deleting the model with Forge deletes the copy

Training data never leaves your browser unless you choose **Copy** (or train with a method that runs on the server).
