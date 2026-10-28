# Forge: predictive models

Forge trains, applies, evaluates and manages machine learning models (classification and regression) on Datagrok
tables. Every function below runs in the browser: call it through the `datagrok_exec` tool with
`grok.functions.call('Forge:<function>', {...})`. They run in the page, not on the server, so a server-side function
call cannot run them. Apart from `shareModel`, which opens the platform's Share dialog, they open no dialog: a mistake
fails the call with one sentence that says what to fix.

## When to use

Any request to predict a column from other columns, or to train, apply, evaluate, compare, find, describe, rename,
tag, share or delete a model. Not for clustering or PCA, and not for the built-in **ML | Models** tool (its menu,
**Browse > Platform > Predictive models** and the table's **Models** pane): that is a different, older tool, and the
two do not share models.

## User intent -> call

| User intent | Call |
|---|---|
| Which methods exist, what settings they take | `Forge:listMethods` |
| Which methods can learn from my table, which one Forge suggests | `Forge:listMethods` with `table` and `target` |
| Which models exist, which fit this table, which predict X | `Forge:listModels` (`table`, `target`, `name`, `task`, `method`, `tags`) |
| Tell me about this model, what does it need, is it any good | `Forge:getModel` (`findings` for what is worth a look) |
| Compare models | `Forge:listModels` (validation metrics per row) or two `Forge:getModel` results |
| Train a model, predict X from Y | `Forge:trainModel`: without `name` it trains and measures only; with `name` it saves, keeping no training data (`storage` `none`) unless the user asks for a copy |
| Apply a model to a table | `Forge:applyModel`, after `getModel` (`applicableTo`) or `listModels` with `table` says the model fits; map a feature whose column has a close name in `columnNamesMap` |
| How good is the model on this data | `Forge:evaluateModel` |
| Rename a model, change its description or tags | `Forge:updateModel` |
| Share a model | `Forge:shareModel` (opens the Share dialog; then fill and confirm the dialog) |
| Delete a model | `Forge:deleteModel` (confirm with the user first) |

`model` in every function is a model id, or a name that only one model you can see has.

## Signatures and results

`?` marks an optional input, `= x` its default. A list is returned as a `DG.DataFrame`, one model as an object.

### listMethods

```
Forge:listMethods(table?: dataframe, target?: string) -> dataframe
```

- `table`: a table to check the methods against, with the view's default features
- `target`: name of the column to predict in the table; by default the table's last column

Returns the table **Methods**, one row per method `trainModel` accepts: **Method**, **Package**, **Method type**,
**Runs on** (`browser`, or `server` when the training data leaves the browser), **Keeps training rows** (the model
file holds training rows, as SVM's does), **Retrains live**, **Roles**, **Hyperparameters** (the setting names),
**Settings** (a json text with one entry per setting: `name`, `type` and, when the method declares them, `default`,
`min`, `max`, `choices` and `description`; these are the names and values `trainModel`'s `hyperparameters` accepts)
and, with a `table`, **Applicable** (the method can learn from the table's default features to predict `target`) and
**Suggested** (the method `trainModel` trains with when given no `method`). Both are decided as `trainModel` decides
them with its default `features` and `options`, after the checks the Predictive model view shows.

Errors: `The table has no column 'Specie'.`, `The table has no columns.`, with a `table`
`No machine learning method is installed. Install the EDA package.` and any problem `trainModel` refuses the table
and `target` for, such as `Training needs at least 10 rows; the table has 3.`

### listModels

```
Forge:listModels(name?: string, target?: string, task?: string, method?: string, tags?: string, table?: dataframe,
  limit?: int = 100) -> dataframe
```

- `name`: a part of the model name, ignoring case
- `target`: the name of the column the model predicts, ignoring case
- `task`: `regression` or `classification`
- `method`: the method as Forge names it, for example `XGBoost`
- `tags`: comma-separated tags; a model is listed when it has every one of them, ignoring case
- `table`: only the models whose features all have a column of the same or a close name in the table; checked in the
  browser, nothing is uploaded
- `limit`: at most this many models

Returns the table **Forge models**, one row per model you can see, newest first: **Id**, **Name**, **Description**,
**Method**, **Task**, **Target**, **Features** (the feature names), **Training rows**, **Data storage** (`none`,
`reference` or `copy`), **Tags**, **Author** (a login), **Created**, then `<Metric> (validation)` for every metric a
listed model has (a regression model's cells are empty in a classifier's metrics).

Errors: `The task must be regression or classification.`, `The limit must be a whole number of at least 1.`

### getModel

```
Forge:getModel(model: string) -> object
```

Returns one object with what the model's **Details**, **Performance** and **Activity** panes show:

| Key | Value |
|---|---|
| `id`, `version`, `name`, `description` | the model's record |
| `tags` | a list |
| `method`, `package` | the method and its package |
| `task`, `target` | `regression` or `classification`, the predicted column |
| `targetType`, `classes` | the target's type, a classifier's classes |
| `features` | the feature names in training order |
| `featureTypes` | `{name: type}` |
| `hyperparameters` | the method's settings |
| `preparation` | the recorded preparation |
| `metrics` | `train`, `validation`, for two classes `positiveClass` and, for a classifier, `classes` (per class: `label`, `rows`, `train` and `validation` values of `precision`, `sensitivity`, `f1`, `auc`, `averagePrecision`); null when none were saved |
| `findings` | `{text, detail}` for a metric below its threshold or a validation value far from the training one |
| `seed`, `splitting` | how the data was split |
| `storage` | `mode`, with a reference its `source`, with a copy its `tableId` |
| `dataset` | `name`, `rows` and `summary` (the data summary kept instead of the data) |
| `hasModelFile`, `keepsTrainingRows` | booleans |
| `author` | a login |
| `created`, `updated` | ISO dates |
| `lastRun` | `when` and `status` of the newest application, null before the first one |
| `applications` | the number of applications |
| `applicableTo` | the names of the open tables the model fits (a column of a close name counts) |

The metric keys inside `train`, `validation` and `evaluation`: `mse`, `rmse`, `mae`, `r2` (regression), `accuracy`,
`f1` (classification), `sensitivity`, `specificity`, `precision`, `npv` (two classes), `auc` (classifiers whose
method returns probabilities, and Predict probability; with more than two classes, the average over the classes of
the one-vs-rest AUCs), `balancedAccuracy` (the average over the classes of the share of each class predicted
correctly), `mcc` (Matthews correlation, -1 to 1), `averagePrecision` (the area under the precision-recall curve; as
`auc`, only with probabilities). `auc` (when the report exists), `averagePrecision` and `classes` come from Crux
through the EDA package's `classificationReport`.

Errors (every function with `model`): `No Forge model 'X' is available to you.`,
`Several models are named 'X'. Use the model id.`

### trainModel

```
Forge:trainModel(table: dataframe, target: string, features?: list<string>, method?: string, hyperparameters?: map,
  options?: map, name?: string, description?: string, tags?: string, storage?: string = 'none',
  showProgress?: bool = true) -> object
```

Trains as the Predictive model view does, with its defaults and 5-fold cross-validation:

- `table`, `target` (a column name): a numerical target makes a regression model, a text or boolean one a classifier
- `features`: column names; by default every numerical column except the target, bigint columns and integer columns
  without empty cells whose values are all different (row numbers, ids)
- `method`: `CatBoost`, `Decision Tree`, `Extra Trees`, `Hist Gradient Boosting`, `k-Nearest Neighbors`, `LDA`,
  `Linear Regression`, `Naive Bayes`, `Neural Network`, `PLS Regression`, `Random Forest`, `RealMLP`, `Softmax`,
  `SVM`, `TabM`, `XGBoost` (the installed methods, `listMethods`); by default the method the view suggests for the data
- `hyperparameters`: the method's settings by name, for example `{iterations: 50}`; settings not named keep their
  defaults
- `options`: the preparation, each key optional:

  | Key | Default | Meaning |
  |---|---|---|
  | `missingValues` | `skip`; `keep` for a method that handles missing values (XGBoost, Random Forest, Extra Trees, Hist Gradient Boosting, CatBoost, Decision Tree, TabM) | `skip` leaves out the rows with empty cells; `impute` fills the features' empty cells; `keep` gives the rows with empty cells to the method as they are, only for a method that handles missing values (rows without a target are always left out) |
  | `neighbors`, `distance` | the EDA imputer's (`4`, `Euclidean`) | the settings of `impute` |
  | `oneHot` | on when every text or yes/no feature has at most 20 categories (while `skipUniqueCategories` is on, features whose values are all different do not count) | One-hot encoding of text and yes/no features |
  | `skipUniqueCategories` | `true` | leaves out text features whose values are all different |
  | `predictProbability` | `false` | for a target with exactly two classes: trains on the probability of the first class |
  | `cutoff` | `0.5` | the probability from which a row is predicted as the first class |
  | `seed` | random | the seed of the cross-validation split, a whole number from 0 to 2147483647 |
  | `folds` | `5` | the number of cross-validation parts, at least 2 |

- `name`: saves the model under this name; without a name the model is trained and measured but not saved.
  `description` and `tags` (comma-separated) go with it
- `storage`: what the saved model keeps of its training data: `none`, `reference` (the table must come from a file,
  a query or a script) or `copy` (uploads the feature and target columns). Nothing is uploaded unless `copy` is asked
  for
- `showProgress`: the task-bar progress `Training <method> model`; cancelling it fails the call with
  `Training was cancelled.`

Returns `modelId` (null when not saved), `runId`, `name`, `method`, `task`, `target`, `features`, `rowCount` (the
rows trained on), `skippedRows`, `seed` and `metrics` (`train`, `validation`, for two classes `positiveClass` and, for
a classifier, `classes`). Nothing is added to the table.

Errors, for example:

- `The table has no column 'Specie'.`, `The table has no columns 'a', 'b'.`
- `No method 'Boost' is installed. The methods are: <methods>.`
- `The method XGBoost has no setting 'speed'. Its settings are: <settings>.`,
  `The setting 'rate' of Softmax needs a number.`, `The setting 'iterations' of XGBoost must be at most 100.`
- `Unknown option 'cutof'. The options are: missingValues, neighbors, distance, oneHot, skipUniqueCategories, predictProbability, cutoff, seed, folds.`
- `Predict probability needs a target with exactly two classes.`
- `The option 'missingValues' needs skip, impute, keep (keep: only for a method that handles missing values).`,
  `SVM needs rows without missing values. Choose Skip rows or Impute.`
- `Enter a name for the model.`,
  `The table has no known source, so the model cannot keep a reference to it. Use storage none or copy.`
- any problem the view shows for the selection, such as `Choose at least one feature.`

### applyModel

```
Forge:applyModel(model: string, table: dataframe, columnNamesMap?: map, showProgress?: bool = true) -> dataframe
```

Adds the prediction column `<target> (predicted)` (`<target> (predicted 2)` and so on when the name is taken) to
`table` and returns the table. Features are matched to columns of the same name (ignoring case); `columnNamesMap`
maps the others (`{'Sepal.Length': 'sepal length'}`, model feature name -> table column name), including the
columns of a close name that `applicableTo` and `listModels` count as fitting. Missing values are handled as at
training: a model trained with Skip rows gives the rows with empty cells no prediction, one trained with Impute fills
them with its settings, one trained with the gaps kept predicts every row; a model that does not record it keeps the
gaps when its method handles missing values and skips those rows otherwise.
`showProgress`: the cancellable task-bar progress `Predicting <target>`. The application
is recorded with the model.

Errors, for example:

- `The table has no column 'Sepal.Width' the model needs. Map it in columnNamesMap.`,
  `The table has no columns 'a', 'b' the model needs. Map them in columnNamesMap.`
- `The method 'X' is not installed. Install the <package> package.`
- `Application was cancelled.`

### evaluateModel

```
Forge:evaluateModel(model: string, table?: dataframe, target?: string, columnNamesMap?: map,
  showProgress?: bool = true) -> object
```

Measures a saved model as **Run Evaluation** does:

- `table`: the table to evaluate on; by default the file or query the model keeps a reference to, or its uploaded
  copy, which does not join the workspace
- `target`: the column with the actual values, by default the column named as the model's target
- `columnNamesMap`: the features in columns named differently, as for `applyModel`
- `showProgress`: the task-bar progress `Evaluating <model>`; cancelling it fails the call with
  `Evaluation was cancelled.`

Returns `id`, `name`, `table` (the evaluated table's name), `target`, `rowCount` (the rows predicted), `skippedRows`
(the rows with a missing feature value, handled as `applyModel` does; 0 for a model trained with the gaps kept),
`positiveClass` (for two
classes), `metrics` (`train` and `validation` as saved, `evaluation` on this data), `classes` (a classifier: the
evaluation's per-class values, each `label`, `rows`, `values`) and `insights` (each a `text`, with the numbers behind it in `detail` or `items`).
Nothing is added to the table and the evaluation is not recorded as an application.

Errors, besides the ones of `applyModel`'s mapping:

- `The model keeps no training data. Pass a table with the actual values.`
- `The table has no column 'Specie' with the actual values. Pass the column's name in target.`

### updateModel

```
Forge:updateModel(model: string, name?: string, description?: string, tags?: string) -> object
```

Changes what **Edit model...** changes: `name` (not blank), `description` and `tags` (comma-separated; they replace
the current tags, and an empty text removes them). A field not given stays as it is; the last change wins. Returns
`id`, `name`, `description`, `tags` (a list) and `version` after the change; an open catalog and a context panel
showing the model show the change by themselves.

Errors: `Enter a name for the model.`, `Give a name, a description or tags to change.`

### shareModel

```
Forge:shareModel(model: string)
```

Opens the platform's Share dialog for the model, the one **Share...** in its **Sharing** pane opens, and returns at
once, without waiting for the dialog. The share itself is made in the dialog. A missing model fails the call with no
dialog.

### deleteModel

```
Forge:deleteModel(model: string) -> string
```

Deletes the model as the catalog's trash icon does, with its model file, its uploaded copy and its applications, and
returns its id; its training runs stay, without the model.

## Recipes

Train and save, then read the model back as the check (the assertion runs in its own scope, so it finds the model by
name):

```js
const result = await grok.functions.call('Forge:trainModel', {table: t, target: 'Species', name: 'Iris species'});
return {modelId: result.modelId, method: result.method, validation: result.metrics.validation};
```

```js
// verify
const model = await grok.functions.call('Forge:getModel', {model: 'Iris species'});
return {id: model.id, method: model.method, validation: model.metrics?.validation};
```

The models that fit the current table, as a grid in the chat:

```js
return DG.Viewer.grid(await grok.functions.call('Forge:listModels', {table: t})).root;
```

Describe a model:

```js
const m = await grok.functions.call('Forge:getModel', {model: 'Iris species'});
return ui.tableFromMap({Name: m.name, Method: m.method, Task: m.task, Target: m.target,
  Features: m.features.join(', '), 'Training rows': m.dataset.rows, Author: m.author, Created: m.created,
  Applications: m.applications, 'Applicable to': m.applicableTo.join(', '),
  Findings: m.findings.map((f) => f.text).join(' ')});
```

Evaluate and report the evaluation next to the validation, with the insights as the explanation:

```js
const r = await grok.functions.call('Forge:evaluateModel', {model: 'Iris species', table: t});
const values = {};
for (const [metric, value] of Object.entries(r.metrics.evaluation))
  values[metric] = `${value.toFixed(3)} (validation ${r.metrics.validation[metric]?.toFixed(3) ?? '-'})`;
return ui.divV([ui.tableFromMap(values), ...r.insights.map((i) => ui.divText(i.text))]);
```

Apply to the current table and name the new column:

```js
const before = t.columns.names();
await grok.functions.call('Forge:applyModel', {model: 'Iris species', table: t});
return {column: t.columns.names().find((name) => !before.includes(name))};
```

```js
// verify
return t.columns.names().filter((name) => name.startsWith('Species (predicted'));
```

## Rules

- Never upload data unless the user asks for it: `trainModel` keeps no training data by default; `storage: 'copy'`
  uploads the feature and target columns.
- Say what the model is fit for from `getModel.findings` and `evaluateModel.insights`, in plain words.
- Prefer the view's defaults (features, method): pass `features` and `method` only when the user names them.
- Apply a model only after `getModel` (`applicableTo`) or `listModels` with `table` says the table fits. Both accept
  a column of a close name, while `applyModel` and `evaluateModel` take only the same name (ignoring case): when they
  fail with `The table has no column '…' the model needs. Map it in columnNamesMap.`, pass `columnNamesMap` mapping
  the feature to its close column.
- Models are visible per user: a model another user did not share is not found
  (`No Forge model 'X' is available to you.`).
- After `shareModel`, complete the dialog it opened with the dialog functions (`list_view_widgets`, then
  `getShareInfo` / `addShareGrantee` / `setShareAccess`, or `getDialogInfo` / `setInput` / `clickButton`).
- Confirm with the user before `deleteModel`.
- Forge's menu is **ML | Forge**; **ML | Models** is the built-in tool, not Forge.
