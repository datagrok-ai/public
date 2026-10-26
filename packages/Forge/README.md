# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version trains models with every installed method, retrains them as you change the inputs, saves them with or
without their training data, lists them in the catalog with their details, performance and
activity, compares them, and applies them to new data, from the menu, the model's context menu, the table's context
panel or from scripts and Datagrok AI, which can also list methods and models, describe, train, evaluate, edit,
share and delete models.

## Train a model

Open a table, then **ML | Forge | Train...**. The **Predictive model** view opens with the inputs on the left, in
groups you can fold with their chevrons: **Data** (Table, Target, Features), **Preparation** (see "Preparation"; shown
only while one of its inputs applies) and **Method** (Method and its settings). They start open, and a group with a
problem opens itself. The results are on the right; drag the splitter between them to resize.

* **Table**: the table to learn from, the current one by default
* **Target**: the column to predict, the last column by default. A numerical target makes a regression model;
  a text or boolean target makes a classifier. Rows without a target value are left out of the training; the
  tooltip of **Target** says how many
* **Features**: the columns to learn from. By default, all numerical columns except the target, bigint columns and
  integer columns without empty cells whose values are all different, such as row numbers and ids. Dates and very large whole numbers (bigint
  columns) cannot be features; text and yes/no features need **One-hot encoding** (see "Preparation"). A column
  removed from the table leaves **Features** at once (and empties **Target** when it was the target)
* **Method**: the machine learning method. The list holds the methods that can learn from the chosen target and
  features (a text target, for example, leaves out the regression methods). Until you choose one, Forge suggests a
  method as the built-in tool did: XGBoost for a classifier, PLS Regression for a regression with five or more
  numerical features, Linear Regression for any other regression (Chemprop for a molecule feature, where it is
  installed). The method you choose stays while it can learn from the data; when a change rules it out, a balloon
  says `<Method> cannot be used with this selection; <suggested> is chosen.` When no method can learn from the
  data, **Method** is empty and red: `No method can learn from this selection. Check the features and the target.`
* the method's settings (hyperparameters), with their default values. Values you change are kept per method while
  the view is open, also when you switch the method or the table

Hover an input to see what it is for. The inputs are checked as you change them. A problem is shown on the input it
concerns: the input turns red and its tooltip says what to fix, for example that the target is also checked as a
feature, or that a setting is out of its range (`Value must be less than 100`). Nothing trains while there is a
problem.

The model trains by itself: 200 ms after each change that leaves a valid selection, a new training starts, and a
change during a training stops it (it is not recorded). A method that says it is too slow for this data (for example
SVM on 20000 rows) does not retrain on every change; **Performance** then reads `<Method> takes a while on <N> rows, so it
does not retrain on every change. Click Train.` and a **Train** button appears under the inputs, at their right edge.
It trains the current selection (its tooltip: `Train the model on the current selection.`). **Train** is never shown
for a method that retrains by itself. For a slow method it is disabled, with the reason as its tooltip, while the
selection has a problem (the target, the features, or a setting out of its range), is checked, a training runs, or
after a training until you change an input. When the problem is a setting out of its range (a hyperparameter or a
missing-values setting), **Performance** reads `Fix the settings.` A progress bar `Training <Method> model` shows each
training in the task bar; you can cancel it there (`Training was cancelled.`). While a training runs, a small loader is
shown next to the **Performance** header (the previous results stay while the model retrains by itself). The
**Performance** panel then shows a grid of metrics with a **Metric**, a **Train** and a **Validation** column:

* **Train**: the quality of the model on the rows it was trained on
* **Validation**: the quality on rows the model has not seen. Forge splits the table into five parts at random,
  trains five times on four parts and predicts the fifth, and measures all these predictions together
  (5-fold cross-validation). Each training draws a new random split; its seed is shown and saved with the model

Under the grid, a bullet list gives `Missing values: kept (the method handles them)` when the method kept the empty
cells of the training features, `Rows: <used> used, <skipped> skipped (missing values)` when rows were left out
because of missing values, `Validation: 5-fold cross-validation on <rows> rows`, `Seed: <seed>` (the seed can be
selected, and the copy icon right after it, **Copy the seed**, copies it) and, for two classes, `Positive class:
<class>`.

Hover a column header or a metric name to see what it means; hover a value to see it in full precision:

| Task | Metrics |
|---|---|
| Regression | **MSE**, **RMSE** (mean squared error and its root; lower is better), **MAE** (mean absolute error; lower is better), **R2** (the share of the target's variation the model explains; 1 is perfect) |
| Classification | **Accuracy** (the share of correct predictions), **F1** (balances precision and sensitivity; for more than two classes, the average over the classes) |
| Two classes | also **Sensitivity**, **Specificity**, **Precision**, **Negative Predicted Value**, computed for the positive class named under the table (the first class in alphabetical order) |
| Predict probability | the classification metrics of the labels the cutoff gives, and **AUC-ROC** (the area under the ROC curve of the probabilities: 0.5 is a coin toss, 1 is perfect) |

Training runs in the browser's main thread: on large tables the page pauses for a few seconds.

### Results

The results area holds rows of panels, each under a grey title line (hover it for what the panel shows), every plot and
grid filling its panel under the title; a list or a grid longer than its panel scrolls inside it. The rows
start with set heights and the panels with set shares of the row's width: the plots wider than **Performance** in the
first row of a two-class model, **Insights & Tips** narrower and **Wrong predictions** wider in the second row, equal
widths elsewhere. Drag the divider between two panels or two rows to resize them, and scroll the area for the rows
below; the sizes are kept when you switch to another view and back. A panel that does not apply is absent (its
neighbours share the row) and an empty row disappears:

| Task | Row 1 (350 px) | Row 2 (300 px) | Row 3 (450 px) |
|---|---|---|---|
| Regression | **Predicted <target> vs Actual**, **Residuals**, **Distributions of predicted and actual <target>** | **Insights & Tips**, **Performance**, **Statistics** | **Features**, **Features correlation** |
| Two classes | **Predicted <target> vs Actual**, **ROC Curve**, **Performance** | **Insights & Tips**, **Confusions**, **Wrong predictions** | **Features**, **Features correlation**, **Class balance** |
| More classes | **Predicted <target> vs Actual**, **Confusions** | **Insights & Tips**, **Performance**, **Wrong predictions** | **Features**, **Features correlation**, **Class balance** |

When the second row is left with fewer than three panels (no **Wrong predictions**, or no **Insights & Tips**),
**Class balance** moves to its right end. Then come the method's own views where the method has them (**PLS
Regression**), at most three in a row of 350 px. Before the first result, and while **Train** runs, the area shows
only the **Performance** panel (its title and a hint) with **Insights & Tips** under it.

**Insights & Tips** lists what is worth a look in the data and in the trained model (the panel is absent when there is
nothing to say). Each line is short; hover it to see the numbers behind it, and the tooltip's last line, in grey
(`See Class balance.`), names the panel that shows what it is about, when that panel is shown. That panel's viewer
shows the line as its description too (not **Performance**, whose grid marks its own values):

* about the data, after every check: the target's classes far from equally frequent (`The classes are unbalanced.`;
  hover: each class with its share, `red: 24.6%`, `white: 75.4%`, the shares adding up to 100%; see **Class
  balance**), the same for a text feature (`The classes of <column> are unbalanced.`), numerical columns with a
  correlation above 0.9, in one line (`Features are highly correlated.`; hover: a list of the pairs such as
  `Petal.Length and Petal.Width: 0.96`; see **Features correlation** when a pair of features is among them), a text
  column whose values are almost all different
* about a classifier: few correct predictions (`Only 15% of the rows are predicted correctly.`; see
  **Performance**), one class found much better than the other (`The model finds M much better than F.`; see
  **Confusions**)
* about a regression model: large errors mostly on one side (`... the model over-predicts.` or `under-predicts.`; see
  **Residuals**), statistics of the predictions far from the actual ones, in one line (`The predictions' statistics
  deviate significantly from the actual ones.`; hover: a list such as `variance: predicted 0.214, actual 0.762`; see
  **Statistics**)
* about the metrics: **R2** or **AUC-ROC** below its threshold (`R2 is below 0.3.`, `AUC-ROC is below 0.5.`) and
  validation values far from the training ones, in one line (`The model may be overfitting.`; hover: a list such as
  `Accuracy: train 0.900, validation 0.400`); see **Performance**, whose grid shows each such value in red, with
  `The validation <metric> is far from the training one: the model may be overfitting.` in its tooltip

The viewers:

* a regression model: **Predicted <target> vs Actual** (a scatter plot with a regression line), **Residuals** (the
  prediction minus the actual value, against the prediction), **Distributions of predicted and actual <target>** (a
  histogram in two colors, `TARGET` and `PREDICTION`), **Statistics** (a grid of the eleven summary statistics, `min`
  ... `q3`, of the actual and of the predicted values; hover a statistic for its meaning)
* a classifier: **ROC Curve** (with **Predict probability**: the rate of true positives against the rate of false
  positives at the cutoffs 0, 0.01, ..., 1, drawn as a line with a dot per cutoff, the dashed `Non-informative
  baseline` and the current **Positive class cutoff** marked; the title shows **AUC-ROC** as **Performance** does,
  `ROC Curve (AUC = 0.912)`, the area under all the scores, which can differ slightly from the area under the 101
  points drawn), **Predicted <target> vs Actual** (the first two
  numerical features colored by the target, when there are two; a right prediction is a circle, a wrong one a cross),
  **Confusions** (a confusion matrix, for fewer than 10 classes), **Wrong predictions** (when the model gets a training
  row wrong: a grid of those rows with their row numbers, every column of the table and, last, **Predicted
  <target>**), **Class balance** (a bar per class of the target with its number of rows, each class in its own color,
  without the chart's column selectors; rows skipped for missing values are not counted)
* both: **Features** (a PC plot of the numerical features and, for a regression, the target and the prediction,
  colored by the target) and **Features correlation**, both for 3 to 20 numerical features

The plots show your table's own columns, without copies: hover a point and the table's grid highlights its row, hover
a row of the grid and the point lights up (rows skipped for missing values shift nothing). While **Table**,
**Target**, **Features** and the rows used stay the same (a hyperparameter or **Method** change, a change of **Positive
class cutoff**), a new result keeps the viewers, their zoom and the sizes you dragged, and only the predictions in them
change; a panel that appears or goes (**Wrong predictions**, **Insights & Tips**, the method's own views) lays out again
only the rows it changes (also row 3 when **Class balance** moves), at their starting sizes. Any other change builds the viewers again. A change of **Positive class cutoff**
re-cuts the scores of the last training, without retraining: **ROC Curve** stays as it is and only its mark moves,
the panels of the labels change; a retrain redraws the curve in the same plot.

Under the rows, the panel **Runs** (hover its title: `The trainings of this view, newest first; click one to set its
inputs again.`) lists every training of this view that was not stopped by a change (completed, failed or
cancelled), newest first, the newest one current: **#** (its number in this view), **Method**,
**Target**, **Features** (how many; hover for the names), **Preparation** (`impute (4 neighbors, Euclidean)`,
`one-hot`, `skip unique categories`, `probability (cutoff 0.5)`; empty when none), **Rows**, **Status** (`completed`,
`failed` or `cancelled`; hover a failed one for the error, a cancelled one reads `Stopped from the task bar.`),
**Seed**, and the validation value of every metric any run has (`Accuracy (validation)`; hover for `Train 0.912,
validation 0.845`). Click a row to set its target, features, preparation, method and settings again when they differ
from the form's. A click on a row whose inputs the form already has, while the results show that training, only makes
it current; when only **Positive class cutoff** differs, the click sets the cutoff back and the results are re-cut
without training; after a cancelled or failed training, a click on its row trains it again. The read-only grid moves
nothing with the arrow keys. After a click that sets other inputs, the view
then trains as after any change (a slow method waits for **Train**), and the new training becomes the first row. When
the table no longer has a column of that run, a warning says `The table no longer has the columns of this run:
<columns>.` and nothing changes. A **Table** change starts an empty list; the runs stay recorded.

### Preparation

The **Preparation** group, between **Data** and **Method**, prepares the data before the method gets it. Each input is
shown only while it applies to the checked features and the target, and the group is hidden while none does; an input
that disappears goes back to its default. Hover an input to see what it does.

* **Missing values**: shown when a checked feature has empty cells; its tooltip lists the columns and how many cells
  are empty. **Skip rows** (the default) leaves those rows out. **Impute** fills the empty cells from the most similar
  rows (k nearest neighbors, from the EDA package) on a copy, with the settings **Neighbors** and **Distance** shown
  under it; your table keeps its empty cells. Not shown for a method that handles missing values (see "Methods"): it
  always trains with the empty cells kept, and when the features had empty cells **Performance** (also the saved
  model's) says `Missing values: kept (the method handles them)` and the Runs cell `keep missing values`
* **One-hot encoding**: shown while a checked feature is text or yes/no, and checked when it appears if every such
  feature has at most 20 categories (ids that **Skip unique categories** leaves out do not count), unchecked
  otherwise; once you toggle it, your choice stays until it hides. Its tooltip: "Lets the model use categories as
  numbers without implying an order: each category gets its own 0/1 column." The methods need numbers, so while it is
  unchecked such a feature marks **Features** red: `<Method> needs numerical features. Check One-hot encoding in
  Preparation, or uncheck: <features>.` Checked, it turns each text or yes/no feature into one 0/1 column per value.
  The values are saved with the model, and applying it builds exactly those columns: a value the
  new table lacks, or one the training never saw, does not stop the application (an unseen value counts as none of
  the known ones)
* **Skip unique categories** (checked when it appears): shown while a checked text feature has a different value in
  every row, such as an id. It leaves such features out, so an id never becomes thousands of 0/1 columns. The left-out
  features are saved with the model, and the **Apply predictive model** dialog does not ask for them: a table without
  them can be applied to. When it leaves no feature, **Features** turns red: `Skip unique categories leaves no
  feature. Check more features.`
* **Predict probability** (unchecked by default): shown while the target is text or yes/no with exactly two classes.
  Checked, the model is trained as a regression on the two classes (1 for the positive class, the first one in
  alphabetical order, 0 for the other) and predicts the probability of the positive class; **Method** then lists the
  methods for a regression. The score is the regression model's prediction on the 0/1 target, not a calibrated
  probability, and it can fall outside 0..1; **AUC-ROC** uses the scores, the other metrics the labels after the
  cutoff. **Positive class cutoff** appears under it, a number from 0 to 1 (0.5 by default) with a slider: type a value
  or move the slider, and a row whose score reaches the cutoff is predicted as the positive class. A value that is not
  a number from 0 to 1 marks the input red, **Performance** says `Fix the settings.` and **Save** is unavailable. Changing
  the cutoff does not retrain: the metrics are computed again from the scores of the last training, **AUC-ROC** and
  **ROC Curve** stay (only its mark moves), and **Save** saves the model with the current cutoff. Applied to a table, the model adds the class labels, not the
  scores

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

Open a table, then **ML | Forge | Apply...** (or use the play icon in the catalog or the run icon in the table's
**Forge models** pane). The **Apply predictive model**
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
* **Missing values**: a choice shown only when a chosen column has empty cells, preset to how the model was trained.
  **Skip rows** gives those rows no prediction; **Impute** fills the empty cells on a copy (**Neighbors**,
  **Distance**, preset to the model's, needs the EDA package) and predicts every row except those it cannot fill (for
  example, a row whose features are all empty). A model of a method that handles missing values, trained
  without an explicit Skip rows or Impute (the default, also on data without empty cells), does not show the choice:
  every row is predicted. A model that does not record how it was trained (an older
  one) is applied with the gaps kept when its method handles missing values, otherwise with Skip rows
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

### Models for a table

Show a table in the context panel (for example, click a cell or an empty area of its grid): the pane **Forge models** lists
the models you can see whose features all have a close column in the table, newest first, as the built-in **Models**
pane does. Each line has the small green run icon (**Run prediction**: the **Apply predictive model** dialog with that model and
this table, showing the table after applying), the model's name as a link and the menu arrow. Hover the name for the
model's summary (Created by, Created, Updated, Trained on, Last run, Tags, Inputs, Outputs, Method and, when an open
table fits, Applicable to); click it to show the model in the context panel; right-click it, or click the arrow, for
the model's menu (the platform's row commands with Forge's **Apply...** and **Download**); drag it onto a Space to
link the model there. Under the list, the button **Train** (`Train a new model on this table.`) opens the **Predictive
model** view on the table. When no model fits, the pane says `No model fits this table.` above **Train**. The pane
reads the models each time it opens; nothing is uploaded for that.

## From scripts and Datagrok AI

Nine functions list the methods and the models, describe, train, apply, evaluate, edit, share and delete models without
the UI, from the console (the `~` key), a script or another package. Apart from `shareModel`, which opens the
platform's Share dialog, they open no dialog: a mistake fails the call with one sentence that says what to fix.

### List the methods

```javascript
// Every installed method with its settings
grok.shell.addTableView(await grok.functions.call('Forge:listMethods', {}));

// The same, with the methods that can learn to predict Species from the table and the one Forge suggests
grok.shell.addTableView(await grok.functions.call('Forge:listMethods', {table: grok.shell.t, target: 'Species'}));
```

`Forge:listMethods` returns the table **Methods**, one row per method `trainModel` accepts, in the order of the
catalog's **Methods** grid: **Method**, **Package**, **Method type**, **Runs on** (`browser`, or `server` when the
training data leaves the browser), **Keeps training rows** (the model file holds training rows, as SVM's does),
**Retrains live**, **Roles**, **Hyperparameters** (the setting names) and **Settings**: a json text with one entry per
setting, `name`, `type` and, when the method declares them, `default`, `min`, `max`, `choices` and `description`. These
are the names and values `trainModel`'s `hyperparameters` accepts.

With a `table`, two more columns: **Applicable** (the method can learn from the table's default features, as the
**Predictive model** view picks them, to predict `target`) and **Suggested** (the method `trainModel` trains with when
given no `method`). Both are decided as `trainModel` decides them with its default `features` and `options`.
`target` is a column name, by default the table's last column; a missing one fails the call with `The table has no
column 'Specie'.`, and a table without columns with `The table has no columns.` A table and target `trainModel`
refuses fail the call with the same error, such as `Training needs at least 10 rows; the table has 3.`, and so does
`No machine learning method is installed. Install the EDA package.`

### List the models

```javascript
// Every model you can see, newest first
grok.shell.addTableView(await grok.functions.call('Forge:listModels', {}));

// The classifiers that fit the current table and are tagged demo
grok.shell.addTableView(await grok.functions.call('Forge:listModels',
  {task: 'classification', table: grok.shell.t, tags: 'demo'}));
```

`Forge:listModels` returns the table **Forge models**, one row per model you can see, newest first: **Id**, **Name**,
**Description**, **Method**, **Task**, **Target**, **Features** (the feature names), **Training rows**, **Data
storage** (`none`, `reference` or `copy`), **Tags**, **Author** (a login), **Created**, then `<Metric> (validation)`
for every metric a listed model has (a regression model's cells are empty in a classifier's metrics). Every input is
optional:

* `name`: a part of the model name, ignoring case
* `target`: the name of the column the model predicts, ignoring case
* `task`: `regression` or `classification`
* `method`: the method as Forge names it, for example `XGBoost`
* `tags`: comma-separated tags; a model is listed when it has every one of them, ignoring case
* `table`: only the models whose features all have a close column in the table, as **Applicable to** in the catalog;
  checked in the browser, nothing is uploaded
* `limit` (`100` by default): at most this many models

A wrong task fails the call with `The task must be regression or classification.`, a limit below 1 with `The limit
must be a whole number of at least 1.`

### Describe a model

```javascript
await grok.functions.call('Forge:getModel', {model: 'Iris species'});
```

`Forge:getModel` returns one object with what the model's **Details**, **Performance** and **Activity** show: `id`,
`version`, `name`, `description`, `tags` (a list), `method`, `package`, `task`, `target`, `targetType`, `classes` (a
classifier's classes), `features` (in training order), `featureTypes`, `hyperparameters`, `preparation` (the recorded
preparation), `metrics` (`train`, `validation` and, for two classes, `positiveClass`; null when none were saved),
`findings` (`{text, detail}` for a metric below its threshold or a validation value far from the training one, as
the metrics grid marks them), `seed`, `splitting`, `storage` (`mode`, with a reference its `source`, with a copy its
`tableId`), `dataset` (`name`, `rows` and `summary`, the data summary kept instead of the data), `hasModelFile`,
`keepsTrainingRows`, `author` (a login), `created` and `updated` (ISO dates), `lastRun` (`when` and `status` of the
newest application, null before the first one), `applications` (their number) and `applicableTo` (the open tables the
model fits). `applicableTo` and the `table` of `listModels` count a column of a close name as fitting, while
`applyModel` and `evaluateModel` match only the same name (ignoring case): map such a column in `columnNamesMap`. `model` is an id or a name that only one model you can see has, as for `applyModel`.

### Train

```javascript
// Train and measure, without saving: the result holds the metrics, result.modelId is null
await grok.functions.call('Forge:trainModel', {table: grok.shell.t, target: 'Species'});

// Train and save, with an uploaded copy of the training columns
await grok.functions.call('Forge:trainModel', {table: grok.shell.t, target: 'Species', name: 'Iris species',
  storage: 'copy'});

// Choose the features, the method and its settings
await grok.functions.call('Forge:trainModel', {table: grok.shell.t, target: 'Petal.Length',
  features: ['Sepal.Length', 'Sepal.Width'], method: 'Linear Regression', hyperparameters: {iterations: 50}});
```

`Forge:trainModel` trains as the **Predictive model** view does, with its defaults and 5-fold cross-validation:

* `table`, `target` (a column name): a numerical target makes a regression model, a text or boolean one a classifier
* `features`: column names; by default every numerical column except the target, bigint columns and integer columns
  without empty cells whose values are all different (row numbers, ids)
* `method`: `CatBoost`, `Decision Tree`, `Extra Trees`, `Hist Gradient Boosting`, `k-Nearest Neighbors`, `LDA`,
  `Linear Regression`, `Naive Bayes`, `Neural Network`, `PLS Regression`, `Random Forest`, `RealMLP`, `Softmax`,
  `SVM`, `TabM`, `XGBoost` (the installed methods); by default the method the view suggests for the data
* `hyperparameters`: the method's settings by name; settings not named keep their defaults
* `options`: the preparation, each key optional:

  | Key | Default | Meaning |
  |---|---|---|
  | `missingValues` | `skip`; `keep` for a method that handles missing values | `skip` leaves out the rows with empty cells; `impute` fills the features' empty cells; `keep` gives the rows with empty cells to the method as they are, only for a method that handles missing values (rows without a target are always left out) |
  | `neighbors`, `distance` | the EDA imputer's (`4`, `Euclidean`) | the settings of `impute` |
  | `oneHot` | on when every text or yes/no feature has at most 20 categories (while `skipUniqueCategories` is on, features whose values are all different do not count) | One-hot encoding of text and yes/no features |
  | `skipUniqueCategories` | `true` | leaves out text features whose values are all different |
  | `predictProbability` | `false` | for a target with exactly two classes: trains on the probability of the first class |
  | `cutoff` | `0.5` | the probability from which a row is predicted as the first class |
  | `seed` | random | the seed of the cross-validation split, a whole number from 0 to 2147483647 |
  | `folds` | `5` | the number of cross-validation parts, at least 2 |

* `name`: saves the model under this name; without a name the model is trained and measured but not saved.
  `description` and `tags` (comma-separated) go with it
* `storage` (`none` by default): what the saved model keeps of its training data, as in "Data storage": `none`,
  `reference` (the table must come from a file, a query or a script) or `copy` (uploads the feature and target
  columns). Nothing is uploaded unless you ask for `copy`
* `showProgress` (`true` by default): the task-bar progress `Training <method> model`; cancel it there to stop the
  training, which fails the call with `Training was cancelled.`

The result: `modelId` (null when not saved), `runId`, `name`, `method`, `task`, `target`, `features`, `rowCount` (the
rows trained on), `skippedRows`, `seed` and `metrics` (`train`, `validation` and, for two classes, `positiveClass`).
Nothing is added to the table. Every training that runs is recorded as a training run with the source `api`, saved or
not, cancelled and failed ones too; a call refused for its inputs, or whose data preparation fails, records nothing. The errors a call gets, for example:

* `The table has no column 'Specie'.`, `The table has no columns 'a', 'b'.`
* `No method 'Boost' is installed. The methods are: <methods>.`
* `The method XGBoost has no setting 'speed'. Its settings are: <settings>.`, `The setting 'rate' of Softmax needs a
  number.`, `The setting 'iterations' of XGBoost must be at most 100.`
* `Unknown option 'cutof'. The options are: missingValues, neighbors, distance, oneHot, skipUniqueCategories,
  predictProbability, cutoff, seed, folds.`
* `Predict probability needs a target with exactly two classes.`
* `The option 'missingValues' needs skip, impute, keep (keep: only for a method that handles missing values).`,
  `SVM needs rows without missing values. Choose Skip rows or Impute.`
* `Enter a name for the model.`, `The table has no known source, so the model cannot keep a reference to it. Use
  storage none or copy.`
* any problem the view shows for the selection, such as `Choose at least one feature.`

### Apply

```javascript
await grok.functions.call('Forge:applyModel',
  {model: 'Iris species', table: grok.shell.t, columnNamesMap: {}, showProgress: true});
```

`model` is a model id or a name that only one model you can see has. Features are matched to columns of the same name
(ignoring case); `columnNamesMap` maps the others (`{'Sepal.Length': 'sepal length'}`). A feature without a column, or
a column of the wrong kind, fails the call with a message and no dialog. Missing values are handled as at training:
a model trained with Skip rows gives the rows with empty cells no prediction, one trained with Impute fills them with
its settings, one trained with the gaps kept predicts every row; a model that does not record it keeps the gaps when
its method handles missing values and skips those rows otherwise.

### Evaluate

```javascript
// On the data the model keeps a reference to or a copy of
await grok.functions.call('Forge:evaluateModel', {model: 'Iris species'});

// On a table with the actual values, a feature in a column named differently
await grok.functions.call('Forge:evaluateModel', {model: 'Iris species', table: grok.shell.t,
  columnNamesMap: {'Sepal.Width': 'sepal width'}});
```

`Forge:evaluateModel` measures a saved model as **Run Evaluation** does and returns the metrics next to the saved ones:

* `table`: the table to evaluate on; by default the file or query the model keeps a reference to, or its uploaded
  copy, which does not join your workspace
* `target`: the column with the actual values, by default the column named as the model's target
* `columnNamesMap`: the features in columns named differently, as for `applyModel`; the other features are matched by
  name
* `showProgress` (`true` by default): the task-bar progress `Evaluating <model>`; cancel it there to stop the
  evaluation, which fails the call with `Evaluation was cancelled.`

The result: `id`, `name`, `table` (the evaluated table's name), `target`, `rowCount` (the rows predicted),
`skippedRows` (the rows with a missing feature value, handled as `applyModel` does; 0 for a model trained with the
gaps kept), `positiveClass`
(for two classes), `metrics` (`train` and `validation` as saved, `evaluation` on this data) and `insights` (the **Insights & Tips** of the evaluation and of the
saved metrics: each a `text`, with the numbers behind it in `detail` or `items`). Nothing is added to the table and the
evaluation is not recorded as an application. `model` is an id or a unique name, as for `applyModel`. The errors,
besides the ones of `applyModel`'s mapping:

* `The model keeps no training data. Pass a table with the actual values.`
* `The table has no column 'Specie' with the actual values. Pass the column's name in target.`

### Edit

```javascript
await grok.functions.call('Forge:updateModel', {model: 'Iris species', name: 'Iris classifier', tags: 'iris, demo'});
```

`Forge:updateModel` changes what **Edit model...** changes: `name`, `description` and `tags` (comma-separated; they
replace the current tags, and an empty text removes them). A field not given stays as it is; the last change wins.
It returns `id`, `name`, `description`, `tags` and `version` after the change, and an open catalog and a context panel
showing the model show the change by themselves. `Enter a name for the model.` refuses a blank name, `Give a name, a
description or tags to change.` a call with none of them.

### Share

```javascript
await grok.functions.call('Forge:shareModel', {model: 'Iris species'});
```

`Forge:shareModel` opens the platform's Share dialog for the model, the one **Share...** in its **Sharing** pane opens,
and returns at once, without waiting for the dialog; the share itself is made in the dialog (see "Share a model").
`model` is an id or a unique name, as for `applyModel`; a missing model fails the call with `No Forge model 'X' is
available to you.` and no dialog.

### Delete

```javascript
await grok.functions.call('Forge:deleteModel', {model: 'Iris species'});
```

`Forge:deleteModel` deletes the model as the catalog's trash icon does, with its model file, its uploaded copy and its
applications, and returns its id; its training runs stay, without the model, and a context panel showing it is
emptied. `model` is an id or a unique name, as
for `applyModel`: `No Forge model 'X' is available to you.`, `Several models are named 'X'. Use the model id.`

### Datagrok AI (Grokky)

On a server where the Grokky package is installed, the Datagrok AI assistant works with Forge through the same nine
functions, which it runs in your browser: it can list the methods and your models, describe, train, apply, evaluate,
rename, tag and delete models, and open the Share dialog for a model. Forge publishes what the assistant needs to know
in its `agents/` folder: `package-knowledge.yaml` (what Forge does and where to find it) and `skills.md` (which
function answers which request, the inputs, the results, the error messages and the rules, such as uploading no data
unless you ask for a copy). The built-in **ML | Models** is another tool with its own models; the assistant's Forge
knowledge does not cover it.

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
* **Performance**: the metrics saved with the model (the **Metric / Train / Validation** grid of the training, a value
  below its threshold or far from the training one in red, as in the results) and how it was validated, then **Run
  Evaluation** (for a model trained and saved in Forge; its tooltip: "Apply the model to its training data, or to an
  open table, and show how well it predicts."). It measures the model again on data with known targets: the
  file or query the model keeps a reference to, or its uploaded copy; for a model that keeps neither, the **Evaluate
  model** dialog asks for **Table** (the open tables the model fits) and **Target** (the column with the actual
  values), and without a fitting open table a warning says to open one. While it runs, the task bar shows
  `Evaluating <model>` with a cancel. The view **Evaluation: <model>** then shows the results in rows as after
  training: the insights, the **Performance** panel with an **Evaluation** column next to **Train** and **Validation** and
  the bullet `Evaluation: <table> (<rows> rows)`, and the plots. Nothing is added to the table, a referenced or copied table does not
  join your workspace, and an evaluation is not recorded as an application
* **Activity**: a grid of the model's applications, newest first (up to 100): **When**, **Who**, **Table**, **Rows**,
  **Prediction column**, **Status**, **Source** (`ui` or `api`) and **Duration (ms)**; hover a header, a status or a
  source for what it means. When the shown model is applied in this session (the Apply dialog, the catalog or
  Forge:applyModel), **Activity** and the **Last run** and **Applications** of **Details** show the new application by
  themselves; an application made by another user or in another browser tab shows once you choose another model and
  come back
* **Sharing**: the groups the model is shared with and **Share...** (see "Share a model")
* **History**: every change of the model's record (who, when, what)

Wherever the platform shows a model (the Domains view, Browse, the title in the context panel), its menu has Forge's
**Apply...** (the Apply dialog) and **Download** (the model file, `<model name>.bin`) next to the platform's own
commands. The platform's **Edit...** changes every field of the record and cannot open a saved model yet (a platform
issue), and its **Delete** removes the record only, leaving the model file behind. Right-clicking a row of the catalog
offers **Apply...**, **Edit model...** (Name, Description and Tags, the fields `Forge:updateModel` writes from a
script), **Download** and **Delete model** (with its file and its applications), and **Compare** for a selection of two
or more models.

With two or more models selected in the catalog, the context panel shows them under the title `Compare N models`
side by side as forms across the panel's width, one per model: the name, method, task, target, training rows,
creation date and the validation metrics. The forms come from the **Forms** viewer of the PowerGrid package; without
PowerGrid, the panel shows the comparison as a table and says `Install the PowerGrid package to see the comparison
as forms.` With fewer selected models the panel shows the chosen model again, or, when no model is chosen, the one
selected model; with neither it is emptied.

**Compare** opens the table view **Compare models**: one row per model with its name, description, method, task,
target, training rows and creation date, and the train and validation value of every metric the models have, with the
Forms viewer next to the grid.

In **Browse > Platform > Domains > forge > model**, the card view shows each model in five lines of one height: the
name (with a table icon when an open table fits the model; hover it for `Applicable to <tables>`), what it predicts,
by which features, the method and the creation date; a long line ends with `...`. Hover a card, a model's **Name**
in the catalog grid (its other cells keep the grid's own tooltips), or a model anywhere the platform shows it, for the
full picture: the name in bold, `Predicts <target> from <N> features with <method> ·
<task> · <rows> rows`, the features (a long list ends with `+N more`), and in grey `Created <date> by <author> · Data
storage: <mode> · Tags: <tags>`. The view's filter panel offers **Method**, **Task**, **Data storage**, **Tags** and
**Created**.

A prediction column's context panel has a **Predicted by** pane with the card of the model that made it, with the same
tooltip.

## Share a model

Every saved model belongs to its author, who can see, apply, edit, delete and share it; nobody else sees it until it
is shared. The author shares a model with **Share...** in its **Sharing** pane, and the pane lists the new groups as
soon as the share is done. `Forge:shareModel` opens the same dialog from a script (see "From scripts and Datagrok
AI"). A user a model is shared with sees it in the catalog and applies it, and the application is recorded with the
model. Every user trains, saves and deletes models of their own the same way.

Training runs are not shared with their model: each user sees only their own. On a Datagrok 1.28.0 server an
administrator, too, sees only the models shared with them, not every user's (a platform issue).

To share a model through a Space, drag it onto the Space in **Browse**: its card (**Browse > Platform > Domains >
forge > model**, card mode), its name in a table's **Forge models** pane, or its name in the title of its context
panel. The platform asks to confirm and links the model to the Space; everyone with access to the Space can then see
and apply it. Dropping it again links nothing twice. The platform has no **Add to** command for models yet.

The built-in **ML | Models** tools keep working next to Forge. The two do not share models. **Browse > Platform >
Predictive models** and the table's **Models** pane are the built-in tool's; Forge's are **Browse > Apps > Forge** and
the pane **Forge models**.

## Methods

Forge does not train models itself. Methods live in other packages and in scripts:

* **EDA**: Linear Regression, Softmax, PLS Regression, XGBoost, SVM, Random Forest, Extra Trees, Hist Gradient
  Boosting, CatBoost, Decision Tree, k-Nearest Neighbors, Naive Bayes, LDA, Neural Network, RealMLP, TabM
* **Chem**: Chemprop, on servers where it is published
* **Scripts** that follow the same contract

A method appears only when its package is installed, and every installed method trains in the **Predictive model**
view. Forge passes an empty cell only to a method that handles missing values (XGBoost, Random Forest, Extra Trees,
Hist Gradient Boosting, CatBoost, Decision Tree, TabM): it gets the rows with gaps as they are, and **Missing values** is
not shown for it. For every other method, **Missing values** skips or fills them first. Methods
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
