# Forge

Forge brings predictive modeling to [Datagrok](https://datagrok.ai): train a model on a table, apply it to
new data, and manage your models in one catalog.

This version trains models with every installed method, retrains them as you change the inputs, saves them with or
without their training data, lists them in the catalog with their details, performance and
activity, compares them, and applies them to new data, from the menu, the model's context menu, the table's context
panel or from scripts and Datagrok AI, which can also list methods and models, describe, train, evaluate, edit,
share and delete models.

## Train a model

Open a table, then **ML | Models | Train Model...**. The **Predictive model** view opens with the inputs on the left, in
groups you can fold with their chevrons: **Data** (Table, Target, Features), **Preparation** (see "Preparation"; shown
only while one of its inputs applies), **Method** (Method, **Settings** and the method's settings; see "Search") and **Validation** (see
"Validation"). They start open but **Validation**, which starts folded, and a group with a problem opens itself. The results are on the right; drag the splitter between them to resize.

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
  The first item, **Auto**, finds the best model instead of training one method (see "Find best model"): its
  tooltip reads `Auto: Find tries every checked method and trains the best one`, nothing trains while it is chosen,
  **Train** is hidden, the results of the model shown stay, and when no model is shown (after a change of the data)
  **Performance** reads `Click Find to try every checked method on this data and train the best one.`; the current
  model card stays the last completed training
* the method's settings (hyperparameters), with their default values. Values you change are kept per method while
  the view is open, also when you switch the method or the table. A setting the method does not use with the other
  values or the task is disabled (SVM: **Gamma** for the Linear kernel, **Degree** for any kernel but Polynomial,
  **Offset** for Linear and RBF, **Epsilon** for a classifier; see "Methods"); it keeps its value, is not checked, and goes to the method, which
  ignores it. **Settings** above them is **Manual** (the default: the settings as inputs) or **Auto-tune** (Search
  options and **Tune**, which searches better values on this data; see "Search"); with **Auto** the group shows
  the status `Not searched`, **Search options** (**Metric**, **Strategy**, **Limit**, **Max time** or **Max trials**,
  **Candidates**) and **Find** instead (see
  "Search")

Hover an input to see what it is for. The inputs are checked as you change them. A problem is shown on the input it
concerns: the input turns red and its tooltip says what to fix, for example that the target is also checked as a
feature, or that a setting is out of its range (`Value must be less than 100`). Nothing trains while there is a
problem.

The model trains by itself: 200 ms after each change that leaves a valid selection, a new training starts, and a
change during a training stops it (it is not recorded). A method that says it is too slow for this data (for example
SVM on 20000 rows) does not retrain on every change; **Performance** then reads `<Method> takes a while on <N> rows, so it
does not retrain on every change. Click Train.` and, in **Manual**, a **Train** button appears at the right, right under
the inputs; when the inputs are taller than the panel and scroll, it stays in view at the panel's bottom edge until
the end of the inputs (see "Compare models").
It trains the current selection (its tooltip: `Train the model on the current selection.`). **Train** is never shown
for a method that retrains by itself. For a slow method it is disabled, with the reason as its tooltip, while the
selection has a problem (the target, the features, or a setting out of its range), is checked, a training runs, or
after a training until you change an input. When the problem is a setting out of its range (a hyperparameter or a
missing-values setting), **Performance** reads `Fix the settings.` A progress bar `Training <Method> model` shows each
training in the task bar; you can cancel it there (`Training was cancelled.`). A training started by **Train** locks
the form, as a search does (see "Search"): until it ends (completed, failed, cancelled or superseded), every input of the
form is disabled (Settings, its reset button, the pins and Search options too, the Find best model block under
**Auto**), **Train**, **Tune** and **Find** wait and a click on
a **Runs** or **Search** row or on a pinned model card does nothing; afterwards each input is as it was. A method that retrains by itself does not lock the form. While a training runs, a small loader is
shown next to the **Performance** header (the previous results stay while the model retrains by itself). The
**Performance** panel then shows a grid of metrics with a **Metric**, a **Train** and a **Validation** column:

* **Train**: the quality of the model on the rows it was trained on
* **Validation**: the quality on rows the model has not seen, as the **Validation** group sets it (see "Validation").
  By default Forge splits the table into five parts at random, trains five times on four parts and predicts the
  fifth, and measures all these predictions together (5-fold cross-validation); for a classifier the split is
  stratified: each part holds each class's share of the rows. The split seed (42 by default) is shown and saved with
  the model, so the same seed gives the same parts. The final model is always trained on all the rows
* **Test** (only after a search that held rows out, see "Search"): the quality of the final model on the rows held out
  before the search, which neither the search nor the training saw

Under the grid, a bullet list gives `Missing values: kept (the method handles them)` when the method kept the empty
cells of the training features, `Rows: <used> used, <skipped> skipped (missing values)` when rows were left out
because of missing values, `Validation: 5-fold cross-validation on <rows> rows` (or `Validation: hold out of 20% (<N>
rows)`; `, stratified` after it when the split was), `Split seed: <seed>` (only when the rows were shuffled; the seed
can be selected, and the copy icon right after it, **Copy the seed**, copies it), after a search `Search: Tune by <metric>,
limit 5 min: <N> trials of 1 method, chosen by the one-standard-error rule` (with a trial limit `limit 100 trials: 1
method, ...`; a search continued with the other kind `limit 1 min + 50 trials: <N> trials of ...`; after Find
best model `Search: Find best model by <metric>, limit 5 min: <N> trials of <M> methods, chosen by the highest <metric>
− SE` (`the lowest <metric> + SE` for an error) or `..., limit 100 trials: <M> methods, ...`, M the methods with a
completed trial, the limit followed by `, Focus on leaders` for that strategy)
with either `Test: <N> rows held out before
the search` or, when nothing was held out, `Validation after a search is optimistic: the best of the trials was chosen
on it`, for two classes, `Positive class: <class>`, and, with **Predict probability** on two classes,
`Classification threshold: <t>` (`Classification threshold: <t>, the best F1 on the validation rows: its validation
values are optimistic` while the threshold is the best value of the metric **Best threshold by** chose, see "Preparation").

Hover a column header or a metric name to see what it means; hover a value to see it in full precision:

| Task | Metrics |
|---|---|
| Regression | **MSE**, **RMSE** (mean squared error and its root; lower is better), **MAE** (mean absolute error; lower is better), **R2** (the share of the target's variation the model explains; 1 is perfect) |
| Classification | **Accuracy** (the share of correct predictions), **F1** (balances precision and recall; for more than two classes, the average over the classes), **Balanced accuracy** (the average over the classes of the share of each class predicted correctly), **MCC** (Matthews correlation, -1 to 1) |
| Two classes | also **Recall (Sensitivity)**, **Specificity**, **Precision**, **Negative Predictive Value**, computed for the positive class named under the table (the first class in alphabetical order) |
| Classifier with class probabilities | **AUC-ROC** (the area under the ROC curve of the probabilities: 0.5 is a coin toss, 1 is perfect), for every classifier whose method returns class probabilities (all EDA classifiers; SVM only with **Predict probability**): for two classes, of the positive class's probability; with more than two classes, the average over the classes of the one-vs-rest AUCs; **Average precision** (the area under the precision-recall curve of the positive class's probability; the average over the classes for more than two); **Log loss** (the average of -ln(probability of the actual class): 0 is perfect, lower is better; a confident wrong prediction costs most) and **Brier score** (the mean squared difference between the class probabilities and the actual class, 0 or 1: 0 is perfect, lower is better; 0.25 is guessing 0.5 for two classes; with more than two classes, summed over the classes), computed by Forge. With **Predict probability**, the label metrics are those of the labels the **Classification threshold** gives; the probability metrics do not depend on it |

Average precision, the per-class values and the curves come from Crux (scikit-learn's metrics) through the EDA
package's `classificationReport` function; without it, or when it fails, Forge shows the other metrics only.
They add a few Crux calls to every training of a classifier (about 0.3 s at 100,000 rows).

Training runs in the browser's main thread: on large tables the page pauses for a few seconds.

### Results

The results area holds rows of panels, each under a grey title line (hover it for what the panel shows), every plot and
grid filling its panel under the title; a list or a grid longer than its panel scrolls inside it. The rows
start with set heights and the panels with set shares of the row's width: **Insights & Tips** narrower in the second
row, equal widths elsewhere. Drag the divider between two panels or two rows to resize them, and scroll the area for
the rows below; the sizes are kept when you switch to another view and back. A panel that does not apply is absent
(its neighbours share the row) and an empty row disappears:

| Task | Row 1 (350 px) | Row 2 (350 px) | Row 3 (450 px) |
|---|---|---|---|
| Regression | **Predicted <target> vs Actual**, **Residuals**, **Distributions of predicted and actual <target>** | **Insights & Tips**, **Performance**, **Statistics** | **Features**, **Features correlation** |
| Two classes | **Predicted <target> vs Actual**, **ROC Curve**, **PR Curve** (equal widths) | **Insights & Tips**, **Confusions**, **Performance** | **Features**, **Features correlation**, **Class balance** |
| More classes | **Predicted <target> vs Actual**, **Confusions**, **Per-class** | **Insights & Tips**, **Performance**, **Class balance** | **Features**, **Features correlation** |

**Insights & Tips** is always in the second row: when there is nothing to say it shows one line that reads the model's
quality (`Balanced accuracy 0.9xx on the validation rows: on average 9x% of each class is predicted correctly.` for a
classifier, `R2 ...: the model explains ...% of the target's variation.` for a regression), so the panels never move
when a setting or the threshold changes. Then come the method's own views where the method has them (**PLS Regression**), at most three in
a row of 350 px, then, when the model gets a training row wrong, **Wrong predictions** as the last row (300 px), right
above **Runs**. Before the first result, and while **Train** runs, the area shows
only the **Performance** panel (its title and a hint) with **Insights & Tips** under it.

**Insights & Tips** lists what is worth a look in the data and in the trained model (with a result that has Balanced
accuracy or R2, as every result of the Train view has, it is never empty: see above). Each line is short; hover it to see the numbers behind it, and the tooltip's last line, in grey
(`See Class balance.`), names the panel that shows what it is about, when that panel is shown. That panel's viewer
shows the line as its description too, above the viewer, under the matrix for **Confusions** (not **Performance**,
whose grid marks its own values):

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
* a classifier: **ROC Curve** (for a two-class classifier whose method returns class probabilities: recall against
  the false positive rate at every distinct probability of the positive class, drawn as a thin line (axes
  `FPR (1 − Specificity)` and `Recall (Sensitivity)`), the dashed `Non-informative baseline` and, with **Predict
  probability** only, the current **Classification threshold** circled; the title shows **AUC-ROC** as
  **Performance** does, `ROC Curve (AUC = 0.912)`: the area under the drawn curve), **PR Curve** (the same two-class
  cases: precision against recall at every cutoff of the probability, axes `Recall (Sensitivity)` and `Precision`; its
  dashed `Non-informative baseline` is the positive class's share of the rows, the precision of a model that guesses;
  titled `PR Curve (AP = 0.812)` with **Average precision**, the threshold circled with **Predict probability**), a
  model saved with the former regression mode of **Predict probability** draws both from its scores in Run Evaluation,
  **Per-class** (three or more classes: a parallel coordinates plot
  with one line per class in the class's color - the same colors as **Features** and **Class balance** - across
  **Precision**, **Recall**, **F1** and **AUC-ROC**, the validation values, the evaluation's in Run Evaluation; the
  legend of the classes is on top when it is short, not shown otherwise; no in-chart filters or color selector),
  **Predicted <target> vs Actual** (the first two
  numerical features colored by the target, when there are two; a right prediction is a circle, a wrong one a cross),
  **Confusions** (a confusion matrix, for fewer than 10 classes), **Wrong predictions** (when the model gets a training
  row wrong: a grid of those rows with their row numbers, every column of the table and, last, **Predicted
  <target>**), **Class balance** (a bar per class of the target with its number of rows, each class in its own color,
  without the chart's column selectors; rows skipped for missing values are not counted)
* both: **Features** (a PC plot of the numerical features and, for a regression, the target and the prediction,
  colored by the target) and **Features correlation**, both for 3 to 20 numerical features

The plots show your table's own columns, without copies: hover a point and the table's grid highlights its row, hover
a row of the grid and the point lights up (rows skipped for missing values shift nothing). While **Table**,
**Target**, **Features** and the rows used stay the same (a hyperparameter or **Method** change, a change of **Predict
probability** or **Classification threshold**), a new result keeps the viewers, their zoom and the sizes you dragged,
and only the predictions in them change; a panel that appears or goes (**Wrong predictions**, the
method's own views) lays out again only the rows it changes, at their starting sizes. Any other change builds the viewers again. A change of **Predict probability** (for a method whose
probabilities are not a setting) or of **Classification threshold** re-cuts the class probabilities of the last
training, without retraining: **ROC Curve** and **PR Curve** stay as they are and only their circles move, the panels
of the labels change; a retrain redraws both curves in the same plots.

Under the rows (and **Wrong predictions**), the panel **Runs** (hover its title: `The trainings of this view, newest first; click one to set its
inputs again.`) lists every training of this view that was not stopped by a change (completed, failed or
cancelled), newest first, the newest one current: **#** (its number in this view), **Method**,
**Target**, **Features** (how many; hover for the names), **Preparation** (`impute (4 neighbors, Euclidean)`,
`one-hot`, `skip unique categories`, `probability (threshold 0.42, best F1)` for a completed
two-class training with the threshold optimized for a metric, `probability (threshold 0.42)` for one edited by hand
or a training that did not complete, `probability` for more classes; empty
when none), **Rows**, **Status** (`completed`,
`failed` or `cancelled`; hover a failed one for the error, a cancelled one reads `Stopped from the task bar.`),
**Split seed** (empty when the rows were not shuffled), and the validation value of every metric any run has
(`Accuracy (validation)`; hover for `Train 0.912, validation 0.845`). Click a row to set its target, features,
preparation, validation (with the split seed), method and settings again when they differ from the form's. A click on a row whose inputs the form already has, while the results show that training, only makes
it current; when only **Predict probability**, **Classification threshold** or **Best threshold by** differs, the click sets
them back (an optimized threshold with its metric, a hand-edited one as edited) and the results are re-cut without
training (for SVM, whose probabilities are a setting, **Predict probability** differs as
that setting, so the click trains); after a cancelled or failed training, a click on its row trains it again. The
**Method** tooltip lists the settings without that hidden one. The read-only grid moves
nothing with the arrow keys. After a click that sets other inputs, the view
then trains as after any change (a slow method waits for **Train**), and the new training becomes the first row. When
the table no longer has a column of that run, a warning says `The table no longer has the columns of this run:
<columns>.` and nothing changes. A **Table** change starts an empty list; the runs stay recorded.

### Compare models

Under the inputs (and **Train** when shown), a card shows the current model: the newest completed training of the
table in this view, on a light blue background. The cards are as wide as the inputs (their right edge is **Tune**'s)
and as tall as their content. A failed or cancelled training, or an input changed before the next
training, leaves it as it is. Its header is the method in bold; under it, for the winner of a search, an italic line
reads `Auto-tuned model` (Find best model) or `Auto-tuned settings` (Tune); other cards show nothing between the
header and the tabs (the features and the preparation are in **Runs**). Two tabs follow, both small tables of one look (the header row darker than the
rows, the names light grey, the values a shade darker and right-aligned, no colored values):

* **Performance** (open at first): the metrics of that training in the columns of the **Performance** grid
  (**Metric**, **Train**, **Validation**, and **Test** after a search with a hold out), with three decimals; three
  rows: first the metric its search optimized (Tune or Find best model), else `Balanced accuracy` for a classifier and
  `R2` for a regressor, then the next metrics in the grid's order. **Show more** under the table (`Show all metrics`)
  lists every metric, the first one still first, and turns into **Show less** (`Show fewer metrics`), on that card
  only (a new current card starts with three rows);
  hovering a value shows it in full with what the column measures and, for a value **Performance** marks, the
  finding; hover a header or a metric for its meaning. The values are the ones recorded with the training: a later
  change of **Predict probability** or **Classification threshold** re-cuts **Performance**, not the card.
* **Settings**: every setting of the method (**Setting**, **Value**) with its value as trained, the ones the method
  ignored included, but not SVM's **Class probabilities**, which **Predict probability** sets (hover a row for its
  meaning; `No settings` for a method without any).

The tabs switch together: a click on a tab of one card shows that tab on every card, and a new card opens on it.

The "+" at the right of the current card's header (hover: `Pin to compare with other models`) pins the model: a light
grey card of the same training appears right under the current card, above the cards pinned before, with a trash icon
(`Remove from the comparison`) that removes it. While a pinned card has the current model's inputs (every input but
**Predict probability**, **Classification threshold** and **Best threshold by**; for SVM, whose probabilities are a
setting, **Predict probability** counts), the current card shows no icon; removing that pinned card brings "+" back. Nothing limits the number of pinned cards. Click a pinned card (not its tab headers or its icon) to set
its inputs again, as a click on its **Runs** row does.

The cards compare models of the same data: a card is shown only while **Target** and the **Validation** inputs (the
scheme, its folds or test fraction, shuffle, and the split seed when shuffled) are the ones it was trained with; a
change hides the other cards at once, and setting the inputs back shows them again before any retraining.
**Features**, the preparation, the method and its settings may differ. The cards are kept for the session of the view
and are not saved; a **Table** change starts an empty list. **Train** stays at the panel's bottom edge while the inputs
scroll; once the cards come into view, it scrolls up with the end of the inputs.

### Search

**Settings** under **Method** (a radio, **Manual** by default; hover: `How the settings are chosen: by hand or by a
search on this data`) chooses how the method's settings are set; switching it runs nothing. **Manual** shows the
settings as inputs. Under the radio a status reads `Default` (the method's defaults; hover: `The method's recommended
settings`), `Auto-tuned` (the values the last tune of this method set in this view, on the table, target and features
it ran on, in any order: `Found by the search`) or `Customized` (`Changed by you`; the tuned values
on other data: `Tuned for <target> on <table>; tune again for this data`); it follows the values and the data, so a
value or a target changed back reads as before. The reset
button at the right of the radio sets the method's defaults at once (`Reset <Method> settings to defaults`; grey with
`<Method> settings are already default` when they are); the mode, the pins, the ranges, the metric and the limit stay.
A method with settings but without a search space (Chemprop, scripts) has its **Auto-tune** disabled (`<Method>
declares no search space.`); a method without settings (Naive Bayes) has no Settings row at all: no radio, reset or
status, and no **Tune**.

**Auto-tune** hides **Train** and shows **Search options** (folded at first, its chevron under the **Method** header's
text; a click on it opens it) and **Tune** at the right under it. **Search options** lists first the search's own rows:

* **Metric**: what the search maximizes (or minimizes for an error), measured by the validation the **Validation**
  group sets: `Balanced accuracy` (the default), `Accuracy`, `F1`, `MCC`, `AUC-ROC`, `Log loss`, `Brier score`
  (these two minimized) for a classifier (`AUC-ROC`, `Log loss` and `Brier score` need a method that returns class
  probabilities: with SVM without **Predict probability** every trial fails), `R2` (the default), `RMSE`, `MAE` for a
  regression
* **Limit**: what ends the search, one limit (a radio, hover: `What ends the search: a time or a trial limit`):
  **Time** or **Trials** (the default); under it the one input it chooses:
* **Max time**: whole minutes, `1` by default, with the unit `min` and a plus/minus clicker (step 1) on one underline
  (hover: `Minutes the search may run`)
* **Max trials**: the search starts no more than this many trials (completed, pruned and failed ones alike; 50 by default;
  hover: `Maximum number of trials, pruned and failed ones included`)

then the settings in the same order, each on one line with a pin button at the right:

* a setting the search tunes: its range, two fields `<from> – <to>` (a choice setting, such as SVM's **Kernel**: a field
  listing the checked choices, all at first, that opens a list of checkboxes); the ranges are the ones the method
  declares
* a pinned setting (a click on its pin, `Pin to exclude from tuning. <Setting>: will be tuned`; pinned: `Pin to exclude
  from tuning. <Setting>: pinned, excluded from tuning`): one field (a list for a choice setting, a checkbox for a
  yes/no setting) with the setting's value, which the search keeps; it is the same value as in **Manual**, both ways. A
  click on the filled pin gives the kept range back
* a setting the method declares no range for (for example **Model seed**): one field with its value and a grey pin,
  `<Setting>: Fixed, can’t be tuned`
* a setting the method ignores with the pinned values or the task (SVM: **Gamma** with **Kernel** pinned to Linear,
  **Epsilon** for a classifier) is disabled and not searched; one whose use depends on a tuned setting is searched

The ranges are checked as you type, and a problem turns the fields red with its text on a line under the row:
`“From” must be less than “to”`, `“From” must be at least <min>`, `“To” must be at most <max>`, `The bounds must be
whole numbers` for a whole setting, `“From” must be above 0` where the method searches on a log scale (`“From” must be
at least 1` for a whole setting), `Enter both
bounds`, `Choose at least one` for a choice setting; a pinned value out of its range shows the setting's own message.
**Max time** (1 to 60) and **Max trials** (1 to 10000) are checked by the platform as any whole-number input with a
range: a wrong value gets a red underline, and hovering the field shows its description with the platform's message
under it (`Value must be greater than 1` below the range, `Value must be less than 60` or `Value must be less than
10000` above it, `'abc' is not an integer.`, `Can't be empty`); only the shown one counts. **Search options** opens on a
problem. The pins and the reset are buttons, so the keyboard reaches them. Each user
keeps the mode, the folded state, the pins, the ranges, the metric and the limit (its kind, the minutes and the trials)
per method (an entry that no longer fits the method falls back to its default; a time stored in seconds by an earlier
version becomes whole minutes, rounded up, 1 to 60).

**Tune** (its tooltip: `Find the best settings`) searches better values of the tuned settings on the current
selection. It is disabled, with the reason as its tooltip, while the selection has a problem or is checked, while a
training runs, while a range has a problem (`Fix the highlighted range`), when nothing is left to tune (`All
settings are fixed`) and while the limit is invalid (the field's first message, captioned: `Max time: Value must be
less than 60`, `Max trials: Value must be greater than 1`); a search never starts with one of these problems.
**Predict probability** stays as the form has it: the trials score the labels its threshold gives (for SVM, every
trial fits with **Class probabilities**, about five times slower); the winner's training then optimizes the threshold
for **Best threshold by** again, unless the threshold was edited by hand. **Missing values** stays as the form has it.

Each trial is a regular training of the method with other settings, validated as the **Validation** group says (its
folds, stratification, shuffle and split seed) on the same folds for every trial; with **Hold-out**, each trial is
scored on that one held out part. The first trial takes the form's current values when each tuned one lies in its
range; otherwise it is drawn too. The trials draw values at random (the split seed) from the ranges, on a log scale
where the method declares one, a fractional value rounded to 3 significant digits (`0.245`, `2.12`, `54.3`) and kept within the range, so the
settings show the value the model trains on; the pinned settings keep their values. A setting the method ignores with the drawn
values (see "Methods") takes its default instead, so trials that differ only there are one trial, and the leaderboard
shows the default. After five trials that passed their first two folds, a trial whose
first two folds score below the median of those earlier trials (the ones started before it that were not pruned there)
stops there (pruned), so parallel fits give the same result as one at a time; a hold out prunes nothing.
EDA's methods run up to four fold fits at a time (half the processor's threads, at most four), each in its own Crux
worker. The time limit is checked between fits: no new trial starts after it, and a trial whose first two folds
forecast an overrun stops (the first trial always completes); one fit of a slow method can run past it. The task bar
shows
the method and what is left of the limit: the trials finished (completed, pruned or failed) of this run's limit,
`Tuning CatBoost (12/50)`, or the time left, `Tuning CatBoost (42 s left)`, `Tuning CatBoost (2 min left)`, and once
the time is up the trials still running, `Tuning CatBoost (finishing 3 trials)`; a continued search the same, counting
this run's trials;
**Save**, **Train** and **Tune** wait, and the loader is next to the **Performance** header. Until the search and the
training of its winner end, every input of the form is disabled and a click on a **Runs** or **Search** row or on a
pinned model card does nothing; afterwards each input is as it was (a setting the method ignores stays disabled). To stop the search, cancel
it in the task bar. When 10 trials have failed and none has completed, the search stops.

The chosen trial follows the one-standard-error rule: among the completed trials whose value is within one standard
error (over the folds) of the best trial's, the simplest settings win (the settings the method marks as simpler: fewer
trees, shallower trees, stronger regularization), then the earlier trial.

When the table has 500 rows or more (after the skipped ones), a fifth of them, rounded and kept between 100 and 2000
rows (each class in its share, within one row, for a classifier), is held out before the search; the search never sees it, and **Performance** shows the final
model's values on it as **Test**. Below 500 rows nothing is held out and the bullets say the validation is optimistic.

When the search stops, the balloon `Optimal settings for <Method> found.` appears (cancelled from
the task bar: `Search cancelled; the best so far is trained.`), the settings are set to the chosen trial's values,
**Settings** switches to **Manual** and reads `Auto-tuned`, and the model trains on the search rows with the same split (`Training <Method> model`), so its **Validation** is the
chosen trial's value. That training is a regular run: **Runs** gets its row (hover its **Method**: `Found by Tune: <N>
trials` above the settings) and **Save** saves it with the search. When no trial completed, a warning says `No trial
completed.`, or `No trial completed: <the first failed trial's error>` when a trial failed, and nothing changes.

**Tune** again continues the last search of the method on this table while nothing it depends on changed (the
target, the features, the preparation with its threshold (the trials score the labels it gives; with **Predict
probability** the winner's training usually moves it, so the next search starts anew), the missing values, the validation and its seed, the metric,
the ranges, the pins and the values of the settings not tuned): the new limit is added on top, the earlier trials
count (none is tried twice), and the leaderboard, the counts and the choice cover all of them; its tooltip then reads
`Continue the search: <N> trials so far` and the task bar `Continuing: Tuning <Method>`. A change of any of these, or
of the table, starts a new search. The search is kept in memory only, for the session of the view.

The panel **Search** (above **Runs**; hover its title: `Each method's best settings, ranked by <Metric> − SE.`) appears
after a search of two methods or more (Find best model); after Tune, or a Find of one method, it is not shown, as
Performance and the run's tooltip describe that search. It shows, in a splitter, a bar chart (two thirds of the width at first) left of the
leaderboard, both over one table, so the current row and the selection are shared: one bar per method with a completed
trial, its **<Metric> − SE** (the grid's column), best first, with the value labels, in one color, with the description
`<Metric> minus its standard error: a cautious estimate` (`plus` for an error) above it, without a title, a line, a
legend or selectors. Draggable splitters separate the last results row (**Wrong predictions**), **Search** and
**Runs**; a block not shown leaves no gap. The leaderboard, in the chart's order (by **<Metric> − SE**, best first, equal
values by the method's name, capitals first; the chosen one is always first and the top bar; the methods
without a completed trial last):
**#** (hover the header: `Rank by <Metric> − SE.`), **Method** (hover for the chosen settings, then `Searched: <setting> <from>–<to>, <setting> (<choices>), ...` and
`Pinned: <settings>`), **<Metric> (CV)**, or **<Metric> (hold-out)** with **Hold-out**
(hover for the fold values; one value with a hold out), **SE** (the standard error of the metric over the folds; 0 with a hold out),
**<Metric> − SE** (the score less its standard error, plus for an error minimized: **<Metric> + SE**; hover the
header, e.g. for F1: `F1 minus its standard error over the folds: a cautious estimate of F1. The method with the
highest value is chosen.`, for RMSE `RMSE plus its standard error ... The method with the lowest value is chosen.`;
with a hold out `(0 on a hold out: one score)` instead of `over the folds`),
**Trials** (hover the header: `Every trial of the method, pruned and failed ones included.`) and **Chosen** (hover
the header: `The method with the highest <Metric> − SE.`, the chosen cell: `Chosen: the highest <Metric> − SE.`). A
click on
its row or its bar sets the method and those settings in
**Manual**, and the view trains as after any change (a click on a **Runs** row switches to **Manual** too). The panel stays until the next search or a **Table** change.

#### Find best model

With **Auto** chosen in **Method**, the **Method** group shows, instead of **Settings** and the method's settings, a
block in the look of **Auto-tune** (the same captions column, fonts, row heights and spacing): under **Method** the
italic status `Not searched` (hover: `Click Find to search`), then the **Search options** toggle (folded at first;
each user keeps it open or folded, and it opens on an invalid limit) over these rows:

* **Metric**: as in Tune's Search options (`Balanced accuracy` or `R2` by default)
* **Strategy**: **Equal**, the equal split (the default; hover: `Every method gets an equal share of the limit: a fair
  comparison`) or **Leaders**, focus on leaders (`Each round the worse half drops out; the rest get more trials`)
* **Limit**, then **Max time** or **Max trials**: as in Tune's Search options (**Trials**, 100 of them, by default;
  **Max time** 1 minute):
  the trials of all methods together, completed, pruned and failed alike, so the methods' **Trials** never add up to
  more. With **Equal** and **Trials**, fewer trials than checked methods are a problem of the field, shown the same
  way (`At least one trial per checked method: <N>`), checked again when the methods or the strategy change
* **Candidates** (hover the caption: `The methods to try, each with its defaults first`,
  followed by one `Not listed: <reason>` line per method left out): a checkbox per method that can learn from the
  selection and run in the browser, cheapest first (Linear Regression, PLS Regression, Softmax, Naive Bayes, LDA,
  Decision Tree, k-Nearest Neighbors, XGBoost, Hist Gradient Boosting, CatBoost, Extra Trees, Random Forest, SVM,
  Neural Network, TabM, RealMLP, then the others). A server method is not listed (`Not listed: <Method>: Runs on the
  server; not searched in the browser`), nor, while **Predict probability** is on, a method without class
  probabilities (`Not listed: <Method> gives no class probabilities`), nor, under **AUC-ROC**, **Log loss** or **Brier
  score**, a method without them with the current settings (SVM with **Predict probability** off: `Not listed: SVM
  gives no class probabilities, which AUC-ROC needs`). Above 1 minute every method is checked; at 1 minute, or with
  **Trials**, only the ones that retrain by themselves on this data (not RealMLP, nor TabM above 200 rows, Neural Network above 2000 rows,
  Random Forest above 20000 cells), and not SVM while **Predict probability** is on (its fits are about five times
  slower). A method you check or uncheck keeps your choice; the others follow the time limit, the metric and
  **Predict probability**

**Find** sits at the right under Search options, where **Tune** does, also while Search options is folded. Each user
keeps the metric, the strategy, the limit and the methods they toggled. Choosing **Auto** runs nothing and keeps the
results of
the model shown.
**Predict probability** stays enabled while a checked method gives class probabilities, or while none is checked (Find
is then disabled); otherwise it is unchecked and disabled: `No checked method gives class probabilities`; toggling it changes the list's defaults and trains nothing.
**Missing values** shows for features with empty cells even when the suggested method could keep them: the search
skips those rows, or imputes them with **Impute**, and never keeps them, as every method gets the same folds; the run,
the card and a saved model record that.

**Find** (its tooltip: `Find the best model`) is disabled, with the reason as its tooltip, as **Tune** is, when no
method is checked (`Check at least one method`) and while the limit is invalid (`Max time: ...`, `Max trials: ...`, `Max
trials: At least one trial per checked method: <N>`). The search holds out the same rows as Tune and validates each trial as
the **Validation** group says. With **Equal** the checked methods run one after another, cheapest first, each
with its defaults first and then its own random draws: each gets an equal share of the trials (rounded up) or the
time left when it starts, so what a method does not use (it finished early, its ranges were used up or its trials kept
failing) goes to the methods after it. **Search**
then explains a short **Trials** count in its tooltip: `Its share of the limit ran out` or `Its space was exhausted`.
With **Leaders** it runs successive halving: first every checked method once with its defaults, cheapest
first (several at a time); then, round by round, the better half of the methods with a search space (by their best value so far) get
twice as many new trials as in the round before (2, 4, ...), each over its whole declared ranges from its own random
draws, and the worse half drops out, but never down to fewer than three; the last three or fewer are then tuned in
turn until the limit ends. A method without a search space (Naive Bayes) is tried with its defaults only (its
**Trials** tooltip in **Search**: `No settings to tune`), and a
method whose trials all failed or whose ranges are used up leaves without counting against the others. One limit
covers all the methods; a trial whose first two folds
forecast an overrun stops (abandoned), so a slow method may be left without a completed trial. The trials score the
labels the form's threshold gives, and the winner's training optimizes the threshold for **Best threshold by** again.
The task bar shows the same progress as Tune's (`Finding best model (79/100)`, `Finding best model (finishing 2
trials)`, `Finding best model (42 s left)`,
`Finding best model (2 min left)`). Each method's best settings follow the one-standard-error rule as Tune's;
across the methods the one with the highest **<Metric> − SE** (the lowest **<Metric> + SE** for an error) wins, a tie
(the same value as the chart shows it) to the method whose name comes first, capitals first, as the chart orders its
bars.

The winner trains as Tune's does: the balloon reads `Best model: <Method> (<Metric> <value>).` (cancelled: `Search
cancelled; the best so far is trained.`), **Method** is set to it with **Settings** in **Manual** (`Auto-tuned`, also
when its defaults won; `Customized` once you change a setting) and the model trains on the search rows (`Training <Method> model`); its **Runs** row
reads `Found by Find best model: <N> trials`. To search again, choose **Auto** again: **Find** then continues the last
search as **Tune** does (the same data, validation, metric, strategy and checked methods; the limit added on top, an
equal split splitting the new limit equally again). **Search**
lists the leaders, then every method without a completed trial with an empty value, its trials and the reason as the
**Method** tooltip (`No completed trial: <error>`, `Abandoned: the budget ran out`, `Not started: the budget ran
out`); a click on such a row sets that method with its defaults. A method the halving dropped reads `Eliminated in
round <N>` in its **Trials** and **Method** tooltips (round 1 is the defaults).

### Preparation

The **Preparation** group, between **Data** and **Method**, prepares the data before the method gets it. Each input is
shown only while it applies to the checked features and the target, and the group is hidden while none does; an input
that disappears goes back to its default (except **Classification threshold** and **Best threshold by**, see below). Hover
an input to see what it does.

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
* **Predict probability** (unchecked by default): shown while the target is text or yes/no with two or more classes.
  It uses the method's own class probabilities; the method stays the classifier it is. Its tooltip: `Adds the
  probability of each class when the model is applied; for two classes, lets you choose the classification threshold
  of <positive class>.` (`Adds the probability of each class when the model is applied.` for more classes). For a
  method that gives no class probabilities (Chemprop, a script without `meta.mlprobabilities`) it is unchecked and
  disabled, its tooltip `<Method> gives no class probabilities`. Checked, the saved model adds the class probabilities
  when it is applied (see "Apply a model"), and the labels come from them: for two classes, the positive class (the
  first one in alphabetical order) where its probability reaches **Classification threshold**, else the other class;
  for more classes, the most probable class. When the method returns no probabilities for this data, nothing is added.
  Toggling it does not retrain: the labels and their metrics are computed again from the class probabilities of the
  last training; **AUC-ROC**, **Average precision**, **Log loss** and **Brier score** do not change. SVM is the
  exception: its probabilities are its **Class probabilities** setting (Platt scaling, an internal cross-validation,
  about five times the training time), which **Predict probability** sets and the settings do not show, so toggling it
  retrains SVM and adds a run
* **Classification threshold** (two classes, under a checked **Predict probability**): a number from 0 to 1 (0.5
  until a best is found) with a slider, shown with 3 significant digits (`0.520`; an exact value as it is: `0.5`, `0.9`) while the full value
  is used, saved and shown everywhere else; its tooltip (on the caption and the box): `Rows whose probability of <positive class> is at
  least this value are predicted as <positive class>, the others as <negative class>.` (the classes in bold, `true` /
  `false` for a yes/no target; the comparison is in single precision, as the probabilities are stored). A
  value that is not a number from 0 to 1 marks the input red, **Performance** says `Fix the settings.` and **Save** is
  unavailable. Editing it by hand does not retrain: the label metrics are computed again, **AUC-ROC**, **Average
  precision**, **Log loss**, **Brier score**, **ROC Curve** and **PR Curve** stay (only their circles move), and
  **Save** saves the model with the current threshold. A hand-edited threshold stays across trainings until **Set**.
  While the threshold is the optimized one, **Performance** adds `Classification threshold: <t>, the best <metric> on
  the validation rows: its validation values are optimistic` (the threshold was chosen on those rows), and the saved
  model records the metric (`thresholdMetric` in its options) and its **Performance** pane reads the same line; after a
  hand edit it reads `Classification threshold: <t>` and records no metric. Unchecking **Predict probability** keeps
  **Best threshold by**, the threshold and its hand edit: checked again, a threshold not edited by hand is optimized
  for the kept metric, and a hand-edited one stays until **Set**. A target with other classes sets them back to `F1`
  and 0.5 and forgets the hand edit.
* **Best threshold by** (two classes, under **Classification threshold**; its tooltip: `The metric the classification
  threshold maximizes on the validation rows.`): `F1` (the default), `Balanced accuracy` (its tooltip says it picks
  the threshold Youden's J picks), `MCC` or `Accuracy` (each item's tooltip describes the metric). The threshold is
  optimized at once: when **Predict probability** is checked and after every new training (a search's winner too),
  **Classification threshold** becomes the value that maximizes the metric on the out-of-fold probabilities (the
  midpoint between the best probability and the next lower one, rounded to 3 decimals when that stays between them; a
  tie goes to the one nearest 0.5). Changing the metric does not change the threshold by itself. The button **Set**,
  at the right end of the **Best threshold by** row (a button as **Tune**, its letters where Tune's end), sets the
  threshold to the chosen metric's best (a re-cut, no training); it is shown only while the threshold differs from
  that best (after a metric change or a hand edit; hover: `Set the threshold that maximizes <metric> on the validation
  rows.`). Under **Classification threshold** an italic status, as the one under the method's **Settings**, says what
  the threshold is: `Best for <metric>` (the chosen metric's best, or after a metric change the best of the metric it
  was set for) or `Customized` (edited by hand). Without a best (no trained model yet, no class probabilities on the
  validation rows, or one class in them; the threshold then stays as it is, 0.5 by default) neither the status nor
  **Set** is shown.

### Validation

The **Validation** group, after **Method** and folded at first, sets how **Validation** is
measured. Its values are part of the selection: a change trains again, they stay when you change the table, a **Runs**
row sets them again, and a search uses them (see "Search"). Hover an input to see what it does.

* **Scheme**: `K-fold, stratified` (the default for a classifier, listed only for one), `K-fold` (the default for a
  regression; a target that becomes numerical turns the stratified choice into it) or `Hold-out`. K-fold predicts each
  part by a model of the other parts; stratified keeps each class's share in every part. Hold-out predicts one part,
  held out (each class's share for a classifier), by a model of the other rows; the out-of-fold predictions then cover
  only the held out rows
* **Folds** (K-fold only): the number of parts, from 2 to 20 (5 by default)
* **Test fraction** (Hold-out only): the share of the rows held out, from 0.05 to 0.5 (0.2 by default)
* **Shuffle** (checked by default): shuffles the rows by the split seed before they are split. Unchecked, the parts are
  blocks of consecutive rows (with stratification, each class is dealt in row order) or the last rows are held out
  (the last of each class for a classifier), and **Split seed** is disabled
* **Split seed** (42 by default): the seed of the shuffle, so the same seed gives the same parts; with a search, also of
  its own hold out and of the settings it draws

In every scheme the final model is trained on all the rows. A training needs at least two rows per part, or two held
out rows.

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

Open a table, then **ML | Models | Apply Model...** (or use the play icon in the catalog or the run icon in the table's
**Models** pane). The **Apply predictive model**
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
the name is taken) to the table, tagged with the model's id (`predictive.model`). A model saved with **Predict
probability** also adds, right after it, the probability columns `P(<target> = <class>)` (`P(<target> = <class>) 2`
when taken), tagged the same way: for two classes the positive class's only (`P(SEX = F)`), for more one per class
(`P(Species = setosa)`, ...); hover a column's header for `Probability that <target> is <class>, predicted by <model
name>` (the column's description). A model saved without it adds none. A progress bar in the task bar shows the
batches; cancel it there to stop before the next batch. A balloon names the new columns (`Added the columns "SEX
(predicted)", "P(SEX = F)" to demog.`) and the
number of skipped rows. Every application, completed, failed or cancelled, is recorded with the model.

### Models for a table

Show a table in the context panel (for example, click a cell or an empty area of its grid): the pane **Models** lists
the models you can see whose features all have a close column in the table, newest first. Each line has the small green run icon (**Run prediction**: the **Apply predictive model** dialog with that model and
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
are the names and values `trainModel`'s `hyperparameters` accepts. A method's probability setting (SVM's
`probability`) is not listed: `predictProbability` sets it.

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
preparation), `metrics` (`train`, `validation`, for two classes `positiveClass` and, for a classifier, `classes`:
each class's `label`, `rows` and `train` / `validation` values of `precision`, `sensitivity`, `f1`, `auc`,
`averagePrecision`; `test` after a search that held rows out; null when none were saved), `search` (the search that
chose the settings, as `trainModel` returns it; null without one),
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

`Forge:trainModel` trains as the **Predictive model** view does, with its defaults (5-fold cross-validation, the split
seed 42):

* `table`, `target` (a column name): a numerical target makes a regression model, a text or boolean one a classifier
* `features`: column names; by default every numerical column except the target, bigint columns and integer columns
  without empty cells whose values are all different (row numbers, ids)
* `method`: `CatBoost`, `Decision Tree`, `Extra Trees`, `Hist Gradient Boosting`, `k-Nearest Neighbors`, `LDA`,
  `Linear Regression`, `Naive Bayes`, `Neural Network`, `PLS Regression`, `Random Forest`, `RealMLP`, `Softmax`,
  `SVM`, `TabM`, `XGBoost` (the installed methods), or `auto`: Find best model (see "Search"), every method that can
  learn from the data and runs in the browser tried with its defaults, then an equal share of the limit each (or,
  with `options.search.strategy: 'halving'`, the better half kept by successive halving down to the best three,
  those tuned), the best one trained (`options.search` defaults to `{budgetSeconds: 60}`; `missingValues` defaults to
  `skip`);
  by default the method the view suggests for the data
* `hyperparameters`: the method's settings by name; settings not named keep their defaults (not with `auto`: every
  method starts from its defaults)
* `options`: the preparation and the validation, each key optional:

  | Key | Default | Meaning |
  |---|---|---|
  | `missingValues` | `skip`; `keep` for a method that handles missing values (never with `auto`) | `skip` leaves out the rows with empty cells; `impute` fills the features' empty cells; `keep` gives the rows with empty cells to the method as they are, only for a method that handles missing values and not with method `auto` (rows without a target are always left out) |
  | `neighbors`, `distance` | the EDA imputer's (`4`, `Euclidean`) | the settings of `impute` |
  | `oneHot` | on when every text or yes/no feature has at most 20 categories (while `skipUniqueCategories` is on, features whose values are all different do not count) | One-hot encoding of text and yes/no features |
  | `skipUniqueCategories` | `true` | leaves out text features whose values are all different |
  | `predictProbability` | `false` | for a text or boolean target: the saved model adds the class probabilities when applied, and the labels come from them (for two classes, at `threshold`; for more, the most probable class). It stays on with `search` (the trials score those labels) and, for SVM, sets its `probability` setting |
  | `threshold` | optimized | two classes with `predictProbability`: the probability of the first class from which a row is predicted as it. Without it, the value that maximizes `thresholdMetric` on the validation rows (after a search, on the winner's), stored with the model as `thresholdMetric`; `0.5` when there is none (the validation rows hold one class). A given `threshold` is used as it is and no metric is stored. `cutoff` is accepted as its old name (`threshold` wins when both are given) |
  | `thresholdMetric` | `f1` | the metric the threshold is optimized for: `f1`, `balancedAccuracy`, `mcc` or `accuracy` (as **Best threshold by**); ignored when `threshold` (or `cutoff`) is given |
  | `scheme` | `kfold` | how Validation is measured: `kfold` (cross-validation) or `holdout` (one held out part, then the final model on all rows), as **Scheme** in "Validation" |
  | `folds` | `5` | `kfold`: the number of cross-validation parts, at least 2 |
  | `testFraction` | `0.2` | `holdout`: the share of the rows held out, from 0.05 to 0.5 |
  | `stratified` | `true` | each class keeps its share in every part or in the hold out; ignored for a regression |
  | `shuffle` | `true` | shuffles the rows by `seed` before the split; `false`: blocks of consecutive rows, or the last rows held out |
  | `seed` | `42` | the split seed, a whole number from 0 to 2147483647 (with `search`, also of the search's hold out and of the trials' settings) |
  | `search` | none (`{budgetSeconds: 60}` with `auto`) | `{metric, budgetSeconds, maxTrials, strategy, methods, topK}`, each key optional, but one limit, `budgetSeconds` or `maxTrials`, never both (neither: `budgetSeconds: 60`): with a method, tunes its settings first, as **Tune** does, with `auto` finds the best model, as **Find** does (see "Search"), then trains the chosen ones. `metric`: `r2` (regression default), `rmse`, `mae`, `balancedAccuracy` (classifier default), `accuracy`, `f1`, `mcc`, `auc`, `logLoss`, `brier` (these two minimized); `budgetSeconds`: a whole number of seconds, at least 1; `maxTrials`: starts no more than that many trials (of all methods; completed, pruned and failed ones alike); only with `auto`: `strategy`, `equal` (the default: each method an equal share of the limit, at least one trial each) or `halving` (successive halving, as **Leaders**), `methods`, the method names to try (by default every browser method that can learn from the data, without the ones lacking class probabilities while `predictProbability` is on or under `auc`, `logLoss`, `brier`, and at 60 seconds or less or with `maxTrials` only the ones that retrain by themselves, not SVM with `predictProbability`), and `topK`, how many methods successive halving keeps to tune (default `3`) |

* `name`: saves the model under this name; without a name the model is trained and measured but not saved.
  `description` and `tags` (comma-separated) go with it
* `storage` (`none` by default): what the saved model keeps of its training data, as in "Data storage": `none`,
  `reference` (the table must come from a file, a query or a script) or `copy` (uploads the feature and target
  columns). Nothing is uploaded unless you ask for `copy`
* `showProgress` (`true` by default): the task-bar progress `Training <method> model` (with `search`, first `Tuning
  <method>`, with `auto` first `Finding best model`); cancel it there to stop the training, which fails the call with `Training was cancelled.` (a search
  cancelled after a completed trial trains the best so far)

```javascript
// Tune XGBoost's settings for at most 30 seconds, then train and save the chosen ones
await grok.functions.call('Forge:trainModel', {table: grok.shell.t, target: 'Species', method: 'XGBoost',
  options: {seed: 7, search: {metric: 'accuracy', budgetSeconds: 30}}, name: 'Iris species tuned'});

// Find the best of three methods within 30 seconds; result.method is the one trained
await grok.functions.call('Forge:trainModel', {table: grok.shell.t, target: 'Species', method: 'auto',
  options: {search: {budgetSeconds: 30, methods: ['XGBoost', 'Decision Tree', 'Naive Bayes']}}});
```

The result: `modelId` (null when not saved), `runId`, `name`, `method`, `task`, `target`, `features`, `rowCount` (the
rows trained on; after a search with a hold out, without the held out rows), `skippedRows`, `seed`, `metrics` (`train`,
`validation`, `test` (after a search with 500 rows or more: the final model on the held out rows), each an object of
metric ids: `mse`, `rmse`, `mae`, `r2` for a regression; `accuracy`, `f1`, `balancedAccuracy`, `mcc`, for two classes
`sensitivity`, `specificity`, `precision`, `npv`, and, for a method with class probabilities, `auc`,
`averagePrecision`, `logLoss`, `brier`), for two classes
`positiveClass` and, for a classifier, `classes`, the per-class values as `getModel` gives them) and, after a search,
`search`: `mode` (`tune`, or `best` with `auto`; `method` is then the chosen leader's), `metric`, `budgetSeconds` or
`maxTrials` (the limit; a search continued in the view sums its runs'), `strategy` (`best`: `equal` or `halving`),
`seed`, `trials`, `completed`, `pruned`, `failed`, `elapsedMs`,
`stopped` (`budget`, `cancelled`, `exhausted`, `trials` or `failed`), `holdOut` (`rows`, `fraction`; absent below 500
rows), `space` (`searched`: the dimensions searched, each `{name, type, low?, high?, log?, step?, choices?, simpler?}` (a categorical dimension has no bounds);
`pinned`: the settings kept at their values, none from the API; Tune only), `methods` (`best` only: every method tried,
in cost order, `{method, trials, outcome, note?, eliminatedRound?}` with `outcome` `complete`, `failed` (`note`: its
first error), `abandoned` (only abandoned or pruned trials) or `unstarted`, and `eliminatedRound` the round of successive halving
after which it dropped out, from 0, absent for a method tuned to the end or one that left otherwise, and under an
equal split) and `leaders` (per method, its trial chosen by the one-standard-error rule: `method`, `hyperparameters`, `score`
(the metric by the chosen validation on the search rows), `se`, `foldScores`, `trials`, `chosen`, marking the one
trained: across methods the highest `score` − `se`, the lowest `score` + `se` for an error, a tie by the method name, capitals first).
Nothing is added to the table. Every training that runs is recorded as a training run with the source `api`, saved or
not, cancelled and failed ones too (after a search, the chosen settings' training, with `search`); a call refused for
its inputs, whose data preparation fails or whose search completes no trial records nothing. The errors a call gets, for example:

* `The table has no column 'Specie'.`, `The table has no columns 'a', 'b'.`
* `No method 'Boost' is installed. The methods are: <methods>.`
* `The method XGBoost has no setting 'speed'. Its settings are: <settings>.`, `The setting 'rate' of Softmax needs a
  number.`, `The setting 'iterations' of XGBoost must be at most 100.`
* `Unknown option 'cutof'. The options are: missingValues, neighbors, distance, oneHot, skipUniqueCategories,
  predictProbability, threshold, thresholdMetric, seed, scheme, folds, testFraction, stratified, shuffle, search.`
* `The option 'thresholdMetric' needs one of: f1, balancedAccuracy, mcc, accuracy.`
* `The option 'scheme' needs kfold or holdout.`, `The option 'testFraction' needs a number from 0.05 to 0.5.`,
  `The option 'shuffle' needs true or false.`, `Training needs at least 10 rows; the table has 8.` (two rows per
  part, or two held out rows)
* `The option 'search' needs an object with metric, budgetSeconds, maxTrials, strategy, methods, topK.`, `The search
  metric 'r3' is not one of: r2, rmse, mae, balancedAccuracy, accuracy, f1, mcc, auc, logLoss, brier.`, `The search
  metric 'r2' does not apply to a classification. Choose one of: balancedAccuracy, accuracy, f1, mcc, auc, logLoss,
  brier.`, `The option 'budgetSeconds' needs a whole number of at least 1.`, `The search options 'budgetSeconds' and
  'maxTrials' cannot be combined: choose one limit.`
* `Naive Bayes declares no search space, so it cannot be tuned.`, `No trial completed within the budget. Give the
  search more time.`, `No trial completed: <the first failed trial's error>` (for example `No trial completed: The
  metric AUC-ROC needs class probabilities; SVM gives none.`), `The search was cancelled before a trial completed.`
* with `auto`: `The search option 'methods' needs a list of method names.`, `The option 'topK' needs a whole number of
  at least 1.`, `The option 'missingValues' keep cannot be used with method auto. Use skip or impute.`, `The
  hyperparameters cannot be given with method auto: every method starts from its defaults.`, `Linear Regression cannot
  take part in the search: it cannot learn from this selection.` (`...: it runs on the server.` for a server method),
  `<Method> gives no class probabilities, so Predict probability cannot be used with it.`, `SVM cannot take part in
  the search: it gives no class probabilities, which AUC-ROC needs.`, `No method can take part in the search. Check
  the features and the target.`, `The search option 'strategy' needs equal or halving.`, `At least one trial per
  method: <N>.` (`equal` with fewer `maxTrials` than methods); with a method: `The search options 'strategy',
  'methods' and 'topK' need method auto.`
* `Predict probability needs a text or boolean target.`, `Chemprop gives no class probabilities, so Predict
  probability cannot be used with it.`, `Enter a number from 0 to 1.` (a two-class `threshold` outside 0..1)
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
a column of the wrong kind, fails the call with a message and no dialog. The call adds `<target> (predicted)` and, for
a model saved with Predict probability, its probability columns `P(<target> = <class>)`, as the dialog does. Missing
values are handled as at training:
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
(for two classes), `metrics` (`train` and `validation` as saved, `evaluation` on this data, with the ids of `trainModel`'s
result; a model saved with the former regression mode of Predict probability is measured on its scores, so its
`evaluation` has `auc` and `averagePrecision` but no `logLoss` or `brier`), `classes` (a classifier:
the evaluation's per-class values, each class's `label`, `rows` and `values`) and `insights` (the **Insights & Tips** of the evaluation and of the
saved metrics: each a `text`, with the numbers behind it in `detail` or `items`; empty when there is nothing to say,
without the panel's line that reads Balanced accuracy or R2). Nothing is added to the table and the
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
unless you ask for a copy). Models of the earlier predictive modeling tool are not Forge's; the assistant's Forge
knowledge does not cover them.

## Model catalog

Open the catalog from **Browse > Apps > Forge** or **Browse > Platform > Predictive models**. The view has two sections:

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
  first one), **Applications**, **Features**, **Target**, **Method**, **Task**, **Search** (for a model whose settings a
  search chose: `Tune (Balanced accuracy, limit 1 min, 120 trials)`, with a trial limit `Find best model (Balanced
  accuracy, limit 100 trials)`), **Applicable to** (the open tables that fit), and
  **Tags**, which you can edit right there
* **Performance**: the metrics saved with the model (the **Metric / Train / Validation** grid of the training, with
  **Test** after a search that held rows out, a value below its threshold or far from the training one in red, as in
  the results) and how it was validated (with the search bullets), then **Run
  Evaluation** (for a model trained and saved in Forge; its tooltip: "Apply the model to its training data, or to an
  open table, and show how well it predicts."). It measures the model again on data with known targets: the
  file or query the model keeps a reference to, or its uploaded copy; for a model that keeps neither, the **Evaluate
  model** dialog asks for **Table** (the open tables the model fits) and **Target** (the column with the actual
  values), and without a fitting open table a warning says to open one. While it runs, the task bar shows
  `Evaluating <model>` with a cancel. The view **Evaluation: <model>** then shows the results in rows as after
  training: the insights, the **Performance** panel with an **Evaluation** column next to **Train** and **Validation** and
  the bullet `Evaluation: <table> (<rows> rows)`, and the plots (**Per-class** shows the evaluation's values; a model
  saved before the per-class values has rows such as **MCC** with only an Evaluation value). Nothing is added to the table, a referenced or copied table does not
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
tooltip. A column predicted by a model of the earlier predictive modeling tool says so there, with the model's name:
`Predicted by "<name>", a model of the earlier predictive modeling tool. Forge cannot show it yet.`

## Share a model

Every saved model belongs to its author, who can see, apply, edit, delete and share it; nobody else sees it until it
is shared. The author shares a model with **Share...** in its **Sharing** pane, and the pane lists the new groups as
soon as the share is done. `Forge:shareModel` opens the same dialog from a script (see "From scripts and Datagrok
AI"). A user a model is shared with sees it in the catalog and applies it, and the application is recorded with the
model. Every user trains, saves and deletes models of their own the same way.

Training runs are not shared with their model: each user sees only their own. On a Datagrok 1.28.0 server an
administrator, too, sees only the models shared with them, not every user's (a platform issue).

To share a model through a Space, drag it onto the Space in **Browse**: its card (**Browse > Platform > Domains >
forge > model**, card mode), its name in a table's **Models** pane, or its name in the title of its context
panel. The platform asks to confirm and links the model to the Space; everyone with access to the Space can then see
and apply it. Dropping it again links nothing twice. The platform has no **Add to** command for models yet.

On a server with Forge, **ML | Models**, **Browse > Platform > Predictive models** and the table's **Models** pane are
Forge's. Models of the earlier predictive modeling tool are not shown; a column predicted by such a model says so in
**Predicted by**. Forge needs Datagrok **<minimum platform version: TBD>** or later: on an older server the earlier
tool stays in **ML | Models** and a table shows two **Models** panes.

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
script, a function marked `meta.mlserver: true`, or Chemprop) receives it, and takes no part in Find best model (see
"Search"). Chemprop and script
methods follow the same contract but were not available to test on this version's server. **Impute** uses the EDA package's `knnImpute` function; without it, **Missing values**
offers **Skip rows** only. Imputing a large table takes a while.

To add a method, register functions that share the same `meta.mlname` (the method name) and set
`meta.mlrole` on each of them:

| Role | Required | Purpose |
|---|---|---|
| `train` | yes | Trains a model on a table and a target column, returns the model |
| `apply` | yes | Applies a trained model to a table, returns the predictions: the first column; `proba_<class>` columns after it give the class probabilities (AUC-ROC, Log loss, Brier score, the curves, Predict probability) |
| `isApplicable` | yes | Tells whether the method can learn from the given data |
| `isInteractive` | no | Tells whether training is fast enough to rerun on every change |
| `visualize` | no | Shows method-specific views of a trained model |

The inputs of the `train` function, except the table and the target column, become the method's hyperparameters;
their default values and descriptions (shown as tooltips) come from the function's annotation. The optional key
`meta.mlenabled` says when a setting is used: a JSON object `{"<input>": <condition> | [<condition>, ...]}`, where a
condition is an object whose clauses all hold and a list holds when any of its conditions does. A clause is
`"<other input>": [<values>]` (the other input's value is one of them), `"<other input>": {"gt": <n>}` or
`{"lt": <n>}` (above / below a number; both may be given), or `"task": ["classification" | "regression"]` (the task
the method trains, the target's). An input whose rule fails is disabled in the view (it keeps its
value, is not checked, and is still passed to `train`, which must ignore it), and a search sets it to its default. A
malformed value counts as no rules (the view logs a warning once). EDA declares, for example, SVM
`{"gamma": {"kernel": ["RBF", "Polynomial", "Sigmoid"]}, "degree": {"kernel": ["Polynomial"]}, "coef0": {"kernel":
["Polynomial", "Sigmoid"]}, "epsilon": {"task": ["regression"]}}` and XGBoost `{"seed": [{"subsample": {"lt": 1}}, {"colsample": {"lt": 1}}]}`. Two optional keys of
the `train` function serve **Tune**:

* `meta.mlspace`: the search space, a JSON object over the `train` inputs:
  `{"<input>": {"type": "int" | "float" | "categorical", "low", "high", "log"?, "step"?, "choices"?, "simpler"?: "low" |
  "high" | "first"}}` (the bounds inside the input's min and max, the choices among its own; `[low, high]` and a list
  of choices are short forms); its ranges are the ones **Search options** starts from. A method without it cannot be
  tuned: its **Auto-tune** is disabled (`<Method> declares no search space.`) (EDA declares one for every method but Naive Bayes, which has no settings)
* `meta.mlparallel: 'true'`: the method accepts concurrent `train` and `apply` calls, so a search runs its fold fits in
  parallel (EDA's methods: a pool of Crux workers); without it they run one at a time

One more optional key of the `train` function serves **Predict probability**: `meta.mlprobabilities`, `'true'` when
`apply` returns the `proba_<class>` columns of a classifier, or the name of a bool `train` input that turns them on
(EDA's SVM: `probability`). Such an input is set from **Predict probability**, hidden from the settings, Search options,
the card's **Settings** tab and `Forge:listMethods`, and refused in `trainModel`'s `hyperparameters`. Without the key, **Predict
probability** is disabled for the method. EDA declares `'true'` on every classifier but SVM.

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
