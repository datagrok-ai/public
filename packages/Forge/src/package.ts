/* eslint-disable max-len, max-params */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {runApplyModel} from './apply/apply-model';
import {ModelDetails, runGetModel} from './catalog/get-model';
import {runListModels} from './catalog/list-models';
import {runUpdateModel, UpdatedModel} from './catalog/model-edit';
import {runDeleteModel} from './storage/model-store';
import {EvaluateModelResult, runEvaluateModel} from './training/evaluate-api';
import {runListMethods, runTrainModel, TrainModelResult} from './training/train-api';
import {openApplyDialog} from './ui/apply-model-dialog';
import {ForgeApp} from './ui/forge-app';
import {openShareDialog, registerModelActions} from './ui/model-actions';
import {ModelComparisonHandler} from './ui/model-comparison';
import {ForgeModelHandler} from './ui/model-handler';
import {isForgePrediction, predictedByPane} from './ui/prediction-column-panel';
import {tableModelsPane} from './ui/table-models-panel';
import {TrainView} from './ui/train-view';
export * from './package.g';

export const _package = new DG.Package();

//name: _initForge
//description: Registers the Forge model handler, the model commands and the model comparison handler
//tags: autostart
//meta.autostartImmediate: true
export function _initForge(): void {
  const registered = DG.ObjectHandler.list().map((h) => h.name);
  if (ForgeModelHandler.registerOnce(registered)) {
    registerModelActions();
    ForgeModelHandler.followChanges();
  }
  ModelComparisonHandler.registerOnce(registered);
}

//name: Forge
//description: Predictive modeling: methods and the model catalog
//tags: app
//output: view result
export async function forgeApp(): Promise<DG.ViewBase> {
  return ForgeApp.create();
}

//name: forgeModels
//description: Opens the Forge app
//top-menu: ML | Forge | Models
export async function forgeModels(): Promise<void> {
  await ForgeApp.open();
}

//name: forgeTrain
//description: Trains a predictive model on the current table
//top-menu: ML | Forge | Train...
export async function forgeTrain(): Promise<void> {
  await TrainView.open();
}

//name: forgeApply
//description: Applies a saved Forge model to the current table
//top-menu: ML | Forge | Apply...
export async function forgeApply(): Promise<void> {
  await openApplyDialog(grok.shell.currentTable);
}

//name: applyModel
//description: Applies a saved Forge model to a table and adds the prediction column
//input: string model {description: Model id or name}
//input: dataframe table
//input: map columnNamesMap {optional: true; description: Model feature name -> table column name}
//input: bool showProgress = true {optional: true}
//output: dataframe result
export async function applyModel(model: string, table: DG.DataFrame, columnNamesMap: {[feature: string]: string} | null,
  showProgress: boolean): Promise<DG.DataFrame> {
  return runApplyModel(model, table, columnNamesMap, showProgress);
}

//name: trainModel
//description: Trains a predictive model on a table with 5-fold cross-validation by default, records the training run and, when a name is given, saves the model to the Forge catalog. Returns the metrics and the saved model's id
//input: dataframe table {description: The table to learn from}
//input: string target {description: Name of the column to predict. A numerical column gives a regression model, a text or boolean column a classifier}
//input: list<string> features {optional: true; description: Names of the columns to learn from. By default every numerical column except the target and id-like whole-number columns}
//input: string method {optional: true; description: The method as Forge names it, any method Forge:listMethods lists. By default the method Forge suggests for the data}
//input: map hyperparameters {optional: true; description: The method's settings by name, for example iterations: 50. Settings not named keep the method's defaults}
//input: map options {optional: true; description: Data preparation. missingValues: skip, impute or keep (default keep for a method that handles missing values, otherwise skip), with neighbors and distance for impute. oneHot (default on when every text feature has at most 20 categories), skipUniqueCategories (default true), predictProbability (default false, two-class targets only), cutoff (default 0.5), seed (default random), folds (default 5)}
//input: string name {optional: true; description: Saves the model under this name. Without a name the model is trained and measured but not saved}
//input: string description {optional: true; description: Description of the saved model}
//input: string tags {optional: true; description: Comma-separated tags of the saved model}
//input: string storage = 'none' {optional: true; choices: ["none", "reference", "copy"]; description: What the saved model keeps of its training data. none keeps a summary, reference a link to the table's source (the table must have one), copy uploads the feature and target columns to the server}
//input: bool showProgress = true {optional: true; description: Shows a cancellable progress in the task bar}
//output: object result {description: modelId (null when not saved), runId, name, method, task, target, features, rowCount, skippedRows, seed and metrics with train, validation, positiveClass and, for a classifier, classes (per-class values)}
export async function trainModel(table: DG.DataFrame, target: string, features: string[] | null, method: string | null,
  hyperparameters: {[name: string]: unknown} | null, options: {[name: string]: unknown} | null, name: string | null,
  description: string | null, tags: string | null, storage: string | null,
  showProgress: boolean | null): Promise<TrainModelResult> {
  return runTrainModel({table, target, features, method, hyperparameters, options, name, description, tags, storage,
    showProgress});
}

//name: deleteModel
//description: Deletes a Forge model with its model file, its uploaded training copy and its application records. The model's training runs keep their history without the model
//input: string model {description: Model id, or a name that only one model you can see has}
//output: string id {description: The id of the deleted model}
export async function deleteModel(model: string): Promise<string> {
  return runDeleteModel(model);
}

//name: listMethods
//description: Lists the installed machine learning methods with their package, type, where they run, roles and settings (hyperparameters) with their defaults and ranges. With a table, marks the methods that can learn from it and the one Forge suggests
//input: dataframe table {optional: true; description: A table to check the methods against, with the view's default features}
//input: string target {optional: true; description: Name of the column to predict in the table. By default the table's last column}
//output: dataframe result {description: One row per method, Method, Package, Method type, Runs on, Keeps training rows, Retrains live, Roles, Hyperparameters, Settings (json text, name, type, default, min, max, choices, description per setting) and, with a table, Applicable and Suggested}
export async function listMethods(table: DG.DataFrame | null, target: string | null): Promise<DG.DataFrame> {
  return runListMethods(table, target);
}

//name: listModels
//description: Lists the saved Forge models you can see, newest first, with their validation metrics. Filters by a part of the name, the target, the task, the method, tags and the table the model can be applied to
//input: string name {optional: true; description: A part of the model name, case-insensitive}
//input: string target {optional: true; description: Name of the column the model predicts}
//input: string task {optional: true; description: regression or classification}
//input: string method {optional: true; description: The method as Forge names it, for example XGBoost}
//input: string tags {optional: true; description: Comma-separated tags. A model is listed when it has every one of them}
//input: dataframe table {optional: true; description: Only the models whose features all have a column of the same or a close name in this table. Checked in the browser, nothing is uploaded}
//input: int limit = 100 {optional: true; description: At most this many models, newest first}
//output: dataframe result {description: One row per model, Id, Name, Description, Method, Task, Target, Features, Training rows, Data storage, Tags, Author, Created and a <Metric> (validation) column for every metric a listed model has}
export async function listModels(name: string | null, target: string | null, task: string | null,
  method: string | null, tags: string | null, table: DG.DataFrame | null, limit: number | null): Promise<DG.DataFrame> {
  return runListModels({name, target, task, method, tags, table, limit});
}

//name: getModel
//description: Everything Forge keeps about a saved model, as one object. Use it to explain a model, to check what a table needs before applying it, and to compare models
//input: string model {description: Model id, or a name that only one model you can see has}
//output: object result {description: id, version, name, description, tags, method, package, task, target, targetType, classes, features, featureTypes, hyperparameters, preparation, metrics (train, validation, positiveClass, classes), findings, seed, splitting, storage (mode, source, tableId), dataset (name, rows, summary), hasModelFile, keepsTrainingRows, author, created, updated, lastRun (when, status) or null, applications and applicableTo (the open tables that fit, a column of a close name counts)}
export async function getModel(model: string): Promise<ModelDetails> {
  return runGetModel(model, grok.shell.tables);
}

//name: updateModel
//description: Renames a saved Forge model or changes its description or tags. A field not given stays as it is
//input: string model {description: Model id, or a name that only one model you can see has}
//input: string name {optional: true; description: The new name, not blank}
//input: string description {optional: true; description: The new description}
//input: string tags {optional: true; description: Comma-separated tags that replace the current ones. An empty text removes them}
//output: object result {description: id, name, description, tags and version after the change}
export async function updateModel(model: string, name: string | null, description: string | null,
  tags: string | null): Promise<UpdatedModel> {
  return runUpdateModel({model, name, description, tags});
}

//name: evaluateModel
//description: Measures a saved model on data with known target values, as Run Evaluation does, and returns its metrics next to the stored train and validation ones. Without a table it uses the training data the model keeps a link to or a copy of. Adds no column and records nothing
//input: string model {description: Model id, or a name that only one model you can see has}
//input: dataframe table {optional: true; description: The table to evaluate on. By default the model's referenced or copied training data}
//input: string target {optional: true; description: The column with the actual values. By default the column named as the model's target}
//input: map columnNamesMap {optional: true; description: Model feature name -> table column name, for columns named differently}
//input: bool showProgress = true {optional: true; description: Shows a cancellable progress Evaluating <model> in the task bar}
//output: object result {description: id, name, table, target, rowCount, skippedRows, positiveClass, metrics with train, validation and evaluation, classes (the evaluation's per-class values) and insights}
export async function evaluateModel(model: string, table: DG.DataFrame | null, target: string | null,
  columnNamesMap: {[feature: string]: string} | null, showProgress: boolean | null): Promise<EvaluateModelResult> {
  return runEvaluateModel({model, table, target, columnNamesMap, showProgress});
}

//name: shareModel
//description: Opens the platform's Share dialog for a saved Forge model, the same as Share... in its Sharing pane. Returns as soon as the dialog is shown
//input: string model {description: Model id, or a name that only one model you can see has}
export async function shareModel(model: string): Promise<void> {
  await openShareDialog(model);
}

//name: isPredictionColumn
//description: Whether the column holds the predictions of a Forge model
//input: column col
//output: bool result
export function isPredictionColumn(col: DG.Column): boolean {
  return isForgePrediction(col);
}

//name: Predicted by
//description: The Forge model that predicted the column
//input: column col
//output: widget result
//meta.role: panel
//condition: Forge:isPredictionColumn(col)
export function predictedByPanel(col: DG.Column): DG.Widget {
  return new DG.Widget(predictedByPane(col));
}

//name: Forge models
//description: The Forge models that fit the table, and training a new one
//input: dataframe table
//output: widget result
//meta.role: panel
export function tableModelsPanel(table: DG.DataFrame): DG.Widget {
  return new DG.Widget(tableModelsPane(table));
}
