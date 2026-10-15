/* eslint-disable max-len, max-params */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {runApplyModel} from './apply/apply-model';
import {runDeleteModel} from './storage/model-store';
import {runTrainModel, TrainModelResult} from './training/train-api';
import {openApplyDialog} from './ui/apply-model-dialog';
import {ForgeApp} from './ui/forge-app';
import {registerModelActions} from './ui/model-actions';
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
  if (ForgeModelHandler.registerOnce(registered))
    registerModelActions();
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
//input: string method {optional: true; description: The method as Forge names it (XGBoost, Linear Regression, PLS Regression, Softmax, SVM). By default the method Forge suggests for the data}
//input: map hyperparameters {optional: true; description: The method's settings by name, for example iterations: 50. Settings not named keep the method's defaults}
//input: map options {optional: true; description: Data preparation. missingValues: skip (default) or impute, with neighbors and distance for impute. oneHot (default on when every text feature has at most 20 categories), skipUniqueCategories (default true), predictProbability (default false, two-class targets only), cutoff (default 0.5), seed (default random), folds (default 5)}
//input: string name {optional: true; description: Saves the model under this name. Without a name the model is trained and measured but not saved}
//input: string description {optional: true; description: Description of the saved model}
//input: string tags {optional: true; description: Comma-separated tags of the saved model}
//input: string storage = 'none' {optional: true; choices: ["none", "reference", "copy"]; description: What the saved model keeps of its training data. none keeps a summary, reference a link to the table's source (the table must have one), copy uploads the feature and target columns to the server}
//input: bool showProgress = true {optional: true; description: Shows a cancellable progress in the task bar}
//output: object result {description: modelId (null when not saved), runId, name, method, task, target, features, rowCount, skippedRows, seed and metrics with train, validation and positiveClass}
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
