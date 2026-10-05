/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {runApplyModel} from './apply/apply-model';
import {openApplyDialog} from './ui/apply-model-dialog';
import {ForgeApp} from './ui/forge-app';
import {TrainView} from './ui/train-view';
export * from './package.g';

export const _package = new DG.Package();

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
