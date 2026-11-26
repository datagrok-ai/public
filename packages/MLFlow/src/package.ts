/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {syncModels} from './sync';
import {ModelReference} from './types';

export const _package = new DG.Package();
export * from './package.g';

export class PackageFunctions {
  @grok.decorators.func()
  static info() : void {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func({
    'meta': {'mlname': 'MLflow', 'mlrole': 'sync'},
    'name': 'Sync MLflow models',
    'description': 'Writes the MLflow models you may use into the Forge catalog; returns what was not refreshed',
  })
  static async syncMlflow(): Promise<string> {
    if (DG.Func.find({package: 'Forge', name: 'registerModel'}).length === 0)
      return '';
    // One sync per user at a time, so two tabs do not write a model twice; a page without Web Locks runs unlocked.
    return navigator.locks === undefined ? syncModels() :
      navigator.locks.request(`mlflow-sync-${DG.User.current().id}`, syncModels);
  }

  @grok.decorators.func({
    'meta': {'mlname': 'MLflow', 'mlrole': 'apply'},
    'name': 'Apply MLflow',
    'description': 'Predicts with a published MLflow model in its container, which the first application starts',
  })
  static async applyMlflow(
    @grok.decorators.param({'options': {'description': 'The model\'s input columns, in its order.'}})
      df: DG.DataFrame,
    @grok.decorators.param({'type': 'dynamic',
      'options': {'description': 'The model reference as UTF-8 json bytes.'}})
      model: Uint8Array): Promise<DG.DataFrame> {
    const {connectionId, name, type}: ModelReference = JSON.parse(new TextDecoder().decode(model));
    const result: DG.DataFrame = await grok.functions.call('mlflowApply', {connectionId, name, table: df});
    // Core types a column by its values: whole numbers of a double output come as int.
    const prediction = result.columns.byIndex(0);
    if (type && prediction.type !== type)
      result.columns.replace(prediction, prediction.convertTo(type));
    return result;
  }
}
