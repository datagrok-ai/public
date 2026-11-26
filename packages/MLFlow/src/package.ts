/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {columnsOf, ListedModel, syncDiff, targetOf} from './sync';

export const _package = new DG.Package();
export * from './package.g';

const ENGINE = 'MLflow';
const MAX_MODELS = 1000000;

/** The rows of [frame] as objects of its [columns]' texts, keyed by [columns]' values. */
function rowsOf<T extends string>(frame: DG.DataFrame, columns: {[column: string]: T}): {[key in T]: string}[] {
  return Array.from({length: frame.rowCount}, (_, i) => Object.fromEntries(Object.entries(columns).map(
    ([column, key]) => [key, String(frame.col(column)?.get(i) ?? '')])) as {[key in T]: string});
}

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
    const sync = async (): Promise<string> => {
      const listed: ListedModel[] = rowsOf(await grok.functions.call('MlflowModels', {}),
        {connection_id: 'connectionId', connection: 'connection', name: 'name', version: 'version',
          description: 'description', inputs: 'inputs', outputs: 'outputs', error: 'error'});
      const own = rowsOf(await grok.functions.call('Forge:listModels', {method: ENGINE, own: true, limit: MAX_MODELS}),
        {Id: 'id', Name: 'name', Description: 'description', Target: 'target', Features: 'features',
          Options: 'options'});
      const {register, remove} = syncDiff(listed, own.map((row) => ({...row, options: JSON.parse(row.options)})));
      // One refused write does not stop the others.
      const failed: string[] = [];
      const write = (name: string, call: Promise<unknown>) => call.catch(() => failed.push(name));
      for (const {model, row} of register) {
        const target = targetOf(model);
        await write(model.name, grok.functions.call('Forge:registerModel', {id: row?.id, name: row?.name ?? model.name,
          description: row?.description ?? model.description, engine: ENGINE,
          task: ['string', 'bool'].includes(target.type) ? 'classification' : 'regression', target,
          features: {columns: columnsOf(model.inputs)},
          options: {connectionId: model.connectionId, name: model.name, version: model.version}}));
      }
      for (const row of remove)
        await write(row.name, grok.functions.call('Forge:deleteModel', {model: row.id}));
      return [...listed.filter((m) => m.error !== '' && m.name === '')
        .map((m) => `${m.connection} was not refreshed: ${m.error}`), ...(failed.length === 0 ? [] :
        [`${failed.length} MLflow ${failed.length === 1 ? 'model was' : 'models were'} not updated (${failed[0]}${
          failed.length === 1 ? '' : ', ...'}).`])].join(' ');
    };
    // One sync per user at a time, so two tabs do not write the same model twice; without Web Locks (a page that is
    // not a secure context) it runs unlocked.
    return navigator.locks === undefined ? sync() :
      navigator.locks.request(`mlflow-sync-${DG.User.current().id}`, sync);
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
    const {connectionId, name} = JSON.parse(new TextDecoder().decode(model));
    return grok.functions.call('MlflowApply', {connectionId, name, table: df});
  }
}
