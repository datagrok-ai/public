/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();
export * from './package.g';

const MAX_MODELS = 1000;
/** MLflow signature types as the platform's column types; the other types have the platform's names. */
const COLUMN_TYPES: {[type: string]: string} = {float: 'double', long: 'int', integer: 'int', boolean: 'bool'};

export class PackageFunctions {
  @grok.decorators.func()
  static info() : void {
    grok.shell.info(_package.webRoot);
  }

  @grok.decorators.func({
    'meta': {'mlname': 'MLflow', 'mlrole': 'models'},
    'name': 'MLflow models',
    'description': 'Lists the MLflow models synced from the MLFlow connections: id, name, target, features, created',
  })
  static async mlflowModels(): Promise<DG.DataFrame> {
    const models = await grok.dapi.models.filter('source = "MLFlow"').list({pageSize: MAX_MODELS});
    const features = (m: DG.Model) =>
      JSON.stringify(m.inputs.map((c) => ({name: c.name, type: COLUMN_TYPES[c.type] ?? c.type})));
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('id', models.map((m) => m.nqName)),
      DG.Column.fromStrings('name', models.map((m) => m.name)),
      DG.Column.fromStrings('target', models.map((m) => m.outputs[0]?.name ?? '')),
      DG.Column.fromStrings('features', models.map(features)),
      DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'created', models.map((m) => m.createdOn)),
    ]);
  }

  @grok.decorators.func({
    'meta': {'mlname': 'MLflow', 'mlrole': 'apply'},
    'name': 'Apply MLflow',
    'description': 'Predicts with an MLflow model in its container, which the first application starts',
  })
  static async applyMlflow(
    @grok.decorators.param({'options': {'description': 'The model\'s input columns, in its order.'}})
      df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'The model id from MLflow models, as UTF-8 bytes.'}})
      model: any): Promise<DG.DataFrame> {
    const table = df.clone();
    const inputs = table.columns.names();
    await grok.ml.applyModel(new TextDecoder().decode(model as Uint8Array), table,
      Object.fromEntries(inputs.map((name) => [name, name])), false);
    for (const name of inputs)
      table.columns.remove(name);
    return table;
  }
}
