import * as DG from 'datagrok-api/dg';
import {ApplyModelRow} from '../apply/apply-model';
import {textColumn} from '../columns';
import {Engine} from '../engines/engine';
import {ForgeError} from '../forge-error';

/** The models the method [engine] provides, as catalog rows. Its `meta.mlrole: models` function returns one row per
 * model: `id` (what its `apply` gets as the model, as UTF-8 bytes), `name`, `target`, `features` (the json list
 * `[{name, type, semType?}]`) and, optionally, the date `created`. */
export async function providedModels(engine: Engine): Promise<DG.DataFrame> {
  const func = engine.functions.models;
  if (func === undefined)
    throw new ForgeError(`The method '${engine.name}' provides no models.`);
  const listed: unknown = await func.apply({});
  if (!(listed instanceof DG.DataFrame))
    throw new ForgeError(`The method '${engine.name}' returned no model list.`);
  const rows = Array.from({length: listed.rowCount}, (_, i) => i);
  const texts = (name: string) => rows.map((i) => String(listed.col(name)?.get(i) ?? ''));
  const created = listed.col('created');
  const frame = DG.DataFrame.fromColumns([
    textColumn('id', texts('id')),
    textColumn('name', texts('name')),
    textColumn('engine_name', rows.map(() => engine.name)),
    textColumn('target_name', texts('target')),
    textColumn('features', texts('features').map(featuresText)),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'created_on',
      rows.map((i) => created?.type === DG.COLUMN_TYPE.DATE_TIME ? created.get(i) : null)),
  ]);
  frame.name = `${engine.name} models`;
  return frame;
}

/** A method's json feature list as a model row's `features`; a list that is not json gives a model without features,
 * which is shown and refused when applied. */
function featuresText(list: string): string {
  try {
    const columns: unknown = JSON.parse(list);
    return JSON.stringify(Array.isArray(columns) ? {columns} : {});
  } catch {
    return '{}';
  }
}

/** Row [i] of {@link providedModels} as a model row; its empty `options` tell the card it has the feature list. */
export function providedRow(models: DG.DataFrame, engine: Engine, i: number): ApplyModelRow {
  return {id: models.get('id', i), name: models.get('name', i), engine_name: engine.name,
    engine_namespace: engine.namespace, target_name: models.get('target_name', i),
    features: JSON.parse(models.get('features', i)), options: {}, created_on: models.get('created_on', i)};
}
