import * as DG from 'datagrok-api/dg';
import {featureSchemasOf, requiredFeaturesOf} from '../apply/apply-model';
import {isSuggested} from '../apply/column-matching';
import {ModelRow} from '../generated/db';

type FeatureRow = Pick<ModelRow, 'features' | 'options'>;

/** The [tables] in which every feature the model needs ({@link requiredFeaturesOf}) has a close column; none for a
 * model without a feature list. */
export function applicableTables(row: FeatureRow, tables: DG.DataFrame[]): DG.DataFrame[] {
  const required = requiredFeaturesOf(row);
  return required === null ? [] : tables.filter((table) => isSuggested(required, table));
}

/** The names of every feature of the model and its {@link applicableTables}. */
export function featureFit(row: FeatureRow, tables: DG.DataFrame[]): {names: string[]; tables: DG.DataFrame[]} {
  return {names: (featureSchemasOf(row.features) ?? []).map((f) => f.name), tables: applicableTables(row, tables)};
}
