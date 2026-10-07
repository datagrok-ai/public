import * as DG from 'datagrok-api/dg';
import {featureSchemasOf} from '../apply/apply-model';
import {isSuggested} from '../apply/column-matching';
import {ModelRow} from '../generated/db';

/** The [tables] in which every feature of the model has a close column; none for a model without a feature list. */
export function applicableTables(row: Pick<ModelRow, 'features'>, tables: DG.DataFrame[]): DG.DataFrame[] {
  return featureFit(row, tables).tables;
}

/** The model's feature names and its {@link applicableTables}, from one read of its feature list. */
export function featureFit(row: Pick<ModelRow, 'features'>, tables: DG.DataFrame[]):
  {names: string[]; tables: DG.DataFrame[]} {
  const features = featureSchemasOf(row.features);
  return features === null ? {names: [], tables: []} :
    {names: features.map((f) => f.name), tables: tables.filter((table) => isSuggested(features, table))};
}
