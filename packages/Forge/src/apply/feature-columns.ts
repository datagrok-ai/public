import * as DG from 'datagrok-api/dg';
import {ColumnSchema} from '../training/train-model';
import {ColumnMapping, mappedColumns} from './column-matching';

/** The mapped columns of [table] in training order, the table's own; a column named unlike its feature is a renamed
 * copy. */
export function featureColumns(table: DG.DataFrame, features: ColumnSchema[], mapping: ColumnMapping): DG.Column[] {
  const columns: DG.Column[] = [];
  for (const feature of features) {
    for (const [name, column] of mappedColumns(feature, mapping)) {
      const col = table.getCol(column);
      columns.push(col.name === name ? col : renamedCopy(col, name));
    }
  }
  return columns;
}

// Renaming a shared column would rename it in the user's table too.
function renamedCopy(col: DG.Column, name: string): DG.Column {
  const copy = col.clone();
  copy.name = name;
  return copy;
}
