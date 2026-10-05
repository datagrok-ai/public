import * as DG from 'datagrok-api/dg';
import {sharedFrame} from '../preparation/shared-frame';
import {ColumnSchema} from '../training/train-model';
import {ColumnMapping, mappedColumns} from './column-matching';

/** The mapped columns in training order, shared with [table]; a column named unlike its feature is a renamed copy.
 * Release the frame with `releaseFrame` after use. */
export function featureFrame(table: DG.DataFrame, features: ColumnSchema[], mapping: ColumnMapping): DG.DataFrame {
  const columns: DG.Column[] = [];
  for (const feature of features) {
    for (const [name, column] of mappedColumns(feature, mapping)) {
      const col = table.getCol(column);
      columns.push(col.name === name ? col : renamedCopy(col, name));
    }
  }
  return sharedFrame(columns);
}

// Renaming a shared column would rename it in the user's table too.
function renamedCopy(col: DG.Column, name: string): DG.Column {
  const copy = col.clone();
  copy.name = name;
  return copy;
}
