import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {DimensionalityReducer} from '@datagrok-libraries/utils/src/reduce_dimensionality';
import {getSequenceMolecularWeight} from './molecular-measure';

export function peptideSimilaritySpace(
  table: DG.DataFrame,
  alignedSequencesColumn: DG.Column,
  method: string,
  measure: string,
  cyclesCount: number,
  activityColumnName: string,
) {
  const axesNames = ['~X', '~Y', 'MW'];

  const reducer = new DimensionalityReducer(
    alignedSequencesColumn.toList(),
    method,
    measure,
    {cycles: cyclesCount},
  );
  const embcols = reducer.transform(true);
  const columns = Array.from(embcols, (v: Float32Array, k) => (DG.Column.fromFloat32Array(axesNames[k], v)));

  const sequences = alignedSequencesColumn.toList();
  const mw: Float32Array = new Float32Array(sequences.length).fill(0);

  let currentSequence;
  for (let i = 0; i < sequences.length; ++i) {
    currentSequence = sequences[i];
    mw[i] = currentSequence == null ? 0 : getSequenceMolecularWeight(currentSequence);
  }

  columns.push(DG.Column.fromFloat32Array('MW', mw));

  const edf = DG.DataFrame.fromColumns(columns);

  // Add new axes.
  for (const axis of axesNames) {
    const col = table.col(axis);

    if (col == null) {
      table.columns.insert(edf.getCol(axis));
    } else {
      table.columns.replace(col, edf.getCol(axis));
    }
  }

  const view = (grok.shell.v as DG.TableView);

  const viewer = view.addViewer(DG.VIEWER.SCATTER_PLOT, {x: '~X', y: '~Y', color: activityColumnName, size: 'MW'});
  // (viewer as DG.ScatterPlotViewer).zoom(
  //   edf.getCol('~X').min,
  //   edf.getCol('~Y').min,
  //   edf.getCol('~X').max,
  //   edf.getCol('~Y').max,
  // );
  return viewer;
}
