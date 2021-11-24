import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {getSequenceMolecularWeight} from './molecular-measure';
import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';
import {Coordinates} from '@datagrok-libraries/utils/src/type_declarations';

/**
 * Creates a worker to perform a dimensionality reduction.
 *
 * @param {any[]} columnData Samples to process.
 * @param {string} method Embedding method.
 * @param {string} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @return {*} Promise.
 */
function createDimensinalityReducingWorker(columnData: any[], method: string, measure: string, cyclesCount: number) {
  return new Promise(function(resolve) {
    const worker = new Worker(new URL('../workers/dimensionality-reducer.ts', import.meta.url));
    worker.postMessage({
      columnData: columnData,
      method: method,
      measure: measure,
      cyclesCount: cyclesCount,
    });
    worker.onmessage = ({data: {embedding}}) => {
      resolve(embedding);
    };
  });
}

/**
 * Creates scatter plot with sequences embeded.
 *
 * @export
 * @param {DG.DataFrame} table The table containing samples.
 * @param {DG.Column} alignedSequencesColumn The samples column.
 * @param {string} method Embedding method.
 * @param {string} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @param {string} activityColumnName A columns containing an activity to assign it to points radius.
 * @return {*} The scatter plot viewer.
 */
export async function peptideSimilaritySpace(
  table: DG.DataFrame,
  alignedSequencesColumn: DG.Column,
  method: string,
  measure: string,
  cyclesCount: number,
  activityColumnName: string,
) {
  const axesNames = ['~X', '~Y', 'MW'];
  let columnData = alignedSequencesColumn.toList();

  columnData = columnData.map((v, _) => AlignedSequenceEncoder.clean(v));

  const embcols = await createDimensinalityReducingWorker(columnData, method, measure, cyclesCount);
  const columns = Array.from(
    embcols as Coordinates,
    (v: Float32Array, k) => (DG.Column.fromFloat32Array(axesNames[k], v)),
  );

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
  return viewer;
}
