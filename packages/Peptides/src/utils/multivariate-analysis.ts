import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';

export async function callMVA(
  tableGrid: DG.Grid,
  view: DG.View,
  currentDf: DG.DataFrame,
  options: {[name: string]: string},
  sequencesCol: DG.Column,
) {
  const activityCol = await _scaleColumn(currentDf.getCol(options['activityColumnName']), options['scaling']);
  const encDf = _encodeSequences(sequencesCol);
  const scaledColName = `${options['activityColumnName']}scaled`;

  _insertColumns(
    currentDf,
    [DG.Column.fromList('double', scaledColName, activityCol.toList())],
  );
  _insertColumns(currentDf, encDf.columns);

  const res = await grok.functions.call('MultivariateAnalysis', {
    table: currentDf,
    features: encDf.columns.names(),
    prediction: scaledColName,
    components: 10,
    showScores: true,
    showRegresCoefs: true,
  });
  console.log(res);
}

/**
 * Encodes a series of sequences into a certain scale.
 *
 * @param {string[]} sequencesCol Column containing the sequences.
 * @return {DG.DataFrame} The data frame with seqences encoded.
 */
function _encodeSequences(sequencesCol: DG.Column): DG.DataFrame {
  const nRows = sequencesCol.length;
  const nCols = AlignedSequenceEncoder.clean(sequencesCol.get(0)).length;
  const enc = new AlignedSequenceEncoder('WimleyWhite');
  const positions = new Array(nCols).fill(0).map((_) => new Float32Array(nRows));

  for (let j = 0; j < nRows; ++j) {
    const s = AlignedSequenceEncoder.clean(sequencesCol.get(j));
    for (let i = 0; i < nCols; ++i) {
      positions[i][j] = enc.encodeLettter(s[i]);
    }
  }
  const df = DG.DataFrame.fromColumns(positions.map(
    (v, i) => DG.Column.fromFloat32Array((i+1).toString(), v),
  ));
  return df;
}

async function _scaleColumn(column: DG.Column, method: string): Promise<DG.Column> {
  if (method == 'none') {
    return column;
  }

  const formula = (method.startsWith('-') ? '0-' : '')+'Log10(${'+column.name+'})';
  const newCol = await column.applyFormula(formula);

  if (newCol == null) {
    throw new Error('Column formula returned unexpected null.');
  }
  return newCol!;
}

function _insertColumns(targetDf: DG.DataFrame, columns: DG.Column[]): DG.DataFrame {
  for (const col of columns) {
    targetDf.columns.add(col);
  }
  return targetDf;
}

