/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
//import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';
import {assert, transposeMatrix} from '@datagrok-libraries/utils/src/operations';
import {Vector, Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {kendallsTau} from '@datagrok-libraries/statistics/src/correlation-coefficient';

/**
 * Converts a Matrix into a DataFrame.
 *
 * @export
 * @param {Matrix} matrix A matrix.
 * @return {DG.DataFrame} The data frame.
 */
export function matrix2DataFrame(matrix: Matrix): DG.DataFrame {
  return DG.DataFrame.fromColumns(matrix.map((v, i) => DG.Column.fromFloat32Array(`${i+1}`, v)));
}

/**
 * Encodes amino acid sequences into a numeric representation.
 *
 * @param {DG.Column} col A column containing the sequences.
 * @return {DG.DataFrame} The resulting data frame.
 */
function calcPositions(col: DG.Column): DG.DataFrame {
  const sequences = col.toList().map((v, _) => AlignedSequenceEncoder.clean(v));
  const enc = new AlignedSequenceEncoder();
  const encSeqs = sequences.map((v) => Vector.from(enc.encode(v)));
  const positions = transposeMatrix(encSeqs);
  return matrix2DataFrame(positions);
}

/**
 * Unfolds a data frame into <category>-<value> format.
 *
 * @param {DG.DataFrame} df A data frame to unfold.
 * @return {DG.DataFrame} The resulting data frame.
 */
function melt(df: DG.DataFrame): DG.DataFrame {
  let keys: string[] = [];
  const values: Float32Array = new Float32Array(df.columns.length*df.rowCount);
  let i = 0;

  for (const c of df.columns.toList()) {
    keys = keys.concat(Array<string>(c.length).fill(c.name));
    values.set(c.getRawData(), i);
    i += df.rowCount;
  }
  assert(keys.length == values.length);
  return DG.DataFrame.fromColumns([DG.Column.fromStrings('keys', keys), DG.Column.fromFloat32Array('values', values)]);
}

/**
 * Calculates Spearman's rho rank correlation coefficient.
 *
 * @param {DG.DataFrame} df A data frame to process.
 * @return {DG.DataFrame} The correlation matrix.
 */
// eslint-disable-next-line no-unused-vars
function calcSpearmanRhoMatrix(df: DG.DataFrame): DG.DataFrame {
  const nItems = df.columns.length;
  const rho = new Array(nItems).fill(0).map((_) => new Float32Array(nItems).fill(0));

  for (let i = 0; i < nItems; ++i) {
    for (let j = i+1; j < nItems; ++j) {
      rho[i][j] = df.columns.byIndex(i).stats.spearmanCorr(df.columns.byIndex(j));
      rho[j][i] = rho[i][j];
    }
  }
  return matrix2DataFrame(rho);
}

/**
 * Calculates Kendall's tau rank correlation coefficient.
 *
 * @param {DG.DataFrame} df A data frame to process.
 * @param {number} [alpha=0.05] The significance threshold.
 * @return {DG.DataFrame} The correlation matrix.
 */
function calcKendallTauMatrix(df: DG.DataFrame, alpha: number = 0.05): DG.DataFrame {
  const nItems = df.columns.length;
  const tau = new Array(nItems).fill(0).map((_) => new Float32Array(nItems).fill(0));

  for (let i = 0; i < nItems; ++i) {
    for (let j = i+1; j < nItems; ++j) {
      const res = kendallsTau(df.columns.byIndex(i).getRawData(), df.columns.byIndex(j).getRawData());
      tau[i][j] = res.prob < alpha ? res.test : 0;
      tau[j][i] = tau[i][j];
    }
  }
  return matrix2DataFrame(tau);
}

/**
 * Creates acorrelation plot and a box plot to perform correlation analysis.
 *
 * @export
 * @param {DG.Column} sequencesColumn A column containing amino acid sequences.
 * @return {[DG.Viewer, DG.Viewer]} These two plots.
 */
export function correlationAnalysisPlots(sequencesColumn: DG.Column): [DG.Viewer, DG.Viewer] {
  const posDF = calcPositions(sequencesColumn);
  const cpviewer = DG.Viewer.fromType(
    DG.VIEWER.CORR_PLOT,
    posDF,
    {
      'xColumnNames': posDF.columns.names(),
      'yColumnNames': posDF.columns.names(),
      'correlationType': 'Spearman',
    });

  const rhoDF = calcKendallTauMatrix(posDF);
  const meltDF = melt(rhoDF);

  const bpviewer = DG.Viewer.fromType(
    DG.VIEWER.BOX_PLOT,
    meltDF, {
      'categoryColumnName': 'keys',
      'valueColumnName': 'values',
      'statistics': ['min', 'max', 'avg', 'med'],
    });
  return [cpviewer, bpviewer];
}
