/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
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
function calcPositions(col: DG.Column): [DG.DataFrame, string[]] {
  const sequences = col.toList().map((v, _) => AlignedSequenceEncoder.clean(v));
  const enc = new AlignedSequenceEncoder('WimleyWhite');
  const encSeqs = sequences.map((v) => Vector.from(enc.encode(v)));
  console.log(encSeqs);
  const positions = transposeMatrix(encSeqs);
  return [matrix2DataFrame(positions), sequences];
}

/**
 * Unfolds a data frame into <category>-<value> format.
 *
 * @param {DG.DataFrame} df A data frame to unfold.
 * @return {DG.DataFrame} The resulting data frame.
 */

/**
 * Unfolds a data frame into <category>-<value> format.
 *
 * @param {DG.DataFrame} df A data frame to unfold.
 * @param {string} [keysName='Keys'] A name for keys column.
 * @param {string} [valuesName='Values'] A name for values column.
 * @return {DG.DataFrame} The resulting data frame.
 */
function melt(df: DG.DataFrame, keysName = 'Keys', valuesName = 'Values'): DG.DataFrame {
  let keys: string[] = [];
  const values: Float32Array = new Float32Array(df.columns.length*df.rowCount);
  let i = 0;

  for (const c of df.columns.toList()) {
    keys = keys.concat(Array<string>(c.length).fill(c.name));
    values.set(c.getRawData(), i);
    i += df.rowCount;
  }
  assert(keys.length == values.length);
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(keysName, keys),
    DG.Column.fromFloat32Array(valuesName, values),
  ]);
}

/**
 * Formats a matrix into <category1>-<category2>-<value> format.
 *
 * @param {DG.DataFrame} adjMatrix A data matrix to deal with.
 * @return {DG.DataFrame} The resulting data frame.
 */
function createNetwork(adjMatrix: DG.DataFrame): DG.DataFrame {
  const nCols = adjMatrix.columns.length;
  const nRows = adjMatrix.rowCount;

  assert(nCols == nRows);

  const pos1: Array<number> = [];
  const pos2: Array<number> = [];
  const weight: Array<number> = [];

  for (let i = 0; i < nCols; ++i) {
    const c = adjMatrix.columns.byIndex(i);

    for (let j = i+1; j < nRows; ++j) {
      const r = c.getRawData()[j];

      if (Math.abs(r) > 0) {
        pos1.push(i+1);
        pos2.push(j+1);
        weight.push(r);
      }
    }
  }

  const pos1Col = DG.Column.fromList('int', 'pos1', pos1);
  const pos2Col = DG.Column.fromList('int', 'pos2', pos2);
  const weightCol = DG.Column.fromList('double', 'weight', weight);

  return DG.DataFrame.fromColumns([pos1Col, pos2Col, weightCol]);
}

/**
 * Calculates Spearman's rho rank correlation coefficient.
 *
 * @param {DG.DataFrame} df A data frame to process.
 * @return {DG.DataFrame} The correlation matrix.
 */
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

type PairGuide = {[key: string]: number};

function makeSequenceTemplate(pos1: number, pos2: number, sequences: string[]): [string, number] {
  const guide: PairGuide = {};
  const nSeqs = sequences.length;

  for (let i = 0; i < nSeqs; ++i) {
    const s = sequences[i];
    const key: string = [s[pos1-1], s[pos2-1]].join('');
    guide[key] = (guide[key] || 0) + 1;
  }
  for (const key of Object.keys(guide)) {
    guide[key] /= nSeqs;
  }
  const sortableArray = Object.entries(guide);
  const sortedArray = sortableArray.sort(([, a], [, b]) => b - a); // Reverse order
  //const sortedObject = Object.fromEntries(sortedArray);
  return sortedArray[0]; // Pair with max frequency.
}

function analyseCorrelation(network: DG.DataFrame, sequences: string[]) {
  assert(network.columns.length == 3);

  const [pos1Col, pos2Col, weightCol] =
    new Array<any>(network.columns.length)
      .fill(0)
      .map((_, i) => network.columns.byIndex(i).getRawData());

  const pairs = [];

  for (let i = 0; i < network.rowCount; ++i) {
    const [pos1, pos2, weight] = [pos1Col[i], pos2Col[i], weightCol[i]];
    const [pair, freq] = makeSequenceTemplate(pos1, pos2, sequences);
    pairs.push([pos1, pos2, weight, pair, freq]);
  }
  return pairs;
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
      tau[i][j] = (res.prob < alpha) && (Math.abs(res.test) >= 0.5) ? res.test : 0;
      tau[j][i] = tau[i][j];
    }
  }
  return matrix2DataFrame(tau);
}

function describeNetwork(positions: DG.DataFrame, network: any[]): DG.DataFrame {
  const cols = Array.from(positions.columns).map(
    (_, i) => DG.Column.fromList('string', positions.columns.byIndex(i).name, []),
  );
  const df = DG.DataFrame.fromColumns(cols);

  for (const [pos1, pos2, weight, pair, freq] of network) {
    //const [pos1, pos2, weight] = Array.from(r.cells).map((c) => c.value);
    const values = new Array(df.columns.length).fill('');
    values[pos1-1] = pair[0];//(weight as number).toString();
    values[pos2-1] = pair[1];//(weight as number).toString();
    df.rows.addNew(values);
    console.log([pos1, pos2, weight, pair, freq]);
  }
  return df;
}

/**
 * Creates acorrelation plot and a box plot to perform correlation analysis.
 *
 * @export
 * @param {DG.Column} sequencesColumn A column containing amino acid sequences.
 * @return {[DG.Viewer, DG.Viewer]} These two plots.
 */
export function correlationAnalysisPlots(sequencesColumn: DG.Column): [DG.Viewer, DG.Viewer, DG.Viewer] {
  const [posDF, sequences] = calcPositions(sequencesColumn);
  /*const cpviewer = DG.Viewer.fromType(
    DG.VIEWER.CORR_PLOT,
    posDF,
    {
      'xColumnNames': posDF.columns.names(),
      'yColumnNames': posDF.columns.names(),
      'correlationType': 'Spearman',
    });*/

  const ccDF = calcKendallTauMatrix(posDF);

  const hmviewer = DG.Viewer.fromType(
    DG.VIEWER.HEAT_MAP,
    posDF,
  );

  /*const meltDF = melt(posDF, 'Position', 'Scale');

  const bpviewer = DG.Viewer.fromType(
    DG.VIEWER.BOX_PLOT,
    meltDF, {
      'categoryColumnName': 'Position',
      'valueColumnName': 'Scale',
      'statistics': ['min', 'max', 'avg', 'med'],
    });*/

  const nwDF = createNetwork(ccDF);

  const nwviewer = DG.Viewer.fromType(
    DG.VIEWER.NETWORK_DIAGRAM,
    nwDF, {
      'node1ColumnName': 'pos1',
      'node2ColumnName': 'pos2',
      'edgeColorColumnName': 'weight',
      'edgeWidthColumnName': 'weight',
      'edgeWidth': 4,
    });

  const pairs = analyseCorrelation(nwDF, sequences);
  const nwdDF = describeNetwork(posDF, pairs);

  const caviewer = DG.Viewer.fromType(
    DG.VIEWER.GRID,
    nwdDF);

  return [caviewer, hmviewer, nwviewer];
}
