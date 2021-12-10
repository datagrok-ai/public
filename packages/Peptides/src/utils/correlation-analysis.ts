/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as DG from 'datagrok-api/dg';

import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';
import {assert} from '@datagrok-libraries/utils/src/operations';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {kendallsTau} from '@datagrok-libraries/statistics/src/correlation-coefficient';

/**
 * Converts a Matrix into a DataFrame.
 *
 * @export
 * @param {Matrix} matrix A matrix.
 * @return {DG.DataFrame} The data frame.
 */
function matrix2DataFrame(matrix: Matrix): DG.DataFrame {
  return DG.DataFrame.fromColumns(matrix.map((v, i) => DG.Column.fromFloat32Array(`${i+1}`, v)));
}

/**
 * Encodes sequence into a certain scale.
 *
 * @param {DG.DataFrame} df A data frame containing the sequences.
 * @param {string[]} [positionColumns] If given instructs which columns to consider as sequences containing.
 * @return {DG.DataFrame} The data frame with seqences encoded.
 */
function encodeSequences(df: DG.DataFrame, positionColumns?: string[]): DG.DataFrame {
  const [nCols, nRows] = [positionColumns ? positionColumns.length : df.columns.length, df.rowCount];
  const enc = new AlignedSequenceEncoder('WimleyWhite');
  const positions = new Array(nCols).fill(0).map((_) => new Float32Array(nRows));

  for (let i = 0; i < nCols; ++i) {
    const col: DG.Column = positionColumns ? df.getCol(positionColumns[i]) : df.columns.byIndex(i);

    for (let j = 0; j < nRows; ++j) {
      const letter = col.get(j);
      positions[i][j] = enc.encodeLettter(letter);
    }
  }
  const posDF = DG.DataFrame.fromColumns(positions.map((v, i) => DG.Column.fromFloat32Array(df.columns.names()[i], v)));
  return posDF;
}

/**
 * Formats an adjacency matrix into <category1>-<category2>-<value> format.
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
 * Calculates Kendall's tau rank correlation matrix.
 *
 * @param {DG.DataFrame} df A data frame to process.
 * @param {number} [alpha=0.05] The significance threshold.
 * @param {number} [rAbsCutoff=0.5] The absolute R cutoff.
 * @return {DG.DataFrame} The correlation matrix.
 */
function calcKendallTauMatrix(df: DG.DataFrame, alpha: number = 0.05, rAbsCutoff = 0.5): DG.DataFrame {
  const nItems = df.columns.length;
  const tau = new Array(nItems).fill(0).map((_) => new Float32Array(nItems).fill(0));

  for (let i = 0; i < nItems; ++i) {
    for (let j = i+1; j < nItems; ++j) {
      const res = kendallsTau(df.columns.byIndex(i).getRawData(), df.columns.byIndex(j).getRawData());
      tau[i][j] = (res.prob < alpha) && (Math.abs(res.test) >= rAbsCutoff) ? res.test : 0;
      tau[j][i] = tau[i][j];
    }
  }
  return matrix2DataFrame(tau);
}

/**
 * Calculates a correlation matrix via method chosen.
 *
 * @param {DG.DataFrame} df A data frame.
 * @return {DG.DataFrame} The correlation matrix.
 */
function calcCorrelationMatrix(df: DG.DataFrame): DG.DataFrame {
  return calcKendallTauMatrix(df);
}

type Weights = {[pos: number]: number};
type Guide = {[pos: number]: Weights};

/**
 * Calculates a dictionary with the keys containing the first correlating positions.
 * Values correspond to a dictionary containing the positions and corresponding R-value
 * which the given position correlating with.
 *
 * @param {DG.DataFrame} network A network to process.
 * @return {Guide} The formatted dictionary.
 */
function calcGuide(network: DG.DataFrame): Guide {
  assert(network.columns.length == 3);

  const guide: Guide = {};
  let [pos1Col, pos2Col, weightCol] = Array.from(network.columns);

  pos1Col = pos1Col.getRawData();
  pos2Col = pos2Col.getRawData();
  weightCol = weightCol.getRawData();

  function _addWeight(pos1: number, pos2: number, weight: number) {
    if (guide[pos1] == undefined) {
      guide[pos1] = {};
    }
    guide[pos1][pos2] = weight;
  }

  for (let i = 0; i < network.rowCount; ++i) {
    const [pos1, pos2, weight] = [pos1Col[i], pos2Col[i], weightCol[i]];
    _addWeight(pos1, pos2, weight);
    _addWeight(pos2, pos1, weight);
  }
  return guide;
}

function calcCorrelations(df: DG.DataFrame, positionColumns?: string[]): Guide {
  const posDF = encodeSequences(df, positionColumns);
  const ccDF = calcCorrelationMatrix(posDF);
  const nwDF = createNetwork(ccDF);
  const guide = calcGuide(nwDF);
  return guide;
}

/**
 * Formats correlating positions to place in the corresponding tooltips.
 * Higlights correlating positions' headers.
 *
 * @export
 * @class CorrelationAnalysisVisualizer
 */
export class CorrelationAnalysisVisualizer {
  protected guide: Guide;
  protected highlightedColumns: number[];

  /**
   * Creates an instance of CorrelationAnalysisVisualizer.
   * @param {DG.DataFrame} df A data frame to take sequences from.
   * @param {string[]} positionColumns Optional columns list to take the sequences from.
   * @memberof CorrelationAnalysisVisualizer
   */
  constructor(df: DG.DataFrame, positionColumns: string[]) {
    if (df) {
      this.guide = calcCorrelations(df, positionColumns);
      this.highlightedColumns = Object.keys(this.guide).map((v) => parseInt(v));
    } else {
      throw new Error('Dataframe was not found in the grid.');
    }
  }

  /**
   * Returns a dictionary with the correlating positions and their R-value.
   *
   * @readonly
   * @type {Guide} The dictionary.
   * @memberof CorrelationAnalysisVisualizer
   */
  get path(): Guide {
    return this.guide;
  }

  /**
   * Checks if the position column name is found among correlelating ones.
   *
   * @param {string} name The name of the column.
   * @return {boolean} True if the position is correlating with any oter.
   * @memberof CorrelationAnalysisVisualizer
   */
  public isPositionCorrelating(name: string): boolean {
    return this.highlightedColumns.includes(parseInt(name));
  }
}
