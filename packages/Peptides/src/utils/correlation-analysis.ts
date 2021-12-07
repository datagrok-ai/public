/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as ui from 'datagrok-api/ui';
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
export function calcPositions(col: DG.Column): [DG.DataFrame, string[]] {
  const sequences = col.toList().map((v, _) => AlignedSequenceEncoder.clean(v));
  const enc = new AlignedSequenceEncoder('WimleyWhite');
  const encSeqs = sequences.map((v) => Vector.from(enc.encode(v)));
  console.log(encSeqs);
  const positions = transposeMatrix(encSeqs);
  return [matrix2DataFrame(positions), sequences];
}

function encodeSequences(df: DG.DataFrame, positionColumns?: string[]): DG.DataFrame {
  const [nCols, nRows] = [positionColumns ? positionColumns.length : df.columns.length, df.rowCount];
  //const sequences: string[] = new Array(nRows).fill('');
  const enc = new AlignedSequenceEncoder('WimleyWhite');
  const positions = new Array(nCols).fill(0).map((_) => new Float32Array(nRows));

  for (let i = 0; i < nCols; ++i) {
    const col: DG.Column = positionColumns ? df.getCol(positionColumns[i]) : df.columns.byIndex(i);
    //let s = '';

    for (let j = 0; j < nRows; ++j) {
      const letter = col.get(j);
      //s += col[j];
      positions[i][j] = enc.encodeLettter(letter);
    }
    //sequences[i] = s;
  }
  const posDF = DG.DataFrame.fromColumns(positions.map((v, i) => DG.Column.fromFloat32Array(df.columns.names()[i], v)));
  return posDF;//, sequences];
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
/*function melt(df: DG.DataFrame, keysName = 'Keys', valuesName = 'Values'): DG.DataFrame {
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
}*/

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
/*function calcSpearmanRhoMatrix(df: DG.DataFrame): DG.DataFrame {
  const nItems = df.columns.length;
  const rho = new Array(nItems).fill(0).map((_) => new Float32Array(nItems).fill(0));

  for (let i = 0; i < nItems; ++i) {
    for (let j = i+1; j < nItems; ++j) {
      rho[i][j] = df.columns.byIndex(i).stats.spearmanCorr(df.columns.byIndex(j));
      rho[j][i] = rho[i][j];
    }
  }
  return matrix2DataFrame(rho);
}*/

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
  console.log(guide);
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

function calcCorrelationMatrix(df: DG.DataFrame) {
  return calcKendallTauMatrix(df);
}

function describeNetwork(positions: DG.DataFrame, network: any[]): DG.DataFrame {
  const cols = Array.from(positions.columns).map(
    (_, i) => DG.Column.fromList('string', positions.columns.byIndex(i).name, []),
  );
  const df = DG.DataFrame.fromColumns(cols);

  for (const [pos1, pos2, weight, pair, freq] of network) {
    //const [pos1, pos2, weight] = Array.from(r.cells).map((c) => c.value);
    const values = new Array(df.columns.length).fill('');
    values[pos1-1] = (weight as number).toFixed(2);//pair[0];//(weight as number).toString();
    values[pos2-1] = (weight as number).toFixed(2);//pair[1];//(weight as number).toString();
    df.rows.addNew(values);
    console.log([pos1, pos2, weight, pair, freq]);
  }
  return df;
}

type Weights = {[pos: number]: number};
export type Guide = {[pos: number]: Weights};

const api = <any>window;

class CorrelationDescriptionViewer extends DG.Grid {
  protected network: DG.DataFrame;
  protected guide: Guide;

  constructor(d: any, network: DG.DataFrame) {
    super(d);
    this.network = network;
    this.guide = this.calcGuide();
    //this.onCellTooltip(this.customCellTooltip);
  }

  static fromNetwork(table: { d: any; }, network: DG.DataFrame): CorrelationDescriptionViewer {
    return new CorrelationDescriptionViewer(api.grok_Grid_Create(table.d), network);
  }

  protected calcGuide(): Guide {
    assert(this.network.columns.length == 3);

    const guide: Guide = {};
    let [pos1Col, pos2Col, weightCol] = Array.from(this.network.columns);

    pos1Col = pos1Col.getRawData();
    pos2Col = pos2Col.getRawData();
    weightCol = weightCol.getRawData();

    function _addWeight(pos1: number, pos2: number, weight: number) {
      if (guide[pos1] == undefined) {
        guide[pos1] = {};
      }
      guide[pos1][pos2] = weight;
    }

    for (let i = 0; i < this.network.rowCount; ++i) {
      const [pos1, pos2, weight] = [pos1Col[i], pos2Col[i], weightCol[i]];
      _addWeight(pos1, pos2, weight);
      _addWeight(pos2, pos1, weight);
    }
    console.log(guide);
    return guide;
  }

  public customCellTooltip(cell: DG.GridCell, x: number, y: number) {
    if (cell.isColHeader && cell.tableColumn != null) {
      const pos1 = parseInt(cell.tableColumn.name);

      const elements: HTMLElement[] = [];
      elements.push(ui.divText(cell.tableColumn.name, {style: {fontWeight: 'bold', fontSize: 10}}));
      elements.push(ui.divText('Found correlations with:\n'));

      for (const [pos2, weight] of Object.entries(this.guide[pos1])) {
        elements.push(ui.divText(`${pos2}: R = ${weight.toFixed(2)}\n`, {style: {color: weight > 0 ? 'red' : 'blue'}}));
      }

      ui.tooltip.show(ui.divV(elements), x, y);
      return true;
    }
  }
}

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
  console.log(guide);
  return guide;
}

export function correlationAnalysis(df: DG.DataFrame, positionColumns?: string[]): Guide {
  const posDF = encodeSequences(df, positionColumns);
  const ccDF = calcCorrelationMatrix(posDF);
  const nwDF = createNetwork(ccDF);
  //const pairs = analyseCorrelation(nwDF, sequences);
  //const nwdDF = describeNetwork(posDF, pairs);
  const guide = calcGuide(nwDF);
  return guide;
}

/**
 * Creates a correlation plot and a box plot to perform correlation analysis.
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

  /*const caviewer = DG.Viewer.fromType(
    DG.VIEWER.GRID,
    nwdDF,
  );*/
  const caviewer = CorrelationDescriptionViewer.fromNetwork(nwdDF, nwDF);
  caviewer.onCellTooltip(caviewer.customCellTooltip.bind(caviewer));

  return [caviewer, hmviewer, nwviewer];
}

export class CorrelationAnalysisVisualizer {
  protected guide: Guide;
  protected highlightedColumns: number[];
  protected grid: DG.Grid;
  protected positions: string[];

  constructor(grid: DG.Grid, positionColumns: string[]) {
    this.grid = grid;
    this.positions = positionColumns;

    const df = this.grid.dataFrame;

    if (df) {
      this.guide = correlationAnalysis(df, this.positions);
      this.highlightedColumns = Object.keys(this.guide).map((v) => parseInt(v));
    } else {
      throw new Error('Dataframe was not found in the grid.');
    }
  }

  public getTooltipElements(column: DG.Column, positions: string[]): HTMLElement[] {
    const name = column.name;
    const pos1 = parseInt(name);

    if (this.highlightedColumns.includes(pos1)) {
      const padLen = Math.round(Math.log10(positions.length))+1;

      const elements: HTMLElement[] = [];
      elements.push(ui.divText(name, {style: {fontWeight: 'bold', fontSize: 10}}));
      elements.push(ui.divText('Found correlations with:\n'));

      for (const [pos2, weight] of Object.entries(this.guide[pos1])) {
        const w = (weight as number);
        const style = {style: {color: w > 0 ? 'red' : 'blue'}};
        elements.push(ui.divText(`${pos2.padStart(padLen, '0')}: R = ${w.toFixed(2)}\n`, style));
      }

      return elements;
    }
    return [];
  }

  public customGridColumnHeader(args: DG.GridCellRenderArgs) {
    const cell = args.cell;

    assert(cell.isColHeader);
    //console.log('customGridColumnHeader');

    if (cell.tableColumn && this.highlightedColumns.includes(parseInt(cell.tableColumn.name))) {
      //console.log(`onCellPrepare: set background for ${cell.tableColumn.name}`);
      cell.style.backColor = DG.Color.lightBlue;
    }
  }
}
