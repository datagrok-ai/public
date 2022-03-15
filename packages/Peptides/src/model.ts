import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject, Observable} from 'rxjs';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {addViewerToHeader, StackedBarChart} from './viewers/stacked-barchart-viewer';
import {PeptidesController} from './peptides';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {ChemPalette} from './utils/chem-palette';
import {MonomerLibrary} from './monomer-library';


export class PeptidesModel {
  private _dataFrame: DG.DataFrame;
  private _activityColumn: string | null;
  private _activityScaling: string | null;
  private _sourceGrid: DG.Grid | null;
  private _twoColorMode: boolean | null;
  private _initialBitset: DG.BitSet | null;
  private _grouping: boolean = false;
  private _isUpdating: boolean = false;
  private _substFlag = false;
  private _statsDataFrameSubject = new Subject<DG.DataFrame>();
  private _sarGridSubject = new Subject<DG.Grid>();
  private _sarVGridSubject = new Subject<DG.Grid>();
  private _groupMappingSubject = new Subject<StringDictionary>();
  private _substFlagSubject = new Subject<boolean>();
  private static _modelName = 'peptidesModel';

  private constructor(dataFrame: DG.DataFrame) {
    this._dataFrame = dataFrame;
    this._activityColumn = null;
    this._activityScaling = null;
    this._sourceGrid = null;
    this._twoColorMode = null;
    this._initialBitset = null;
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName];
  }

  get dataFrame(): DG.DataFrame {
    return this._dataFrame;
  }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {
    return this._statsDataFrameSubject.asObservable();
  }

  get onSARGridChanged(): Observable<DG.Grid> {
    return this._sarGridSubject.asObservable();
  }

  get onSARVGridChanged(): Observable<DG.Grid> {
    return this._sarVGridSubject.asObservable();
  }

  get onGroupMappingChanged(): Observable<StringDictionary> {
    return this._groupMappingSubject.asObservable();
  }

  get onSubstFlagChanged(): Observable<boolean> {
    return this._substFlagSubject.asObservable();
  }

  async updateData(activityCol: string | null, activityScaling: string | null, sourceGrid: DG.Grid | null,
    twoColorMode: boolean | null, initialBitset: DG.BitSet | null, grouping: boolean | null) {
    this._activityColumn = activityCol ?? this._activityColumn;
    this._activityScaling = activityScaling ?? this._activityScaling;
    this._sourceGrid = sourceGrid ?? this._sourceGrid;
    this._twoColorMode = twoColorMode ?? this._twoColorMode;
    this._initialBitset = initialBitset ?? this._initialBitset;
    this._grouping = grouping ?? this._grouping;
    await this.updateDefault();
  }

  async updateDefault() {
    if (
      this._activityColumn && this._activityScaling && this._sourceGrid && this._twoColorMode !== null &&
      !this._isUpdating
    ) {
      this._isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, groupMapping] = await this.describe(
        this._activityColumn, this._activityScaling, this._twoColorMode, this._grouping);
      this._statsDataFrameSubject.next(statsDf);
      this._groupMappingSubject.next(groupMapping);
      this._sarGridSubject.next(viewerGrid);
      this._sarVGridSubject.next(viewerVGrid);
      this._substFlag = !this._substFlag;
      this._substFlagSubject.next(this._substFlag);

      this._sourceGrid.invalidate();

      this._isUpdating = false;
    }

    await this.updateBarchart();
  }

  async updateBarchart() {
    const stackedBarchart = await this._dataFrame?.plot.fromType('StackedBarChartAA') as StackedBarChart;
    if (stackedBarchart && this._sourceGrid)
      addViewerToHeader(this._sourceGrid, stackedBarchart);
  }

  static get modelName() {
    return PeptidesModel._modelName;
  }

  async describe(
    activityColumn: string, activityScaling: string, twoColorMode: boolean, grouping: boolean,
  ): Promise<[DG.Grid, DG.Grid, DG.DataFrame, StringDictionary]> {
    if (this._sourceGrid === null)
      throw new Error(`Source grid is not initialized`);

    //Split the aligned sequence into separate AARs
    let splitSeqDf: DG.DataFrame | undefined;
    let invalidIndexes: number[];
    const col: DG.Column = (this._dataFrame.columns as DG.ColumnList).bySemType('alignedSequence')!;
    [splitSeqDf, invalidIndexes] = PeptidesController.splitAlignedPeptides(col);
    splitSeqDf.name = 'Split sequence';

    const positionColumns = (splitSeqDf.columns as DG.ColumnList).names();
    const activityColumnScaled = `${activityColumn}Scaled`;
    const renderColNames: string[] = (splitSeqDf.columns as DG.ColumnList).names();
    const positionColName = 'Pos';
    const aminoAcidResidue = 'AAR';

    (splitSeqDf.columns as DG.ColumnList).add(this._dataFrame.getCol(activityColumn));

    joinDataFrames(this._dataFrame, positionColumns, splitSeqDf, activityColumn);

    for (const dfCol of (this._dataFrame.columns as DG.ColumnList)) {
      if (splitSeqDf.col(dfCol.name) && dfCol.name != activityColumn)
        PeptidesController.setAARRenderer(dfCol, this._sourceGrid);
    }

    sortSourceGrid(this._sourceGrid);

    const [scaledDf, newColName] = await PeptidesController.scaleActivity(
      activityScaling, activityColumn, activityColumnScaled, this._dataFrame);
    //TODO: make another func
    const scaledCol = scaledDf.getCol(activityColumnScaled);
    const oldScaledCol = this._dataFrame.getCol(activityColumnScaled);
    const oldScaledColGridName = oldScaledCol.temp['gridName'];
    const oldScaledGridCol = this._sourceGrid.col(oldScaledColGridName);

    (splitSeqDf.columns as DG.ColumnList).add(scaledCol);
    (this._dataFrame.columns as DG.ColumnList).replace(oldScaledCol, scaledCol);
    if (newColName === activityColumn)
      this._sourceGrid.col(activityColumn)!.name = `~${activityColumn}`;
    if (oldScaledGridCol !== null) {
      oldScaledGridCol.name = newColName;
      oldScaledGridCol.visible = true;
    }
    this._sourceGrid.columns.setOrder([newColName]);

    splitSeqDf = splitSeqDf.clone(this._initialBitset);

    //unpivot a table and handle duplicates
    splitSeqDf = splitSeqDf.groupBy(positionColumns)
      .add('med', activityColumnScaled, activityColumnScaled)
      .aggregate();

    const peptidesCount = splitSeqDf.getCol(activityColumnScaled).length;

    let matrixDf = splitSeqDf.unpivot([activityColumnScaled], positionColumns, positionColName, aminoAcidResidue);

    //TODO: move to chem palette
    let groupMapping: StringDictionary = {};
    if (grouping) {
      groupMapping = aarGroups;
      const aarCol = matrixDf.getCol(aminoAcidResidue);
      aarCol.init((index) => groupMapping[aarCol.get(index)[0]] ?? '-');
      aarCol.compact();
    } else
      Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);


    //statistics for specific AAR at a specific position
    const statsDf = await calculateStatistics(
      matrixDf, positionColName, aminoAcidResidue, activityColumnScaled, peptidesCount, splitSeqDf, groupMapping,
    );

    // SAR matrix table
    //pivot a table to make it matrix-like
    matrixDf = statsDf.groupBy([aminoAcidResidue])
      .pivot(positionColName)
      .add('first', 'Mean difference', '')
      .aggregate();
    matrixDf.name = 'SAR';

    // Setting category order
    await setCategoryOrder(twoColorMode, statsDf, aminoAcidResidue, matrixDf);

    // SAR vertical table (naive, choose best Mean difference from pVals <= 0.01)
    const sequenceDf = createVerticalTable(statsDf, aminoAcidResidue, positionColName, twoColorMode);
    renderColNames.push('Mean difference');

    const [sarGrid, sarVGrid] = createGrids(
      matrixDf, aminoAcidResidue, positionColumns, sequenceDf, positionColName, grouping,
    );

    setCellRendererFunc(
      renderColNames, positionColName, aminoAcidResidue, statsDf, twoColorMode, sarGrid, sarVGrid,
    );

    // show all the statistics in a tooltip over cell
    setTooltipFunc(
      renderColNames, statsDf, aminoAcidResidue, positionColName, peptidesCount, grouping, sarGrid, sarVGrid,
      this._dataFrame,
    );

    postProcessGrids(this._sourceGrid, invalidIndexes, grouping, aminoAcidResidue, sarGrid, sarVGrid);

    //TODO: return class instead
    return [sarGrid, sarVGrid, statsDf, groupMapping];
  }
}

export const aarGroups = {
  'R': 'PC',
  'H': 'PC',
  'K': 'PC',
  'D': 'NC',
  'E': 'NC',
  'S': 'U',
  'T': 'U',
  'N': 'U',
  'Q': 'U',
  'C': 'SC',
  'U': 'SC',
  'G': 'SC',
  'P': 'SC',
  'A': 'H',
  'V': 'H',
  'I': 'H',
  'L': 'H',
  'M': 'H',
  'F': 'H',
  'Y': 'H',
  'W': 'H',
  '-': '-',
};

const groupDescription: {[key: string]: {'description': string, 'aminoAcids': string[]}} = {
  'PC': {'description': 'Positive Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['R', 'H', 'K']},
  'NC': {'description': 'Negative Amino Acids, with Electrically Charged Side Chains', 'aminoAcids': ['D', 'E']},
  'U': {'description': 'Amino Acids with Polar Uncharged Side Chains', 'aminoAcids': ['S', 'T', 'N', 'Q']},
  'SC': {'description': 'Special Cases', 'aminoAcids': ['C', 'U', 'G', 'P']},
  'H': {
    'description': 'Amino Acids with Hydrophobic Side Chain',
    'aminoAcids': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
  },
  '-': {'description': 'Unknown Amino Acid', 'aminoAcids': ['-']},
};

function joinDataFrames(df: DG.DataFrame, positionColumns: string[], splitSeqDf: DG.DataFrame, activityColumn: string) {
  // if (df.col(activityColumnScaled))
  //   (df.columns as DG.ColumnList).remove(activityColumnScaled);


  //FIXME: this column usually duplicates, so remove it then
  // if (df.col(`${activityColumnScaled} (2)`))
  //   (df.columns as DG.ColumnList).remove(`${activityColumnScaled} (2)`);


  // append splitSeqDf columns to source table and make sure columns are not added more than once
  const dfColsSet = new Set((df.columns as DG.ColumnList).names());
  if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
    df.join(
      splitSeqDf, [activityColumn], [activityColumn], (df.columns as DG.ColumnList).names(), positionColumns, 'inner',
      true);
  }
}

function sortSourceGrid(sourceGrid: DG.Grid) {
  if (sourceGrid) {
    const colNames: DG.GridColumn[] = [];
    for (let i = 1; i < sourceGrid.columns.length; i++)
      colNames.push(sourceGrid.columns.byIndex(i)!);

    colNames.sort((a, b)=>{
      if (a.column!.semType == 'aminoAcids') {
        if (b.column!.semType == 'aminoAcids')
          return 0;
        return -1;
      }
      if (b.column!.semType == 'aminoAcids')
        return 1;
      return 0;
    });
    sourceGrid.columns.setOrder(colNames.map((v) => v.name));
  }
}

async function calculateStatistics(
  matrixDf: DG.DataFrame, positionColName: string, aminoAcidResidue: string, activityColumnScaled: string,
  peptidesCount: number, splitSeqDf: DG.DataFrame, groupMapping: StringDictionary,
) {
  matrixDf = matrixDf.groupBy([positionColName, aminoAcidResidue])
    .add('count', activityColumnScaled, 'Count')
    .aggregate();

  const countThreshold = 4;
  //@ts-ignore: never gets old
  matrixDf.rows.filter((row) => row.Count >= countThreshold && row.Count <= peptidesCount - countThreshold);
  matrixDf = matrixDf.clone(matrixDf.filter);

  // calculate additional stats
  await (matrixDf.columns as DG.ColumnList).addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

  //calculate p-values based on t-test
  let pvalues: Float32Array = new Float32Array(matrixDf.rowCount).fill(1);
  const mdCol: DG.Column = (matrixDf.columns as DG.ColumnList).addNewFloat('Mean difference');
  const pValCol: DG.Column = (matrixDf.columns as DG.ColumnList).addNewFloat('pValue');
  for (let i = 0; i < matrixDf.rowCount; i++) {
    const position = matrixDf.get(positionColName, i);
    const aar = matrixDf.get(aminoAcidResidue, i);

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] === aar);
    const currentActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] !== aar);
    const otherActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [activityColumnScaled])
      .getCol(activityColumnScaled)
      .toList();

    const testResult = tTest(currentActivity, otherActivity);
    // testResult = uTest(currentActivity, otherActivity);
    const currentMeanDiff = testResult['Mean difference']!;
    const pvalue = testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'];

    mdCol.set(i, currentMeanDiff);
    pvalues[i] = pvalue;
  }

  pvalues = fdrcorrection(pvalues)[1];

  for (let i = 0; i < pvalues.length; ++i)
    pValCol.set(i, pvalues[i]);

  return matrixDf.clone();
}

async function setCategoryOrder(
  twoColorMode: boolean, statsDf: DG.DataFrame, aminoAcidResidue: string, matrixDf: DG.DataFrame,
) {
  const sortArgument = twoColorMode ? 'Absolute Mean difference' : 'Mean difference';
  if (twoColorMode)
    await (statsDf.columns as DG.ColumnList).addNewCalculated('Absolute Mean difference', 'Abs(${Mean difference})');

  const aarWeightsDf = statsDf.groupBy([aminoAcidResidue]).sum(sortArgument, 'weight').aggregate();
  const aarList = aarWeightsDf.getCol(aminoAcidResidue).toList();
  const getWeight = (aar: string) => aarWeightsDf
    .groupBy(['weight'])
    .where(`${aminoAcidResidue} = ${aar}`)
    .aggregate()
    .get('weight', 0);
  aarList.sort((first, second) => getWeight(second) - getWeight(first));

  matrixDf.getCol(aminoAcidResidue).setCategoryOrder(aarList);
}

function createVerticalTable(
  statsDf: DG.DataFrame, aminoAcidResidue: string, positionColName: string, twoColorMode: boolean,
) {
  // TODO: aquire ALL of the positions
  let sequenceDf = statsDf.groupBy(['Mean difference', aminoAcidResidue, positionColName, 'Count', 'Ratio', 'pValue'])
    .where('pValue <= 0.1')
    .aggregate();

  let tempStats: DG.Stats;
  const maxAtPos: {[index: string]: number} = {};
  for (const pos of sequenceDf.getCol(positionColName).categories) {
    tempStats = DG.Stats.fromColumn(
      sequenceDf.getCol('Mean difference'),
      DG.BitSet.create(sequenceDf.rowCount, (i) => sequenceDf.get(positionColName, i) === pos),
    );
    maxAtPos[pos] = twoColorMode ?
      (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) : tempStats.max;
  }
  sequenceDf = sequenceDf.clone(DG.BitSet.create(sequenceDf.rowCount, (i) => {
    return sequenceDf.get('Mean difference', i) === maxAtPos[sequenceDf.get(positionColName, i)];
  }));

  return sequenceDf;
}

function createGrids(
  matrixDf: DG.DataFrame, aminoAcidResidue: string, positionColumns: string[], sequenceDf: DG.DataFrame,
  positionColName: string, grouping: boolean,
) {
  const sarGrid = matrixDf.plot.grid();
  sarGrid.sort([aminoAcidResidue]);
  sarGrid.columns.setOrder([aminoAcidResidue].concat(positionColumns));

  const sarVGrid = sequenceDf.plot.grid();
  sarVGrid.sort([positionColName]);
  sarVGrid.col('pValue')!.format = 'four digits after comma';
  sarVGrid.col('pValue')!.name = 'P-Value';

  if (!grouping) {
    let tempCol = (matrixDf.columns as DG.ColumnList).byName(aminoAcidResidue);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);

    tempCol = (sequenceDf.columns as DG.ColumnList).byName(aminoAcidResidue);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);
  }

  return [sarGrid, sarVGrid];
}

function setCellRendererFunc(
  renderColNames: string[], positionColName: string, aminoAcidResidue: string, statsDf: DG.DataFrame,
  twoColorMode: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid,
) {
  const mdCol = statsDf.getCol('Mean difference');
  const cellRendererFunc = function(args: DG.GridCellRenderArgs) {
    args.g.save();
    args.g.beginPath();
    args.g.rect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
    args.g.clip();

    if (args.cell.isRowHeader && args.cell.gridColumn.visible) {
      args.cell.gridColumn.visible = false;
      args.preventDefault();
      return;
    }

    if (
      args.cell.isTableCell &&
      args.cell.tableRowIndex !== null &&
      args.cell.tableColumn !== null &&
      args.cell.cell.value !== null
    ) {
      if (renderColNames.indexOf(args.cell.tableColumn.name) !== -1) {
        const currentPosition = args.cell.tableColumn.name !== 'Mean difference' ?
          args.cell.tableColumn.name : args.cell.grid.table.get(positionColName, args.cell.tableRowIndex);
        const query =
          `${aminoAcidResidue} = ${args.cell.grid.table.get(aminoAcidResidue, args.cell.tableRowIndex)} ` +
          `and ${positionColName} = ${currentPosition}`;

        const pVal: number = statsDf.groupBy(['pValue']).where(query).aggregate().get('pValue', 0);

        let coef;
        const variant = args.cell.cell.value < 0;
        if (pVal < 0.01)
          coef = variant && twoColorMode ? '#FF7900' : '#299617';
        else if (pVal < 0.05)
          coef = variant && twoColorMode ? '#FFA500' : '#32CD32';
        else if (pVal < 0.1)
          coef = variant && twoColorMode ? '#FBCEB1' : '#98FF98';
        else
          coef = DG.Color.toHtml(DG.Color.lightLightGray);


        const chooseMin = () => twoColorMode ? 0 : mdCol.min;
        const chooseMax = () => twoColorMode ? Math.max(Math.abs(mdCol.min), mdCol.max) : mdCol.max;
        const chooseCurrent = () => twoColorMode ? Math.abs(args.cell.cell.value) : args.cell.cell.value;

        const rCoef = (chooseCurrent() - chooseMin()) / (chooseMax() - chooseMin());

        const maxRadius = 0.9 * (args.bounds.width > args.bounds.height ? args.bounds.height : args.bounds.width) / 2;
        const radius = Math.floor(maxRadius * rCoef);

        args.g.beginPath();
        args.g.fillStyle = coef;
        args.g.arc(
          args.bounds.x + args.bounds.width / 2, args.bounds.y + args.bounds.height / 2, radius < 3 ? 3 : radius, 0,
          Math.PI * 2, true,
        );
        args.g.closePath();

        args.g.fill();
        args.preventDefault();
      }
    }
    args.g.restore();
  };
  sarGrid.onCellRender.subscribe(cellRendererFunc);
  sarVGrid.onCellRender.subscribe(cellRendererFunc);
}

function setTooltipFunc(
  renderColNames: string[], statsDf: DG.DataFrame, aminoAcidResidue: string, positionColName: string,
  peptidesCount: number, grouping: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid, sourceDf: DG.DataFrame,
) {
  const onCellTooltipFunc = async function(cell: DG.GridCell, x: number, y: number) {
    if (
      !cell.isRowHeader && !cell.isColHeader && cell.tableColumn !== null && cell.cell.value !== null &&
        cell.tableRowIndex !== null && renderColNames.indexOf(cell.tableColumn.name) !== -1) {
      const tooltipMap: { [index: string]: string } = {};

      for (const col of (statsDf.columns as DG.ColumnList).names()) {
        if (col !== aminoAcidResidue && col !== positionColName) {
          const currentPosition = cell.tableColumn.name !== 'Mean difference' ?
            cell.tableColumn.name : cell.grid.table.get(positionColName, cell.tableRowIndex);
          const query =
            `${aminoAcidResidue} = ${cell.grid.table.get(aminoAcidResidue, cell.tableRowIndex)} ` +
            `and ${positionColName} = ${currentPosition}`;
          const textNum = statsDf.groupBy([col]).where(query).aggregate().get(col, 0);
          let text = `${col === 'Count' ? textNum : textNum.toFixed(5)}`;

          if (col === 'Count')
            text += ` / ${peptidesCount}`;
          else if (col === 'pValue')
            text = parseFloat(text) !== 0 ? text : '<0.01';


          tooltipMap[col === 'pValue' ? 'p-value' : col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
    }
    if (
      !cell.isColHeader &&
      cell.tableColumn !== null &&
      cell.tableColumn.name == aminoAcidResidue &&
      cell.cell.value !== null &&
      cell.tableRowIndex !== null
    ) {
      if (grouping) {
        const currentGroup = groupDescription[cell.cell.value];
        const divText = ui.divText('Amino Acids in this group: ' + currentGroup['aminoAcids'].join(', '));
        ui.tooltip.show(ui.divV([ui.h3(currentGroup['description']), divText]), x, y);
      } else {
        const monomerLib = sourceDf.temp[MonomerLibrary.id];
        ChemPalette.showTooltip(cell, x, y, monomerLib);
      }
    }
    return true;
  };
  sarGrid.onCellTooltip(onCellTooltipFunc);
  sarVGrid.onCellTooltip(onCellTooltipFunc);
}

function postProcessGrids(
  sourceGrid: DG.Grid, invalidIndexes: number[], grouping: boolean, aminoAcidResidue: string, sarGrid: DG.Grid,
  sarVGrid: DG.Grid,
) {
  sourceGrid.onCellPrepare((cell: DG.GridCell) => {
    const currentRowIndex = cell.tableRowIndex;
    if (currentRowIndex && invalidIndexes.includes(currentRowIndex) && !cell.isRowHeader)
      cell.style.backColor = DG.Color.lightLightGray;
  });

  const mdCol: DG.GridColumn = sarVGrid.col('Mean difference')!;
  mdCol.name = 'Diff';

  for (const grid of [sarGrid, sarVGrid]) {
    grid.props.rowHeight = 20;
    grid.columns.rowHeader!.width = 20;
    for (let i = 0; i < grid.columns.length; ++i) {
      const col = grid.columns.byIndex(i)!;
      if (grid == sarVGrid && col.name !== 'Diff' && col.name !== 'AAR')
        col.width = 45;
      else
        col.width = grid.props.rowHeight;
    }
  }

  if (grouping) {
    sarGrid.col(aminoAcidResidue)!.name = 'Groups';
    sarVGrid.col(aminoAcidResidue)!.name = 'Groups';
  }

  sarGrid.props.allowEdit = false;
  sarVGrid.props.allowEdit = false;
}
