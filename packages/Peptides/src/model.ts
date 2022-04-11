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
import * as C from './utils/constants';


export class PeptidesModel {
  _dataFrame: DG.DataFrame;
  _activityScaling: string | null = null;
  _sourceGrid: DG.Grid | null = null;
  _twoColorMode: boolean | null = null;
  _initialBitset: DG.BitSet | null = null;
  _grouping: boolean = false;
  _isUpdating: boolean = false;
  _substitutionTable: DG.DataFrame | null = null;
  _statsDataFrameSubject = new Subject<DG.DataFrame>();
  _sarGridSubject = new Subject<DG.Grid>();
  _sarVGridSubject = new Subject<DG.Grid>();
  _groupMappingSubject = new Subject<StringDictionary>();
  _substitutionTableSubject = new Subject<DG.DataFrame>();
  static _modelName = 'peptidesModel';
  _activityLimit: number | null = null;
  _maxSubstitutions: number | null = null;
  _isSubstitutionOn: boolean | null = null;
  _isSubstInitialized = false;
  _substTableTooltipData: { [aar: string]: number[][][]; } | null = null;

  private constructor(dataFrame: DG.DataFrame) {
    this._dataFrame = dataFrame;
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName];
  }

  get dataFrame(): DG.DataFrame { return this._dataFrame; }

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> { return this._statsDataFrameSubject.asObservable(); }

  get onSARGridChanged(): Observable<DG.Grid> { return this._sarGridSubject.asObservable(); }

  get onSARVGridChanged(): Observable<DG.Grid> { return this._sarVGridSubject.asObservable(); }

  get onGroupMappingChanged(): Observable<StringDictionary> { return this._groupMappingSubject.asObservable(); }

  get onSubstTableChanged(): Observable<DG.DataFrame> { return this._substitutionTableSubject.asObservable(); }

  get substTooltipData(): { [aar: string]: number[][][]; } { return this._substTableTooltipData!; }

  async updateData(
    activityScaling?: string, sourceGrid?: DG.Grid, twoColorMode?: boolean, initialBitset?: DG.BitSet,
    grouping?: boolean, activityLimit?: number, maxSubstitutions?: number, isSubstitutionOn?: boolean) {
    this._activityScaling = activityScaling ?? this._activityScaling;
    this._sourceGrid = sourceGrid ?? this._sourceGrid;
    this._twoColorMode = twoColorMode ?? this._twoColorMode;
    this._initialBitset = initialBitset ?? this._initialBitset;
    this._grouping = grouping ?? this._grouping;
    this._activityLimit = activityLimit ?? this._activityLimit;
    this._maxSubstitutions = maxSubstitutions ?? this._maxSubstitutions;
    this._isSubstitutionOn = isSubstitutionOn ?? this._isSubstitutionOn;

    await this.updateDefault();
  }

  async updateDefault() {
    if (this._activityScaling && this._sourceGrid && this._twoColorMode !== null && !this._isUpdating) {
      this._isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf, substTable, groupMapping] = await this.initializeViewersComponents();
      this._statsDataFrameSubject.next(statsDf);
      this._groupMappingSubject.next(groupMapping);
      this._sarGridSubject.next(viewerGrid);
      this._sarVGridSubject.next(viewerVGrid);
      if (this._isSubstitutionOn) {
        this._substitutionTableSubject.next(substTable);
        this._isSubstInitialized = true;
      }

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

  async initializeViewersComponents(): Promise<[DG.Grid, DG.Grid, DG.DataFrame, DG.DataFrame, StringDictionary]> {
    if (this._sourceGrid === null)
      throw new Error(`Source grid is not initialized`);

    //Split the aligned sequence into separate AARs
    let splitSeqDf: DG.DataFrame | undefined;
    let invalidIndexes: number[];
    const col: DG.Column = (this._dataFrame.columns as DG.ColumnList).bySemType(C.SEM_TYPES.ALIGNED_SEQUENCE)!;
    [splitSeqDf, invalidIndexes] = PeptidesController.splitAlignedPeptides(col);
    splitSeqDf.name = 'Split sequence';

    const positionColumns = (splitSeqDf.columns as DG.ColumnList).names();
    const renderColNames: string[] = (splitSeqDf.columns as DG.ColumnList).names();

    (splitSeqDf.columns as DG.ColumnList).add(this._dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY));

    joinDataFrames(this._dataFrame, positionColumns, splitSeqDf);

    for (const dfCol of (this._dataFrame.columns as DG.ColumnList)) {
      if (splitSeqDf.col(dfCol.name) && dfCol.name != C.COLUMNS_NAMES.ACTIVITY)
        PeptidesController.setAARRenderer(dfCol, this._sourceGrid);
    }

    sortSourceGrid(this._sourceGrid);

    await createScaledCol(this._activityScaling!, this._dataFrame, this._sourceGrid, splitSeqDf);
    
    splitSeqDf = splitSeqDf.clone(this._initialBitset);

    //unpivot a table and handle duplicates
    splitSeqDf = splitSeqDf.groupBy(positionColumns)
      .add('med', C.COLUMNS_NAMES.ACTIVITY_SCALED, C.COLUMNS_NAMES.ACTIVITY_SCALED)
      .aggregate();

    const peptidesCount = splitSeqDf.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).length;

    let matrixDf = splitSeqDf.unpivot(
      [C.COLUMNS_NAMES.ACTIVITY_SCALED], positionColumns, C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);

    //TODO: move to chem palette
    let groupMapping: StringDictionary = {};
    if (this._grouping) {
      groupMapping = aarGroups;
      const aarCol = matrixDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
      aarCol.init((index) => groupMapping[aarCol.get(index)[0]] ?? '-');
      aarCol.compact();
    } else
      Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);


    //statistics for specific AAR at a specific position
    const statsDf = await calculateStatistics(matrixDf, peptidesCount, splitSeqDf, groupMapping);

    // SAR matrix table
    //pivot a table to make it matrix-like
    matrixDf = statsDf.groupBy([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.MEAN_DIFFERENCE, '')
      .aggregate();
    matrixDf.name = 'SAR';

    // Setting category order
    await setCategoryOrder(this._twoColorMode!, statsDf, matrixDf);

    // SAR vertical table (naive, choose best Mean difference from pVals <= 0.01)
    const sequenceDf = createVerticalTable(statsDf, this._twoColorMode!);
    renderColNames.push(C.COLUMNS_NAMES.MEAN_DIFFERENCE);

    let substTable: DG.DataFrame | null = null;
    if (this._isSubstitutionOn || !this._isSubstInitialized)
      substTable = this.calcSubstitutions();

    const [sarGrid, sarVGrid] = createGrids(matrixDf, positionColumns, sequenceDf, this._grouping);

    setCellRendererFunc(
      renderColNames, statsDf, this._twoColorMode!, sarGrid, sarVGrid, this._isSubstitutionOn!, substTable!);

    // show all the statistics in a tooltip over cell
    setTooltipFunc(renderColNames, statsDf, peptidesCount, this._grouping, sarGrid, sarVGrid, this._dataFrame);

    postProcessGrids(this._sourceGrid, invalidIndexes, this._grouping, sarGrid, sarVGrid);

    //TODO: return class instead
    return [sarGrid, sarVGrid, statsDf, substTable!, groupMapping];
  }

  calcSubstitutions() {
    const col: DG.Column = this.dataFrame.columns.bySemType(C.SEM_TYPES.ALIGNED_SEQUENCE);
    // let values: number[] = df.columns.byName('IC50').toList();
    const values: number[] = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).toList();
    // values = values;
    const splitedMatrix = split(col);

    const tableValues: { [aar: string]: number[] } = {};
    const tableTooltips: { [aar: string]: {}[][] } = {};
    const tableCases: { [aar: string]: number[][][] } = {};

    const nRows = splitedMatrix.length;
    const nCols = splitedMatrix[0].length;
    const nColsArray = Array(nCols);

    for (let i = 0; i < nRows - 1; i++) {
      for (let j = i + 1; j < nRows; j++) {
        let substCounter = 0;
        const subst1: { [pos: number]: [string, {}] } = {};
        const subst2: { [pos: number]: [string, {}] } = {};
        const delta = values[i] - values[j];
        
        if (Math.abs(delta) < this._activityLimit!)
          continue;
        
        for (let k = 0; k < nCols; k++) {
          const smik = splitedMatrix[i][k];
          const smjk = splitedMatrix[j][k];
          if (smik != smjk) {
            const vi = values[i].toFixed(2);
            const vj = values[j].toFixed(2);
            substCounter++;
            subst1[k] = [
              smik,
              {
                key: `${smik === '-' ? 'Empty' : smik} → ${smjk === '-' ? 'Empty' : smjk}`,
                value: `${vi} → ${vj}`,
                diff: values[j] - values[i],
              },
            ];
            subst2[k] = [
              smjk,
              {
                key: `${smjk === '-' ? 'Empty' : smjk} → ${smik === '-' ? 'Empty' : smik}`,
                value: `${vj} → ${vi}`,
                diff: values[i] - values[j],
              },
            ];
          }
        }

        if (substCounter > this._maxSubstitutions! || substCounter === 0)
          continue;

        for (const subst of [subst1, subst2]) {
          Object.keys(subst).forEach((pos) => {
            const posInt = parseInt(pos);
            const aar = subst[posInt][0];
            if (!Object.keys(tableValues).includes(aar)) {
              tableValues[aar] = Array(...nColsArray).map(() => DG.INT_NULL);
              tableTooltips[aar] = Array(...nColsArray).map(() => []);
              tableCases[aar] = Array(...nColsArray).map(() => []);
            }

            tableValues[aar][posInt] = tableValues[aar][posInt] === DG.INT_NULL ? 1 : tableValues[aar][posInt] + 1;
            tableTooltips[aar][posInt] = !tableTooltips[aar][posInt].length ?
              [{key: 'Substitution', value: 'Values'}] : tableTooltips[aar][posInt];
            tableTooltips[aar][posInt].push(subst[posInt][1]);
            tableCases[aar][posInt].push([i, j, subst == subst1 ? delta : -delta]);
          });
        }
      }
    }

    const tableValuesKeys = Object.keys(tableValues);
    const dfLength = tableValuesKeys.length;
    const cols = [...nColsArray.keys()].map((v) => DG.Column.int(`${v < 10 ? 0 : ''}${v}`, dfLength));
    cols.forEach((currentCol) => currentCol.semType = 'Substitution');
    const aarCol = DG.Column.string(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, dfLength);
    cols.splice(0, 1, aarCol);
    const table = DG.DataFrame.fromColumns(cols);

    for (let i = 0; i < dfLength; i++) {
      const aar = tableValuesKeys[i];
      tableValues[aar].splice(0, 1);
      table.rows.setValues(i, [aar, ...tableValues[aar]]);
    }

    // let groupMapping: { [key: string]: string } = {};

    //TODO: enable grouping
    // Object.keys(aarGroups).forEach((value) => groupMapping[value] = value);
    this._substTableTooltipData = tableCases;

    return table;

    // this.viewerGrid = table.plot.grid();

    // setAARRenderer(aarCol, this.viewerGrid);

    // this.viewerGrid.onCellTooltip(
    //   (gCell, x, y) => {
    //     if (gCell.cell.value !== DG.INT_NULL && gCell.tableColumn !== null && gCell.tableRowIndex !== null) {
    //       const colName = gCell.tableColumn.name;
    //       if (colName !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE) {
    //         const aar = this.viewerGrid!.table.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, gCell.tableRowIndex);
    //         const pos = parseInt(colName);
    //         const lengthTableTooltip = tableTooltips[aar][pos].length;
    //         const sortedTableTooltips = [];
    //         const resTooltip: {[index: string]: string}[] = [];
    //         let tooltipText: any = ui.divText('No substitutions');
    //         let haveEllipsis = false;

    //         if (lengthTableTooltip) {
    //           const mn = Math.min(5, lengthTableTooltip);
    //           for (let i = 0; i < lengthTableTooltip; ++i) {
    //             const val: {[key: string]: any} = tableTooltips[aar][pos][i];
    //             sortedTableTooltips.push([i, val['diff'], val]);
    //           }
    //           sortedTableTooltips.sort(function(a, b) {
    //             return b[1] - a[1];
    //           });
    //           for (let i = 0; i < mn; ++i) {
    //             const idx = sortedTableTooltips[i][0];
    //             resTooltip.push(tableTooltips[aar][pos][idx]);
    //           }
    //           if (lengthTableTooltip > mn) {
    //             for (let i = Math.max(lengthTableTooltip - mn, mn); i < lengthTableTooltip; ++i) {
    //               const idx = sortedTableTooltips[i][0];
    //               if (lengthTableTooltip > 2 * mn && !haveEllipsis) {
    //                 haveEllipsis = true;
    //                 resTooltip.push({key: '...', value: '...'});
    //               }
    //               resTooltip.push(tableTooltips[aar][pos][idx]);
    //             }
    //           }
    //           tooltipText = DG.HtmlTable.create(
    //             resTooltip, (item: {[index: string]: string}, idx: number) => [item.key, item.value],
    //           ).root;
    //         }
    //         ui.tooltip.show(tooltipText, x, y);
    //       }
    //     }
    //     return true;
    //   },
    // );

    // this.viewerGrid.columns.rowHeader!.width = 30;
    // this.viewerGrid.props.rowHeight = 20;
    // for (const col of table.columns.names()) {
    //   this.viewerGrid.col(col)!.width = this.viewerGrid.props.rowHeight;
    //   this.viewerGrid.col(col)!.width = 30;
    // }

    // this.viewerGrid.onCellRender.subscribe((args) => {
    //   if (args.cell.isRowHeader && args.cell.gridColumn.visible) {
    //     args.cell.gridColumn.visible = false;
    //     args.preventDefault();
    //   }
    // });

    // this.viewerGrid.props.allowEdit = false;

    // table.onCurrentCellChanged.subscribe((_) => {
    //   if (table.currentCol?.name !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE && table.currentCell.value !== null) {
    //     const aar = table.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, table.currentRowIdx);
    //     const pos = parseInt(table.currentCol.name);
    //     const currentCase = tableCases[aar][pos];
    //     const tempDfLength = currentCase.length;
    //     const initCol = DG.Column.string('Initial', tempDfLength);
    //     const subsCol = DG.Column.string('Substituted', tempDfLength);

    //     const tempDf = DG.DataFrame.fromColumns([
    //       initCol,
    //       subsCol,
    //       DG.Column.float('Difference', tempDfLength),
    //     ]);

    //     for (let i = 0; i < tempDfLength; i++) {
    //       const row = currentCase[i];
    //       tempDf.rows.setValues(i, [col.get(row[0]), col.get(row[1]), row[2]]);
    //     }

    //     tempDf.temp['isReal'] = true;

    //     initCol.semType = 'alignedSequence';
    //     initCol.temp['isAnalysisApplicable'] = false;
    //     subsCol.semType = 'alignedSequence';
    //     subsCol.temp['isAnalysisApplicable'] = false;

    //     grok.shell.o = DG.SemanticValue.fromValueType(tempDf, 'Substitution');
    //   }
    // });
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

const groupDescription: {[key: string]: {'description': string, aminoAcids: string[]}} = {
  'PC': {'description': 'Positive Amino Acids, with Electrically Charged Side Chains', aminoAcids: ['R', 'H', 'K']},
  'NC': {'description': 'Negative Amino Acids, with Electrically Charged Side Chains', aminoAcids: ['D', 'E']},
  'U': {'description': 'Amino Acids with Polar Uncharged Side Chains', aminoAcids: ['S', 'T', 'N', 'Q']},
  'SC': {'description': 'Special Cases', aminoAcids: ['C', 'U', 'G', 'P']},
  'H': {
    'description': 'Amino Acids with Hydrophobic Side Chain',
    aminoAcids: ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'],
  },
  '-': {'description': 'Unknown Amino Acid', aminoAcids: ['-']},
};

function joinDataFrames(df: DG.DataFrame, positionColumns: string[], splitSeqDf: DG.DataFrame) {
  // append splitSeqDf columns to source table and make sure columns are not added more than once
  const dfColsSet = new Set((df.columns as DG.ColumnList).names());
  if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
    df.join(
      splitSeqDf, [C.COLUMNS_NAMES.ACTIVITY], [C.COLUMNS_NAMES.ACTIVITY], (df.columns as DG.ColumnList).names(),
      positionColumns, 'inner', true);
  }
}

function sortSourceGrid(sourceGrid: DG.Grid) {
  if (sourceGrid) {
    const colNames: DG.GridColumn[] = [];
    for (let i = 1; i < sourceGrid.columns.length; i++)
      colNames.push(sourceGrid.columns.byIndex(i)!);

    colNames.sort((a, b)=>{
      if (a.column!.semType == C.SEM_TYPES.AMINO_ACIDS) {
        if (b.column!.semType == C.SEM_TYPES.AMINO_ACIDS)
          return 0;
        return -1;
      }
      if (b.column!.semType == C.SEM_TYPES.AMINO_ACIDS)
        return 1;
      return 0;
    });
    sourceGrid.columns.setOrder(colNames.map((v) => v.name));
  }
}

export async function createScaledCol(
  activityScaling: string, df: DG.DataFrame, sourceGrid: DG.Grid, splitSeqDf: DG.DataFrame,
) {
  const [scaledDf, newColName] = await PeptidesController.scaleActivity(
    activityScaling, df, df.temp[C.COLUMNS_NAMES.ACTIVITY]);
  //TODO: make another func
  const scaledCol = scaledDf.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  (splitSeqDf.columns as DG.ColumnList).add(scaledCol);
  const oldScaledCol = df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  (df.columns as DG.ColumnList).replace(oldScaledCol, scaledCol);
  const gridCol = sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  if (gridCol !== null) {
    gridCol.name = newColName;
    df.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED] = newColName;
  }

  sourceGrid.columns.setOrder([newColName]);
}

async function calculateStatistics(
  matrixDf: DG.DataFrame, peptidesCount: number, splitSeqDf: DG.DataFrame, groupMapping: StringDictionary,
) {
  matrixDf = matrixDf.groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE])
    .add('count', C.COLUMNS_NAMES.ACTIVITY_SCALED, 'Count')
    .aggregate();

  const countThreshold = 4;
  //@ts-ignore: never gets old
  matrixDf.rows.filter((row) => row.Count >= countThreshold && row.Count <= peptidesCount - countThreshold);
  matrixDf = matrixDf.clone(matrixDf.filter);

  // calculate additional stats
  await (matrixDf.columns as DG.ColumnList).addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

  //calculate p-values based on t-test
  let pvalues: Float32Array = new Float32Array(matrixDf.rowCount).fill(1);
  const mdCol: DG.Column = (matrixDf.columns as DG.ColumnList).addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
  const pValCol: DG.Column = (matrixDf.columns as DG.ColumnList).addNewFloat(C.COLUMNS_NAMES.P_VALUE);
  for (let i = 0; i < matrixDf.rowCount; i++) {
    const position = matrixDf.get(C.COLUMNS_NAMES.POSITION, i);
    const aar = matrixDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i);

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] === aar);
    const currentActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [C.COLUMNS_NAMES.ACTIVITY_SCALED])
      .getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED)
      .toList();

    //@ts-ignore
    splitSeqDf.rows.select((row) => groupMapping[row[position]] !== aar);
    const otherActivity: number[] = splitSeqDf
      .clone(splitSeqDf.selection, [C.COLUMNS_NAMES.ACTIVITY_SCALED])
      .getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED)
      .toList();

    const testResult = tTest(currentActivity, otherActivity);
    // testResult = uTest(currentActivity, otherActivity);
    const currentMeanDiff = testResult[C.COLUMNS_NAMES.MEAN_DIFFERENCE]!;
    const pvalue = testResult[currentMeanDiff >= 0 ? 'p-value more' : 'p-value less'];

    mdCol.set(i, currentMeanDiff);
    pvalues[i] = pvalue;
  }

  pvalues = fdrcorrection(pvalues)[1];

  for (let i = 0; i < pvalues.length; ++i)
    pValCol.set(i, pvalues[i]);

  return matrixDf.clone();
}

async function setCategoryOrder(twoColorMode: boolean, statsDf: DG.DataFrame, matrixDf: DG.DataFrame) {
  const absMD = 'Absolute Mean difference';
  const sortArgument = twoColorMode ? absMD : C.COLUMNS_NAMES.MEAN_DIFFERENCE;
  if (twoColorMode)
    await (statsDf.columns as DG.ColumnList).addNewCalculated(absMD, 'Abs(${Mean difference})');

  const aarWeightsDf = statsDf.groupBy([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE]).sum(sortArgument, 'weight').aggregate();
  const aarList = aarWeightsDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE).toList();
  const getWeight = (aar: string) => aarWeightsDf
    .groupBy(['weight'])
    .where(`${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ${aar}`)
    .aggregate()
    .get('weight', 0);
  aarList.sort((first, second) => getWeight(second) - getWeight(first));

  matrixDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE).setCategoryOrder(aarList);
}

function createVerticalTable(statsDf: DG.DataFrame, twoColorMode: boolean) {
  // TODO: aquire ALL of the positions
  let sequenceDf = statsDf.groupBy([C.COLUMNS_NAMES.MEAN_DIFFERENCE, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, C.COLUMNS_NAMES.POSITION, 'Count', 'Ratio', C.COLUMNS_NAMES.P_VALUE])
    .where('pValue <= 0.1')
    .aggregate();

  let tempStats: DG.Stats;
  const maxAtPos: {[index: string]: number} = {};
  for (const pos of sequenceDf.getCol(C.COLUMNS_NAMES.POSITION).categories) {
    tempStats = DG.Stats.fromColumn(
      sequenceDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE),
      DG.BitSet.create(sequenceDf.rowCount, (i) => sequenceDf.get(C.COLUMNS_NAMES.POSITION, i) === pos),
    );
    maxAtPos[pos] = twoColorMode ?
      (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) : tempStats.max;
  }
  sequenceDf = sequenceDf.clone(DG.BitSet.create(sequenceDf.rowCount, (i) => {
    return sequenceDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, i) === maxAtPos[sequenceDf.get(C.COLUMNS_NAMES.POSITION, i)];
  }));

  return sequenceDf;
}

function createGrids(
  matrixDf: DG.DataFrame, positionColumns: string[], sequenceDf: DG.DataFrame, grouping: boolean,
) {
  const sarGrid = matrixDf.plot.grid();
  sarGrid.sort([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE]);
  sarGrid.columns.setOrder([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE].concat(positionColumns as C.COLUMNS_NAMES[]));

  const sarVGrid = sequenceDf.plot.grid();
  sarVGrid.sort([C.COLUMNS_NAMES.POSITION]);
  sarVGrid.col(C.COLUMNS_NAMES.P_VALUE)!.format = 'four digits after comma';
  sarVGrid.col(C.COLUMNS_NAMES.P_VALUE)!.name = 'P-Value';

  if (!grouping) {
    let tempCol = (matrixDf.columns as DG.ColumnList).byName(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);

    tempCol = (sequenceDf.columns as DG.ColumnList).byName(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);
  }

  return [sarGrid, sarVGrid];
}

function setCellRendererFunc(
  renderColNames: string[], statsDf: DG.DataFrame, twoColorMode: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid,
  isSubstitutionOn: boolean, substTable?: DG.DataFrame,
) {
  const mdCol = statsDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
  const cellRendererFunc = function(args: DG.GridCellRenderArgs) {
    const canvasContext = args.g;
    const bound = args.bounds;
    const cell = args.cell;
    const tableColName = cell.tableColumn?.name;
    const tableRowIndex = cell.tableRowIndex!;
    const cellValue = cell.cell.value;
    const midX = bound.x + bound.width / 2;
    const midY = bound.y + bound.height / 2;

    canvasContext.save();
    canvasContext.beginPath();
    canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
    canvasContext.clip();

    if (cell.isRowHeader && cell.gridColumn.visible) {
      cell.gridColumn.visible = false;
      args.preventDefault();
      return;
    }

    if (tableColName === C.COLUMNS_NAMES.MEAN_DIFFERENCE)
      console.log('md col');

    if (cell.isTableCell && tableColName && tableRowIndex !== null && renderColNames.indexOf(tableColName) !== -1) {
      const gridTable = cell.grid.table;
      const currentPosition = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
        tableColName : gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
      const currentAAR = gridTable.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, tableRowIndex);
      const queryAAR = `${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ${currentAAR}`;
      if (cellValue) {
        const query = `${queryAAR} and ${C.COLUMNS_NAMES.POSITION} = ${currentPosition}`;
        const pVal: number = statsDf
          .groupBy([C.COLUMNS_NAMES.P_VALUE])
          .where(query)
          .aggregate()
          .get(C.COLUMNS_NAMES.P_VALUE, 0);

        let coef: string;
        const variant = cellValue < 0;
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
        const chooseCurrent = () => twoColorMode ? Math.abs(cellValue) : cellValue;

        const rCoef = (chooseCurrent() - chooseMin()) / (chooseMax() - chooseMin());

        const maxRadius = 0.9 * (bound.width > bound.height ? bound.height : bound.width) / 2;
        const radius = Math.floor(maxRadius * rCoef);

        canvasContext.beginPath();
        canvasContext.fillStyle = coef;
        canvasContext.arc(midX, midY, radius < 3 ? 3 : radius, 0, Math.PI * 2, true);
        canvasContext.closePath();

        canvasContext.fill();
      }
      if (isSubstitutionOn) {
        canvasContext.textBaseline = 'middle';
        canvasContext.textAlign = 'center';
        canvasContext.fillStyle = DG.Color.toHtml(DG.Color.black);
        canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
        const substValue = substTable?.groupBy([currentPosition])
          .where(queryAAR)
          .aggregate()
          .get(currentPosition, 0);
        if (substValue && substValue !== DG.INT_NULL)
          canvasContext.fillText(substValue, midX, midY);
      }
      args.preventDefault();
    }
    canvasContext.restore();
  };
  sarGrid.onCellRender.subscribe(cellRendererFunc);
  sarVGrid.onCellRender.subscribe(cellRendererFunc);
}

function setTooltipFunc(
  renderColNames: string[], statsDf: DG.DataFrame, peptidesCount: number, grouping: boolean, sarGrid: DG.Grid,
  sarVGrid: DG.Grid, sourceDf: DG.DataFrame,
) {
  const onCellTooltipFunc = async function(cell: DG.GridCell, x: number, y: number) {
    if (
      !cell.isRowHeader && !cell.isColHeader && cell.tableColumn !== null && cell.cell.value !== null &&
        cell.tableRowIndex !== null && renderColNames.indexOf(cell.tableColumn.name) !== -1) {
      const tooltipMap: { [index: string]: string } = {};

      for (const col of (statsDf.columns as DG.ColumnList).names()) {
        if (col !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE && col !== C.COLUMNS_NAMES.POSITION) {
          const currentPosition = cell.tableColumn.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
            cell.tableColumn.name : cell.grid.table.get(C.COLUMNS_NAMES.POSITION, cell.tableRowIndex);
          const query =
            `${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ${cell.grid.table.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, cell.tableRowIndex)} ` +
            `and ${C.COLUMNS_NAMES.POSITION} = ${currentPosition}`;
          const textNum = statsDf.groupBy([col]).where(query).aggregate().get(col, 0);
          let text = `${col === 'Count' ? textNum : textNum.toFixed(5)}`;

          if (col === 'Count')
            text += ` / ${peptidesCount}`;
          else if (col === C.COLUMNS_NAMES.P_VALUE)
            text = parseFloat(text) !== 0 ? text : '<0.01';


          tooltipMap[col === C.COLUMNS_NAMES.P_VALUE ? 'p-value' : col] = text;
        }
      }

      ui.tooltip.show(ui.tableFromMap(tooltipMap), x, y);
    }
    if (!cell.isColHeader && cell.tableColumn?.name == C.COLUMNS_NAMES.AMINO_ACID_RESIDUE) {
      if (grouping) {
        const currentGroup = groupDescription[cell.cell.value];
        const divText = ui.divText('Amino Acids in this group: ' + currentGroup[C.SEM_TYPES.AMINO_ACIDS].join(', '));
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
  sourceGrid: DG.Grid, invalidIndexes: number[], grouping: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid,
) {
  sourceGrid.onCellPrepare((cell: DG.GridCell) => {
    const currentRowIndex = cell.tableRowIndex;
    if (currentRowIndex && invalidIndexes.includes(currentRowIndex) && !cell.isRowHeader)
      cell.style.backColor = DG.Color.lightLightGray;
  });

  const mdCol: DG.GridColumn = sarVGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
  mdCol.name = 'Diff';

  for (const grid of [sarGrid, sarVGrid]) {
    grid.props.rowHeight = 20;
    grid.columns.rowHeader!.width = 20;
    for (let i = 0; i < grid.columns.length; ++i) {
      const col = grid.columns.byIndex(i)!;
      if (grid == sarVGrid && col.name !== 'Diff' && col.name !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)
        col.width = 45;
      else
        col.width = grid.props.rowHeight;
    }
  }

  if (grouping) {
    sarGrid.col(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)!.name = 'Groups';
    sarVGrid.col(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)!.name = 'Groups';
  }

  sarGrid.props.allowEdit = false;
  sarVGrid.props.allowEdit = false;
}

export function split(peptideColumn: DG.Column, filter: boolean = true): string[][] {
  const splitPeptidesArray: string[][] = [];
  let currentSplitPeptide: string[];
  let modeMonomerCount = 0;
  let currentLength;
  const colLength = peptideColumn.length;

  // splitting data
  const monomerLengths: { [index: string]: number } = {};
  for (let i = 0; i < colLength; i++) {
    currentSplitPeptide = peptideColumn.get(i).split('-').map((value: string) => value ? value : '-');
    splitPeptidesArray.push(currentSplitPeptide);
    currentLength = currentSplitPeptide.length;
    monomerLengths[currentLength + ''] =
      monomerLengths[currentLength + ''] ? monomerLengths[currentLength + ''] + 1 : 1;
  }
  //@ts-ignore: what I do here is converting string to number the most effective way I could find. parseInt is slow
  modeMonomerCount = 1 * Object.keys(monomerLengths).reduce((a, b) => monomerLengths[a] > monomerLengths[b] ? a : b);

  // making sure all of the sequences are of the same size
  // and marking invalid sequences
  let nTerminal: string;
  const invalidIndexes: number[] = [];
  let splitColumns: string[][] = Array.from({length: modeMonomerCount}, (_) => []);
  modeMonomerCount--; // minus N-terminal
  for (let i = 0; i < colLength; i++) {
    currentSplitPeptide = splitPeptidesArray[i];
    nTerminal = currentSplitPeptide.pop()!; // it is guaranteed that there will be at least one element
    currentLength = currentSplitPeptide.length;
    if (currentLength !== modeMonomerCount)
      invalidIndexes.push(i);

    for (let j = 0; j < modeMonomerCount; j++)
      splitColumns[j].push(j < currentLength ? currentSplitPeptide[j] : '-');

    splitColumns[modeMonomerCount].push(nTerminal);
  }
  modeMonomerCount--; // minus C-terminal

  //create column names list
  const columnNames = Array.from({length: modeMonomerCount}, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1}`);
  columnNames.splice(0, 0, 'N-terminal');
  columnNames.push('C-terminal');

  // filter out the columns with the same values
  if (filter) {
    splitColumns = splitColumns.filter((positionArray, index) => {
      const isRetained = new Set(positionArray).size > 1;
      if (!isRetained)
        columnNames.splice(index, 1);

      return isRetained;
    });
  }

  return splitPeptidesArray;
}
