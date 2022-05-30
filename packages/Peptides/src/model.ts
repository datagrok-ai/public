import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Subject, Observable} from 'rxjs';
import {addViewerToHeader, StackedBarChart} from './viewers/stacked-barchart-viewer';
import {PeptidesController} from './peptides';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {ChemPalette} from './utils/chem-palette';
import {MonomerLibrary} from './monomer-library';
import * as C from './utils/constants';
import * as type from './utils/types';
import {FilteringStatistics} from './utils/filtering-statistics';
import {getTypedArrayConstructor, stringToBool} from './utils/misc';

export class PeptidesModel {
  static _modelName = 'peptidesModel';

  _statsDataFrameSubject = new Subject<DG.DataFrame>();
  _sarGridSubject = new Subject<DG.Grid>();
  _sarVGridSubject = new Subject<DG.Grid>();
  _substitutionTableSubject = new Subject<type.SubstitutionsInfo>();

  _isUpdating: boolean = false;
  _isSubstInitialized = false;
  isBitsetChangedInitialized = false;

  //viewer properties
  _filterMode!: boolean;
  _twoColorMode!: boolean;
  _activityScaling!: string;
  _isSubstitutionOn!: boolean;
  _activityLimit!: number;
  _maxSubstitutions!: number;

  _sarGrid!: DG.Grid;
  _sarVGrid!: DG.Grid;
  _sourceGrid!: DG.Grid;
  _dataFrame: DG.DataFrame;
  splitCol!: DG.Column;
  stackedBarchart!: StackedBarChart;

  substitutionsInfo: type.SubstitutionsInfo = new Map();

  private constructor(dataFrame: DG.DataFrame) {
    this._dataFrame = dataFrame;
    this._dataFrame.temp[C.PEPTIDES_ANALYSIS] = true;

    this.updateProperties();
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  static get modelName(): string {return PeptidesModel._modelName;}

  get dataFrame(): DG.DataFrame {return this._dataFrame;}

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {return this._statsDataFrameSubject.asObservable();}

  get onSARGridChanged(): Observable<DG.Grid> {return this._sarGridSubject.asObservable();}

  get onSARVGridChanged(): Observable<DG.Grid> {return this._sarVGridSubject.asObservable();}

  get onSubstTableChanged(): Observable<type.SubstitutionsInfo> {return this._substitutionTableSubject.asObservable();}

  updateProperties(): void {
    this._activityScaling = this._dataFrame.tags['scaling'];
    this._filterMode = stringToBool(this._dataFrame.tags['filterMode']);
    this._twoColorMode = stringToBool(this._dataFrame.tags['bidirectionalAnalysis']);
    this._isSubstitutionOn = stringToBool(this._dataFrame.tags['showSubstitution']);
    this._maxSubstitutions = parseInt(this._dataFrame.tags['maxSubstitutions']);
    this._activityLimit = parseFloat(this._dataFrame.tags['activityLimit']);
  }

  setProperties(
    activityScaling: string, filterMode: boolean, twoColorMode: boolean, isSubstitutionOn: boolean,
    maxSubstitutions: number, activityLimit: number, forceUpdate = false,
  ): void {
    const chooseAction =
      (value: string, defaultValue: string | boolean | number): string | boolean | number =>
        forceUpdate ? value : defaultValue ?? value;
    this._dataFrame.tags['scaling'] = chooseAction(`${activityScaling}`, this._dataFrame.tags['scaling']);
    this._dataFrame.tags['filterMode'] = chooseAction(`${filterMode}`, this._dataFrame.tags['filterMode']);
    this._dataFrame.tags['bidirectionalAnalysis'] =
      chooseAction(`${twoColorMode}`, this._dataFrame.tags['bidirectionalAnalysis']);
    this._dataFrame.tags['showSubstitution'] =
      chooseAction(`${isSubstitutionOn}`, this._dataFrame.tags['showSubstitution']);
    this._dataFrame.tags['maxSubstitutions'] =
      chooseAction(`${maxSubstitutions}`, this._dataFrame.tags['maxSubstitutions']);
    this._dataFrame.tags['activityLimit'] = chooseAction(`${activityLimit}`, this._dataFrame.tags['activityLimit']);

    this.updateProperties();
  }

  async updateData(
    activityScaling?: string, sourceGrid?: DG.Grid, twoColorMode?: boolean, activityLimit?: number,
    maxSubstitutions?: number, isSubstitutionOn?: boolean, filterMode?: boolean,
  ): Promise<void> {
    //FIXME: threre are too many assignments, some are duplicating
    this._activityScaling = activityScaling ?? this._activityScaling;
    this._sourceGrid = sourceGrid ?? this._sourceGrid;
    this._twoColorMode = twoColorMode ?? this._twoColorMode;
    this._activityLimit = activityLimit ?? this._activityLimit;
    this._maxSubstitutions = maxSubstitutions ?? this._maxSubstitutions;
    this._isSubstitutionOn = isSubstitutionOn ?? this._isSubstitutionOn;
    this._filterMode = filterMode ?? this._filterMode;
    this.setProperties(this._activityScaling, this._filterMode, this._twoColorMode, this._isSubstitutionOn,
      this._maxSubstitutions, this._activityLimit, true);

    await this.updateDefault();
  }

  async updateDefault(): Promise<void> {
    if (this._activityScaling && this._sourceGrid && this._twoColorMode !== null && !this._isUpdating) {
      this._isUpdating = true;
      const [viewerGrid, viewerVGrid, statsDf] = await this.initializeViewersComponents();
      //FIXME: modify during the initializeViewersComponents stages
      this._statsDataFrameSubject.next(statsDf);
      this._sarGridSubject.next(viewerGrid);
      this._sarVGridSubject.next(viewerVGrid);
      if (this._isSubstitutionOn) {
        this._substitutionTableSubject.next(this.substitutionsInfo);
        this._isSubstInitialized = true;
      }
    }
    await this.updateBarchart();
    this.invalidateGrids();

    this._isUpdating = false;
  }

  async updateBarchart(): Promise<void> {
    this.stackedBarchart ??= await this._dataFrame?.plot.fromType('StackedBarChartAA') as StackedBarChart;
    if (this.stackedBarchart && this._sourceGrid)
      addViewerToHeader(this._sourceGrid, this.stackedBarchart);
  }

  async initializeViewersComponents(): Promise<[DG.Grid, DG.Grid, DG.DataFrame]> {
    if (this._sourceGrid === null)
      throw new Error(`Source grid is not initialized`);

    //Split the aligned sequence into separate AARs
    let splitSeqDf: DG.DataFrame | undefined;
    let invalidIndexes: number[];
    const col: DG.Column = this._dataFrame.columns.bySemType(C.SEM_TYPES.ALIGNED_SEQUENCE)!;
    [splitSeqDf, invalidIndexes] = PeptidesController.splitAlignedPeptides(col);

    const positionColumns = splitSeqDf.columns.names();
    const renderColNames: string[] = splitSeqDf.columns.names();

    splitSeqDf.columns.add(this._dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY));

    this.joinDataFrames(positionColumns, splitSeqDf);

    for (const dfCol of this._dataFrame.columns) {
      if (splitSeqDf.col(dfCol.name) && dfCol.name != C.COLUMNS_NAMES.ACTIVITY)
        PeptidesController.setAARRenderer(dfCol, this._sourceGrid);
    }

    this.sortSourceGrid(this._sourceGrid);

    await this.createScaledCol(this._activityScaling, this._dataFrame, this._sourceGrid, splitSeqDf);

    //unpivot a table and handle duplicates
    splitSeqDf = splitSeqDf.groupBy(positionColumns)
      .add('med', C.COLUMNS_NAMES.ACTIVITY_SCALED, C.COLUMNS_NAMES.ACTIVITY_SCALED)
      .aggregate();

    const peptidesCount = splitSeqDf.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).length;

    let matrixDf = splitSeqDf.unpivot(
      [C.COLUMNS_NAMES.ACTIVITY_SCALED], positionColumns, C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);

    //FIXME: for some reason Mean difference is not calculated for all the AARs
    //statistics for specific AAR at a specific position
    const statsDf = await this.calculateStatistics(matrixDf, peptidesCount, splitSeqDf);

    // SAR matrix table
    //pivot a table to make it matrix-like
    matrixDf = statsDf.groupBy([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.MEAN_DIFFERENCE, '')
      .aggregate();
    matrixDf.name = 'SAR';

    // Setting category order
    await this.setCategoryOrder(this._twoColorMode, statsDf, matrixDf);

    // SAR vertical table (naive, choose best Mean difference from pVals <= 0.01)
    const sequenceDf = this.createVerticalTable(statsDf, this._twoColorMode);
    renderColNames.push(C.COLUMNS_NAMES.MEAN_DIFFERENCE);


    if (this._isSubstitutionOn || !this._isSubstInitialized)
      this.calcSubstitutions();

    //TODO: move everything below out to controller
    const [sarGrid, sarVGrid] = this.createGrids(matrixDf, positionColumns, sequenceDf);

    this._sarGrid = sarGrid;
    this._sarVGrid = sarVGrid;

    this.setCellRenderers(
      renderColNames, statsDf, this._twoColorMode, sarGrid, sarVGrid, this._isSubstitutionOn);

    // show all the statistics in a tooltip over cell
    this.setTooltips(renderColNames, statsDf, peptidesCount, sarGrid, sarVGrid, this._dataFrame);

    this.setInteractionCallback();

    this.modifyOrCreateSplitCol(C.CATEGORIES.ALL, C.CATEGORIES.ALL);

    this.setBitsetCallback();

    this.postProcessGrids(this._sourceGrid, invalidIndexes, sarGrid, sarVGrid);

    if (this.dataFrame.tags[C.TAGS.AAR] !== '' && this.dataFrame.tags[C.TAGS.POSITION] !== '') {
      const sarDf = sarGrid.dataFrame;
      const rowCount = sarDf.rowCount;
      let index = -1;
      for (let i = 0; i < rowCount; i++) {
        if (sarDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i) === this.dataFrame.tags[C.TAGS.AAR]) {
          index = i;
          break;
        }
      }
      sarDf.currentCell = sarDf.cell(index, this.dataFrame.tags[C.TAGS.POSITION]);
    }

    //TODO: return class instead
    return [sarGrid, sarVGrid, statsDf];
  }

  //TODO: move to controller?
  calcSubstitutions(): void {
    const activityValues: DG.Column<number> = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    const columnList: DG.Column<string>[] = this.dataFrame.columns.bySemTypeAll(C.SEM_TYPES.AMINO_ACIDS);
    const nCols = columnList.length;
    if (nCols == 0)
      throw new Error(`Couldn't find any column of semType '${C.SEM_TYPES.AMINO_ACIDS}'`);

    this.substitutionsInfo = new Map();
    const nRows = this.dataFrame.rowCount;
    for (let seq1Idx = 0; seq1Idx < nRows - 1; seq1Idx++) {
      for (let seq2Idx = seq1Idx + 1; seq2Idx < nRows; seq2Idx++) {
        let substCounter = 0;
        const activityValSeq1 = activityValues.get(seq1Idx);
        const activityValSeq2 = activityValues.get(seq2Idx);
        const delta = activityValSeq1 - activityValSeq2;
        if (Math.abs(delta) < this._activityLimit)
          continue;

        let substCounterFlag = false;
        const tempData: {pos: string, seq1monomer: string, seq2monomer: string, seq1Idx: number, seq2Idx: number}[] =
          [];
        for (const currentPosCol of columnList) {
          const seq1monomer = currentPosCol.get(seq1Idx);
          const seq2monomer = currentPosCol.get(seq2Idx);
          if (seq1monomer == seq2monomer)
            continue;

          substCounter++;
          substCounterFlag = substCounter > this._maxSubstitutions;
          if (substCounterFlag)
            break;

          tempData.push({
            pos: currentPosCol.name,
            seq1monomer: seq1monomer,
            seq2monomer: seq2monomer,
            seq1Idx: seq1Idx,
            seq2Idx: seq2Idx,
          });
        }

        if (substCounterFlag || substCounter == 0)
          continue;

        for (const tempDataElement of tempData) {
          const position = tempDataElement.pos;

          //Working with seq1monomer
          const seq1monomer = tempDataElement.seq1monomer;
          if (!this.substitutionsInfo.has(seq1monomer))
            this.substitutionsInfo.set(seq1monomer, new Map());

          let positionsMap = this.substitutionsInfo.get(seq1monomer)!;
          if (!positionsMap.has(position))
            positionsMap.set(position, new Map());

          let indexes = positionsMap.get(position)!;

          !indexes.has(seq1Idx) ? indexes.set(seq1Idx, [seq2Idx]) : (indexes.get(seq1Idx)! as number[]).push(seq2Idx);

          //Working with seq2monomer
          const seq2monomer = tempDataElement.seq2monomer;
          if (!this.substitutionsInfo.has(seq2monomer))
            this.substitutionsInfo.set(seq2monomer, new Map());

          positionsMap = this.substitutionsInfo.get(seq2monomer)!;
          if (!positionsMap.has(position))
            positionsMap.set(position, new Map());

          indexes = positionsMap.get(position)!;
          !indexes.has(seq2Idx) ? indexes.set(seq2Idx, [seq1Idx]) : (indexes.get(seq2Idx)! as number[]).push(seq1Idx);
        }
      }
    }

    const TypedArray = getTypedArrayConstructor(nRows);
    for (const positionMap of this.substitutionsInfo.values()) {
      for (const indexMap of positionMap.values()) {
        for (const [index, indexArray] of indexMap.entries())
          indexMap.set(index, new TypedArray(indexArray));
      }
    }
  }

  joinDataFrames(positionColumns: string[], splitSeqDf: DG.DataFrame): void {
    // append splitSeqDf columns to source table and make sure columns are not added more than once
    const dfColsSet = new Set(this._dataFrame.columns.names());
    if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
      this._dataFrame.join(splitSeqDf, [C.COLUMNS_NAMES.ACTIVITY], [C.COLUMNS_NAMES.ACTIVITY],
        this._dataFrame.columns.names(), positionColumns, 'inner', true);
    }
  }

  sortSourceGrid(sourceGrid: DG.Grid): void {
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

  //TODO: move out, merge with controller's scaleActivity
  async createScaledCol(
    activityScaling: string, df: DG.DataFrame, sourceGrid: DG.Grid, splitSeqDf: DG.DataFrame,
  ): Promise<void> {
    const [scaledDf, newColName] = await PeptidesController.scaleActivity(
      activityScaling, df, df.temp[C.COLUMNS_NAMES.ACTIVITY]);
    //TODO: make another func
    const scaledCol = scaledDf.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    splitSeqDf.columns.add(scaledCol);
    const oldScaledCol = df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    df.columns.replace(oldScaledCol, scaledCol);
    const gridCol = sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    if (gridCol !== null) {
      gridCol.name = newColName;
      df.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED] = newColName;
    }

    sourceGrid.columns.setOrder([newColName]);
  }

  //TODO: move out
  async calculateStatistics(
    matrixDf: DG.DataFrame, peptidesCount: number, splitSeqDf: DG.DataFrame): Promise<DG.DataFrame> {
    matrixDf = matrixDf.groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE])
      .add('count', C.COLUMNS_NAMES.ACTIVITY_SCALED, 'Count')
      .aggregate();

    const countThreshold = 4;
    matrixDf.rows.filter((row) => row.Count >= countThreshold && row.Count <= peptidesCount - countThreshold);
    matrixDf = matrixDf.clone(matrixDf.filter);

    // calculate additional stats
    await matrixDf.columns.addNewCalculated('Ratio', '${count}/'.concat(`${peptidesCount}`));

    //calculate p-values based on t-test
    let pvalues: Float32Array = new Float32Array(matrixDf.rowCount).fill(1);
    const mdCol: DG.Column = matrixDf.columns.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const pValCol: DG.Column = matrixDf.columns.addNewFloat(C.COLUMNS_NAMES.P_VALUE);
    for (let i = 0; i < matrixDf.rowCount; i++) {
      const position = matrixDf.get(C.COLUMNS_NAMES.POSITION, i);
      const aar = matrixDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i);

      splitSeqDf.rows.select((row) => row[position] === aar);
      const currentActivity: number[] = splitSeqDf
        .clone(splitSeqDf.selection, [C.COLUMNS_NAMES.ACTIVITY_SCALED])
        .getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED)
        .toList();

      splitSeqDf.rows.select((row) => row[position] !== aar);
      const otherActivity: number[] = splitSeqDf
        .clone(splitSeqDf.selection, [C.COLUMNS_NAMES.ACTIVITY_SCALED])
        .getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED)
        .toList();

      const testResult = tTest(currentActivity, otherActivity);
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

  async setCategoryOrder(twoColorMode: boolean, statsDf: DG.DataFrame, matrixDf: DG.DataFrame): Promise<void> {
    const absMD = 'Absolute Mean difference';
    const sortArgument = twoColorMode ? absMD : C.COLUMNS_NAMES.MEAN_DIFFERENCE;
    if (twoColorMode)
      await statsDf.columns.addNewCalculated(absMD, 'Abs(${Mean difference})');

    const aarWeightsDf = statsDf.groupBy([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE]).sum(sortArgument, 'weight').aggregate();
    const aarList = aarWeightsDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE).toList();
    const getWeight = (aar: string): number => aarWeightsDf
      .groupBy(['weight'])
      .where(`${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ${aar}`)
      .aggregate()
      .get('weight', 0) as number;
    aarList.sort((first, second) => getWeight(second) - getWeight(first));

    matrixDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE).setCategoryOrder(aarList);
  }

  createVerticalTable(statsDf: DG.DataFrame, twoColorMode: boolean): DG.DataFrame {
    // TODO: aquire ALL of the positions
    const columns = [C.COLUMNS_NAMES.MEAN_DIFFERENCE, C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, C.COLUMNS_NAMES.POSITION,
      'Count', 'Ratio', C.COLUMNS_NAMES.P_VALUE];
    let sequenceDf = statsDf.groupBy(columns)
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
        (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) :
        tempStats.max;
    }
    sequenceDf = sequenceDf.clone(DG.BitSet.create(sequenceDf.rowCount, (i) =>
      sequenceDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, i) === maxAtPos[sequenceDf.get(C.COLUMNS_NAMES.POSITION, i)]));

    return sequenceDf;
  }

  createGrids(matrixDf: DG.DataFrame, positionColumns: string[], sequenceDf: DG.DataFrame): DG.Grid[] {
    const sarGrid = matrixDf.plot.grid();
    sarGrid.sort([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE]);
    sarGrid.columns.setOrder([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE].concat(positionColumns as C.COLUMNS_NAMES[]));

    const sarVGrid = sequenceDf.plot.grid();
    sarVGrid.sort([C.COLUMNS_NAMES.POSITION]);
    sarVGrid.col(C.COLUMNS_NAMES.P_VALUE)!.format = 'four digits after comma';
    sarVGrid.col(C.COLUMNS_NAMES.P_VALUE)!.name = 'P-Value';

    let tempCol = matrixDf.columns.byName(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);

    tempCol = sequenceDf.columns.byName(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      PeptidesController.setAARRenderer(tempCol, sarGrid);

    return [sarGrid, sarVGrid];
  }

  //TODO: move out
  setCellRenderers(
    renderColNames: string[], statsDf: DG.DataFrame, twoColorMode: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid,
    isSubstitutionOn: boolean,
  ): void {
    const mdCol = statsDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const cellRendererAction = (args: DG.GridCellRenderArgs): void => {
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

      if (cell.isTableCell && tableColName && tableRowIndex !== null && renderColNames.indexOf(tableColName) !== -1) {
        const gridTable = cell.grid.table;
        const currentPosition = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
          tableColName : gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
        const currentAAR = gridTable.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, tableRowIndex);
        if (currentAAR === 'Aib' && currentPosition === '02')
          console.log('stop');

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


          const chooseMin = (): number => twoColorMode ? 0 : mdCol.min;
          const chooseMax = (): number => twoColorMode ? Math.max(Math.abs(mdCol.min), mdCol.max) : mdCol.max;
          const chooseCurrent = (): any => twoColorMode ? Math.abs(cellValue) : cellValue;

          const rCoef = (chooseCurrent() - chooseMin()) / (chooseMax() - chooseMin());

          const maxRadius = 0.9 * (bound.width > bound.height ? bound.height : bound.width) / 2;
          const radius = Math.floor(maxRadius * rCoef);

          canvasContext.beginPath();
          canvasContext.fillStyle = coef;
          canvasContext.arc(midX, midY, radius < 3 ? 3 : radius, 0, Math.PI * 2, true);
          canvasContext.closePath();

          canvasContext.fill();
          if (isSubstitutionOn) {
            canvasContext.textBaseline = 'middle';
            canvasContext.textAlign = 'center';
            canvasContext.fillStyle = DG.Color.toHtml(DG.Color.getContrastColor(DG.Color.fromHtml(coef)));
            canvasContext.font = '13px Roboto, Roboto Local, sans-serif';
            let substValue = 0;
            this.substitutionsInfo.get(currentAAR)?.get(currentPosition)?.forEach((idxs) => substValue += idxs.length);
            if (substValue && substValue != 0)
              canvasContext.fillText(substValue.toString(), midX, midY);
          }
        }
        args.preventDefault();
      }
      canvasContext.restore();
    };
    sarGrid.onCellRender.subscribe(cellRendererAction);
    sarVGrid.onCellRender.subscribe(cellRendererAction);
  }

  //TODO: move out
  setTooltips(
    renderColNames: string[], statsDf: DG.DataFrame, peptidesCount: number, sarGrid: DG.Grid, sarVGrid: DG.Grid,
    sourceDf: DG.DataFrame,
  ): void {
    const onCellTooltipAction = async (cell: DG.GridCell, x: number, y: number): Promise<boolean> => {
      if (
        !cell.isRowHeader && !cell.isColHeader && cell.tableColumn !== null && cell.cell.value !== null &&
          cell.tableRowIndex !== null && renderColNames.indexOf(cell.tableColumn.name) !== -1) {
        const tooltipMap: { [index: string]: string } = {};

        for (const col of statsDf.columns.names()) {
          if (col !== C.COLUMNS_NAMES.AMINO_ACID_RESIDUE && col !== C.COLUMNS_NAMES.POSITION) {
            const currentPosition = cell.tableColumn.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
              cell.tableColumn.name : cell.grid.table.get(C.COLUMNS_NAMES.POSITION, cell.tableRowIndex);
            const query =
              `${C.COLUMNS_NAMES.AMINO_ACID_RESIDUE} = ` +
              `${cell.grid.table.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, cell.tableRowIndex)} ` +
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
        const monomerLib = sourceDf.temp[MonomerLibrary.id];
        ChemPalette.showTooltip(cell, x, y, monomerLib);
      }
      return true;
    };
    sarGrid.onCellTooltip(onCellTooltipAction);
    sarVGrid.onCellTooltip(onCellTooltipAction);
  }

  //TODO: think about it, move out?
  setInteractionCallback(): void {
    const sarDf = this._sarGrid.dataFrame;
    const sarVDf = this._sarVGrid.dataFrame;

    const getAARandPosition = (isVertical = false): [string, string] => {
      let aar : string;
      let position: string;
      if (isVertical) {
        const currentRowIdx = sarVDf.currentRowIdx;
        aar = sarVDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, currentRowIdx);
        position = sarVDf.get(C.COLUMNS_NAMES.POSITION, currentRowIdx);
      } else {
        aar = sarDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, sarDf.currentRowIdx);
        position = sarDf.currentCol.name;
      }
      return [aar, position];
    };

    this._sarGrid.onCurrentCellChanged.subscribe((gc) => {
      const isNegativeRowIndex = sarDf.currentRowIdx === -1;
      if (!sarDf.currentCol || (!sarDf.currentCell.value && !isNegativeRowIndex))
        return;
      this.syncGrids(false, sarDf, sarVDf);
      let aar: string = C.CATEGORIES.ALL;
      let position: string = C.CATEGORIES.ALL;
      if (!isNegativeRowIndex) {
        [aar, position] = getAARandPosition();
        this.dataFrame.tags[C.TAGS.AAR] = aar;
        this.dataFrame.tags[C.TAGS.POSITION] = position;
      } else
        this.dataFrame.tags[C.TAGS.AAR] = this.dataFrame.tags[C.TAGS.POSITION] = '';

      // this.dataFrame.temp['substTable'] = this.getSubstitutionTable();
      this.modifyOrCreateSplitCol(aar, position);
      this.fireBitsetChanged();
      this.invalidateGrids();
      grok.shell.o = this.dataFrame;
    });

    this._sarVGrid.onCurrentCellChanged.subscribe((gc) => {
      if (!sarVDf.currentCol || sarVDf.currentRowIdx === -1)
        return;
      this.syncGrids(true, sarDf, sarVDf);
    });
  }

  invalidateGrids(): void {
    this.stackedBarchart?.computeData();
    this._sarGrid.invalidate();
    this._sarVGrid.invalidate();
    this._sourceGrid?.invalidate();
    //TODO: this.peptideSpaceGrid.invalidate();
  }

  setBitsetCallback(): void {
    if (this.isBitsetChangedInitialized)
      return;
    const filter = this.dataFrame.filter;
    const selection = this.dataFrame.selection;

    const changeBitset = (currentBitset: DG.BitSet, previousBitset: DG.BitSet): void => {
      previousBitset.setAll(!this._filterMode, false);
      currentBitset.init((i) => {
        const currentCategory = this.splitCol.get(i);
        return currentCategory !== C.CATEGORIES.OTHER && currentCategory !== C.CATEGORIES.ALL;
      }, false);
    };

    const recalculateStatistics =
      (bitset: DG.BitSet): void => (this.dataFrame.temp[C.STATS] as FilteringStatistics).setMask(bitset);

    filter.onChanged.subscribe(() => {
      changeBitset(filter, selection);
      recalculateStatistics(filter);
    });
    selection.onChanged.subscribe(() => {
      changeBitset(selection, filter);
      recalculateStatistics(selection);
    });
    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(): void {(this._filterMode ? this._dataFrame.filter : this._dataFrame.selection).fireChanged();}

  //TODO: move out
  postProcessGrids(sourceGrid: DG.Grid, invalidIndexes: number[], sarGrid: DG.Grid, sarVGrid: DG.Grid): void {
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

    sarGrid.props.allowEdit = false;
    sarVGrid.props.allowEdit = false;
  }

  //TODO: refactor, use this.sarDf and accept aar & position as parameters
  syncGrids(sourceVertical: boolean, sarDf: DG.DataFrame, sarVDf: DG.DataFrame): void {
    let otherColName: string;
    let otherRowIndex: number;
    const otherDf = sourceVertical ? sarDf : sarVDf;

    if (otherDf.temp[C.FLAGS.CELL_CHANGING])
      return;

    //on vertical SAR viewer click
    if (sourceVertical) {
      const currentRowIdx = sarVDf.currentRowIdx;
      const currentColName = sarVDf.currentCol.name;
      if (currentColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      otherColName = sarVDf.get(C.COLUMNS_NAMES.POSITION, currentRowIdx);
      const otherRowName: string = sarVDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, currentRowIdx);
      otherRowIndex = -1;
      const rows = otherDf.rowCount;
      for (let i = 0; i < rows; i++) {
        if (otherDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i) === otherRowName) {
          otherRowIndex = i;
          break;
        }
      }
    //on SAR viewer click
    } else {
      otherColName = C.COLUMNS_NAMES.MEAN_DIFFERENCE;
      const otherPos: string = sarDf.currentCol.name;
      if (otherPos === C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)
        return;

      const otherAAR: string =
        sarDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, sarDf.currentRowIdx);
      otherRowIndex = -1;
      for (let i = 0; i < sarVDf.rowCount; i++) {
        if (
          sarVDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, i) === otherAAR &&
          sarVDf.get(C.COLUMNS_NAMES.POSITION, i) === otherPos
        ) {
          otherRowIndex = i;
          break;
        }
      }
    }
    otherDf.temp[C.FLAGS.CELL_CHANGING] = true;
    otherDf.currentCell = otherDf.cell(otherRowIndex, otherColName);
    otherDf.temp[C.FLAGS.CELL_CHANGING] = false;
  }

  getSplitColValueAt(index: number, aar: string, position: string, aarLabel: string): string {
    const currentAAR = this.dataFrame.get(position, index) as string;
    return currentAAR === aar ? aarLabel : C.CATEGORIES.OTHER;
  }

  modifyOrCreateSplitCol(aar: string, position: string): void {
    const df = this.dataFrame;
    this.splitCol = df.col(C.COLUMNS_NAMES.SPLIT_COL) ??df.columns.addNew(C.COLUMNS_NAMES.SPLIT_COL, 'string');

    if (aar === C.CATEGORIES.ALL && position === C.CATEGORIES.ALL) {
      this.splitCol.init(() => C.CATEGORIES.ALL);
      return;
    }

    const aarLabel = `${aar === '-' ? 'Gap' : aar} : ${position}`;
    this.splitCol.init((i) => this.getSplitColValueAt(i, aar, position, aarLabel));

    this.splitCol.setCategoryOrder([aarLabel]);
    this.splitCol.compact();

    const colorMap: {[index: string]: string | number} = {};

    colorMap[C.CATEGORIES.OTHER] = DG.Color.blue;
    colorMap[aarLabel] = DG.Color.orange;
    df.getCol(C.COLUMNS_NAMES.SPLIT_COL).colors.setCategorical(colorMap);
  }
}
