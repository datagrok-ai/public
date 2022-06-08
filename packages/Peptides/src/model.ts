import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {Subject, Observable} from 'rxjs';
import {addViewerToHeader, StackedBarChart} from './viewers/stacked-barchart-viewer';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {ChemPalette} from './utils/chem-palette';
import {MonomerLibrary} from './monomer-library';
import * as C from './utils/constants';
import * as type from './utils/types';
import {FilteringStatistics} from './utils/filtering-statistics';
import {getSeparator, getTypedArrayConstructor, stringToBool} from './utils/misc';
import {_package} from './package';
import {SARViewer, SARViewerVertical} from './viewers/sar-viewer';
import {PeptideSpaceViewer} from './viewers/peptide-space-viewer';
import {setAARRenderer} from './utils/cell-renderer';

export class PeptidesModel {
  static _modelName = 'peptidesModel';

  _statsDataFrameSubject = new Subject<DG.DataFrame>();
  _sarGridSubject = new Subject<DG.Grid>();
  _sarVGridSubject = new Subject<DG.Grid>();
  _substitutionTableSubject = new Subject<type.SubstitutionsInfo>();

  _isUpdating: boolean = false;
  _isSubstInitialized = false;
  isBitsetChangedInitialized = false;
  isCellChanging = false;

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
  splitCol!: DG.Column<boolean>;
  stackedBarchart!: StackedBarChart;

  substitutionsInfo: type.SubstitutionsInfo = new Map();
  isInitialized: boolean = false;
  currentView!: DG.TableView;

  _currentSelection!: type.SelectionObject;

  private constructor(dataFrame: DG.DataFrame) {
    this._dataFrame = dataFrame;
    this._dataFrame.temp[C.PEPTIDES_ANALYSIS] = true;

    this.updateProperties();
  }

  static async getInstance(dataFrame: DG.DataFrame, dgPackage?: DG.Package): Promise<PeptidesModel> {
    if (dataFrame.temp[MonomerLibrary.id] === null) {
      const sdf = await (dgPackage ?? _package).files.readAsText('HELMMonomers_June10.sdf');
      dataFrame.temp[MonomerLibrary.id] ??= new MonomerLibrary(sdf);
    }
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    await (dataFrame.temp[PeptidesModel.modelName] as PeptidesModel).init();
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  static get modelName(): string {return PeptidesModel._modelName;}

  get onStatsDataFrameChanged(): Observable<DG.DataFrame> {return this._statsDataFrameSubject.asObservable();}

  get onSARGridChanged(): Observable<DG.Grid> {return this._sarGridSubject.asObservable();}

  get onSARVGridChanged(): Observable<DG.Grid> {return this._sarVGridSubject.asObservable();}

  get onSubstTableChanged(): Observable<type.SubstitutionsInfo> {return this._substitutionTableSubject.asObservable();}

  get currentSelection(): type.SelectionObject {
    this._currentSelection ??= JSON.parse(this._dataFrame.tags[C.TAGS.SELECTION] || '{}');
    return this._currentSelection;
  }
  set currentSelection(selection: type.SelectionObject) {
    this._currentSelection = selection;
    this._dataFrame.tags[C.TAGS.SELECTION] = JSON.stringify(selection);

    this.modifyOrCreateSplitCol();
    this.fireBitsetChanged();
    this.invalidateGrids();
    grok.shell.o = this._dataFrame;
  }

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
    [splitSeqDf, invalidIndexes] = PeptidesModel.splitAlignedPeptides(col);

    const positionColumns = splitSeqDf.columns.names();
    const renderColNames: string[] = splitSeqDf.columns.names();

    const activityCol = this._dataFrame.columns.bySemType(C.SEM_TYPES.ACTIVITY)!;
    splitSeqDf.columns.add(activityCol);

    this.joinDataFrames(positionColumns, splitSeqDf);

    for (const dfCol of this._dataFrame.columns) {
      if (positionColumns.includes(dfCol.name))
        setAARRenderer(dfCol, this._sourceGrid);
    }

    this.sortSourceGrid(this._sourceGrid);

    await this.createScaledCol(this._activityScaling, this._dataFrame, this._sourceGrid, splitSeqDf);

    //unpivot a table and handle duplicates
    splitSeqDf = splitSeqDf.groupBy(positionColumns)
      .add('med', C.COLUMNS_NAMES.ACTIVITY_SCALED, C.COLUMNS_NAMES.ACTIVITY_SCALED)
      .aggregate();

    const peptidesCount = splitSeqDf.rowCount;

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

    this.modifyOrCreateSplitCol();

    this.setBitsetCallback();

    this.postProcessGrids(this._sourceGrid, invalidIndexes, sarGrid, sarVGrid);

    //TODO: return class instead
    return [sarGrid, sarVGrid, statsDf];
  }

  //TODO: move to controller?
  calcSubstitutions(): void {
    const activityValues: DG.Column<number> = this._dataFrame.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
    const columnList: DG.Column<string>[] = this._dataFrame.columns.bySemTypeAll(C.SEM_TYPES.AMINO_ACIDS);
    const nCols = columnList.length;
    if (nCols == 0)
      throw new Error(`Couldn't find any column of semType '${C.SEM_TYPES.AMINO_ACIDS}'`);

    this.substitutionsInfo = new Map();
    const nRows = this._dataFrame.rowCount;
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
    const name = this._dataFrame.name;
    const dfColsSet = new Set(this._dataFrame.columns.names());
    if (!positionColumns.every((col: string) => dfColsSet.has(col))) {
      this._dataFrame.join(splitSeqDf, [C.COLUMNS_NAMES.ACTIVITY], [C.COLUMNS_NAMES.ACTIVITY],
        this._dataFrame.columns.names(), positionColumns, 'inner', true);
    }
    this._dataFrame.name = name;
    this.currentView.name = name;
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
    const [scaledDf, newColName] = await PeptidesModel.scaleActivity(
      activityScaling, df, df.temp[C.COLUMNS_NAMES.ACTIVITY]);
    //TODO: make another func
    const scaledCol = scaledDf.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    scaledCol.semType = C.SEM_TYPES.ACTIVITY_SCALED;
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
    const aarCol = matrixDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    const posCol = matrixDf.getCol(C.COLUMNS_NAMES.POSITION);
    for (let i = 0; i < matrixDf.rowCount; i++) {
      const position = posCol.get(i);
      const aar = aarCol.get(i);

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
    const posColCategories = sequenceDf.getCol(C.COLUMNS_NAMES.POSITION).categories;
    const mdCol = sequenceDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const posCol = sequenceDf.getCol(C.COLUMNS_NAMES.POSITION);
    const rowCount = sequenceDf.rowCount;
    for (const pos of posColCategories) {
      tempStats = DG.Stats.fromColumn(mdCol, DG.BitSet.create(rowCount, (i) => posCol.get(i) === pos));
      maxAtPos[pos] = twoColorMode ?
        (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) :
        tempStats.max;
    }
    sequenceDf = sequenceDf.clone(DG.BitSet.create(rowCount, (i) => mdCol.get(i) === maxAtPos[posCol.get(i)]));

    return sequenceDf;
  }

  createGrids(matrixDf: DG.DataFrame, positionColumns: string[], sequenceDf: DG.DataFrame): DG.Grid[] {
    const sarGrid = matrixDf.plot.grid();
    sarGrid.sort([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE]);
    sarGrid.columns.setOrder([C.COLUMNS_NAMES.AMINO_ACID_RESIDUE].concat(positionColumns as C.COLUMNS_NAMES[]));

    const sarVGrid = sequenceDf.plot.grid();
    sarVGrid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = sarVGrid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = 'four digits after comma';
    pValGridCol.name = 'P-Value';

    let tempCol = matrixDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      setAARRenderer(tempCol, sarGrid);

    tempCol = sequenceDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    if (tempCol)
      setAARRenderer(tempCol, sarGrid);

    return [sarGrid, sarVGrid];
  }

  //TODO: move out
  setCellRenderers(
    renderColNames: string[], statsDf: DG.DataFrame, twoColorMode: boolean, sarGrid: DG.Grid, sarVGrid: DG.Grid,
    isSubstitutionOn: boolean,
  ): void {
    const mdCol = statsDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    //decompose into two different renering funcs
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
        const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
          tableColName : gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
        const currentAAR: string = gridTable.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, tableRowIndex);

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

          //TODO: frame based on currentSelection
          // if (currentAAR == 'D' && currentPosition == '11') {
          //   canvasContext.strokeStyle = '#000';
          //   canvasContext.lineWidth = 1;
          //   canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
          // }
          const aarSelection = this.currentSelection[currentPosition];
          if (aarSelection && aarSelection.includes(currentAAR)) {
            canvasContext.strokeStyle = '#000';
            canvasContext.lineWidth = 1;
            canvasContext.strokeRect(bound.x + 1, bound.y + 1, bound.width - 1, bound.height - 1);
          }
        }
        args.preventDefault();
      }
      canvasContext.restore();
    };
    sarGrid.onCellRender.subscribe(cellRendererAction);
    sarVGrid.onCellRender.subscribe(cellRendererAction);
  }

  //FIXME: doesn't work at all
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

    const chooseAction = (aar: string, position: string, isAltPressed: boolean) => {
      isAltPressed ? this.modifyCurrentSelection(aar, position) : this.initCurrentSelection(aar, position);
    };

    const gridCellValidation = (gc: DG.GridCell | null) => !gc || !gc.tableColumn || gc.tableRowIndex == null ||
      gc.tableRowIndex == -1;
    this._sarGrid.root.addEventListener('click', (ev) => {
      const gridCell = this._sarGrid.hitTest(ev.offsetX, ev.offsetY);
      if (gridCellValidation(gridCell) || gridCell!.tableColumn!.name == C.COLUMNS_NAMES.AMINO_ACID_RESIDUE)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = sarDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, gridCell!.tableRowIndex!);
      chooseAction(aar, position, ev.altKey);
    });

    this._sarVGrid.root.addEventListener('click', (ev) => {
      const gridCell = this._sarVGrid.hitTest(ev.offsetX, ev.offsetY);
      if (gridCellValidation(gridCell) || gridCell!.tableColumn!.name != C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const tableRowIdx = gridCell!.tableRowIndex!;
      const position = sarVDf.get(C.COLUMNS_NAMES.POSITION, tableRowIdx);
      const aar = sarVDf.get(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE, tableRowIdx);
      chooseAction(aar, position, ev.altKey);
    });

    const cellChanged = (table: DG.DataFrame) => {
      if (this.isCellChanging)
        return;
      this.isCellChanging = true;
      table.currentRowIdx = -1;
      this.isCellChanging = false;
    };
    this._sarGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(sarDf));
    this._sarVGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(sarVDf));
  }

  modifyCurrentSelection(aar: string, position: string): void {
    const tempSelection = this.currentSelection;
    if (!tempSelection.hasOwnProperty(position))
      tempSelection[position] = [aar];
    else {
      const tempSelectionAt = tempSelection[position];
      const aarIndex = tempSelectionAt.indexOf(aar);
      aarIndex == -1 ? tempSelectionAt.push(aar) :
        tempSelectionAt.length == 1 ? delete tempSelection[position] :
          tempSelectionAt.splice(aarIndex, 1);
    }

    this.currentSelection = tempSelection;
  }

  initCurrentSelection(aar: string, position: string): void {
    const tempSelection: type.SelectionObject = {};
    tempSelection[position] = [aar];
    this.currentSelection = tempSelection;
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
    const filter = this._dataFrame.filter;
    const selection = this._dataFrame.selection;

    const changeBitset = (currentBitset: DG.BitSet, previousBitset: DG.BitSet): void => {
      previousBitset.setAll(!this._filterMode, false);
      // currentBitset.init((i) => {
      //   const currentCategory = this.splitCol.get(i);
      //   return currentCategory !== C.CATEGORIES.OTHER && currentCategory !== C.CATEGORIES.ALL;
      // }, false);
      currentBitset.init((i) => this.splitCol.get(i), false);
    };

    const recalculateStatistics =
      (bitset: DG.BitSet): void => (this._dataFrame.temp[C.STATS] as FilteringStatistics).setMask(bitset);

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

  getSplitColValueAt(index: number, aar: string, position: string, aarLabel: string): string {
    const currentAAR = this._dataFrame.get(position, index) as string;
    return currentAAR === aar ? aarLabel : C.CATEGORIES.OTHER;
  }

  // modifyOrCreateSplitCol(aar: string, position: string): void {
  modifyOrCreateSplitCol(): void {
    const df = this._dataFrame;
    this.splitCol = df.col(C.COLUMNS_NAMES.SPLIT_COL) ?? df.columns.addNewBool(C.COLUMNS_NAMES.SPLIT_COL);

    const positionList = Object.keys(this.currentSelection);
    if (positionList.length == 0) {
      this.splitCol.init(() => true);
      return;
    }

    //TODO: move out
    const splitColInit = (i: number) => {
      for (const position of positionList) {
        const positionCol: DG.Column<string> = this._dataFrame.getCol(position);
        if (this._currentSelection[position].includes(positionCol.get(i)))
          return true;
      }
      return false;
    };
    this.splitCol.init((i) => splitColInit(i));

    this.splitCol.compact();
  }

  static async scaleActivity(
    activityScaling: string, df: DG.DataFrame, originalActivityName?: string, cloneBitset = false,
  ): Promise<[DG.DataFrame, string]> {
    let currentActivityColName = originalActivityName ?? C.COLUMNS_NAMES.ACTIVITY;
    const flag = df.columns.names().includes(currentActivityColName) &&
      currentActivityColName === originalActivityName;
    currentActivityColName = flag ? currentActivityColName : C.COLUMNS_NAMES.ACTIVITY;
    const tempDf = df.clone(cloneBitset ? df.filter : null, [currentActivityColName]);

    let formula = '${' + currentActivityColName + '}';
    let newColName = 'activity'; //originalActivityName ?? df.temp[C.COLUMNS_NAMES.ACTIVITY] ?? currentActivityColName;
    switch (activityScaling) {
    case 'none':
      break;
    case 'lg':
      formula = `Log10(${formula})`;
      newColName = `Log10(${newColName})`;
      break;
    case '-lg':
      formula = `-1*Log10(${formula})`;
      newColName = `-Log10(${newColName})`;
      break;
    default:
      throw new Error(`ScalingError: method \`${activityScaling}\` is not available.`);
    }

    await tempDf.columns.addNewCalculated(C.COLUMNS_NAMES.ACTIVITY_SCALED, formula);
    df.tags['scaling'] = activityScaling;

    return [tempDf, newColName];
  }

  static splitAlignedPeptides(peptideColumn: DG.Column<string>, filter: boolean = true): [DG.DataFrame, number[]] {
    const separator = peptideColumn.tags[C.TAGS.SEPARATOR] ?? getSeparator(peptideColumn);
    const splitPeptidesArray: string[][] = [];
    let currentSplitPeptide: string[];
    let modeMonomerCount = 0;
    let currentLength;
    const colLength = peptideColumn.length;

    // splitting data
    const monomerLengths: {[index: string]: number} = {};
    for (let i = 0; i < colLength; i++) {
      currentSplitPeptide = peptideColumn.get(i).split(separator).map((value: string) => value ? value : '-');
      splitPeptidesArray.push(currentSplitPeptide);
      currentLength = currentSplitPeptide.length;
      monomerLengths[currentLength + ''] =
        monomerLengths[currentLength + ''] ? monomerLengths[currentLength + ''] + 1 : 1;
    }
    modeMonomerCount =
      parseInt(Object.keys(monomerLengths).reduce((a, b) => monomerLengths[a] > monomerLengths[b] ? a : b));

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
    const columnNames = Array.from({length: modeMonomerCount}, (_, index) => `${index + 1 < 10 ? 0 : ''}${index + 1 }`);
    columnNames.splice(0, 0, 'N');
    columnNames.push('C');

    // filter out the columns with the same values
    if (filter) {
      splitColumns = splitColumns.filter((positionArray, index) => {
        const isRetained = new Set(positionArray).size > 1;
        if (!isRetained)
          columnNames.splice(index, 1);

        return isRetained;
      });
    }

    return [
      DG.DataFrame.fromColumns(splitColumns.map((positionArray, index) => {
        return DG.Column.fromList('string', columnNames[index], positionArray);
      })),
      invalidIndexes,
    ];
  }

  syncProperties(isSourceSAR = true): void {
    const sarViewer = this._dataFrame.temp['sarViewer'];
    const sarViewerVertical = this._dataFrame.temp['sarViewerVertical'];
    const sourceViewer = isSourceSAR ? sarViewer : sarViewerVertical;
    const targetViewer = isSourceSAR ? sarViewerVertical : sarViewer;
    const properties = sourceViewer.props.getProperties();
    for (const property of properties)
      targetViewer.props.set(property.name, property.get(sourceViewer));
  }

  //TODO: move to viewer
  setSARGridCellAt(aar: string, position: string): void {
    const sarDf = this._sarGrid.dataFrame;
    const aarCol = sarDf.getCol(C.COLUMNS_NAMES.AMINO_ACID_RESIDUE);
    const aarColLen = aarCol.length;
    let index = -1;
    for (let i = 0; i < aarColLen; i++) {
      if (aarCol.get(i) === aar) {
        index = i;
        break;
      }
    }
    position = position === C.CATEGORIES.ALL ? C.COLUMNS_NAMES.AMINO_ACID_RESIDUE : position;
    sarDf.currentCell = sarDf.cell(index, position);
  }

  static get chemPalette(): typeof ChemPalette {return ChemPalette;}

  /** Class initializer */
  async init(): Promise<void> {
    if (this.isInitialized)
      return;
    this.isInitialized = true;
    //calculate initial stats
    const stats = new FilteringStatistics();
    const activityScaledCol = this._dataFrame.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
    stats.setData(activityScaledCol.getRawData() as Float32Array);
    stats.setMask(this._dataFrame.selection);
    this._dataFrame.temp[C.STATS] = stats;

    this.currentView = grok.shell.addTableView(this._dataFrame);
    const sourceGrid = this.currentView.grid;
    sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED)!.name = this._dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED];
    sourceGrid.columns.setOrder([this._dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]]);

    this._dataFrame.temp[C.EMBEDDING_STATUS] = false;
    const adjustCellSize = (grid: DG.Grid): void => {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    };

    for (let i = 0; i < sourceGrid.columns.length; i++) {
      const aarCol = sourceGrid.columns.byIndex(i);
      if (aarCol && aarCol.name && aarCol.column?.semType !== C.SEM_TYPES.AMINO_ACIDS &&
        aarCol.name !== this._dataFrame.temp[C.COLUMNS_NAMES.ACTIVITY_SCALED]
      )
        sourceGrid.columns.byIndex(i)!.visible = false;
    }

    const options = {scaling: this._dataFrame.tags['scaling']};
    await this.updateData(this._dataFrame.tags['scaling'], sourceGrid, false, 1, 2, false, false);

    const dockManager = this.currentView.dockManager;

    const sarViewer =
      await this._dataFrame.plot.fromType('peptide-sar-viewer', options) as SARViewer;
    this._dataFrame.temp['sarViewer'] = sarViewer;

    const sarViewerVertical =
      await this._dataFrame.plot.fromType('peptide-sar-viewer-vertical', options) as SARViewerVertical;
    this._dataFrame.temp['sarViewerVertical'] = sarViewer;

    const sarViewersGroup: viewerTypes[] = [sarViewer, sarViewerVertical];

    const peptideSpaceViewerOptions = {method: 't-SNE', measure: 'Levenshtein', cyclesCount: 100};
    const peptideSpaceViewer =
      await this._dataFrame.plot.fromType('peptide-space-viewer', peptideSpaceViewerOptions) as PeptideSpaceViewer;
    dockManager.dock(peptideSpaceViewer, DG.DOCK_TYPE.RIGHT, null, 'Peptide Space Viewer');

    dockViewers(sarViewersGroup, DG.DOCK_TYPE.RIGHT, dockManager, DG.DOCK_TYPE.DOWN);

    sourceGrid.props.allowEdit = false;
    adjustCellSize(sourceGrid);

    this.invalidateGrids();
  }

  invalidateSourceGrid(): void {this._sourceGrid.invalidate();}
}

type viewerTypes = SARViewer | SARViewerVertical;

function dockViewers(
  viewerList: viewerTypes[], attachDirection: DG.DockType, dockManager: DG.DockManager,
  initialAttachDirection?: DG.DockType): DG.DockNode[] | null {
  const viewerListLength = viewerList.length;
  if (viewerListLength === 0)
    return null;

  let currentViewer = viewerList[0];
  const nodeList = [dockManager.dock(currentViewer, initialAttachDirection, null, currentViewer.name ?? '')];
  const ratio = 1 / viewerListLength;

  for (let i = 1; i < viewerListLength; i++) {
    currentViewer = viewerList[i];
    nodeList.push(dockManager.dock(currentViewer, attachDirection, nodeList[i - 1], currentViewer.name ?? '', ratio));
  }
  return nodeList;
}
