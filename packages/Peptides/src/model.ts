import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';

import wu from 'wu';
import * as rxjs from 'rxjs';

import * as C from './utils/constants';
import * as type from './utils/types';
import {calculateBarsData, getTypedArrayConstructor, isGridCellInvalid, scaleActivity} from './utils/misc';
import {MutationCliffsViewer, SARViewerBase, MostPotentResiduesViewer} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getDistributionAndStats, getDistributionWidget} from './widgets/distribution';
import {getStats, Stats} from './utils/statistics';
import {LogoSummary} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';

export class PeptidesModel {
  static modelName = 'peptidesModel';

  mutationCliffsGridSubject = new rxjs.Subject<DG.Grid>();
  mostPotentResiduesGridSubject = new rxjs.Subject<DG.Grid>();
  logoSummaryGridSubject = new rxjs.Subject<DG.Grid>();

  _isUpdating: boolean = false;
  isBitsetChangedInitialized = false;
  isCellChanging = false;

  mutationCliffsGrid!: DG.Grid;
  mostPotentResiduesGrid!: DG.Grid;
  logoSummaryGrid!: DG.Grid;
  sourceGrid!: DG.Grid;
  df: DG.DataFrame;
  splitCol!: DG.Column<boolean>;
  edf: DG.DataFrame | null = null;
  monomerPositionStatsDf!: DG.DataFrame;
  clusterStatsDf!: DG.DataFrame;
  _mutationCliffsSelection!: type.PositionToAARList;
  _invariantMapSelection!: type.PositionToAARList;
  _logoSummarySelection!: number[];
  substitutionsInfo: type.SubstitutionsInfo = new Map();
  isInitialized = false;
  currentView!: DG.TableView;

  isPeptideSpaceChangingBitset = false;
  isChangingEdfBitset = false;

  monomerMap: { [key: string]: { molfile: string, fullName: string } } = {};
  barData: type.MonomerDfStats = {};
  barsBounds: { [position: string]: type.BarCoordinates } = {};
  cachedBarchartTooltip: { bar: string, tooltip: null | HTMLDivElement } = {bar: '', tooltip: null};
  monomerLib: bio.IMonomerLib | null = null; // To get monomers from lib(s)
  monomerWorks: bio.MonomerWorks | null = null; // To get processed monomers

  _settings!: type.PeptidesSettings;
  isRibbonSet = false;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
  }

  static async getInstance(dataFrame: DG.DataFrame): Promise<PeptidesModel> {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    await (dataFrame.temp[PeptidesModel.modelName] as PeptidesModel).init();
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  get onMutationCliffsGridChanged(): rxjs.Observable<DG.Grid> {
    return this.mutationCliffsGridSubject.asObservable();
  }

  get onMostPotentResiduesGridChanged(): rxjs.Observable<DG.Grid> {
    return this.mostPotentResiduesGridSubject.asObservable();
  }

  get onLogoSummaryGridChanged(): rxjs.Observable<DG.Grid> {
    return this.logoSummaryGridSubject.asObservable();
  }

  get mutationCliffsSelection(): type.PositionToAARList {
    this._mutationCliffsSelection ??= JSON.parse(this.df.tags[C.TAGS.SELECTION] || '{}');
    return this._mutationCliffsSelection;
  }

  set mutationCliffsSelection(selection: type.PositionToAARList) {
    this._mutationCliffsSelection = selection;
    this.df.tags[C.TAGS.SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.invalidateGrids();
  }

  get invariantMapSelection(): type.PositionToAARList {
    this._invariantMapSelection ??= JSON.parse(this.df.tags[C.TAGS.FILTER] || '{}');
    return this._invariantMapSelection;
  }

  set invariantMapSelection(selection: type.PositionToAARList) {
    this._invariantMapSelection = selection;
    this.df.tags[C.TAGS.FILTER] = JSON.stringify(selection);
    this.df.filter.fireChanged();
    this.invalidateGrids();
  }

  get logoSummarySelection(): number[] {
    this._logoSummarySelection ??= JSON.parse(this.df.tags[C.TAGS.CLUSTER_SELECTION] || '[]');
    return this._logoSummarySelection;
  }

  set logoSummarySelection(selection: number[]) {
    this._logoSummarySelection = selection;
    this.df.tags[C.TAGS.CLUSTER_SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.invalidateGrids();
  }

  get splitByPos(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] ?? '00')[0];
    return splitByPosFlag == '1' ? true : false;
  }

  set splitByPos(flag: boolean) {
    const splitByAARFlag = (this.df.tags['distributionSplit'] ?? '00')[1];
    this.df.tags['distributionSplit'] = `${flag ? 1 : 0}${splitByAARFlag}`;
  }

  get splitByAAR(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] ?? '00')[1];
    return splitByPosFlag == '1' ? true : false;
  }

  set splitByAAR(flag: boolean) {
    const splitByAARFlag = (this.df.tags['distributionSplit'] ?? '00')[0];
    this.df.tags['distributionSplit'] = `${splitByAARFlag}${flag ? 1 : 0}`;
  }

  get isInvariantMap(): boolean {
    return this.df.getTag('isInvariantMap') === '1';
  }

  set isInvariantMap(x: boolean) {
    this.df.setTag('isInvariantMap', x ? '1' : '0');
  }

  get isMutationCliffSelectionEmpty(): boolean {
    for (const aarList of Object.values(this.mutationCliffsSelection)) {
      if (aarList.length !== 0)
        return false;
    }
    return true;
  }

  get isLogoSummarySelectionEmpty(): boolean {
    return this.logoSummarySelection.length === 0;
  }

  get settings(): type.PeptidesSettings {
    this._settings ??= JSON.parse(this.df.getTag('settings') ?? '{}');
    return this._settings;
  }
  set settings(s: type.PeptidesSettings) {
    for (const [key, value] of Object.entries(s))
      this._settings[key as keyof type.PeptidesSettings] = value as any;
    this.df.setTag('settings', JSON.stringify(this._settings));
    //TODO: update only needed components
    this.updateDefault();
  }

  createAccordion(): DG.Accordion {
    const acc = ui.accordion();
    acc.root.style.width = '100%';
    acc.addTitle(ui.h1(`${this.df.selection.trueCount} selected rows`));
    acc.addPane('Mutation Cliff pairs', () => mutationCliffsWidget(this.df, this).root, true);
    acc.addPane('Distribution', () => getDistributionWidget(this.df, this).root, true);

    return acc;
  }

  updateDefault(): void {
    if ((this.sourceGrid && !this._isUpdating) || !this.isInitialized) {
      this.isInitialized = true;
      this._isUpdating = true;
      this.initializeViewersComponents();
      //FIXME: modify during the initializeViewersComponents stages
      this.mutationCliffsGridSubject.next(this.mutationCliffsGrid);
      this.mostPotentResiduesGridSubject.next(this.mostPotentResiduesGrid);
      if (this.df.getTag(C.TAGS.CLUSTERS))
        this.logoSummaryGridSubject.next(this.logoSummaryGrid);

      this.fireBitsetChanged();
      this.invalidateGrids();
      this._isUpdating = false;
    }
  }

  initializeViewersComponents(): void {
    if (this.sourceGrid === null)
      throw new Error(`Source grid is not initialized`);

    //Split the aligned sequence into separate AARs
    const col = this.df.getCol(C.COLUMNS_NAMES.MACROMOLECULE);
    const alphabet = col.tags['alphabet'];
    const splitSeqDf = splitAlignedSequences(col);

    this.barData = calculateBarsData(splitSeqDf.columns.toList(), this.df.selection);

    const positionColumns = splitSeqDf.columns.names();

    const activityCol = this.df.columns.bySemType(C.SEM_TYPES.ACTIVITY)!;
    splitSeqDf.columns.add(activityCol);

    this.joinDataFrames(positionColumns, splitSeqDf, alphabet);

    this.sortSourceGrid();

    this.createScaledCol(splitSeqDf);

    //unpivot a table and handle duplicates
    let matrixDf = splitSeqDf.groupBy(positionColumns).aggregate();

    matrixDf = matrixDf.unpivot([], positionColumns, C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER);

    //statistics for specific AAR at a specific position
    this.monomerPositionStatsDf = this.calculateMonomerPositionStatistics(matrixDf);

    // SAR matrix table
    //pivot a table to make it matrix-like
    matrixDf = this.monomerPositionStatsDf.groupBy([C.COLUMNS_NAMES.MONOMER])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.MEAN_DIFFERENCE, '')
      .aggregate();
    matrixDf.name = 'SAR';

    // Setting category order
    this.setCategoryOrder(matrixDf);

    // SAR vertical table (naive, choose best Mean difference from pVals <= 0.01)
    const sequenceDf = this.createVerticalTable();

    this.calcSubstitutions();

    [this.mutationCliffsGrid, this.mostPotentResiduesGrid] =
      this.createGrids(matrixDf, sequenceDf, positionColumns, alphabet);

    if (this.df.getTag(C.TAGS.CLUSTERS)) {
      this.clusterStatsDf = this.calculateClusterStatistics();
      this.logoSummaryGrid = this.createLogoSummaryGrid();
    }

    // init invariant map & mutation cliffs selections
    this.initSelections(positionColumns);

    positionColumns.push(C.COLUMNS_NAMES.MEAN_DIFFERENCE);

    this.setBarChartInteraction();

    this.setCellRenderers(positionColumns);

    // show all the statistics in a tooltip over cell
    this.setTooltips(positionColumns);

    this.setInteractionCallback();

    this.setBitsetCallback();

    this.postProcessGrids();
  }

  //TODO: move out
  calcSubstitutions(): void {
    const activityValues: DG.Column<number> = this.df.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
    const columnList: DG.Column<string>[] = this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER);
    const nCols = columnList.length;
    if (nCols == 0)
      throw new Error(`Couldn't find any column of semType '${C.SEM_TYPES.MONOMER}'`);

    this.substitutionsInfo = new Map();
    const nRows = this.df.rowCount;
    for (let seq1Idx = 0; seq1Idx < nRows - 1; seq1Idx++) {
      for (let seq2Idx = seq1Idx + 1; seq2Idx < nRows; seq2Idx++) {
        let substCounter = 0;
        const activityValSeq1 = activityValues.get(seq1Idx)!;
        const activityValSeq2 = activityValues.get(seq2Idx)!;
        const delta = activityValSeq1 - activityValSeq2;
        if (Math.abs(delta) < (this.settings.minActivityDelta ?? 0))
          continue;

        let substCounterFlag = false;
        const tempData: { pos: string, seq1monomer: string, seq2monomer: string, seq1Idx: number, seq2Idx: number }[] =
          [];
        for (const currentPosCol of columnList) {
          const seq1monomer = currentPosCol.get(seq1Idx)!;
          const seq2monomer = currentPosCol.get(seq2Idx)!;
          if (seq1monomer == seq2monomer)
            continue;

          substCounter++;
          substCounterFlag = substCounter > (this.settings.maxMutations ?? 1);
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

  initSelections(positionColumns: string[]): void {
    const tempInvariantMapSelection: type.PositionToAARList = this.invariantMapSelection;
    const mutationCliffsSelection: type.PositionToAARList = this.mutationCliffsSelection;
    for (const pos of positionColumns) {
      tempInvariantMapSelection[pos] ??= [];
      mutationCliffsSelection[pos] ??= [];
    }
    this.invariantMapSelection = tempInvariantMapSelection;
    this.mutationCliffsSelection = mutationCliffsSelection;
    this.barData = calculateBarsData(this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER), this.df.selection);
  }

  joinDataFrames(positionColumns: string[], splitSeqDf: DG.DataFrame, alphabet: string): void {
    // append splitSeqDf columns to source table and make sure columns are not added more than once
    const name = this.df.name;
    const cols = this.df.columns;
    for (const colName of positionColumns) {
      const col = this.df.col(colName);
      const newCol = splitSeqDf.getCol(colName);
      if (col === null)
        cols.add(newCol);
      else {
        cols.remove(colName);
        cols.add(newCol);
      }
      CR.setAARRenderer(newCol, alphabet, this.sourceGrid);
    }
    this.df.name = name;
    this.currentView.name = name;
  }

  sortSourceGrid(): void {
    const colNames: DG.GridColumn[] = [];
    for (let i = 1; i < this.sourceGrid.columns.length; i++)
      colNames.push(this.sourceGrid.columns.byIndex(i)!);

    colNames.sort((a, b) => {
      if (a.column!.semType == C.SEM_TYPES.MONOMER) {
        if (b.column!.semType == C.SEM_TYPES.MONOMER)
          return 0;
        return -1;
      }
      if (b.column!.semType == C.SEM_TYPES.MONOMER)
        return 1;
      return 0;
    });
    this.sourceGrid.columns.setOrder(colNames.map((v) => v.name));
  }

  createScaledCol(splitSeqDf: DG.DataFrame): void {
    const scaledCol = scaleActivity(this.df.getCol(C.COLUMNS_NAMES.ACTIVITY), this.settings.scaling);
    //TODO: make another func
    splitSeqDf.columns.add(scaledCol);
    this.df.columns.replace(C.COLUMNS_NAMES.ACTIVITY_SCALED, scaledCol);
    const gridCol = this.sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    if (gridCol)
      gridCol.name = scaledCol.getTag('gridName');

    this.sourceGrid.columns.setOrder([scaledCol.getTag('gridName')]);
  }

  calculateMonomerPositionStatistics(matrixDf: DG.DataFrame): DG.DataFrame {
    matrixDf = matrixDf.groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER]).aggregate();

    //calculate p-values based on t-test
    const matrixCols = matrixDf.columns;
    const mdCol = matrixCols.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const pValCol = matrixCols.addNewFloat(C.COLUMNS_NAMES.P_VALUE);
    const countCol = matrixCols.addNewInt(C.COLUMNS_NAMES.COUNT);
    const ratioCol = matrixCols.addNewFloat(C.COLUMNS_NAMES.RATIO);
    const aarCol = matrixDf.getCol(C.COLUMNS_NAMES.MONOMER);
    const posCol = matrixDf.getCol(C.COLUMNS_NAMES.POSITION);
    const activityCol: number[] = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).toList();
    const sourceDfLen = activityCol.length;

    for (let i = 0; i < matrixDf.rowCount; i++) {
      const position: string = posCol.get(i);
      const aar: string = aarCol.get(i);
      const mask = DG.BitSet.create(sourceDfLen, (j) => this.df.get(position, j) === aar);
      const stats = getStats(activityCol, mask);

      mdCol.set(i, stats.meanDifference);
      pValCol.set(i, stats.pValue);
      countCol.set(i, stats.count);
      ratioCol.set(i, stats.ratio);
    }

    const monomerCol = matrixDf.getCol(C.COLUMNS_NAMES.MONOMER);
    matrixDf = matrixDf.clone(DG.BitSet.create(matrixDf.rowCount, (i) => monomerCol.get(i) ? true : false));
    return matrixDf as DG.DataFrame;
  }

  calculateClusterStatistics(): DG.DataFrame {
    const originalClustersCol = this.df.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const statsDf = this.df.groupBy([C.COLUMNS_NAMES.CLUSTERS]).aggregate();
    const clustersCol = statsDf.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const statsDfCols = statsDf.columns;
    const mdCol = statsDfCols.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const pValCol = statsDfCols.addNewFloat(C.COLUMNS_NAMES.P_VALUE);
    const countCol = statsDfCols.addNewInt(C.COLUMNS_NAMES.COUNT);
    const ratioCol = statsDfCols.addNewFloat(C.COLUMNS_NAMES.RATIO);
    const activityList: number[] = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).toList();

    for (let rowIdx = 0; rowIdx < clustersCol.length; ++rowIdx) {
      const cluster = clustersCol.get(rowIdx);
      const mask = DG.BitSet.create(activityList.length, (bitIdx) => originalClustersCol.get(bitIdx) === cluster);
      const stats = getStats(activityList, mask);

      mdCol.set(rowIdx, stats.meanDifference);
      pValCol.set(rowIdx, stats.pValue);
      countCol.set(rowIdx, stats.count);
      ratioCol.set(rowIdx, stats.ratio);
    }
    return statsDf;
  }

  setCategoryOrder(matrixDf: DG.DataFrame): void {
    let sortArgument: string = C.COLUMNS_NAMES.MEAN_DIFFERENCE;
    if (this.settings.isBidirectional) {
      const mdCol = this.monomerPositionStatsDf.getCol(sortArgument);
      sortArgument = 'Absolute Mean difference';
      const absMDCol = this.monomerPositionStatsDf.columns.addNewFloat(sortArgument);
      absMDCol.init((i) => Math.abs(mdCol.get(i)));
    }

    const aarWeightsDf = this.monomerPositionStatsDf.groupBy([C.COLUMNS_NAMES.MONOMER]).sum(sortArgument, 'weight')
      .aggregate();
    const aarList = aarWeightsDf.getCol(C.COLUMNS_NAMES.MONOMER).toList();
    const getWeight = (aar: string): number => aarWeightsDf
      .groupBy(['weight'])
      .where(`${C.COLUMNS_NAMES.MONOMER} = ${aar}`)
      .aggregate()
      .get('weight', 0) as number;
    aarList.sort((first, second) => getWeight(second) - getWeight(first));

    matrixDf.getCol(C.COLUMNS_NAMES.MONOMER).setCategoryOrder(aarList);
  }

  createVerticalTable(): DG.DataFrame {
    // TODO: aquire ALL of the positions
    const columns = [C.COLUMNS_NAMES.MEAN_DIFFERENCE, C.COLUMNS_NAMES.MONOMER, C.COLUMNS_NAMES.POSITION,
      'Count', 'Ratio', C.COLUMNS_NAMES.P_VALUE];
    let sequenceDf = this.monomerPositionStatsDf.groupBy(columns)
      .where('pValue <= 0.1')
      .aggregate();

    let tempStats: DG.Stats;
    const maxAtPos: { [index: string]: number } = {};
    const posColCategories = sequenceDf.getCol(C.COLUMNS_NAMES.POSITION).categories;
    const mdCol = sequenceDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const posCol = sequenceDf.getCol(C.COLUMNS_NAMES.POSITION);
    const rowCount = sequenceDf.rowCount;
    for (const pos of posColCategories) {
      tempStats = DG.Stats.fromColumn(mdCol, DG.BitSet.create(rowCount, (i) => posCol.get(i) === pos));
      maxAtPos[pos] = this.settings.isBidirectional ?
        (tempStats.max > Math.abs(tempStats.min) ? tempStats.max : tempStats.min) :
        tempStats.max;
    }
    sequenceDf = sequenceDf.clone(DG.BitSet.create(rowCount, (i) => mdCol.get(i) === maxAtPos[posCol.get(i)]));

    return sequenceDf;
  }

  createGrids(mutationCliffsDf: DG.DataFrame, mostPotentResiduesDf: DG.DataFrame, positionColumns: string[],
    alphabet: string): [DG.Grid, DG.Grid] {
    // Creating Mutation Cliffs grid and sorting columns
    const mutationCliffsGrid = mutationCliffsDf.plot.grid();
    mutationCliffsGrid.sort([C.COLUMNS_NAMES.MONOMER]);
    mutationCliffsGrid.columns.setOrder([C.COLUMNS_NAMES.MONOMER].concat(positionColumns as C.COLUMNS_NAMES[]));

    // Creating Monomer-Position grid, sorting and setting column format
    const mostPotentResiduesGrid = mostPotentResiduesDf.plot.grid();
    mostPotentResiduesGrid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = mostPotentResiduesGrid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = '#.000';
    pValGridCol.name = 'P-value';

    // Setting Monomer column renderer
    CR.setAARRenderer(mutationCliffsDf.getCol(C.COLUMNS_NAMES.MONOMER), alphabet, mutationCliffsGrid);
    CR.setAARRenderer(mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.MONOMER), alphabet, mostPotentResiduesGrid);

    return [mutationCliffsGrid, mostPotentResiduesGrid];
  }

  createLogoSummaryGrid(): DG.Grid {
    const summaryTable = this.df.groupBy([C.COLUMNS_NAMES.CLUSTERS]).aggregate();
    const summaryTableLength = summaryTable.rowCount;
    const clustersCol: DG.Column<number> = summaryTable.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const membersCol: DG.Column<number> = summaryTable.columns.addNewInt('Members');
    const webLogoCol: DG.Column<string> = summaryTable.columns.addNew('WebLogo', DG.COLUMN_TYPE.STRING);
    const tempDfList: DG.DataFrame[] = new Array(summaryTableLength);
    const originalClustersCol = this.df.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const peptideCol: DG.Column<string> = this.df.getCol(C.COLUMNS_NAMES.MACROMOLECULE);

    for (let index = 0; index < summaryTableLength; ++index) {
      const indexes: number[] = [];
      for (let j = 0; j < originalClustersCol.length; ++j) {
        if (originalClustersCol.get(j) === clustersCol.get(index))
          indexes.push(j);
      }
      const tCol = DG.Column.string('peptides', indexes.length);
      tCol.init((i) => peptideCol.get(indexes[i]));

      for (const tag of peptideCol.tags)
        tCol.setTag(tag[0], tag[1]);

      const dfSlice = DG.DataFrame.fromColumns([tCol]);
      tempDfList[index] = dfSlice;
      webLogoCol.set(index, index.toString());
      membersCol.set(index, dfSlice.rowCount);
    }
    webLogoCol.setTag(DG.TAGS.CELL_RENDERER, 'html');

    const grid = summaryTable.plot.grid();
    const gridClustersCol = grid.col(C.COLUMNS_NAMES.CLUSTERS)!;
    gridClustersCol.name = 'Clusters';
    gridClustersCol.visible = true;
    grid.columns.rowHeader!.visible = false;
    grid.props.rowHeight = 55;
    grid.onCellPrepare((cell) => {
      if (cell.isTableCell && cell.tableColumn?.name === 'WebLogo') {
        tempDfList[parseInt(cell.cell.value)].plot.fromType('WebLogo', {maxHeight: 50})
          .then((viewer) => cell.element = viewer.root);
      }
    });
    grid.root.addEventListener('click', (ev) => {
      const cell = grid.hitTest(ev.offsetX, ev.offsetY);
      if (!cell || !cell.isTableCell)
        return;

      const cluster = clustersCol.get(cell.tableRowIndex!)!;
      summaryTable.currentRowIdx = -1;
      if (ev.shiftKey)
        this.modifyClusterSelection(cluster);
      else
        this.initClusterSelection(cluster);
      this.barData = calculateBarsData(this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER), this.df.selection);
    });
    grid.onCellRender.subscribe((gridCellArgs) => {
      const gc = gridCellArgs.cell;
      if (gc.tableColumn?.name !== C.COLUMNS_NAMES.CLUSTERS || gc.isColHeader)
        return;
      const canvasContext = gridCellArgs.g;
      const bound = gridCellArgs.bounds;
      canvasContext.save();
      canvasContext.beginPath();
      canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
      canvasContext.clip();
      CR.renderLogoSummaryCell(canvasContext, gc.cell.value, this.logoSummarySelection, bound);
      gridCellArgs.preventDefault();
      canvasContext.restore();
    });
    grid.onCellTooltip((cell, x, y) => {
      if (!cell.isColHeader && cell.tableColumn?.name === C.COLUMNS_NAMES.CLUSTERS)
        this.showTooltipCluster(cell.cell.value, x, y);
      return true;
    });
    const webLogoGridCol = grid.columns.byName('WebLogo')!;
    webLogoGridCol.cellType = 'html';
    webLogoGridCol.width = 350;

    return grid;
  }

  modifyClusterSelection(cluster: number): void {
    const tempSelection = this.logoSummarySelection;
    const idx = tempSelection.indexOf(cluster);
    if (idx !== -1)
      tempSelection.splice(idx, 1);
    else
      tempSelection.push(cluster);

    this.logoSummarySelection = tempSelection;
  }

  initClusterSelection(cluster: number): void {
    this.logoSummarySelection = [cluster];
  }

  setBarChartInteraction(): void {
    const eventAction = (ev: MouseEvent): void => {
      const cell = this.sourceGrid.hitTest(ev.offsetX, ev.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType == C.SEM_TYPES.MONOMER) {
        const newBarPart = this.findAARandPosition(cell, ev);
        this.requestBarchartAction(ev, newBarPart);
      }
    };

    // The following events makes the barchart interactive
    rxjs.fromEvent<MouseEvent>(this.sourceGrid.overlay, 'mousemove')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
    rxjs.fromEvent<MouseEvent>(this.sourceGrid.overlay, 'click')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  }

  findAARandPosition(cell: DG.GridCell, ev: MouseEvent): { monomer: string, position: string } | null {
    const barCoords = this.barsBounds[cell.tableColumn!.name];
    for (const [monomer, coords] of Object.entries(barCoords)) {
      const isIntersectingX = ev.offsetX >= coords.x && ev.offsetX <= coords.x + coords.width;
      const isIntersectingY = ev.offsetY >= coords.y && ev.offsetY <= coords.y + coords.height;
      if (isIntersectingX && isIntersectingY)
        return {monomer: monomer, position: cell.tableColumn!.name};
    }

    return null;
  }

  requestBarchartAction(ev: MouseEvent, barPart: { position: string, monomer: string } | null): void {
    if (!barPart)
      return;
    const monomer = barPart.monomer;
    const position = barPart.position;
    if (ev.type === 'click') {
      ev.shiftKey ? this.modifyMonomerPositionSelection(monomer, position, true) :
        this.initMonomerPositionSelection(monomer, position, true);
      this.barData = calculateBarsData(this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER), this.df.selection);
    } else {
      const bar = `${monomer}:${position}`;
      if (this.cachedBarchartTooltip.bar == bar)
        ui.tooltip.show(this.cachedBarchartTooltip.tooltip!, ev.clientX, ev.clientY);
      else
        this.cachedBarchartTooltip = {bar: bar, tooltip: this.showTooltipAt(monomer, position, ev.clientX, ev.clientY)};
    }
  }

  setCellRenderers(renderColNames: string[]): void {
    const mdCol = this.monomerPositionStatsDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    //decompose into two different renering funcs
    const renderCell = (args: DG.GridCellRenderArgs): void => {
      const canvasContext = args.g;
      const bound = args.bounds;

      canvasContext.save();
      canvasContext.beginPath();
      canvasContext.rect(bound.x, bound.y, bound.width, bound.height);
      canvasContext.clip();

      // Hide row column
      const cell = args.cell;
      if (cell.isRowHeader && cell.gridColumn.visible) {
        cell.gridColumn.visible = false;
        args.preventDefault();
        return;
      }

      const tableColName = cell.tableColumn?.name;
      const tableRowIndex = cell.tableRowIndex!;
      if (cell.isTableCell && tableColName && tableRowIndex !== null && renderColNames.indexOf(tableColName) !== -1) {
        const cellValue: number | null = cell.cell.value;

        if (cellValue && cellValue !== DG.INT_NULL && cellValue !== DG.FLOAT_NULL) {
          const gridTable = cell.grid.table;
          const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
            tableColName : gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
          const currentAAR: string = gridTable.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);

          if (this.isInvariantMap) {
            const value: number = this.monomerPositionStatsDf
              .groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER, C.COLUMNS_NAMES.COUNT])
              .where(`${C.COLUMNS_NAMES.POSITION} = ${currentPosition} and ${C.COLUMNS_NAMES.MONOMER} = ${currentAAR}`)
              .aggregate().get(C.COLUMNS_NAMES.COUNT, 0);
            CR.renderInvaraintMapCell(
              canvasContext, currentAAR, currentPosition, this.invariantMapSelection, value, bound);
          } else {
            CR.renderMutationCliffCell(canvasContext, currentAAR, currentPosition, this.monomerPositionStatsDf,
              mdCol, bound, cellValue, this.mutationCliffsSelection, this.substitutionsInfo,
              this.settings.isBidirectional);
          }
        }
        args.preventDefault();
      }
      canvasContext.restore();
    };
    this.mutationCliffsGrid.onCellRender.subscribe(renderCell);
    this.mostPotentResiduesGrid.onCellRender.subscribe(renderCell);

    this.sourceGrid.setOptions({'colHeaderHeight': 130});
    this.sourceGrid.onCellRender.subscribe((gcArgs) => {
      const context = gcArgs.g;
      const bounds = gcArgs.bounds;
      const col = gcArgs.cell.tableColumn;

      context.save();
      context.beginPath();
      context.rect(bounds.x, bounds.y, bounds.width, bounds.height);
      context.clip();

      if (gcArgs.cell.isColHeader && col?.semType == C.SEM_TYPES.MONOMER) {
        const barBounds = CR.renderBarchart(context, col, this.barData[col.name], bounds, this.df.filter.trueCount);
        this.barsBounds[col.name] = barBounds;
        gcArgs.preventDefault();
      }

      context.restore();
    });
  }

  setTooltips(renderColNames: string[]): void {
    const showTooltip = (cell: DG.GridCell, x: number, y: number): boolean => {
      const tableCol = cell.tableColumn;
      const tableColName = tableCol?.name;
      const tableRowIndex = cell.tableRowIndex;

      if (!cell.isRowHeader && !cell.isColHeader && tableCol && tableRowIndex != null) {
        const table = cell.grid.table;
        const currentAAR = table.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);

        if (tableCol.semType == C.SEM_TYPES.MONOMER)
          this.showMonomerTooltip(currentAAR, x, y);
        else if (cell.cell.value && renderColNames.includes(tableColName!)) {
          const currentPosition = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
            table.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);

          this.showTooltipAt(currentAAR, currentPosition, x, y);
        }
      }
      return true;
    };

    this.mutationCliffsGrid.onCellTooltip(showTooltip);
    this.mostPotentResiduesGrid.onCellTooltip(showTooltip);
    this.sourceGrid.onCellTooltip((cell, x, y) => {
      const col = cell.tableColumn;
      const cellValue = cell.cell.value;
      if (cellValue && col && col.semType === C.SEM_TYPES.MONOMER)
        this.showMonomerTooltip(cellValue, x, y);
      return true;
    });
  }

  showMonomerTooltip(aar: string, x: number, y: number): void {
    const tooltipElements: HTMLDivElement[] = [];
    const monomerName = aar.toLowerCase();

    const monomer: bio.Monomer | null = wu(['HELM_AA', 'HELM_CHEM'])
      .map((monomerType) => this.monomerWorks!.getCappedMonomer(monomerType, monomerName))
      .find((m) => m != null) ?? null;

    if (monomer) {
      tooltipElements.push(ui.div(monomer.n));
      const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      tooltipElements.push(grok.chem.svgMol(monomer.m, undefined, undefined, options));
    } else
      tooltipElements.push(ui.div(aar));

    ui.tooltip.show(ui.divV(tooltipElements), x, y);
  }

  showTooltipAt(aar: string, position: string, x: number, y: number): HTMLDivElement | null {
    const currentStatsDf = this.monomerPositionStatsDf.rows.match({Pos: position, AAR: aar}).toDataFrame();
    const activityCol = this.df.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
    //TODO: use bitset instead of splitCol
    const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, activityCol.length);
    const currentPosCol = this.df.getCol(position);
    splitCol.init((i) => currentPosCol.get(i) == aar);
    const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
    const stats: Stats = {
      count: currentStatsDf.get(C.COLUMNS_NAMES.COUNT, 0),
      ratio: currentStatsDf.get(C.COLUMNS_NAMES.RATIO, 0),
      pValue: currentStatsDf.get(C.COLUMNS_NAMES.P_VALUE, 0),
      meanDifference: currentStatsDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, 0),
    };
    if (!stats.count)
      return null;

    const tooltip = getDistributionAndStats(distributionTable, stats, `${position} : ${aar}`, 'Other', true);

    ui.tooltip.show(tooltip, x, y);

    return tooltip;
  }

  showTooltipCluster(cluster: number, x: number, y: number): HTMLDivElement | null {
    const currentStatsDf = this.clusterStatsDf.rows.match({clusters: cluster}).toDataFrame();
    const activityCol = this.df.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
    //TODO: use bitset instead of splitCol
    const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, activityCol.length);
    const currentClusterCol = this.df.getCol(C.COLUMNS_NAMES.CLUSTERS);
    splitCol.init((i) => currentClusterCol.get(i) == cluster);
    const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
    const stats: Stats = {
      count: currentStatsDf.get(C.COLUMNS_NAMES.COUNT, 0),
      ratio: currentStatsDf.get(C.COLUMNS_NAMES.RATIO, 0),
      pValue: currentStatsDf.get(C.COLUMNS_NAMES.P_VALUE, 0),
      meanDifference: currentStatsDf.get(C.COLUMNS_NAMES.MEAN_DIFFERENCE, 0),
    };
    if (!stats.count)
      return null;

    const tooltip = getDistributionAndStats(distributionTable, stats, `Cluster: ${cluster}`, 'Other', true);

    ui.tooltip.show(tooltip, x, y);

    return tooltip;
  }

  setInteractionCallback(): void {
    const mutationCliffsDf = this.mutationCliffsGrid.dataFrame;
    const mostPotentResiduesDf = this.mostPotentResiduesGrid.dataFrame;

    const chooseAction =
      (aar: string, position: string, isShiftPressed: boolean, isInvariantMapSelection: boolean = true): void => {
        isShiftPressed ? this.modifyMonomerPositionSelection(aar, position, isInvariantMapSelection) :
          this.initMonomerPositionSelection(aar, position, isInvariantMapSelection);
        this.barData = calculateBarsData(this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER), this.df.selection);
      };

    this.mutationCliffsGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.mutationCliffsGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name == C.COLUMNS_NAMES.MONOMER)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = mutationCliffsDf.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!);
      chooseAction(aar, position, ev.shiftKey, this.isInvariantMap);
    });

    this.mostPotentResiduesGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.mostPotentResiduesGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name != C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const tableRowIdx = gridCell!.tableRowIndex!;
      const position = mostPotentResiduesDf.get(C.COLUMNS_NAMES.POSITION, tableRowIdx);
      const aar = mostPotentResiduesDf.get(C.COLUMNS_NAMES.MONOMER, tableRowIdx);
      chooseAction(aar, position, ev.shiftKey, false);
    });

    const cellChanged = (table: DG.DataFrame): void => {
      if (this.isCellChanging)
        return;
      this.isCellChanging = true;
      table.currentRowIdx = -1;
      this.isCellChanging = false;
    };
    this.mutationCliffsGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(mutationCliffsDf));
    this.mostPotentResiduesGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(mostPotentResiduesDf));
  }

  modifyMonomerPositionSelection(aar: string, position: string, isInvariantMapSelection: boolean): void {
    const tempSelection = isInvariantMapSelection ? this.invariantMapSelection : this.mutationCliffsSelection;
    const tempSelectionAt = tempSelection[position];
    const aarIndex = tempSelectionAt.indexOf(aar);
    if (aarIndex === -1)
      tempSelectionAt.push(aar);
    else
      tempSelectionAt.splice(aarIndex, 1);

    if (isInvariantMapSelection)
      this.invariantMapSelection = tempSelection;
    else
      this.mutationCliffsSelection = tempSelection;
  }

  initMonomerPositionSelection(aar: string, position: string, isInvariantMapSelection: boolean): void {
    const tempSelection = isInvariantMapSelection ? this.invariantMapSelection : this.mutationCliffsSelection;
    for (const key of Object.keys(tempSelection))
      tempSelection[key] = [];
    tempSelection[position] = [aar];

    if (isInvariantMapSelection)
      this.invariantMapSelection = tempSelection;
    else
      this.mutationCliffsSelection = tempSelection;
  }

  invalidateGrids(): void {
    this.mutationCliffsGrid.invalidate();
    this.mostPotentResiduesGrid.invalidate();
    this.logoSummaryGrid?.invalidate();
    this.sourceGrid?.invalidate();
  }

  setBitsetCallback(): void {
    if (this.isBitsetChangedInitialized)
      return;
    const selection = this.df.selection;
    const filter = this.df.filter;
    const clusterCol = this.df.col(C.COLUMNS_NAMES.CLUSTERS);

    const changeSelectionBitset = (currentBitset: DG.BitSet): void => {
      const edfSelection = this.edf?.selection;
      if (this.isPeptideSpaceChangingBitset) {
        if (edfSelection == null)
          return;

        currentBitset.init((i) => edfSelection.get(i) ?? false, false);
        return;
      }

      const updateEdfSelection = (): void => {
        this.isChangingEdfBitset = true;
        edfSelection?.copyFrom(currentBitset);
        this.isChangingEdfBitset = false;
      };

      const positionList = Object.keys(this.mutationCliffsSelection);

      //TODO: move out
      const getBitAt = (i: number): boolean => {
        for (const position of positionList) {
          const positionCol: DG.Column<string> = this.df.getCol(position);
          if (this.mutationCliffsSelection[position].includes(positionCol.get(i)!))
            return true;
        }
        if (this.logoSummarySelection.includes(clusterCol?.get(i)!))
          return true;
        return false;
      };
      currentBitset.init(getBitAt, false);

      updateEdfSelection();
    };

    selection.onChanged.subscribe(() => changeSelectionBitset(selection));
    filter.onChanged.subscribe(() => {
      const positionList = Object.keys(this.invariantMapSelection);
      for (let index = 0; index < this.df.rowCount; ++index) {
        let result = true;
        for (const position of positionList) {
          const aarList = this.invariantMapSelection[position];
          result &&= aarList.length === 0 || aarList.includes(this.df.get(position, index));
          if (!result)
            break;
        }
        filter.set(index, filter.get(index) && result, false);
      }
    });
    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(isPeptideSpaceSource: boolean = false): void {
    this.isPeptideSpaceChangingBitset = isPeptideSpaceSource;
    this.df.selection.fireChanged();
    this.modifyOrCreateSplitCol();
    grok.shell.o = this.createAccordion().root;
    this.isPeptideSpaceChangingBitset = false;
  }

  postProcessGrids(): void {
    const mdCol: DG.GridColumn = this.mostPotentResiduesGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';

    for (const grid of [this.mutationCliffsGrid, this.mostPotentResiduesGrid]) {
      const gridProps = grid.props;
      gridProps.rowHeight = 20;
      const girdCols = grid.columns;
      const colNum = girdCols.length;
      for (let i = 0; i < colNum; ++i) {
        const col = girdCols.byIndex(i)!;
        const colName = col.name;
        col.width =
          grid == this.mostPotentResiduesGrid && colName !== 'Diff' && colName !== C.COLUMNS_NAMES.MONOMER ? 50 :
            gridProps.rowHeight + 10;
      }
    }

    const setViewerGridProps = (grid: DG.Grid): void => {
      const gridProps = grid.props;
      gridProps.allowEdit = false;
      gridProps.allowRowSelection = false;
      gridProps.allowBlockSelection = false;
      gridProps.allowColSelection = false;
    };

    setViewerGridProps(this.mutationCliffsGrid);
    setViewerGridProps(this.mostPotentResiduesGrid);
    if (this.df.getTag(C.TAGS.CLUSTERS))
      setViewerGridProps(this.logoSummaryGrid);
  }

  getSplitColValueAt(index: number, aar: string, position: string, aarLabel: string): string {
    const currentAAR = this.df.get(position, index) as string;
    return currentAAR === aar ? aarLabel : C.CATEGORIES.OTHER;
  }

  modifyOrCreateSplitCol(): void {
    const bs = this.df.selection;
    this.splitCol = this.df.col(C.COLUMNS_NAMES.SPLIT_COL) ??
      this.df.columns.addNewBool(C.COLUMNS_NAMES.SPLIT_COL);
    this.splitCol.init((i) => bs.get(i));
    this.splitCol.compact();
  }

  /** Class initializer */
  async init(): Promise<void> {
    if (this.isInitialized)
      return;

    // Get monomer library through bio library
    this.monomerLib = await bio.getMonomerLib();
    this.monomerLib.onChanged.subscribe(() => {
      this.sourceGrid.invalidate();
    });
    this.monomerWorks = new bio.MonomerWorks(this.monomerLib);

    this.currentView = wu(grok.shell.tableViews).find(({dataFrame}) => dataFrame.tags[C.PEPTIDES_ANALYSIS] === 'true') ??
      grok.shell.addTableView(this.df);
    if (!this.isRibbonSet) {
      this.currentView.setRibbonPanels([[ui.icons.settings(() => getSettingsDialog(this))]], false);
      this.isRibbonSet = true;
    }
    grok.shell.v = this.currentView;
    this.sourceGrid = this.currentView.grid;
    if (this.df.tags[C.PEPTIDES_ANALYSIS] === 'true')
      return;

    this.df.tags[C.PEPTIDES_ANALYSIS] = 'true';
    const scaledGridCol = this.sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED)!;
    scaledGridCol.name = scaledGridCol.column!.getTag('gridName');
    scaledGridCol.format = '#.000';
    this.sourceGrid.columns.setOrder([scaledGridCol.name]);
    this.sourceGrid.props.allowColSelection = false;

    this.df.temp[C.EMBEDDING_STATUS] = false;
    const adjustCellSize = (grid: DG.Grid): void => {
      const colNum = grid.columns.length;
      for (let i = 0; i < colNum; ++i) {
        const iCol = grid.columns.byIndex(i)!;
        iCol.width = isNaN(parseInt(iCol.name)) ? 50 : 40;
      }
      grid.props.rowHeight = 20;
    };

    for (let i = 0; i < this.sourceGrid.columns.length; i++) {
      const currentCol = this.sourceGrid.columns.byIndex(i);
      if (currentCol?.column?.getTag(C.TAGS.VISIBLE) === '0')
        currentCol.visible = false;
    }

    const options = {scaling: this.df.tags['scaling']};

    const dockManager = this.currentView.dockManager;

    const mutationCliffsViewer = await this.df.plot.fromType('peptide-sar-viewer', options) as MutationCliffsViewer;

    const mostPotentResiduesViewer =
      await this.df.plot.fromType('peptide-sar-viewer-vertical', options) as MostPotentResiduesViewer;

    if (this.df.getTag(C.TAGS.CLUSTERS)) {
      const logoSummary = await this.df.plot.fromType('logo-summary-viewer') as LogoSummary;
      dockManager.dock(logoSummary, DG.DOCK_TYPE.RIGHT, null, 'Logo Summary Table');
    }

    this.updateDefault();

    const mcNode =
      dockManager.dock(mutationCliffsViewer, DG.DOCK_TYPE.DOWN, null, mutationCliffsViewer.name);

    dockManager.dock(mostPotentResiduesViewer, DG.DOCK_TYPE.RIGHT, mcNode, mostPotentResiduesViewer.name, 0.3);


    this.sourceGrid.props.allowEdit = false;
    adjustCellSize(this.sourceGrid);

    this.invalidateGrids();
  }
}
