import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';

import wu from 'wu';
import * as rxjs from 'rxjs';

import * as C from './utils/constants';
import * as type from './utils/types';
import {calculateSelected, extractMonomerInfo, scaleActivity} from './utils/misc';
import {MonomerPosition, MostPotentResiduesViewer} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getDistributionAndStats, getDistributionWidget} from './widgets/distribution';
import {getStats, Stats} from './utils/statistics';
import {LogoSummary} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {getMonomerWorks} from './package';
import * as bio from '@datagrok-libraries/bio';
import {findMutations} from './utils/algorithms';

export class PeptidesModel {
  static modelName = 'peptidesModel';

  settingsSubject: rxjs.Subject<type.PeptidesSettings> = new rxjs.Subject();
  _mutatinCliffsSelectionSubject: rxjs.Subject<undefined> = new rxjs.Subject();

  _isUpdating: boolean = false;
  isBitsetChangedInitialized = false;
  isCellChanging = false;

  df: DG.DataFrame;
  splitCol!: DG.Column<boolean>;
  edf: DG.DataFrame | null = null;
  _monomerPositionStatsDf?: DG.DataFrame;
  _clusterStatsDf?: DG.DataFrame;
  _mutationCliffsSelection!: type.PositionToAARList;
  _invariantMapSelection!: type.PositionToAARList;
  _logoSummarySelection!: number[];
  _substitutionsInfo?: type.SubstitutionsInfo;
  isInitialized = false;
  _analysisView?: DG.TableView;

  isPeptideSpaceChangingBitset = false;
  isChangingEdfBitset = false;

  monomerMap: { [key: string]: { molfile: string, fullName: string } } = {};
  monomerLib: bio.IMonomerLib | null = null; // To get monomers from lib(s)
  monomerWorks: bio.MonomerWorks | null = null; // To get processed monomers

  _settings!: type.PeptidesSettings;
  isRibbonSet = false;

  _cp?: bio.SeqPalette;
  initBitset: DG.BitSet;
  isInvariantMapTrigger: boolean = false;
  headerSelectedMonomers: type.MonomerSelectionStats = {};
  webLogoBounds: {[positon: string]: {[monomer: string]: DG.Rect}} = {};
  cachedWebLogoTooltip: {bar: string; tooltip: HTMLDivElement | null;} = {bar: '', tooltip: null};
  _monomerPositionDf?: DG.DataFrame;
  _alphabet?: string;
  _mostPotentResiduesDf?: DG.DataFrame;
  _matrixDf?: DG.DataFrame;
  _splitSeqDf?: DG.DataFrame;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
    this.initBitset = this.df.filter.clone();
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    (dataFrame.temp[PeptidesModel.modelName] as PeptidesModel).init();
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  get monomerPositionDf(): DG.DataFrame {
    this._monomerPositionDf ??= this.createMonomerPositionDf();
    return this._monomerPositionDf;
  }
  set monomerPositionDf(df: DG.DataFrame) {
    this._monomerPositionDf = df;
  }

  get monomerPositionStatsDf(): DG.DataFrame {
    this._monomerPositionStatsDf ??= this.calculateMonomerPositionStatistics();
    return this._monomerPositionStatsDf;
  }
  set monomerPositionStatsDf(df: DG.DataFrame) {
    this._monomerPositionStatsDf = df;
  }

  get matrixDf(): DG.DataFrame {
    this._matrixDf ??= this.buildMatrixDf();
    return this._matrixDf;
  }
  set matrixDf(df: DG.DataFrame) {
    this._matrixDf = df;
  }

  get splitSeqDf(): DG.DataFrame {
    this._splitSeqDf ??= this.buildSplitSeqDf();
    return this._splitSeqDf;
  }
  set splitSeqDf(df: DG.DataFrame) {
    this._splitSeqDf = df;
  }

  get mostPotentResiduesDf(): DG.DataFrame {
    this._mostPotentResiduesDf ??= this.createVerticalTable();
    return this._mostPotentResiduesDf;
  }
  set mostPtentResiduesDf(df: DG.DataFrame) {
    this._mostPotentResiduesDf = df;
  }

  get alphabet(): string {
    const col = this.settings.sequenceColumnName ? this.df.getCol(this.settings.sequenceColumnName) :
      this.df.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)!;
    return col.getTag(bio.TAGS.alphabet);
  }

  get substitutionsInfo(): type.SubstitutionsInfo {
    if (this._substitutionsInfo)
      return this._substitutionsInfo;

    const scaledActivityCol: DG.Column<number> = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    //TODO: set categories ordering the same to share compare indexes instead of strings
    const monomerColumns: type.RawColumn[] = this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER).map(extractMonomerInfo);
    this._substitutionsInfo = findMutations(scaledActivityCol.getRawData(), monomerColumns, this.settings);
    return this._substitutionsInfo;
  }
  set subsitutionsInfo(si: type.SubstitutionsInfo) {
    this._substitutionsInfo = si;
  }

  get clusterStatsDf(): DG.DataFrame {
    this._clusterStatsDf ??= this.calculateClusterStatistics();
    return this._clusterStatsDf;
  }
  set clusterStatsDf(df: DG.DataFrame) {
    this._clusterStatsDf = df;
  }

  get cp(): bio.SeqPalette {
    this._cp ??= bio.pickUpPalette(this.df.getCol(C.COLUMNS_NAMES.MACROMOLECULE));
    return this._cp;
  }
  set cp(_cp: bio.SeqPalette) {
    this._cp = _cp;
  }

  get analysisView(): DG.TableView {
    const shell = grok.shell;
    if (this.df.getTag('newAnalysis') !== '1') {
      this._analysisView = wu(shell.tableViews).find(({dataFrame}) => dataFrame.tags[C.PEPTIDES_ANALYSIS] === '1')!;
      grok.shell.v = this._analysisView;
    }

    this._analysisView ??= shell.addTableView(this.df);
    this.df.setTag('newAnalysis', '');
    return this._analysisView;
  }

  get onMutationCliffsSelectionChanged(): rxjs.Observable<undefined> {
    return this._mutatinCliffsSelectionSubject.asObservable();
  }

  get onSettingsChanged(): rxjs.Observable<type.PeptidesSettings> {
    return this.settingsSubject.asObservable();
  }

  get mutationCliffsSelection(): type.PositionToAARList {
    this._mutationCliffsSelection ??= JSON.parse(this.df.tags[C.TAGS.SELECTION] || '{}');
    return this._mutationCliffsSelection;
  }

  set mutationCliffsSelection(selection: type.PositionToAARList) {
    this._mutationCliffsSelection = selection;
    this.df.tags[C.TAGS.SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this._mutatinCliffsSelectionSubject.next();
  }

  get invariantMapSelection(): type.PositionToAARList {
    this._invariantMapSelection ??= JSON.parse(this.df.tags[C.TAGS.FILTER] || '{}');
    return this._invariantMapSelection;
  }

  set invariantMapSelection(selection: type.PositionToAARList) {
    this._invariantMapSelection = selection;
    this.df.tags[C.TAGS.FILTER] = JSON.stringify(selection);
    this.isInvariantMapTrigger = true;
    this.df.filter.fireChanged();
    this.isInvariantMapTrigger = false;
    this.analysisView.grid.invalidate();
  }

  get logoSummarySelection(): number[] {
    this._logoSummarySelection ??= JSON.parse(this.df.tags[C.TAGS.CLUSTER_SELECTION] || '[]');
    return this._logoSummarySelection;
  }

  set logoSummarySelection(selection: number[]) {
    this._logoSummarySelection = selection;
    this.df.tags[C.TAGS.CLUSTER_SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.analysisView.grid.invalidate();
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
    this.updateDefault();
    this.settingsSubject.next(this.settings);
  }

  createMonomerPositionDf(): DG.DataFrame {
    const matrixDf = this.monomerPositionStatsDf.groupBy([C.COLUMNS_NAMES.MONOMER])
      .pivot(C.COLUMNS_NAMES.POSITION)
      .add('first', C.COLUMNS_NAMES.MEAN_DIFFERENCE, '')
      .aggregate();
    matrixDf.name = 'SAR';

    return matrixDf;
  }

  buildMatrixDf(): DG.DataFrame {
    const splitSeqDfColumns = this.splitSeqDf.columns;
    const positionColumns = splitSeqDfColumns.names();
    return this.splitSeqDf
      .groupBy(positionColumns)
      .aggregate()
      .unpivot([], positionColumns, C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER);
  }

  buildSplitSeqDf(): DG.DataFrame {
    const sequenceCol = this.df.getCol(C.COLUMNS_NAMES.MACROMOLECULE);
    const splitSeqDf = splitAlignedSequences(sequenceCol);
    const splitSeqDfColumns = splitSeqDf.columns;

    for (const splitSeqCol of splitSeqDfColumns)
      splitSeqCol.name = `p${splitSeqCol.name}`;

    return splitSeqDf;
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
    if (!this._isUpdating || !this.isInitialized) {
      this._isUpdating = true;
      this.initializeViewersComponents();

      this.analysisView.grid.invalidate();
      this._isUpdating = false;
    }
  }

  initializeViewersComponents(): void {
    this.joinDataFrames();

    this.sortSourceGrid();

    this.createScaledCol();

    this.initSelections();

    this.setWebLogoInteraction();
    this.webLogoBounds = {};

    this.setCellRenderers();

    this.setTooltips();

    this.setBitsetCallback();

    this.postProcessGrids();
  }

  initSelections(): void {
    const tempInvariantMapSelection: type.PositionToAARList = this.invariantMapSelection;
    const mutationCliffsSelection: type.PositionToAARList = this.mutationCliffsSelection;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns) {
      tempInvariantMapSelection[pos] ??= [];
      mutationCliffsSelection[pos] ??= [];
    }
    this.invariantMapSelection = tempInvariantMapSelection;
    this.mutationCliffsSelection = mutationCliffsSelection;
  }

  joinDataFrames(): void {
    // append splitSeqDf columns to source table and make sure columns are not added more than once
    const name = this.df.name;
    const cols = this.df.columns;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const colName of positionColumns) {
      const col = this.df.col(colName);
      const newCol = this.splitSeqDf.getCol(colName);
      if (col === null)
        cols.add(newCol);
      else {
        cols.remove(colName);
        cols.add(newCol);
      }
      CR.setAARRenderer(newCol, this.alphabet, this.analysisView.grid);
    }
    this.df.name = name;
  }

  sortSourceGrid(): void {
    const colNames: DG.GridColumn[] = [];
    const sourceGridCols = this.analysisView.grid.columns;
    const sourceGridColsCount = sourceGridCols.length;
    for (let i = 1; i < sourceGridColsCount; i++)
      colNames.push(sourceGridCols.byIndex(i)!);

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
    sourceGridCols.setOrder(colNames.map((v) => v.name));
  }

  createScaledCol(): void {
    const sourceGrid = this.analysisView.grid;
    const scaledCol = scaleActivity(this.df.getCol(C.COLUMNS_NAMES.ACTIVITY), this.settings.scaling);
    //TODO: make another func
    this.df.columns.replace(C.COLUMNS_NAMES.ACTIVITY_SCALED, scaledCol);
    const gridCol = sourceGrid.col(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    if (gridCol)
      gridCol.name = scaledCol.getTag('gridName');

    sourceGrid.columns.setOrder([scaledCol.getTag('gridName')]);
  }

  calculateMonomerPositionStatistics(): DG.DataFrame {
    const matrixDf = this.matrixDf.groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER]).aggregate();
    const matrixLen = matrixDf.rowCount;

    const posRawColumns: type.RawColumn[] = [];

    const posCol: DG.Column<string> = matrixDf.getCol(C.COLUMNS_NAMES.POSITION);
    const posColData = posCol.getRawData();
    const posColCategories = posCol.categories;
    for (const position of posColCategories) {
      const currentCol = this.df.getCol(position);
      posRawColumns.push({
        name: position,
        rawData: currentCol.getRawData(),
        cat: currentCol.categories,
      });
    }

    const monomerCol: DG.Column<string> = matrixDf.getCol(C.COLUMNS_NAMES.MONOMER);
    const monomerColData = monomerCol.getRawData();
    const monomerColCategories = monomerCol.categories;

    //calculate p-values based on t-test
    const matrixCols = matrixDf.columns;
    const mdColData = matrixCols.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE).getRawData();
    const pValColData = matrixCols.addNewFloat(C.COLUMNS_NAMES.P_VALUE).getRawData();
    const countColData = matrixCols.addNewInt(C.COLUMNS_NAMES.COUNT).getRawData();
    const ratioColData = matrixCols.addNewFloat(C.COLUMNS_NAMES.RATIO).getRawData();
    const activityColData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();
    const sourceDfLen = activityColData.length;

    for (let i = 0; i < matrixLen; i++) {
      const positionRawIdx = posColData[i];
      const currentPosRawCol = posRawColumns[positionRawIdx];
      const monomerRawIdx = monomerColData[i];
      const mask: boolean[] = new Array(sourceDfLen);

      let trueCount = 0;
      for (let j = 0; j < sourceDfLen; ++j) {
        mask[j] = currentPosRawCol.cat![currentPosRawCol.rawData[j]] == monomerColCategories[monomerRawIdx];

        if (mask[j])
          ++trueCount;
      }

      const maskInfo = {
        trueCount: trueCount,
        falseCount: sourceDfLen - trueCount,
        mask: mask,
      };

      const stats = getStats(activityColData, maskInfo);

      mdColData[i] = stats.meanDifference;
      pValColData[i] = stats.pValue;
      countColData[i] = stats.count;
      ratioColData[i] = stats.ratio;
    }

    return matrixDf as DG.DataFrame;
  }

  calculateClusterStatistics(): DG.DataFrame {
    const originalClustersCol = this.df.getCol(C.COLUMNS_NAMES.CLUSTERS);
    const originalClustersColData = originalClustersCol.getRawData();
    const originalClustersColCategories = originalClustersCol.categories;

    const statsDf = this.df.groupBy([C.COLUMNS_NAMES.CLUSTERS]).aggregate();
    const clustersCol = statsDf.getCol(C.COLUMNS_NAMES.CLUSTERS);
    clustersCol.setCategoryOrder(originalClustersColCategories);
    const clustersColData = clustersCol.getRawData();

    const statsDfCols = statsDf.columns;
    const mdColData = statsDfCols.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE).getRawData();
    const pValColData = statsDfCols.addNewFloat(C.COLUMNS_NAMES.P_VALUE).getRawData();
    const countColData = statsDfCols.addNewInt(C.COLUMNS_NAMES.COUNT).getRawData();
    const ratioColData = statsDfCols.addNewFloat(C.COLUMNS_NAMES.RATIO).getRawData();
    const activityColData: type.RawData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();
    const activityColLen = activityColData.length;

    for (let rowIdx = 0; rowIdx < clustersColData.length; ++rowIdx) {
      const clusterIdx = clustersColData[rowIdx];
      const mask = new Array(activityColLen);
      let trueCount = 0;
      for (let maskIdx = 0; maskIdx < activityColLen; ++maskIdx) {
        mask[maskIdx] = clusterIdx == originalClustersColData[maskIdx];

        if (mask[maskIdx])
          ++trueCount;
      }

      const maskInfo = {
        trueCount: trueCount,
        falseCount: activityColLen - trueCount,
        mask: mask,
      };

      const stats = getStats(activityColData, maskInfo);

      mdColData[rowIdx] = stats.meanDifference;
      pValColData[rowIdx] = stats.pValue;
      countColData[rowIdx] = stats.count;
      ratioColData[rowIdx] = stats.ratio;
    }
    return statsDf;
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

  setWebLogoInteraction(): void {
    const sourceView = this.analysisView.grid;
    const eventAction = (ev: MouseEvent): void => {
      const cell = sourceView.hitTest(ev.offsetX, ev.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType == C.SEM_TYPES.MONOMER) {
        const newBarPart = this.findAARandPosition(cell, ev);
        this.requestBarchartAction(ev, newBarPart);
      }
    };

    // The following events makes the barchart interactive
    rxjs.fromEvent<MouseEvent>(sourceView.overlay, 'mousemove')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
    rxjs.fromEvent<MouseEvent>(sourceView.overlay, 'click')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  }

  findAARandPosition(cell: DG.GridCell, ev: MouseEvent): { monomer: string, position: string } | null {
    const barCoords = this.webLogoBounds[cell.tableColumn!.name];
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
      ev.shiftKey ? this.modifyMonomerPositionSelection(monomer, position, false) :
        this.initMonomerPositionSelection(monomer, position, false);
    } else {
      const bar = `${position} = ${monomer}`;
      if (this.cachedWebLogoTooltip.bar == bar)
        ui.tooltip.show(this.cachedWebLogoTooltip.tooltip!, ev.clientX, ev.clientY);
      else
        this.cachedWebLogoTooltip = {bar: bar, tooltip: this.showTooltipAt(monomer, position, ev.clientX, ev.clientY)};

      //TODO: how to unghighlight?
      // this.df.rows.match(bar).highlight();
    }
  }

  setCellRenderers(): void {
    const sourceGrid = this.analysisView.grid;
    sourceGrid.setOptions({'colHeaderHeight': 130});
    sourceGrid.onCellRender.subscribe((gcArgs) => {
      const ctx = gcArgs.g;
      const bounds = gcArgs.bounds;
      const col = gcArgs.cell.tableColumn;

      ctx.save();
      ctx.beginPath();
      ctx.rect(bounds.x, bounds.y, bounds.width, bounds.height);
      ctx.clip();

      //TODO: optimize
      if (gcArgs.cell.isColHeader && col?.semType == C.SEM_TYPES.MONOMER) {
        const monomerStatsCol: DG.Column<string> = this.monomerPositionStatsDf.getCol(C.COLUMNS_NAMES.MONOMER);
        const positionStatsCol: DG.Column<string> = this.monomerPositionStatsDf.getCol(C.COLUMNS_NAMES.POSITION);
        const rowMask = DG.BitSet.create(this.monomerPositionStatsDf.rowCount,
          (i) => positionStatsCol.get(i) === col.name);
        //TODO: precalc on stats creation
        const sortedStatsOrder = this.monomerPositionStatsDf.getSortedOrder([C.COLUMNS_NAMES.COUNT], [false], rowMask)
          .sort((a, b) => {
            if (monomerStatsCol.get(a) === '-' || monomerStatsCol.get(a) === '')
              return -1;
            else if (monomerStatsCol.get(b) === '-' || monomerStatsCol.get(b) === '')
              return +1;
            return 0;
          });
        const statsInfo: type.StatsInfo = {
          countCol: this.monomerPositionStatsDf.getCol(C.COLUMNS_NAMES.COUNT),
          monomerCol: monomerStatsCol,
          orderedIndexes: sortedStatsOrder,
        };

        this.webLogoBounds[col.name] =
          CR.drawLogoInBounds(ctx, bounds, statsInfo, this.df.rowCount, this.cp, this.headerSelectedMonomers[col.name]);
        gcArgs.preventDefault();
      }

      ctx.restore();
    });
  }

  setTooltips(): void {
    this.analysisView.grid.onCellTooltip((cell, x, y) => {
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

    const mw = getMonomerWorks();
    const mol = mw?.getCappedRotatedMonomer('PEPTIDE', aar);

    if (mol) {
      tooltipElements.push(ui.div(monomerName));
      const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      tooltipElements.push(grok.chem.svgMol(mol, undefined, undefined, options));
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

  setBitsetCallback(): void {
    if (this.isBitsetChangedInitialized)
      return;
    const selection = this.df.selection;
    const filter = this.df.filter;
    const clusterCol = this.df.col(C.COLUMNS_NAMES.CLUSTERS);
    const clusterColData = clusterCol?.getRawData();

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
        if (clusterColData && this.logoSummarySelection.includes(clusterColData[i]))
          return true;
        return false;
      };
      currentBitset.init(getBitAt, false);

      updateEdfSelection();
    };

    selection.onChanged.subscribe(() => changeSelectionBitset(selection));

    filter.onChanged.subscribe(() => {
      const positionList = Object.keys(this.invariantMapSelection);
      const invariantMapBitset = DG.BitSet.create(filter.length, (index) => {
        let result = true;
        for (const position of positionList) {
          const aarList = this.invariantMapSelection[position];
          result &&= aarList.length === 0 || aarList.includes(this.df.get(position, index));
          if (!result)
            return result;
        }
        return result;
      });

      if (!this.isInvariantMapTrigger)
        this.initBitset = filter.clone();

      const temp = invariantMapBitset.and(this.initBitset);
      filter.init((i) => temp.get(i), false);
    });
    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(isPeptideSpaceSource: boolean = false): void {
    this.isPeptideSpaceChangingBitset = isPeptideSpaceSource;
    this.df.selection.fireChanged();
    this.modifyOrCreateSplitCol();
    this.headerSelectedMonomers = calculateSelected(this.df);
    grok.shell.o = this.createAccordion().root;
    this.isPeptideSpaceChangingBitset = false;
  }

  postProcessGrids(): void {
    const posCols = this.splitSeqDf.columns.names();
    const sourceGrid = this.analysisView.grid;
    const sourceGridCols = sourceGrid.columns;
    const sourceGridColsLen = sourceGridCols.length;
    const visibleColumns = Object.keys(this.settings.columns ?? {});
    for (let gcIndex = 0; gcIndex < sourceGridColsLen; ++gcIndex) {
      const gridCol = sourceGridCols.byIndex(gcIndex)!;
      const tableCol = gridCol.column;
      if (!tableCol)
        continue;

      const gridColName = gridCol.name;
      if (posCols.includes(gridColName))
        gridCol.name = gridColName.substring(1);

      const tableColName = tableCol.name;
      gridCol.visible =
        tableCol.semType === C.SEM_TYPES.MONOMER ||
        tableColName === C.COLUMNS_NAMES.ACTIVITY_SCALED ||
        visibleColumns.includes(tableColName ?? '');
    }

    const sourceGridProps = sourceGrid.props;
    sourceGridProps.allowColSelection = false;
    sourceGridProps.allowEdit = false;
    setTimeout(() => sourceGridProps.rowHeight = 20, 10);
    sourceGridProps.allowRowResizing = false;
    sourceGridProps.showCurrentRowIndicator = false;
    this.df.temp[C.EMBEDDING_STATUS] = false;

    for (let i = 0; i < sourceGridColsLen; i++) {
      const currentCol = sourceGridCols.byIndex(i);
      if (currentCol) {
        if (currentCol.column?.getTag(C.TAGS.VISIBLE) === '0')
          currentCol.visible = false;

        currentCol.width = isNaN(parseInt(currentCol.name)) ? 50 : 40;
      }
    }
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
  init(): void {
    if (this.isInitialized)
      return;
    this.isInitialized = true;

    if (!this.isRibbonSet) {
      //TODO: don't pass model, pass parameters instead
      const settingsButton = ui.bigButton('Settings', () => getSettingsDialog(this), 'Peptides analysis settings');
      this.analysisView.setRibbonPanels([[settingsButton]], false);
      this.isRibbonSet = true;
    }

    this.df.tags[C.PEPTIDES_ANALYSIS] = '1';

    this.updateDefault();

    this.analysisView.grid.invalidate();
  }

  async addViewers(): Promise<void> {
    const dockManager = this.analysisView.dockManager;
    const dfPlt = this.df.plot;

    const mutationCliffsViewer = await dfPlt.fromType('peptide-sar-viewer') as MonomerPosition;
    const mostPotentResiduesViewer = await dfPlt.fromType('peptide-sar-viewer-vertical') as MostPotentResiduesViewer;
    if (this.df.getTag(C.TAGS.CLUSTERS))
      await this.addLogoSummaryTableViewer();

    const mcNode = dockManager.dock(mutationCliffsViewer, DG.DOCK_TYPE.DOWN, null, mutationCliffsViewer.name);

    dockManager.dock(mostPotentResiduesViewer, DG.DOCK_TYPE.RIGHT, mcNode, mostPotentResiduesViewer.name, 0.3);
  }

  async addLogoSummaryTableViewer(): Promise<void> {
    const logoSummary = await this.df.plot.fromType('logo-summary-viewer') as LogoSummary;
    this.analysisView.dockManager.dock(logoSummary, DG.DOCK_TYPE.RIGHT, null, 'Logo Summary Table');
  }
}
