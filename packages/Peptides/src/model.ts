import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import * as rxjs from 'rxjs';
import * as uuid from 'uuid';

import * as C from './utils/constants';
import * as type from './utils/types';
import {calculateSelected, extractMonomerInfo, scaleActivity, wrapDistroAndStatsDefault} from './utils/misc';
import {MonomerPosition, MostPotentResiduesViewer} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getDistributionAndStats, getDistributionWidget} from './widgets/distribution';
import {getStats, Stats} from './utils/statistics';
import {LogoSummary} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {_package, getMonomerWorksInstance, getTreeHelperInstance} from './package';
import {findMutations} from './utils/algorithms';
import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {MonomerWorks} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {pickUpPalette, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';
import {createDistanceMatrixWorker} from './utils/worker-creator';

export type SummaryStats = {
  minCount: number, maxCount: number,
  minMeanDifference: number, maxMeanDifference: number,
  minPValue: number, maxPValue: number,
  minRatio: number, maxRatio: number,
};
export type PositionStats = { [monomer: string]: Stats } & { general: SummaryStats };
export type MonomerPositionStats = { [position: string]: PositionStats } & { general: SummaryStats };
export type ClusterStats = {[cluster: string]: Stats};
export enum CLUSTER_TYPE {
  ORIGINAL = 'original',
  CUSTOM = 'custom',
};
export type ClusterType = `${CLUSTER_TYPE}`;
export type ClusterTypeStats = {[clusterType in ClusterType]: ClusterStats};

export class PeptidesModel {
  static modelName = 'peptidesModel';

  settingsSubject: rxjs.Subject<type.PeptidesSettings> = new rxjs.Subject();
  _mutatinCliffsSelectionSubject: rxjs.Subject<undefined> = new rxjs.Subject();
  _newClusterSubject: rxjs.Subject<undefined> = new rxjs.Subject();
  _removeClusterSubject: rxjs.Subject<undefined> = new rxjs.Subject();
  _filterChangedSubject: rxjs.Subject<undefined> = new rxjs.Subject();

  _isUpdating: boolean = false;
  isBitsetChangedInitialized = false;
  isCellChanging = false;
  isUserChangedSelection = true;

  df: DG.DataFrame;
  splitCol!: DG.Column<boolean>;
  // edf: DG.DataFrame | null = null;
  _monomerPositionStats?: MonomerPositionStats;
  _clusterStats?: ClusterTypeStats;
  _mutationCliffsSelection!: type.PositionToAARList;
  _invariantMapSelection!: type.PositionToAARList;
  _logoSummarySelection!: string[];
  _substitutionsInfo?: type.SubstitutionsInfo;
  isInitialized = false;
  _analysisView?: DG.TableView;

  // isPeptideSpaceChangingBitset = false;
  // isChangingEdfBitset = false;

  monomerMap: { [key: string]: { molfile: string, fullName: string } } = {};
  monomerLib: IMonomerLib | null = null; // To get monomers from lib(s)
  monomerWorks: MonomerWorks | null = null; // To get processed monomers

  _settings!: type.PeptidesSettings;
  isRibbonSet = false;

  _cp?: SeqPalette;
  initBitset: DG.BitSet;
  isInvariantMapTrigger: boolean = false;
  headerSelectedMonomers: type.MonomerSelectionStats = {};
  webLogoBounds: { [positon: string]: { [monomer: string]: DG.Rect } } = {};
  cachedWebLogoTooltip: { bar: string; tooltip: HTMLDivElement | null; } = {bar: '', tooltip: null};
  _monomerPositionDf?: DG.DataFrame;
  _alphabet?: string;
  _mostPotentResiduesDf?: DG.DataFrame;
  _matrixDf?: DG.DataFrame;
  _splitSeqDf?: DG.DataFrame;
  _distanceMatrix!: DistanceMatrix;
  _treeHelper!: ITreeHelper;
  _dm!: DistanceMatrix;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
    this.initBitset = this.df.filter.clone();
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    (dataFrame.temp[PeptidesModel.modelName] as PeptidesModel).init();
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  // get distanceMatrix(): DistanceMatrix {
  //   this._distanceMatrix ??= this.calculateDistanceMatrix();
  //   return this._distanceMatrix;
  // }

  get treeHelper(): ITreeHelper {
    this._treeHelper ??= getTreeHelperInstance();
    return this._treeHelper;
  }

  get monomerPositionDf(): DG.DataFrame {
    this._monomerPositionDf ??= this.createMonomerPositionDf();
    return this._monomerPositionDf;
  }

  set monomerPositionDf(df: DG.DataFrame) {
    this._monomerPositionDf = df;
  }

  get monomerPositionStats(): MonomerPositionStats {
    this._monomerPositionStats ??= this.calculateMonomerPositionStatistics();
    return this._monomerPositionStats;
  }

  set monomerPositionStats(mps: MonomerPositionStats) {
    this._monomerPositionStats = mps;
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
    this._mostPotentResiduesDf ??= this.createMostPotentResiduesDf();
    return this._mostPotentResiduesDf;
  }

  set mostPotentResiduesDf(df: DG.DataFrame) {
    this._mostPotentResiduesDf = df;
  }

  get alphabet(): string {
    const col = this.df.getCol(this.settings.sequenceColumnName!);
    return col.getTag(bioTAGS.alphabet);
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

  set substitutionsInfo(si: type.SubstitutionsInfo) {
    this._substitutionsInfo = si;
  }

  get clusterStats(): ClusterTypeStats {
    this._clusterStats ??= this.calculateClusterStatistics();
    return this._clusterStats;
  }

  set clusterStats(clusterStats: ClusterTypeStats) {
    this._clusterStats = clusterStats;
  }

  get cp(): SeqPalette {
    this._cp ??= pickUpPalette(this.df.getCol(this.settings.sequenceColumnName!));
    return this._cp;
  }

  set cp(_cp: SeqPalette) {
    this._cp = _cp;
  }

  get analysisView(): DG.TableView {
    this._analysisView ??=
      wu(grok.shell.tableViews).find(({dataFrame}) => dataFrame.getTag(C.TAGS.UUID) == this.df.getTag(C.TAGS.UUID)) ??
        grok.shell.addTableView(this.df);
    if (this.df.getTag(C.MULTIPLE_VIEWS) != '1')
      grok.shell.v = this._analysisView;

    return this._analysisView;
  }

  get onMutationCliffsSelectionChanged(): rxjs.Observable<undefined> {
    return this._mutatinCliffsSelectionSubject.asObservable();
  }

  get onSettingsChanged(): rxjs.Observable<type.PeptidesSettings> {
    return this.settingsSubject.asObservable();
  }

  get onNewCluster(): rxjs.Observable<undefined> {
    return this._newClusterSubject.asObservable();
  }

  get onRemoveCluster(): rxjs.Observable<undefined> {
    return this._removeClusterSubject.asObservable();
  }

  get onFilterChanged(): rxjs.Observable<undefined> {
    return this._filterChangedSubject.asObservable();
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
    this.analysisView.grid.invalidate();
  }

  get invariantMapSelection(): type.PositionToAARList {
    this._invariantMapSelection ??= JSON.parse(this.df.tags[C.TAGS.FILTER] || '{}');
    return this._invariantMapSelection;
  }

  set invariantMapSelection(selection: type.PositionToAARList) {
    this._invariantMapSelection = selection;
    this.df.tags[C.TAGS.FILTER] = JSON.stringify(selection);
    this.isInvariantMapTrigger = true;
    this.fireBitsetChanged(false, true);
    this.isInvariantMapTrigger = false;
    this.analysisView.grid.invalidate();
  }

  get logoSummarySelection(): string[] {
    this._logoSummarySelection ??= JSON.parse(this.df.tags[C.TAGS.CLUSTER_SELECTION] || '[]');
    return this._logoSummarySelection;
  }

  set logoSummarySelection(selection: string[]) {
    this._logoSummarySelection = selection;
    this.df.tags[C.TAGS.CLUSTER_SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.analysisView.grid.invalidate();
  }

  get splitByPos(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[0];
    return splitByPosFlag == '1' ? true : false;
  }

  set splitByPos(flag: boolean) {
    const splitByAARFlag = (this.df.tags['distributionSplit'] || '00')[1];
    this.df.tags['distributionSplit'] = `${flag ? 1 : 0}${splitByAARFlag}`;
  }

  get splitByAAR(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[1];
    return splitByPosFlag == '1' ? true : false;
  }

  set splitByAAR(flag: boolean) {
    const splitByAARFlag = (this.df.tags['distributionSplit'] || '00')[0];
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

  get customClusters(): Iterable<DG.Column<boolean>> {
    const query: { [key: string]: string } = {};
    query[C.TAGS.CUSTOM_CLUSTER] = '1';
    return this.df.columns.byTags(query);
  }

  get settings(): type.PeptidesSettings {
    this._settings ??= JSON.parse(this.df.getTag('settings') || '{}');
    return this._settings;
  }

  set settings(s: type.PeptidesSettings) {
    const newSettingsEntries = Object.entries(s);
    const updateVars: Set<string> = new Set();
    for (const [key, value] of newSettingsEntries) {
      this._settings[key as keyof type.PeptidesSettings] = value as any;
      switch (key) {
      case 'scaling':
        updateVars.add('activity');
        updateVars.add('mutationCliffs');
        updateVars.add('stats');
        break;
      // case 'columns':
      //   updateVars.add('grid');
      //   break;
      case 'maxMutations':
      case 'minActivityDelta':
        updateVars.add('mutationCliffs');
        break;
      case 'showDendrogram':
        updateVars.add('dendrogram');
        break;
      }
    }
    this.df.setTag('settings', JSON.stringify(this._settings));
    for (const variable of updateVars) {
      switch (variable) {
      case 'activity':
        this.createScaledCol();
        break;
      case 'mutationCliffs':
        const scaledActivityCol: DG.Column<number> = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
        //TODO: set categories ordering the same to share compare indexes instead of strings
        const monomerCols: type.RawColumn[] = this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER).map(extractMonomerInfo);
        this.substitutionsInfo = findMutations(scaledActivityCol.getRawData(), monomerCols, this.settings);
        break;
      case 'stats':
        this.monomerPositionStats = this.calculateMonomerPositionStatistics();
        this.monomerPositionDf = this.createMonomerPositionDf();
        this.mostPotentResiduesDf = this.createMostPotentResiduesDf();
        this.clusterStats = this.calculateClusterStatistics();
        break;
      case 'grid':
        this.updateGrid();
        break;
      case 'dendrogram':
        this.settings.showDendrogram ? this.addDendrogram() : this.closeDendrogram();
        break;
      }
    }

    //TODO: handle settings change
    this.settingsSubject.next(this.settings);
  }

  createMonomerPositionDf(): DG.DataFrame {
    const positions = this.splitSeqDf.columns.names();
    const matrixDf = this.matrixDf
      .groupBy([C.COLUMNS_NAMES.MONOMER])
      .aggregate();
    for (const pos of positions)
      matrixDf.columns.addNewString(pos);

    const monomerCol = matrixDf.getCol(C.COLUMNS_NAMES.MONOMER);
    for (let i = 0; i < monomerCol.length; ++i) {
      if (monomerCol.get(i) == '') {
        matrixDf.rows.removeAt(i);
        break;
      }
    }
    matrixDf.name = 'SAR';

    return matrixDf;
  }

  buildMatrixDf(): DG.DataFrame {
    const splitSeqDfColumns = this.splitSeqDf.columns;
    const positionColumns = splitSeqDfColumns.names();
    return this.splitSeqDf
      .groupBy(positionColumns)
      .aggregate()
      .unpivot([], positionColumns, C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER)
      .groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER])
      .aggregate();
  }

  buildSplitSeqDf(): DG.DataFrame {
    const sequenceCol = this.df.getCol(this.settings.sequenceColumnName!);
    const splitSeqDf = splitAlignedSequences(sequenceCol);

    return splitSeqDf;
  }

  getCompoundBitest(): DG.BitSet {
    return this.df.selection.clone().and(this.df.filter);
  }

  createAccordion(): DG.Accordion | null {
    const trueModel: PeptidesModel | undefined = grok.shell.t.temp[PeptidesModel.modelName];
    if (!trueModel)
      return null;

    const acc = ui.accordion();
    acc.root.style.width = '100%';
    const filterAndSelectionBs = trueModel.getCompoundBitest();
    const filteredTitlePart = trueModel.df.filter.anyFalse ? ` among ${trueModel.df.filter.trueCount} filtered` : '';
    acc.addTitle(ui.h1(`${filterAndSelectionBs.trueCount} selected rows${filteredTitlePart}`));
    if (filterAndSelectionBs.anyTrue) {
      acc.addPane('Actions', () => {
        const newViewButton = ui.button('New view', async () => await trueModel.createNewView(),
          'Creates a new view from current selection');
        const newCluster = ui.button('New cluster', () => trueModel._newClusterSubject.next(),
          'Creates a new cluster from selection');
        const removeCluster = ui.button('Remove cluster', () => trueModel._removeClusterSubject.next(),
          'Removes currently selected custom cluster');
        return ui.divV([newViewButton, newCluster, removeCluster]);
      });
    }
    const table = trueModel.df.filter.anyFalse ? trueModel.df.clone(trueModel.df.filter, null, true) : trueModel.df;
    acc.addPane('Mutation Cliff pairs', () => mutationCliffsWidget(trueModel.df, trueModel).root);
    acc.addPane('Distribution', () => getDistributionWidget(table, trueModel).root);

    return acc;
  }

  updateGrid(): void {
    this.joinDataFrames();

    this.sortSourceGrid();

    this.createScaledCol();

    this.initInvariantMapSelection();
    this.initMutationCliffsSelection();

    this.setWebLogoInteraction();
    this.webLogoBounds = {};

    this.setCellRenderers();

    this.setTooltips();

    this.setBitsetCallback();

    this.postProcessGrids();
  }

  initInvariantMapSelection(cleanInit = false): void {
    const tempInvariantMapSelection: type.PositionToAARList = this.invariantMapSelection;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns) {
      if (cleanInit || !tempInvariantMapSelection.hasOwnProperty(pos))
        tempInvariantMapSelection[pos] = [];
    }

    this.invariantMapSelection = tempInvariantMapSelection;
  }

  initMutationCliffsSelection(cleanInit = false): void {
    const mutationCliffsSelection: type.PositionToAARList = this.mutationCliffsSelection;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns) {
      if (cleanInit || !mutationCliffsSelection.hasOwnProperty(pos))
        mutationCliffsSelection[pos] = [];
    }

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
      CR.setAARRenderer(newCol, this.alphabet);
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
    const scaledCol = scaleActivity(this.df.getCol(this.settings.activityColumnName!), this.settings.scaling);
    //TODO: make another func
    this.df.columns.replace(C.COLUMNS_NAMES.ACTIVITY_SCALED, scaledCol);

    sourceGrid.columns.setOrder([scaledCol.name]);
  }

  calculateMonomerPositionStatistics(): MonomerPositionStats {
    const positionColumns = this.splitSeqDf.columns.toList();
    const monomerPositionObject = {general: {}} as MonomerPositionStats & { general: SummaryStats };
    const activityColData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();
    const sourceDfLen = activityColData.length;

    for (const posCol of positionColumns) {
      const posColData = posCol.getRawData();
      const currentMonomerSet = posCol.categories;
      const currentPositionObject = {general: {}} as PositionStats & { general: SummaryStats };

      for (const [categoryIndex, monomer] of currentMonomerSet.entries()) {
        if (monomer == '')
          continue;

        const mask: boolean[] = new Array(sourceDfLen);
        let trueCount = 0;
        for (let j = 0; j < sourceDfLen; ++j) {
          mask[j] = posColData[j] == categoryIndex;

          if (mask[j])
            ++trueCount;
        }

        const maskInfo = {
          trueCount: trueCount,
          falseCount: sourceDfLen - trueCount,
          mask: mask,
        };

        const stats = getStats(activityColData, maskInfo);
        currentPositionObject[monomer] = stats;

        this.getSummaryStats(currentPositionObject.general, stats);
      }
      monomerPositionObject[posCol.name] = currentPositionObject;
      this.getSummaryStats(monomerPositionObject.general, null, currentPositionObject.general);
    }
    return monomerPositionObject;
  }

  getSummaryStats(genObj: SummaryStats, stats: Stats | null = null, summaryStats: SummaryStats | null = null): void {
    if (stats == null && summaryStats == null)
      throw new Error(`MonomerPositionStatsError: either stats or summaryStats must be present`);

    const possibleMaxCount = stats?.count ?? summaryStats!.maxCount;
    genObj.maxCount ??= possibleMaxCount;
    if (genObj.maxCount < possibleMaxCount)
      genObj.maxCount = possibleMaxCount;

    const possibleMinCount = stats?.count ?? summaryStats!.minCount;
    genObj.minCount ??= possibleMinCount;
    if (genObj.minCount > possibleMinCount)
      genObj.minCount = possibleMinCount;

    const possibleMaxMeanDifference = stats?.meanDifference ?? summaryStats!.maxMeanDifference;
    genObj.maxMeanDifference ??= possibleMaxMeanDifference;
    if (genObj.maxMeanDifference < possibleMaxMeanDifference)
      genObj.maxMeanDifference = possibleMaxMeanDifference;

    const possibleMinMeanDifference = stats?.meanDifference ?? summaryStats!.minMeanDifference;
    genObj.minMeanDifference ??= possibleMinMeanDifference;
    if (genObj.minMeanDifference > possibleMinMeanDifference)
      genObj.minMeanDifference = possibleMinMeanDifference;

    const possibleMaxPValue = stats?.pValue ?? summaryStats!.maxPValue;
    genObj.maxPValue ??= possibleMaxPValue;
    if (genObj.maxPValue < possibleMaxPValue)
      genObj.maxPValue = possibleMaxPValue;

    const possibleMinPValue = stats?.pValue ?? summaryStats!.minPValue;
    genObj.minPValue ??= possibleMinPValue;
    if (genObj.minPValue > possibleMinPValue)
      genObj.minPValue = possibleMinPValue;

    const possibleMaxRatio = stats?.ratio ?? summaryStats!.maxRatio;
    genObj.maxRatio ??= possibleMaxRatio;
    if (genObj.maxRatio < possibleMaxRatio)
      genObj.maxRatio = possibleMaxRatio;

    const possibleMinRatio = stats?.ratio ?? summaryStats!.minRatio;
    genObj.minRatio ??= possibleMinRatio;
    if (genObj.minRatio > possibleMinRatio)
      genObj.minRatio = possibleMinRatio;
  }

  calculateClusterStatistics(): ClusterTypeStats {
    const origClustCol = this.df.getCol(this.settings.clustersColumnName!);
    const origClustColData = origClustCol.getRawData();
    const origClustColCat = origClustCol.categories;

    const customClustColList = wu(this.customClusters).toArray();
    const customClustColDataList = customClustColList.map((v) => v.toList() as boolean[]);
    const customClustColNamesList = customClustColList.map((v) => v.name);

    const rowCount = this.df.rowCount;

    const origClustMasks: boolean[][] = Array.from({length: origClustColCat.length},
      () => new Array(rowCount) as boolean[]);
    const customClustMasks: boolean[][] = Array.from({length: customClustColList.length},
      () => new Array(rowCount) as boolean[]);

    // get original cluster masks in one pass
    // complexity is O(N * (M + 1)) where N is the number of rows and M is the number of custom clusters
    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx) {
      origClustMasks[origClustColData[rowIdx]][rowIdx] = true;
      for (let customClustIdx = 0; customClustIdx < customClustColList.length; ++customClustIdx)
        customClustMasks[customClustIdx][rowIdx] = customClustColDataList[customClustIdx][rowIdx];
    }

    const activityColData: type.RawData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();

    const origClustStats: ClusterStats = {};
    const customClustStats: ClusterStats = {};

    for (let clustType = 0; clustType < 2; ++clustType) {
      const masks = clustType == 0 ? origClustMasks : customClustMasks;
      const clustNames = clustType == 0 ? origClustColCat : customClustColNamesList;
      const resultStats = clustType == 0 ? origClustStats : customClustStats;
      for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
        const mask = masks[maskIdx];
        const trueCount = mask.filter((v) => v).length;
        const maskInfo = {trueCount: trueCount, falseCount: rowCount - trueCount, mask: mask};
        const stats = getStats(activityColData, maskInfo);
        resultStats[clustNames[maskIdx]] = stats;
      }
    }

    const resultStats = {} as ClusterTypeStats;
    resultStats[CLUSTER_TYPE.ORIGINAL] = origClustStats;
    resultStats[CLUSTER_TYPE.CUSTOM] = customClustStats;
    return resultStats;
  }

  createMostPotentResiduesDf(): DG.DataFrame {
    const monomerPositionStatsEntries = Object.entries(this.monomerPositionStats) as [string, PositionStats][];
    const mprDf = DG.DataFrame.create(monomerPositionStatsEntries.length - 1); // Subtract 'general' entry from mp-stats
    const mprDfCols = mprDf.columns;
    const posCol = mprDfCols.addNewInt(C.COLUMNS_NAMES.POSITION);
    const monomerCol = mprDfCols.addNewString(C.COLUMNS_NAMES.MONOMER);
    const mdCol = mprDfCols.addNewFloat(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
    const pValCol = mprDfCols.addNewFloat(C.COLUMNS_NAMES.P_VALUE);
    const countCol = mprDfCols.addNewInt(C.COLUMNS_NAMES.COUNT);
    const ratioCol = mprDfCols.addNewFloat(C.COLUMNS_NAMES.RATIO);

    let i = 0;
    for (const [position, positionStats] of monomerPositionStatsEntries) {
      const generalPositionStats = positionStats.general;
      if (!generalPositionStats)
        continue;

      const filteredMonomerStats = Object.entries(positionStats).filter((v) => {
        const key = v[0];
        if (key == 'general')
          return false;

        return (v[1] as Stats).pValue == generalPositionStats.minPValue;
      }) as [string, Stats][];

      let maxEntry: [string, Stats];
      for (const [monomer, monomerStats] of filteredMonomerStats) {
        if (typeof maxEntry! == 'undefined' || maxEntry[1].meanDifference < monomerStats.meanDifference)
          maxEntry = [monomer, monomerStats];
      }

      posCol.set(i, parseInt(position));
      monomerCol.set(i, maxEntry![0]);
      mdCol.set(i, maxEntry![1].meanDifference);
      pValCol.set(i, maxEntry![1].pValue);
      countCol.set(i, maxEntry![1].count);
      ratioCol.set(i, maxEntry![1].ratio);
      ++i;
    }
    return mprDf;
  }

  modifyClusterSelection(cluster: string): void {
    const tempSelection = this.logoSummarySelection;
    const idx = tempSelection.indexOf(cluster);
    if (idx !== -1)
      tempSelection.splice(idx, 1);
    else
      tempSelection.push(cluster);

    this.logoSummarySelection = tempSelection;
  }

  initClusterSelection(cluster: string): void {
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
      try {
        ctx.beginPath();
        ctx.rect(bounds.x, bounds.y, bounds.width, bounds.height);
        ctx.clip();

        //TODO: optimize
        if (gcArgs.cell.isColHeader && col?.semType == C.SEM_TYPES.MONOMER) {
          const stats = this.monomerPositionStats[col.name];
          //TODO: precalc on stats creation
          const sortedStatsOrder = Object.keys(stats).sort((a, b) => {
            if (a == '' || a == '-')
              return -1;
            else if (b == '' || b == '-')
              return +1;
            return 0;
          }).filter((v) => v != 'general');

          this.webLogoBounds[col.name] = CR.drawLogoInBounds(ctx, bounds, stats, sortedStatsOrder, this.df.rowCount,
            this.cp, this.headerSelectedMonomers[col.name]);
          gcArgs.preventDefault();
        }
      } catch (e) {
        console.warn(`PeptidesHeaderLogoError: couldn't render WebLogo for column \`${col!.name}\`. ` +
          `See original error below.`);
        console.warn(e);
      } finally {
        ctx.restore();
      }
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

    const mw = getMonomerWorksInstance();
    const mol = mw?.getCappedRotatedMonomer('PEPTIDE', aar);

    if (mol) {
      tooltipElements.push(ui.div(monomerName));
      const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      tooltipElements.push(grok.chem.svgMol(mol, undefined, undefined, options));
    } else
      tooltipElements.push(ui.div(aar));

    ui.tooltip.show(ui.divV(tooltipElements), x, y);
  }

  //TODO: move out to viewer code
  showTooltipAt(aar: string, position: string, x: number, y: number): HTMLDivElement | null {
    const stats = this.monomerPositionStats[position][aar];
    if (!stats?.count)
      return null;

    const activityCol = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    //TODO: use bitset instead of splitCol
    const splitCol = DG.Column.bool(C.COLUMNS_NAMES.SPLIT_COL, activityCol.length);
    const currentPosCol = this.df.getCol(position);
    const indexes: number[] = [];
    splitCol.init((i) => {
      const sameMonomer = currentPosCol.get(i) == aar;
      if (sameMonomer)
        indexes.push(i);

      return sameMonomer;
    });
    const colResults: { [colName: string]: number } = {};
    for (const [col, agg] of Object.entries(this.settings.columns || {})) {
      const currentCol = this.df.getCol(col);
      const currentColData = currentCol.getRawData();
      const tempCol = DG.Column.float('', indexes.length);
      tempCol.init((i) => currentColData[indexes[i]]);
      colResults[`${agg}(${col})`] = tempCol.stats[agg as keyof DG.Stats] as number;
    }

    const distributionTable = DG.DataFrame.fromColumns([activityCol, splitCol]);
    const das = getDistributionAndStats(distributionTable, stats, `${position} : ${aar}`, 'Other', true);
    const resultMap: { [key: string]: any } = {...das.tableMap, ...colResults};
    const distroStatsElem = wrapDistroAndStatsDefault(das.labels, das.histRoot, resultMap, true);

    ui.tooltip.show(distroStatsElem, x, y);

    return distroStatsElem;
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
    const clusterCol = this.df.col(this.settings.clustersColumnName!);

    const changeSelectionBitset = (currentBitset: DG.BitSet, posList: type.RawColumn[], clustColCat: string[],
      clustColData: type.RawData, customClust: {[key: string]: boolean[]}): void => {
      const getBitAt = (i: number): boolean => {
        for (const posRawCol of posList) {
          if (this.mutationCliffsSelection[posRawCol.name].includes(posRawCol.cat![posRawCol.rawData[i]]))
            return true;
        }

        const currentOrigClust = clustColCat[clustColData[i]];
        if (typeof currentOrigClust === undefined)
          return false;

        for (const clust of this.logoSummarySelection) {
          if (clust === currentOrigClust)
            return true;

          if (Object.hasOwn(customClust, clust) && customClust[clust][i] === true)
            return true;
        }

        return false;
      };
      currentBitset.init((i) => getBitAt(i), false);
    };

    selection.onChanged.subscribe(() => {
      if (this.isUserChangedSelection)
        return;

      const positionList: type.RawColumn[] = Object.keys(this.mutationCliffsSelection).map((pos) => {
        const posCol = this.df.getCol(pos);
        return {name: pos, cat: posCol.categories, rawData: posCol.getRawData()};
      });

      const clustColCat = clusterCol?.categories ?? [];
      const clustColData = clusterCol?.getRawData() ?? new Int32Array(0);
      const customClust: {[key: string]: boolean[]} = {};
      for (const clust of this.customClusters)
        customClust[clust.name] = clust.toList();

      changeSelectionBitset(selection, positionList, clustColCat, clustColData, customClust);
    });

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

      this._filterChangedSubject.next();
    });
    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(isPeptideSpaceSource: boolean = false, fireFilterChanged: boolean = false): void {
    this.isUserChangedSelection = false;
    // this.isPeptideSpaceChangingBitset = isPeptideSpaceSource;
    this.df.selection.fireChanged();
    if (fireFilterChanged)
      this.df.filter.fireChanged();
    this.modifyOrCreateSplitCol();
    this.headerSelectedMonomers = calculateSelected(this.df);

    const acc = this.createAccordion();
    if (acc != null) {
      grok.shell.o = acc.root;
      for (const pane of acc.panes)
        pane.expanded = true;
    }
    this.isUserChangedSelection = true;
    // this.isPeptideSpaceChangingBitset = false;
  }

  postProcessGrids(): void {
    const posCols = this.splitSeqDf.columns.names();
    const sourceGrid = this.analysisView.grid;
    const sourceGridCols = sourceGrid.columns;
    const sourceGridColsLen = sourceGridCols.length;
    const visibleColumns = Object.keys(this.settings.columns || {});
    const sourceGridProps = sourceGrid.props;
    sourceGridProps.allowColSelection = false;
    sourceGridProps.allowEdit = false;
    sourceGridProps.allowRowResizing = false;
    sourceGridProps.showCurrentRowIndicator = false;
    this.df.temp[C.EMBEDDING_STATUS] = false;
    for (let colIdx = 1; colIdx < sourceGridColsLen; ++colIdx) {
      const gridCol = sourceGridCols.byIndex(colIdx)!;
      const tableColName = gridCol.column!.name;
      gridCol.visible = posCols.includes(tableColName) || (tableColName === C.COLUMNS_NAMES.ACTIVITY_SCALED) ||
        visibleColumns.includes(tableColName);
      gridCol.width = 60;
    }
  }

  closeDendrogram(): void {
    for (const node of this.analysisView.dockManager.rootNode.children) {
      if (node.container.containerElement.innerHTML.includes('Dendrogram')) {
        this.analysisView.dockManager.close(node);
        break;
      }
    }
    const viewer = wu(this.analysisView.viewers).find((v) => v.type === 'Dendrogram');
    viewer?.detach();
    viewer?.close();
  }

  async addDendrogram(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Calculating distance matrix...');
    try {
      const pepColValues: string[] = this.df.getCol(this.settings.sequenceColumnName!).toList();
      this._dm ??= new DistanceMatrix(await createDistanceMatrixWorker(pepColValues, StringMetricsNames.Levenshtein));
      const leafCol = this.df.col('~leaf-id') ?? this.df.columns.addNewString('~leaf-id').init((i) => i.toString());
      const treeNode = await this.treeHelper.hierarchicalClusteringByDistance(this._dm, 'ward');

      this.df.setTag(treeTAGS.NEWICK, this.treeHelper.toNewick(treeNode));
      const leafOrdering = this.treeHelper.getLeafList(treeNode).map((leaf) => parseInt(leaf.name));
      this.analysisView.grid.setRowOrder(leafOrdering);
      const dendrogramViewer = await this.df.plot.fromType('Dendrogram', {nodeColumnName: leafCol.name}) as DG.JsViewer;

      this.analysisView.dockManager.dock(dendrogramViewer, DG.DOCK_TYPE.LEFT, null, 'Dendrogram', 0.25);
    } catch (e) {
      _package.logger.error(e as string);
    } finally {
      pi.close();
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

    if (!this.isRibbonSet && this.df.getTag(C.MULTIPLE_VIEWS) != '1') {
      //TODO: don't pass model, pass parameters instead
      const settingsButton = ui.iconFA('wrench', () => getSettingsDialog(this), 'Peptides analysis settings');
      this.analysisView.setRibbonPanels([[settingsButton]], false);
      this.isRibbonSet = true;
      grok.events.onResetFilterRequest.subscribe(() => {
        this.isInvariantMapTrigger = true;
        this.initInvariantMapSelection(true);
        this.isInvariantMapTrigger = false;
      });
    }

    this.updateGrid();
    this.fireBitsetChanged(false, true);
    this.analysisView.grid.invalidate();
  }

  async addViewers(): Promise<void> {
    const dockManager = this.analysisView.dockManager;
    const dfPlt = this.df.plot;

    const mutationCliffsViewer = await dfPlt.fromType('peptide-sar-viewer') as MonomerPosition;
    const mostPotentResiduesViewer = await dfPlt.fromType('peptide-sar-viewer-vertical') as MostPotentResiduesViewer;
    if (this.settings.clustersColumnName)
      await this.addLogoSummaryTableViewer();

    const mcNode = dockManager.dock(mutationCliffsViewer, DG.DOCK_TYPE.DOWN, null, mutationCliffsViewer.name);

    dockManager.dock(mostPotentResiduesViewer, DG.DOCK_TYPE.RIGHT, mcNode, mostPotentResiduesViewer.name, 0.3);
  }

  async addLogoSummaryTableViewer(): Promise<void> {
    const logoSummary = await this.df.plot.fromType('logo-summary-viewer') as LogoSummary;
    this.analysisView.dockManager.dock(logoSummary, DG.DOCK_TYPE.RIGHT, null, 'Logo Summary Table');
  }

  addNewCluster(clusterName: string): void {
    const newClusterCol = DG.Column.fromBitSet(clusterName, this.getCompoundBitest());
    newClusterCol.setTag(C.TAGS.CUSTOM_CLUSTER, '1');
    this.df.columns.add(newClusterCol);
    this.analysisView.grid.col(newClusterCol.name)!.visible = false;
  }

  async createNewView(): Promise<void> {
    const rowMask = this.getCompoundBitest();
    if (!rowMask.anyTrue)
      return grok.shell.warning('Cannot create a new view, there are no visible selected rows in your dataset');

    const newDf = this.df.clone(rowMask);
    for (const [tag, value] of newDf.tags)
      newDf.setTag(tag, tag == C.TAGS.SETTINGS ? value : '');
    newDf.name = 'Peptides Multiple Views';
    newDf.setTag(C.MULTIPLE_VIEWS, '1');
    newDf.setTag(C.TAGS.UUID, uuid.v4());
    const view = grok.shell.addTableView(newDf);
    view.addViewer('logo-summary-viewer');
  }
}
