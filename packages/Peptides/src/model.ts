import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {pickUpPalette, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {calculateScores, SCORE} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

import wu from 'wu';
import * as rxjs from 'rxjs';
import * as uuid from 'uuid';
import $ from 'cash-dom';

import * as C from './utils/constants';
import * as type from './utils/types';
import {calculateSelected, extractColInfo, scaleActivity, getStatsSummary, prepareTableForHistogram} from './utils/misc';
import {MONOMER_POSITION_PROPERTIES, MonomerPosition, MostPotentResidues, SELECTION_MODE} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getActivityDistribution, getDistributionLegend, getDistributionWidget, getStatsTableMap} from './widgets/distribution';
import {getAggregatedValue, getStats, Stats} from './utils/statistics';
import {LogoSummaryTable} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {_package, getMonomerWorksInstance, getTreeHelperInstance} from './package';
import {findMutations} from './utils/algorithms';
import {createDistanceMatrixWorker} from './utils/worker-creator';
import {getSelectionWidget} from './widgets/selection';

export type SummaryStats = {
  minCount: number, maxCount: number,
  minMeanDifference: number, maxMeanDifference: number,
  minPValue: number, maxPValue: number,
  minRatio: number, maxRatio: number,
};
export type PositionStats = {[monomer: string]: Stats | undefined} & {general: SummaryStats};
export type MonomerPositionStats = {[position: string]: PositionStats | undefined} & {general: SummaryStats};
export type ClusterStats = {[cluster: string]: Stats};
export enum CLUSTER_TYPE {
  ORIGINAL = 'original',
  CUSTOM = 'custom',
};
export type ClusterType = `${CLUSTER_TYPE}`;
export type ClusterTypeStats = {[clusterType in ClusterType]: ClusterStats};
export enum VIEWER_TYPE {
  MONOMER_POSITION = 'Monomer-Position',
  MOST_POTENT_RESIDUES = 'Most Potent Residues',
  LOGO_SUMMARY_TABLE = 'Logo Summary Table',
  DENDROGRAM = 'Dendrogram',
};

export const getAggregatedColName = (aggF: string, colName: string): string => `${aggF}(${colName})`;

export class PeptidesModel {
  static modelName = 'peptidesModel';

  isBitsetChangedInitialized = false;
  isUserChangedSelection = true;

  df: DG.DataFrame;
  _monomerPositionStats?: MonomerPositionStats;
  _clusterStats?: ClusterTypeStats;
  _mutationCliffsSelection!: type.Selection;
  _invariantMapSelection!: type.Selection;
  _clusterSelection!: type.Selection;
  _mutationCliffs: type.MutationCliffs | null = null;
  isInitialized = false;
  _analysisView?: DG.TableView;

  _settings!: type.PeptidesSettings;
  isRibbonSet = false;

  _cp?: SeqPalette;
  headerSelectedMonomers: type.SelectionStats = {};
  webLogoBounds: {[positon: string]: {[monomer: string]: DG.Rect}} = {};
  cachedWebLogoTooltip: {bar: string, tooltip: HTMLDivElement | null} = {bar: '', tooltip: null};
  _alphabet?: string;
  _splitSeqDf?: DG.DataFrame;
  _dm!: DistanceMatrix;
  _layoutEventInitialized = false;

  subs: rxjs.Subscription[] = [];
  isHighlighting: boolean = false;
  latestSelectionItem: (type.SelectionItem & {kind: SELECTION_MODE | 'Cluster'}) | null = null;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    (dataFrame.temp[PeptidesModel.modelName] as PeptidesModel).init();
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  get id(): string {
    const id = this.df.getTag(C.TAGS.UUID);
    if (id === null || id === '')
      throw new Error('PeptidesError: UUID is not defined');

    return id;
  }

  get monomerPositionStats(): MonomerPositionStats {
    this._monomerPositionStats ??= this.calculateMonomerPositionStatistics();
    return this._monomerPositionStats;
  }

  set monomerPositionStats(mps: MonomerPositionStats) {
    this._monomerPositionStats = mps;
  }

  get splitSeqDf(): DG.DataFrame {
    this._splitSeqDf ??= this.buildSplitSeqDf();
    return this._splitSeqDf;
  }

  set splitSeqDf(df: DG.DataFrame) {
    this._splitSeqDf = df;
  }

  get alphabet(): string {
    const col = this.df.getCol(this.settings.sequenceColumnName!);
    return col.getTag(bioTAGS.alphabet);
  }

  get mutationCliffs(): type.MutationCliffs | null {
    return this._mutationCliffs!;
  }

  set mutationCliffs(si: type.MutationCliffs | null) {
    this._mutationCliffs = si;
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
    if (this._analysisView === undefined) {
      this._analysisView = wu(grok.shell.tableViews).find(({dataFrame}) => dataFrame?.getTag(C.TAGS.UUID) === this.id);
      if (this._analysisView === undefined) {
        this._analysisView = grok.shell.addTableView(this.df);
        const posCols = this.splitSeqDf.columns.names();

        for (let colIdx = 1; colIdx < this._analysisView.grid.columns.length; ++colIdx) {
          const gridCol = this._analysisView.grid.columns.byIndex(colIdx)!;
          gridCol.visible =
            posCols.includes(gridCol.column!.name) || (gridCol.column!.name === C.COLUMNS_NAMES.ACTIVITY_SCALED);
        }
      }
    }

    if (this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1' && !this._layoutEventInitialized)
      grok.shell.v = this._analysisView;

    this._analysisView.grid.invalidate();
    return this._analysisView;
  }

  get mutationCliffsSelection(): type.Selection {
    const tagSelection = this.df.getTag(C.TAGS.MUTATION_CLIFFS_SELECTION) ?? this.df.getTag(C.TAGS.SELECTION);
    if (tagSelection === null)
      this.initMutationCliffsSelection({notify: false});
    this._mutationCliffsSelection ??= JSON.parse(tagSelection ?? this.df.getTag(C.TAGS.MUTATION_CLIFFS_SELECTION) ?? this.df.getTag(C.TAGS.SELECTION)!);
    return this._mutationCliffsSelection;
  }

  set mutationCliffsSelection(selection: type.Selection) {
    this._mutationCliffsSelection = selection;
    // TODO: Remove in 1.14.0
    this.df.setTag(C.TAGS.SELECTION, JSON.stringify(selection));
    this.df.setTag(C.TAGS.MUTATION_CLIFFS_SELECTION, JSON.stringify(selection));
    this.fireBitsetChanged();

    const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    mpViewer?.viewerGrid.invalidate();
    const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    mprViewer?.viewerGrid.invalidate();

    this.analysisView.grid.invalidate();
  }

  get invariantMapSelection(): type.Selection {
    const tagSelection = this.df.getTag(C.TAGS.INVARIANT_MAP_SELECTION) ?? this.df.getTag(C.TAGS.FILTER);
    if (tagSelection === null)
      this.initInvariantMapSelection({notify: false});
    this._invariantMapSelection ??= JSON.parse(tagSelection ?? this.df.getTag(C.TAGS.INVARIANT_MAP_SELECTION) ?? this.df.getTag(C.TAGS.FILTER)!);
    return this._invariantMapSelection;
  }

  set invariantMapSelection(selection: type.Selection) {
    this._invariantMapSelection = selection;
    this.df.setTag(C.TAGS.INVARIANT_MAP_SELECTION, JSON.stringify(selection));
    // TODO: Remove in 1.14.0
    this.df.setTag(C.TAGS.FILTER, JSON.stringify(selection));
    this.fireBitsetChanged();
    this.analysisView.grid.invalidate();
  }

  get clusterSelection(): type.Selection {
    const tagSelection = this.df.getTag(C.TAGS.CLUSTER_SELECTION);
    if (tagSelection === null)
      this.initClusterSelection({notify: false});
    this._clusterSelection ??= JSON.parse(tagSelection ?? this.df.getTag(C.TAGS.CLUSTER_SELECTION)!);
    // TODO: Remove in 1.14.0
    if (Array.isArray(this._clusterSelection)) {
      const newSelection: type.Selection = {};
      newSelection[CLUSTER_TYPE.ORIGINAL] = [];
      newSelection[CLUSTER_TYPE.CUSTOM] = [];
      for (const cluster of this._clusterSelection) {
        if (wu(this.customClusters).some((col) => col.name === cluster))
          newSelection[CLUSTER_TYPE.CUSTOM].push(cluster);
        else
          newSelection[CLUSTER_TYPE.ORIGINAL].push(cluster);
      }
      this._clusterSelection = newSelection;
    }
    return this._clusterSelection;
  }

  set clusterSelection(selection: type.Selection) {
    this._clusterSelection = selection;
    this.df.tags[C.TAGS.CLUSTER_SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.analysisView.grid.invalidate();
  }

  get splitByPos(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[0];
    return splitByPosFlag === '1' ? true : false;
  }

  set splitByPos(flag: boolean) {
    const splitByMonomerFlag = (this.df.tags['distributionSplit'] || '00')[1];
    this.df.tags['distributionSplit'] = `${flag ? 1 : 0}${splitByMonomerFlag}`;
  }

  get splitByMonomer(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[1];
    return splitByPosFlag === '1' ? true : false;
  }

  set splitByMonomer(flag: boolean) {
    const splitByMonomerFlag = (this.df.tags['distributionSplit'] || '00')[0];
    this.df.tags['distributionSplit'] = `${splitByMonomerFlag}${flag ? 1 : 0}`;
  }

  get isMutationCliffsSelectionEmpty(): boolean {
    for (const monomerList of Object.values(this.mutationCliffsSelection)) {
      if (monomerList.length !== 0)
        return false;
    }
    return true;
  }

  get isInvariantMapSelectionEmpty(): boolean {
    for (const monomerList of Object.values(this.invariantMapSelection)) {
      if (monomerList.length !== 0)
        return false;
    }
    return true;
  }

  get isClusterSelectionEmpty(): boolean {
    return (this.clusterSelection[CLUSTER_TYPE.ORIGINAL].length + this.clusterSelection[CLUSTER_TYPE.CUSTOM].length) === 0;
  }

  get customClusters(): Iterable<DG.Column<boolean>> {
    const query: { [key: string]: string } = {};
    query[C.TAGS.CUSTOM_CLUSTER] = '1';
    return this.df.columns.byTags(query);
  }

  get settings(): type.PeptidesSettings {
    this._settings ??= JSON.parse(this.df.getTag('settings')!);
    return this._settings;
  }

  set settings(s: type.PeptidesSettings) {
    const newSettingsEntries = Object.entries(s);
    const updateVars: Set<string> = new Set();
    for (const [key, value] of newSettingsEntries) {
      this._settings[key as keyof type.PeptidesSettings] = value as any;
      switch (key) {
      case 'activityColumnName':
      case 'scaling':
        updateVars.add('activity');
        updateVars.add('mutationCliffs');
        updateVars.add('stats');
        break;
      case 'columns':
        updateVars.add('grid');
        break;
      case 'maxMutations':
      case 'minActivityDelta':
        updateVars.add('mutationCliffs');
        break;
      case 'showDendrogram':
        updateVars.add('dendrogram');
        break;
      case 'showLogoSummaryTable':
        updateVars.add('logoSummaryTable');
        break;
      case 'showMonomerPosition':
        updateVars.add('monomerPosition');
        break;
      case 'showMostPotentResidues':
        updateVars.add('mostPotentResidues');
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
        this.updateMutationCliffs();
        break;
      case 'stats':
        this.monomerPositionStats = this.calculateMonomerPositionStatistics();
        this.clusterStats = this.calculateClusterStatistics();
        const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition;
        mpViewer.createMonomerPositionGrid();
        mpViewer.render();
        const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues;
        mprViewer.createMostPotentResiduesGrid();
        mprViewer.render();
        break;
      case 'grid':
        this.setGridProperties();
        break;
      case 'dendrogram':
        this.settings.showDendrogram ? this.addDendrogram() : this.closeViewer(VIEWER_TYPE.DENDROGRAM);
        break;
      case 'logoSummaryTable':
        this.settings.showLogoSummaryTable ? this.addLogoSummaryTable() :
          this.closeViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE);
        break;
      case 'monomerPosition':
        this.settings.showMonomerPosition ? this.addMonomerPosition() :
          this.closeViewer(VIEWER_TYPE.MONOMER_POSITION);
        break;
      case 'mostPotentResidues':
        this.settings.showMostPotentResidues ? this.addMostPotentResidues() :
          this.closeViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES);
        break;
      }
    }

    //TODO: handle settings change
    const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    mpViewer?.createMonomerPositionGrid();
    mpViewer?.render();
    const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    mprViewer?.createMostPotentResiduesGrid();
    mprViewer?.render();
    const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    lstViewer?.createLogoSummaryTableGrid();
    lstViewer?.render();
  }

  get identityTemplate(): string {
    return this.df.getTag(C.TAGS.IDENTITY_TEMPLATE) ?? '';
  }

  set identityTemplate(template: string) {
    this.df.setTag(C.TAGS.IDENTITY_TEMPLATE, template);
  }

  updateMutationCliffs(notify: boolean = true): void {
    const scaledActivityCol: DG.Column<number> = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    //TODO: set categories ordering the same to share compare indexes instead of strings
    const monomerCols: type.RawColumn[] = this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER).map(extractColInfo);
    const targetCol = typeof this.settings.targetColumnName !== 'undefined' ?
      extractColInfo(this.df.getCol(this.settings.targetColumnName)) : null;
    const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    const currentTarget = mpViewer?.getProperty(MONOMER_POSITION_PROPERTIES.TARGET)?.get(mpViewer);
    const targetOptions = {targetCol: targetCol, currentTarget: currentTarget};
    const mutationCliffs = findMutations(scaledActivityCol.getRawData(), monomerCols, this.settings, targetOptions);
    if (notify)
      this.mutationCliffs = mutationCliffs;
    else
      this._mutationCliffs = mutationCliffs;
  }

  buildSplitSeqDf(): DG.DataFrame {
    const sequenceCol = this.df.getCol(this.settings.sequenceColumnName!);
    const splitSeqDf = splitAlignedSequences(sequenceCol);

    return splitSeqDf;
  }

  getCompoundBitset(): DG.BitSet {
    return this.df.selection.clone().and(this.df.filter);
  }

  createAccordion(): DG.Accordion | null {
    const trueModel: PeptidesModel | undefined = grok.shell.t.temp[PeptidesModel.modelName];
    if (!trueModel)
      return null;

    const acc = ui.accordion();
    acc.root.style.width = '100%';
    const filterAndSelectionBs = trueModel.getCompoundBitset();
    const filteredTitlePart = trueModel.df.filter.anyFalse ? ` among ${trueModel.df.filter.trueCount} filtered` : '';
    acc.addTitle(ui.h1(`${filterAndSelectionBs.trueCount} selected rows${filteredTitlePart}`));
    if (filterAndSelectionBs.anyTrue) {
      acc.addPane('Actions', () => {
        const newView = ui.label('New view');
        $(newView).addClass('d4-link-action');
        newView.onclick = (): string => trueModel.createNewView();
        newView.onmouseover =
          (ev): void => ui.tooltip.show('Creates a new view from current selection', ev.clientX + 5, ev.clientY + 5);
        const newCluster = ui.label('New cluster');
        $(newCluster).addClass('d4-link-action');
        newCluster.onclick = (): void => {
          const lstViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
          if (lstViewer === null)
            throw new Error('Logo summary table viewer is not found');
          lstViewer.clusterFromSelection();
        };
        newCluster.onmouseover =
          (ev): void => ui.tooltip.show('Creates a new cluster from selection', ev.clientX + 5, ev.clientY + 5);
        const removeCluster = ui.label('Remove cluster');
        $(removeCluster).addClass('d4-link-action');
        removeCluster.onclick = (): void => {
          const lstViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
          if (lstViewer === null)
            throw new Error('Logo summary table viewer is not found');
          lstViewer.removeCluster();
        };
        removeCluster.onmouseover =
          (ev): void => ui.tooltip.show('Removes currently selected custom cluster', ev.clientX + 5, ev.clientY + 5);
        removeCluster.style.visibility = trueModel.clusterSelection[CLUSTER_TYPE.CUSTOM].length === 0 ? 'hidden' : 'visible';
        return ui.divV([newView, newCluster, removeCluster]);
      });
    }
    const table = trueModel.df.filter.anyFalse ? trueModel.df.clone(trueModel.df.filter, null, true) : trueModel.df;
    acc.addPane('Mutation Cliffs pairs', () => mutationCliffsWidget(trueModel.df, trueModel).root);
    acc.addPane('Distribution', () => getDistributionWidget(table, trueModel).root);
    acc.addPane('Selection', () => getSelectionWidget(trueModel.df, trueModel).root);

    return acc;
  }

  updateGrid(): void {
    this.joinDataFrames();

    this.createScaledCol();

    this.initInvariantMapSelection({notify: false});
    this.initMutationCliffsSelection({notify: false});
    this.initClusterSelection({notify: false});

    this.setWebLogoInteraction();
    this.webLogoBounds = {};

    this.setCellRenderers();

    this.setTooltips();

    this.setBitsetCallback();

    this.setGridProperties();
  }

  initInvariantMapSelection(options: {notify?: boolean} = {}): void {
    options.notify ??= true;

    const tempSelection: type.Selection = {};
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns)
      tempSelection[pos] = [];

    if (options.notify)
      this.invariantMapSelection = tempSelection;
    else
      this._invariantMapSelection = tempSelection;
  }

  initMutationCliffsSelection(options: {notify?: boolean} = {}): void {
    options.notify ??= true;

    const tempSelection: type.Selection = {};
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns)
      tempSelection[pos] = [];

    if (options.notify)
      this.mutationCliffsSelection = tempSelection;
    else
      this._mutationCliffsSelection = tempSelection;
  }

  joinDataFrames(): void {
    // append splitSeqDf columns to source table and make sure columns are not added more than once
    const name = this.df.name;
    const cols = this.df.columns;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const colName of positionColumns) {
      let col = this.df.col(colName);
      const newCol = this.splitSeqDf.getCol(colName);
      if (col !== null)
        cols.remove(colName);

      const newColCat = newCol.categories;
      const newColData = newCol.getRawData();
      col = cols.addNew(newCol.name, newCol.type).init((i) => newColCat[newColData[i]]);
      CR.setMonomerRenderer(col, this.alphabet);
    }
    this.df.name = name;
  }

  createScaledCol(): void {
    const sourceGrid = this.analysisView.grid;
    const scaledCol = scaleActivity(this.df.getCol(this.settings.activityColumnName!), this.settings.scaling);
    //TODO: make another func
    this.df.columns.replace(C.COLUMNS_NAMES.ACTIVITY_SCALED, scaledCol);

    sourceGrid.columns.setOrder([scaledCol.name]);
  }

  calculateMonomerPositionStatistics(options: {isFiltered?: boolean, columns?: string[]} = {}): MonomerPositionStats {
    options.isFiltered ??= false;
    const monomerPositionObject = {general: {}} as MonomerPositionStats & { general: SummaryStats };
    let activityColData: Float64Array = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData() as Float64Array;
    let positionColumns = this.splitSeqDf.columns.toList();
    let sourceDfLen = this.df.rowCount;

    if (options.isFiltered) {
      sourceDfLen = this.df.filter.trueCount;
      const tempActivityData = new Float64Array(sourceDfLen);
      const selectedIndexes = this.df.filter.getSelectedIndexes();
      for (let i = 0; i < sourceDfLen; ++i)
        tempActivityData[i] = activityColData[selectedIndexes[i]];
      activityColData = tempActivityData;
      positionColumns = this.splitSeqDf.clone(this.df.filter).columns.toList();
    }
    options.columns ??= positionColumns.map((col) => col.name);

    for (const posCol of positionColumns) {
      if (!options.columns.includes(posCol.name))
        continue;
      const posColData = posCol.getRawData();
      const posColCateogries = posCol.categories;
      const currentPositionObject = {general: {}} as PositionStats & {general: SummaryStats};

      for (let categoryIndex = 0; categoryIndex < posColCateogries.length; ++categoryIndex) {
        const monomer = posColCateogries[categoryIndex];
        if (monomer === '')
          continue;

        const boolArray: boolean[] = new Array(sourceDfLen).fill(false);
        for (let i = 0; i < sourceDfLen; ++i) {
          if (posColData[i] === categoryIndex)
            boolArray[i] = true;
        }
        const bitArray = BitArray.fromValues(boolArray);
        const stats = bitArray.allFalse || bitArray.allTrue ?
          {count: sourceDfLen, meanDifference: 0, ratio: 1.0, pValue: null, mask: bitArray} :
          getStats(activityColData, bitArray);
        currentPositionObject[monomer] = stats;
        this.getSummaryStats(currentPositionObject.general, stats);
      }
      monomerPositionObject[posCol.name] = currentPositionObject;
      this.getSummaryStats(monomerPositionObject.general, null, currentPositionObject.general);
    }
    return monomerPositionObject;
  }

  getSummaryStats(genObj: SummaryStats, stats: Stats | null = null, summaryStats: SummaryStats | null = null): void {
    if (stats === null && summaryStats === null)
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

    if (!isNaN(stats?.pValue ?? NaN)) {
      const possibleMaxPValue = stats?.pValue ?? summaryStats!.maxPValue;
      genObj.maxPValue ??= possibleMaxPValue;
      if (genObj.maxPValue < possibleMaxPValue)
        genObj.maxPValue = possibleMaxPValue;

      const possibleMinPValue = stats?.pValue ?? summaryStats!.minPValue;
      genObj.minPValue ??= possibleMinPValue;
      if (genObj.minPValue > possibleMinPValue)
        genObj.minPValue = possibleMinPValue;
    }

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
    const rowCount = this.df.rowCount;
    const origClustCol = this.df.getCol(this.settings.clustersColumnName!);
    const origClustColData = origClustCol.getRawData();
    const origClustColCat = origClustCol.categories;
    const origClustMasks: BitArray[] = Array.from({length: origClustColCat.length},
      () => new BitArray(rowCount, false));
    for (let rowIdx = 0; rowIdx < rowCount; ++rowIdx)
      origClustMasks[origClustColData[rowIdx]].setTrue(rowIdx);


    const customClustColList = wu(this.customClusters).toArray();
    const customClustMasks = customClustColList.map(
      (v) => BitArray.fromUint32Array(rowCount, v.getRawData() as Uint32Array));
    const customClustColNamesList = customClustColList.map((v) => v.name);

    const activityColData: type.RawData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();

    const origClustStats: ClusterStats = {};
    const customClustStats: ClusterStats = {};

    for (const clustType of Object.values(CLUSTER_TYPE)) {
      const masks = clustType === CLUSTER_TYPE.ORIGINAL ? origClustMasks : customClustMasks;
      const clustNames = clustType === CLUSTER_TYPE.ORIGINAL ? origClustColCat : customClustColNamesList;
      const resultStats = clustType === CLUSTER_TYPE.ORIGINAL ? origClustStats : customClustStats;
      for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
        const mask = masks[maskIdx];
        const stats = mask.allTrue || mask.allFalse ? {count: mask.length, meanDifference: 0, ratio: 1.0, pValue: null, mask: mask} :
          getStats(activityColData, mask);
        resultStats[clustNames[maskIdx]] = stats;
      }
    }

    const resultStats = {} as ClusterTypeStats;
    resultStats[CLUSTER_TYPE.ORIGINAL] = origClustStats;
    resultStats[CLUSTER_TYPE.CUSTOM] = customClustStats;
    return resultStats;
  }

  initClusterSelection(options: {notify?: boolean} = {}): void {
    options.notify ??= true;

    const newClusterSelection = {} as type.Selection;
    newClusterSelection[CLUSTER_TYPE.ORIGINAL] = [];
    newClusterSelection[CLUSTER_TYPE.CUSTOM] = [];
    if (options.notify)
      this.clusterSelection = newClusterSelection;
    else
      this._clusterSelection = newClusterSelection;
  }

  setWebLogoInteraction(): void {
    const sourceView = this.analysisView.grid;
    const eventAction = (ev: MouseEvent): void => {
      const cell = sourceView.hitTest(ev.offsetX, ev.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType === C.SEM_TYPES.MONOMER) {
        const monomerPosition = this.findWebLogoMonomerPosition(cell, ev);
        if (monomerPosition === null) {
          this.unhighlight();
          return;
        }
        this.requestBarchartAction(ev, monomerPosition);
        this.highlightMonomerPosition(monomerPosition);
      }
    };

    // The following events makes the barchart interactive
    rxjs.fromEvent<MouseEvent>(sourceView.overlay, 'mousemove')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
    rxjs.fromEvent<MouseEvent>(sourceView.overlay, 'click')
      .subscribe((mouseMove: MouseEvent) => eventAction(mouseMove));
  }

  highlightMonomerPosition(monomerPosition: type.SelectionItem): void {
    const bitArray = new BitArray(this.df.rowCount);
    if (monomerPosition.positionOrClusterType === C.COLUMNS_NAMES.MONOMER) {
      const positionStats = Object.values(this.monomerPositionStats);
      for (const posStat of positionStats) {
        const monomerPositionStats = (posStat as PositionStats)[monomerPosition.monomerOrCluster];
        if (monomerPositionStats ?? false)
          bitArray.or(monomerPositionStats!.mask);
      }
    } else {
      const monomerPositionStats = this.monomerPositionStats[monomerPosition.positionOrClusterType]![monomerPosition.monomerOrCluster];
      if (monomerPositionStats ?? false)
        bitArray.or(monomerPositionStats!.mask);
    }

    this.df.rows.highlight((i) => bitArray.getBit(i));
    this.isHighlighting = true;
  }

  highlightCluster(cluster: type.SelectionItem): void {
    const bitArray = this.clusterStats[cluster.positionOrClusterType as ClusterType][cluster.monomerOrCluster].mask;
    this.df.rows.highlight((i) => bitArray.getBit(i));
    this.isHighlighting = true;
  }

  unhighlight(): void {
    if (!this.isHighlighting)
      return;
    this.df.rows.highlight(null);
    this.isHighlighting = false;
  }

  findWebLogoMonomerPosition(cell: DG.GridCell, ev: MouseEvent): type.SelectionItem | null {
    const barCoords = this.webLogoBounds[cell.tableColumn!.name];
    for (const [monomer, coords] of Object.entries(barCoords)) {
      const isIntersectingX = ev.offsetX >= coords.x && ev.offsetX <= coords.x + coords.width;
      const isIntersectingY = ev.offsetY >= coords.y && ev.offsetY <= coords.y + coords.height;
      if (isIntersectingX && isIntersectingY)
        return {monomerOrCluster: monomer, positionOrClusterType: cell.tableColumn!.name};
    }

    return null;
  }

  requestBarchartAction(ev: MouseEvent, monomerPosition: type.SelectionItem): void {
    if (ev.type === 'click')
      this.modifyInvariantMapSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
    else {
      const bar = `${monomerPosition.positionOrClusterType} = ${monomerPosition.monomerOrCluster}`;
      if (this.cachedWebLogoTooltip.bar === bar)
        ui.tooltip.show(this.cachedWebLogoTooltip.tooltip!, ev.clientX, ev.clientY);
      else
        this.cachedWebLogoTooltip = {bar: bar, tooltip: this.showTooltipAt(monomerPosition, ev.clientX, ev.clientY)};
    }
  }

  setCellRenderers(): void {
    const sourceGrid = this.analysisView.grid;
    sourceGrid.setOptions({'colHeaderHeight': 130});
    const headerRenderer = (gcArgs: DG.GridCellRenderArgs): void => {
      const ctx = gcArgs.g;
      const bounds = gcArgs.bounds;
      const col = gcArgs.cell.tableColumn;

      ctx.save();
      try {
        ctx.beginPath();
        ctx.rect(bounds.x, bounds.y, bounds.width, bounds.height);
        ctx.clip();

        //TODO: optimize
        if (gcArgs.cell.isColHeader && col?.semType === C.SEM_TYPES.MONOMER) {
          const isDfFiltered = this.df.filter.anyFalse;
          const stats = (isDfFiltered ? this.calculateMonomerPositionStatistics({isFiltered: true, columns: [col.name]}) :
            this.monomerPositionStats)[col.name];
          if (!stats)
            return;
          //TODO: precalc on stats creation
          const sortedStatsOrder = Object.keys(stats).sort((a, b) => {
            if (a === '' || a === '-')
              return +1;
            else if (b === '' || b === '-')
              return -1;
            return 0;
          }).filter((v) => v !== 'general');

          this.webLogoBounds[col.name] = CR.drawLogoInBounds(ctx, bounds, stats, col.name, sortedStatsOrder,
            this.df.filter.trueCount, this.cp, this.headerSelectedMonomers[col.name]);
          gcArgs.preventDefault();
        }
      } catch (e) {
        console.warn(`PeptidesHeaderLogoError: couldn't render WebLogo for column \`${col!.name}\`. ` +
          `See original error below.`);
        console.warn(e);
      } finally {
        ctx.restore();
      }
    };
    sourceGrid.onCellRender.subscribe((gcArgs) => headerRenderer(gcArgs));

    if (!this._layoutEventInitialized) {
      grok.events.onViewLayoutApplied.subscribe((layout) => {
        if (layout.view.id === this.analysisView.id)
          this.updateGrid();
      });
      this._layoutEventInitialized = true;
    }
  }

  setTooltips(): void {
    this.analysisView.grid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader && cell.tableColumn!.semType === C.SEM_TYPES.MONOMER)
        return true;
      if (!(cell.isTableCell && cell.tableColumn!.semType === C.SEM_TYPES.MONOMER))
        return false;

      this.showMonomerTooltip(cell.cell.value, x, y);
      return true;
    });
  }

  showMonomerTooltip(monomer: string, x: number, y: number): boolean {
    const tooltipElements: HTMLDivElement[] = [];
    const monomerName = monomer.toLowerCase();

    const mw = getMonomerWorksInstance();
    const mol = mw?.getCappedRotatedMonomer('PEPTIDE', monomer);

    if (mol) {
      tooltipElements.push(ui.div(monomerName));
      const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      tooltipElements.push(grok.chem.svgMol(mol, undefined, undefined, options));
    } else if (monomer !== '')
      tooltipElements.push(ui.div(monomer));
    else
      return true;


    ui.tooltip.show(ui.divV(tooltipElements), x, y);

    return true;
  }

  showTooltip(monomerPosition: type.SelectionItem, x: number, y: number, fromViewer: boolean = false): boolean {
    if (monomerPosition.positionOrClusterType === C.COLUMNS_NAMES.MONOMER)
      this.showMonomerTooltip(monomerPosition.monomerOrCluster, x, y);
    else
      this.showTooltipAt(monomerPosition, x, y, fromViewer);
    return true;
  }
  //TODO: move out to viewer code
  showTooltipAt(monomerPosition: type.SelectionItem, x: number, y: number, fromViewer: boolean = false): HTMLDivElement | null {
    const stats = this.monomerPositionStats[monomerPosition.positionOrClusterType]![monomerPosition.monomerOrCluster];
    if (!stats?.count)
      return null;

    const activityCol = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    const mask = DG.BitSet.fromBytes(stats.mask.buffer.buffer, activityCol.length);
    const distributionTable = DG.DataFrame.fromColumns(
      [activityCol, DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask)]);
    const hist = getActivityDistribution(prepareTableForHistogram(distributionTable), true);

    const tableMap = getStatsTableMap(stats);
    if (fromViewer) {
      tableMap['Mean difference'] = `${tableMap['Mean difference']} (size)`;
      if (tableMap['p-value'])
        tableMap['p-value'] = `${tableMap['p-value']} (color)`;
    }
    const aggregatedColMap = this.getAggregatedColumnValues({mask: mask});
    const resultMap = {...tableMap, ...aggregatedColMap};

    const labels = getDistributionLegend(`${monomerPosition.positionOrClusterType} : ${monomerPosition.monomerOrCluster}`, 'Other');
    const distroStatsElem = getStatsSummary(labels, hist, resultMap);

    ui.tooltip.show(distroStatsElem, x, y);

    return distroStatsElem;
  }

  getAggregatedColumnValues(options: {filterDf?: boolean, mask?: DG.BitSet, fractionDigits?: number} = {},
  ): StringDictionary {
    options.filterDf ??= false;
    options.fractionDigits ??= 3;

    const filteredDf = options.filterDf && this.df.filter.anyFalse ? this.df.clone(this.df.filter) : this.df;

    const colResults: StringDictionary = {};
    for (const [colName, aggFn] of Object.entries(this.settings.columns!)) {
      const newColName = getAggregatedColName(aggFn, colName);
      const value = getAggregatedValue(filteredDf.getCol(colName), aggFn, options.mask);
      colResults[newColName] = value.toFixed(options.fractionDigits);
    }
    return colResults;
  }

  setBitsetCallback(): void {
    if (this.isBitsetChangedInitialized)
      return;
    const selection = this.df.selection;
    const filter = this.df.filter;

    const getCombinedSelection = (): DG.BitSet => {
      const combinedSelection = new BitArray(this.df.rowCount, false);
      // Invariant map selection
      for (const [position, monomerList] of Object.entries(this.invariantMapSelection)) {
        for (const monomer of monomerList) {
          const monomerPositionStats = this.monomerPositionStats[position]![monomer]!;
          combinedSelection.or(monomerPositionStats.mask);
        }
      }

      // Mutation cliffs selection
      for (const [position, monomerList] of Object.entries(this.mutationCliffsSelection)) {
        for (const monomer of monomerList) {
          const substitutions = this.mutationCliffs?.get(monomer)?.get(position) ?? null;
          if (substitutions === null)
            continue;
          for (const [key, value] of substitutions.entries()) {
            combinedSelection.setTrue(key);
            for (const v of value)
              combinedSelection.setTrue(v);
          }
        }
      }

      // Cluster selection
      for (const clustType of Object.keys(this.clusterSelection)) {
        for (const clust of this.clusterSelection[clustType]) {
          const clusterStats = this.clusterStats[clustType as CLUSTER_TYPE][clust]!;
          combinedSelection.or(clusterStats.mask);
        }
      }

      return DG.BitSet.fromBytes(combinedSelection.buffer.buffer, combinedSelection.length);
    };

    const getLatestSelection = (): DG.BitSet => {
      if (this.latestSelectionItem === null)
        return getCombinedSelection();
      if (this.latestSelectionItem.kind === SELECTION_MODE.INVARIANT_MAP) {
        const monomerPositionStats = this.monomerPositionStats[this.latestSelectionItem.positionOrClusterType]![this.latestSelectionItem.monomerOrCluster]!;
        return DG.BitSet.fromBytes(monomerPositionStats.mask.buffer.buffer, monomerPositionStats.mask.length);
      } else if (this.latestSelectionItem.kind === SELECTION_MODE.MUTATION_CLIFFS) {
        const substitutions = this.mutationCliffs?.get(this.latestSelectionItem.monomerOrCluster)?.get(this.latestSelectionItem.positionOrClusterType) ?? null;
        if (substitutions === null)
          throw new Error(`Couldn't find substitutions for ${this.latestSelectionItem.monomerOrCluster} at ${this.latestSelectionItem.positionOrClusterType}`);
        const latestSelection = new BitArray(this.df.rowCount, false);
        for (const [key, value] of substitutions.entries()) {
          latestSelection.setTrue(key);
          for (const v of value)
            latestSelection.setTrue(v);
        }
        return DG.BitSet.fromBytes(latestSelection.buffer.buffer, latestSelection.length);
      } else if (this.latestSelectionItem.kind === 'Cluster') {
        const clusterStats = this.clusterStats[this.latestSelectionItem.positionOrClusterType as CLUSTER_TYPE][this.latestSelectionItem.monomerOrCluster]!;
        return DG.BitSet.fromBytes(clusterStats.mask.buffer.buffer, clusterStats.mask.length);
      }
      throw new Error(`Unknown selection kind: ${this.latestSelectionItem.kind}`);
    };

    const showAccordion = (): void => {
      const acc = this.createAccordion();
      if (acc === null)
        return;
      grok.shell.o = acc.root;
      for (const pane of acc.panes)
        pane.expanded = true;
    };

    DG.debounce(selection.onChanged, 500).subscribe(() => {
      if (!this.isUserChangedSelection)
        selection.copyFrom(getLatestSelection(), false);
      showAccordion();
      this.isUserChangedSelection = true;
    });

    DG.debounce(filter.onChanged, 500).subscribe(() => {
      const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
      if (lstViewer !== null && typeof lstViewer.model !== 'undefined') {
        lstViewer.createLogoSummaryTableGrid();
        lstViewer.render();
      }
      showAccordion();
    });

    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(fireFilterChanged: boolean = false): void {
    this.isUserChangedSelection = false;
    this.df.selection.fireChanged();
    if (fireFilterChanged)
      this.df.filter.fireChanged();
    this.headerSelectedMonomers = calculateSelected(this.df);
  }

  setGridProperties(props?: DG.IGridLookSettings): void {
    const sourceGrid = this.analysisView.grid;
    const sourceGridProps = sourceGrid.props;
    sourceGridProps.allowColSelection = props?.allowColSelection ?? false;
    sourceGridProps.allowEdit = props?.allowEdit ?? false;
    sourceGridProps.showCurrentRowIndicator = props?.showCurrentRowIndicator ?? false;
    this.df.temp[C.EMBEDDING_STATUS] = false;
  }

  closeViewer(viewerType: VIEWER_TYPE): void {
    const viewer = this.findViewer(viewerType);
    viewer?.detach();
    viewer?.close();
  }

  findViewerNode(viewerType: VIEWER_TYPE): DG.DockNode | null {
    for (const node of this.analysisView.dockManager.rootNode.children) {
      if (node.container.containerElement.innerHTML.includes(viewerType))
        return node;
    }
    return null;
  }

  async addDendrogram(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Calculating distance matrix...');
    try {
      const pepColValues: string[] = this.df.getCol(this.settings.sequenceColumnName!).toList();
      this._dm ??= new DistanceMatrix(await createDistanceMatrixWorker(pepColValues, StringMetricsNames.Levenshtein));
      const leafCol = this.df.col('~leaf-id') ?? this.df.columns.addNewString('~leaf-id').init((i) => i.toString());
      const treeHelper: ITreeHelper = getTreeHelperInstance()!;
      const treeNode = await treeHelper.hierarchicalClusteringByDistance(this._dm, 'ward');

      this.df.setTag(treeTAGS.NEWICK, treeHelper.toNewick(treeNode));
      const leafOrdering = treeHelper.getLeafList(treeNode).map((leaf) => parseInt(leaf.name));
      this.analysisView.grid.setRowOrder(leafOrdering);
      const dendrogramViewer = await this.df.plot.fromType('Dendrogram', {nodeColumnName: leafCol.name}) as DG.JsViewer;

      this.analysisView.dockManager.dock(dendrogramViewer, DG.DOCK_TYPE.LEFT, null, 'Dendrogram', 0.25);
    } catch (e) {
      _package.logger.error(e as string);
    } finally {
      pi.close();
    }
  }

  /** Class initializer */
  init(): void {
    if (this.isInitialized)
      return;
    this.isInitialized = true;

    if (!this.isRibbonSet && this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1') {
      //TODO: don't pass model, pass parameters instead
      const settingsButton = ui.iconFA('wrench', () => getSettingsDialog(this), 'Peptides analysis settings');
      this.analysisView.setRibbonPanels([[settingsButton]], false);
      this.isRibbonSet = true;
      this.updateGrid();
    }

    this.subs.push(grok.events.onAccordionConstructed.subscribe((acc) => {
      if (!(grok.shell.o instanceof DG.SemanticValue || (grok.shell.o instanceof DG.Column && this.df.columns.toList().includes(grok.shell.o))))
        return;

      const actionsPane = acc.getPane('Actions');

      const actionsHost = $(actionsPane.root).find('.d4-flex-col');
      const calculateIdentity = ui.label('Calculate identity');
      calculateIdentity.classList.add('d4-link-action');
      ui.tooltip.bind(calculateIdentity, 'Adds a column with fractions of matching monomers against sequence in the current row');
      calculateIdentity.onclick = (): void => {
        const seqCol = this.df.getCol(this.settings.sequenceColumnName!);
        calculateScores(this.df, seqCol, seqCol.get(this.df.currentRowIdx), SCORE.IDENTITY);
      };
      actionsHost.append(ui.span([calculateIdentity], 'd4-markdown-row'));

      const calculateSimilarity = ui.label('Calculate similarity');
      calculateSimilarity.classList.add('d4-link-action');
      ui.tooltip.bind(calculateSimilarity, 'Adds a column with sequence similarity scores against sequence in the current row');
      calculateSimilarity.onclick = (): void => {
        const seqCol = this.df.getCol(this.settings.sequenceColumnName!);
        calculateScores(this.df, seqCol, seqCol.get(this.df.currentRowIdx), SCORE.SIMILARITY);
      };
      actionsHost.append(ui.span([calculateSimilarity], 'd4-markdown-row'));
    }));

    this.subs.push(grok.events.onViewRemoved.subscribe((view) => {
      if (view.id === this.analysisView.id)
        this.subs.forEach((v) => v.unsubscribe());
      grok.log.debug(`Peptides: view ${view.name} removed`);
    }));
    this.subs.push(grok.events.onTableRemoved.subscribe((table: DG.DataFrame) => {
      if (table.id === this.df.id)
        this.subs.forEach((v) => v.unsubscribe());
      grok.log.debug(`Peptides: table ${table.name} removed`);
    }));
    this.subs.push(grok.events.onProjectClosed.subscribe((project: DG.Project) => {
      if (project.id === grok.shell.project.id)
        this.subs.forEach((v) => v.unsubscribe());
      grok.log.debug(`Peptides: project ${project.name} closed`);
    }));

    this.fireBitsetChanged(true);
    if (typeof this.settings.targetColumnName === 'undefined')
      this.updateMutationCliffs();
    this.analysisView.grid.invalidate();
  }

  findViewer(viewerType: VIEWER_TYPE): DG.Viewer | null {
    return wu(this.analysisView.viewers).find((v) => v.type === viewerType) || null;
  }

  async addLogoSummaryTable(): Promise<void> {
    this.closeViewer(VIEWER_TYPE.MONOMER_POSITION);
    this.closeViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES);
    const logoSummaryTable = await this.df.plot.fromType(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;
    this.analysisView.dockManager.dock(logoSummaryTable, DG.DOCK_TYPE.RIGHT, null, VIEWER_TYPE.LOGO_SUMMARY_TABLE);
    if (this.settings.showMonomerPosition)
      await this.addMonomerPosition();
    if (this.settings.showMostPotentResidues)
      await this.addMostPotentResidues();
    logoSummaryTable.viewerGrid.invalidate();
  }

  async addMonomerPosition(): Promise<void> {
    const monomerPosition = await this.df.plot.fromType(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition;
    const mostPotentResidues = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    const dm = this.analysisView.dockManager;
    const [dockType, refNode, ratio] = mostPotentResidues === null ? [DG.DOCK_TYPE.DOWN, null, undefined] :
      [DG.DOCK_TYPE.LEFT, this.findViewerNode(VIEWER_TYPE.MOST_POTENT_RESIDUES), 0.7];
    dm.dock(monomerPosition, dockType, refNode, VIEWER_TYPE.MONOMER_POSITION, ratio);
    if (typeof this.settings.targetColumnName !== 'undefined') {
      const target = monomerPosition.getProperty(MONOMER_POSITION_PROPERTIES.TARGET)!;
      const choices = this.df.getCol(this.settings.targetColumnName!).categories;
      target.choices = choices;
      target.set(monomerPosition, choices[0]);
    }
  }

  async addMostPotentResidues(): Promise<void> {
    const mostPotentResidues =
      await this.df.plot.fromType(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues;
    const monomerPosition = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    const dm = this.analysisView.dockManager;
    const [dockType, refNode, ratio] = monomerPosition === null ? [DG.DOCK_TYPE.DOWN, null, undefined] :
      [DG.DOCK_TYPE.RIGHT, this.findViewerNode(VIEWER_TYPE.MONOMER_POSITION), 0.3];
    dm.dock(mostPotentResidues, dockType, refNode, VIEWER_TYPE.MOST_POTENT_RESIDUES, ratio);
  }

  addNewCluster(clusterName: string): void {
    const newClusterCol = DG.Column.fromBitSet(clusterName, this.getCompoundBitset());
    newClusterCol.setTag(C.TAGS.CUSTOM_CLUSTER, '1');
    this.df.columns.add(newClusterCol);
    this.analysisView.grid.col(newClusterCol.name)!.visible = false;
  }

  createNewView(): string {
    const rowMask = this.getCompoundBitset();
    const newDfId = uuid.v4();

    const newDf = this.df.clone(rowMask);
    for (const [tag, value] of newDf.tags)
      newDf.setTag(tag, tag === C.TAGS.SETTINGS ? value : '');

    newDf.name = 'Peptides Multiple Views';
    newDf.setTag(C.TAGS.MULTIPLE_VIEWS, '1');
    newDf.setTag(C.TAGS.UUID, newDfId);

    const view = grok.shell.addTableView(newDf);
    view.addViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE);

    return newDfId;
  }

  modifyInvariantMapSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {shiftPressed: false, ctrlPressed: false}, notify: boolean = true): void {
    if (notify)
      this.invariantMapSelection = this.modifySelection(this.invariantMapSelection, monomerPosition, options);
    else
      this._invariantMapSelection = this.modifySelection(this._invariantMapSelection, monomerPosition, options);
  }

  modifyMutationCliffsSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {shiftPressed: false, ctrlPressed: false}, notify: boolean = true): void {
    if (notify)
      this.mutationCliffsSelection = this.modifySelection(this.mutationCliffsSelection, monomerPosition, options);
    else
      this._mutationCliffsSelection = this.modifySelection(this._mutationCliffsSelection, monomerPosition, options);
  }

  modifyClusterSelection(cluster: type.SelectionItem, options: type.SelectionOptions = {shiftPressed: false, ctrlPressed: false}, notify: boolean = true): void {
    if (notify)
      this.clusterSelection = this.modifySelection(this.clusterSelection, cluster, options);
    else
      this._clusterSelection = this.modifySelection(this._clusterSelection, cluster, options);
  }

  modifySelection(selection: type.Selection, clusterOrMonomerPosition: type.SelectionItem, options: type.SelectionOptions): type.Selection {
    const monomerList = selection[clusterOrMonomerPosition.positionOrClusterType];
    const monomerIndex = monomerList.indexOf(clusterOrMonomerPosition.monomerOrCluster);
    if (options.shiftPressed && options.ctrlPressed) {
      if (monomerIndex !== -1)
        monomerList.splice(monomerIndex, 1);
    } else if (options.ctrlPressed) {
      if (monomerIndex === -1)
        monomerList.push(clusterOrMonomerPosition.monomerOrCluster);
      else
        monomerList.splice(monomerIndex, 1);
    } else if (options.shiftPressed) {
      if (monomerIndex === -1)
        monomerList.push(clusterOrMonomerPosition.monomerOrCluster);
    } else {
      const selectionKeys = Object.keys(selection);
      selection = {};
      for (const posOrClustType of selectionKeys) {
        selection[posOrClustType] = [];
        if (posOrClustType === clusterOrMonomerPosition.positionOrClusterType)
          selection[posOrClustType].push(clusterOrMonomerPosition.monomerOrCluster);
      }
    }
    return selection;
  }
}
