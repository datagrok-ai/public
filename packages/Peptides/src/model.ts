import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {MonomerWorks} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {pickUpPalette, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
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
import {calculateSelected, extractColInfo, scaleActivity, getStatsSummary, prepareTableForHistogram, getTemplate} from './utils/misc';
import {MONOMER_POSITION_PROPERTIES, MonomerPosition, MostPotentResidues} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getActivityDistribution, getDistributionLegend, getDistributionWidget, getStatsTableMap,
} from './widgets/distribution';
import {getAggregatedValue, getStats, Stats} from './utils/statistics';
import {LogoSummaryTable} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {_package, getMonomerWorksInstance, getTreeHelperInstance} from './package';
import {findMutations} from './utils/algorithms';
import {createDistanceMatrixWorker} from './utils/worker-creator';
import {calculateIdentity, calculateSimilarity} from './widgets/similarity';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

export type SummaryStats = {
  minCount: number, maxCount: number,
  minMeanDifference: number, maxMeanDifference: number,
  minPValue: number, maxPValue: number,
  minRatio: number, maxRatio: number,
};
export type PositionStats = {[monomer: string]: Stats} & {general: SummaryStats};
export type MonomerPositionStats = {[position: string]: PositionStats} & {general: SummaryStats};
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

  _isUpdating: boolean = false;
  isBitsetChangedInitialized = false;
  isCellChanging = false;
  isUserChangedSelection = true;

  df: DG.DataFrame;
  splitCol!: DG.Column<boolean>;
  _monomerPositionStats?: MonomerPositionStats;
  _clusterStats?: ClusterTypeStats;
  _mutationCliffsSelection!: type.PositionToAARList;
  _invariantMapSelection!: type.PositionToAARList;
  _clusterSelection!: string[];
  _mutationCliffs: type.MutationCliffs | null = null;
  isInitialized = false;
  _analysisView?: DG.TableView;

  monomerMap: {[key: string]: {molfile: string, fullName: string}} = {};
  monomerLib: IMonomerLib | null = null; // To get monomers from lib(s)
  monomerWorks: MonomerWorks | null = null; // To get processed monomers

  _settings!: type.PeptidesSettings;
  isRibbonSet = false;

  _cp?: SeqPalette;
  headerSelectedMonomers: type.MonomerSelectionStats = {};
  webLogoBounds: {[positon: string]: {[monomer: string]: DG.Rect}} = {};
  cachedWebLogoTooltip: {bar: string, tooltip: HTMLDivElement | null} = {bar: '', tooltip: null};
  _monomerPositionDf?: DG.DataFrame;
  _alphabet?: string;
  _mostPotentResiduesDf?: DG.DataFrame;
  _matrixDf?: DG.DataFrame;
  _splitSeqDf?: DG.DataFrame;
  _distanceMatrix!: DistanceMatrix;
  _dm!: DistanceMatrix;
  _layoutEventInitialized = false;
  isToolboxSet: boolean = false;

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

  get mutationCliffsSelection(): type.PositionToAARList {
    this._mutationCliffsSelection ??= JSON.parse(this.df.tags[C.TAGS.SELECTION] || '{}');
    return this._mutationCliffsSelection;
  }

  set mutationCliffsSelection(selection: type.PositionToAARList) {
    this._mutationCliffsSelection = selection;
    this.df.tags[C.TAGS.SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();

    const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    mpViewer?.viewerGrid.invalidate();
    const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    mprViewer?.viewerGrid.invalidate();

    this.analysisView.grid.invalidate();
  }

  get invariantMapSelection(): type.PositionToAARList {
    this._invariantMapSelection ??=
      JSON.parse(this.df.tags[C.TAGS.INVARIANT_MAP_SELECTION] || this.df.tags[C.TAGS.FILTER] || '{}');
    return this._invariantMapSelection;
  }

  set invariantMapSelection(selection: type.PositionToAARList) {
    this._invariantMapSelection = selection;
    this.df.tags[C.TAGS.INVARIANT_MAP_SELECTION] = JSON.stringify(selection);
    this.fireBitsetChanged();
    this.analysisView.grid.invalidate();
  }

  get clusterSelection(): string[] {
    this._clusterSelection ??= JSON.parse(this.df.tags[C.TAGS.CLUSTER_SELECTION] || '[]');
    return this._clusterSelection;
  }

  set clusterSelection(selection: string[]) {
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
    const splitByAARFlag = (this.df.tags['distributionSplit'] || '00')[1];
    this.df.tags['distributionSplit'] = `${flag ? 1 : 0}${splitByAARFlag}`;
  }

  get splitByAAR(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[1];
    return splitByPosFlag === '1' ? true : false;
  }

  set splitByAAR(flag: boolean) {
    const splitByAARFlag = (this.df.tags['distributionSplit'] || '00')[0];
    this.df.tags['distributionSplit'] = `${splitByAARFlag}${flag ? 1 : 0}`;
  }

  get isMonomerPositionSelectionEmpty(): boolean {
    for (const aarList of Object.values(this.mutationCliffsSelection)) {
      if (aarList.length !== 0)
        return false;
    }
    return true;
  }

  get isClusterSelectionEmpty(): boolean {
    return this.clusterSelection.length === 0;
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
        this.mostPotentResiduesDf = this.createMostPotentResiduesDf();
        this.clusterStats = this.calculateClusterStatistics();
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

  createMonomerPositionDf(): DG.DataFrame {
    const uniqueMonomers = new Set<string>();
    const splitSeqCols = this.splitSeqDf.columns;
    for (const col of splitSeqCols) {
      const colCat = col.categories;
      for (const cat of colCat) {
        if (cat !== '')
          uniqueMonomers.add(cat);
      }
    }

    const monomerCol = DG.Column.fromStrings(C.COLUMNS_NAMES.MONOMER, Array.from(uniqueMonomers));
    const monomerPositionDf = DG.DataFrame.fromColumns([monomerCol]);
    monomerPositionDf.name = 'SAR';
    for (const col of splitSeqCols)
      monomerPositionDf.columns.addNewBool(col.name);

    return monomerPositionDf;
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
        newView.onclick = () => trueModel.createNewView();
        newView.onmouseover = (ev) => ui.tooltip.show('Creates a new view from current selection', ev.clientX + 5, ev.clientY + 5);
        const newCluster = ui.label('New cluster');
        $(newCluster).addClass('d4-link-action');
        newCluster.onclick = () => {
          const lstViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
          if (lstViewer === null)
            throw new Error('Logo summary table viewer is not found');
          lstViewer.clusterFromSelection();
        };
        newCluster.onmouseover = (ev) => ui.tooltip.show('Creates a new cluster from selection', ev.clientX + 5, ev.clientY + 5);
        const removeCluster = ui.label('Remove cluster');
        $(removeCluster).addClass('d4-link-action');
        removeCluster.onclick = () => {
          const lstViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
          if (lstViewer === null)
            throw new Error('Logo summary table viewer is not found');
          lstViewer.removeCluster();
        };
        removeCluster.onmouseover = (ev) => ui.tooltip.show('Removes currently selected custom cluster', ev.clientX + 5, ev.clientY + 5);
        removeCluster.style.visibility = trueModel.clusterSelection.length === 0 ||
          !wu(this.customClusters).some((c) => trueModel.clusterSelection.includes(c.name)) ? 'hidden' : 'visible';
        return ui.divV([newView, newCluster, removeCluster]);
      });
    }
    const table = trueModel.df.filter.anyFalse ? trueModel.df.clone(trueModel.df.filter, null, true) : trueModel.df;
    acc.addPane('Mutation Cliffs pairs', () => mutationCliffsWidget(trueModel.df, trueModel).root);
    acc.addPane('Distribution', () => getDistributionWidget(table, trueModel).root);

    return acc;
  }

  updateGrid(): void {
    this.joinDataFrames();

    this.createScaledCol();

    this.initInvariantMapSelection({notify: false});
    this.initMutationCliffsSelection({notify: false});

    this.setWebLogoInteraction();
    this.webLogoBounds = {};

    this.setCellRenderers();

    this.setTooltips();

    this.setBitsetCallback();

    this.setGridProperties();
  }

  initInvariantMapSelection(options: {cleanInit?: boolean, notify?: boolean} = {}): void {
    options.cleanInit ??= false;
    options.notify ??= true;

    const tempFilter: type.PositionToAARList = this.invariantMapSelection;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns) {
      if (options.cleanInit || !tempFilter.hasOwnProperty(pos))
        tempFilter[pos] = [];
    }

    if (options.notify)
      this.invariantMapSelection = tempFilter;
    else
      this._invariantMapSelection = tempFilter;
  }

  initMutationCliffsSelection(options: {cleanInit?: boolean, notify?: boolean} = {}): void {
    options.cleanInit ??= false;
    options.notify ??= true;

    const tempSelection: type.PositionToAARList = this.mutationCliffsSelection;
    const positionColumns = this.splitSeqDf.columns.names();
    for (const pos of positionColumns) {
      if (options.cleanInit || !tempSelection.hasOwnProperty(pos))
        tempSelection[pos] = [];
    }

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
      CR.setAARRenderer(col, this.alphabet);
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

  calculateMonomerPositionStatistics(): MonomerPositionStats {
    const positionColumns = this.splitSeqDf.columns.toList();
    const sourceDfLen = this.df.rowCount;
    const monomerPositionObject = {general: {}} as MonomerPositionStats & { general: SummaryStats };
    const activityColData = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).getRawData();

    for (const posCol of positionColumns) {
      const posColData = posCol.getRawData();
      const posColCateogries = posCol.categories;
      const currentPositionObject = {general: {}} as PositionStats & { general: SummaryStats };

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
        const stats = getStats(activityColData, bitArray);
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

    for (let clustType = 0; clustType < 2; ++clustType) {
      const masks = clustType === 0 ? origClustMasks : customClustMasks;
      const clustNames = clustType === 0 ? origClustColCat : customClustColNamesList;
      const resultStats = clustType === 0 ? origClustStats : customClustStats;
      for (let maskIdx = 0; maskIdx < masks.length; ++maskIdx) {
        const mask = masks[maskIdx];
        const stats = getStats(activityColData, mask);
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
        if (key === 'general')
          return false;

        return (v[1] as Stats).pValue === generalPositionStats.minPValue;
      }) as [string, Stats][];

      let maxEntry: [string, Stats];
      for (const [monomer, monomerStats] of filteredMonomerStats) {
        if (typeof maxEntry! === 'undefined' || maxEntry[1].meanDifference < monomerStats.meanDifference)
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
    const tempSelection = this.clusterSelection;
    const idx = tempSelection.indexOf(cluster);
    if (idx !== -1)
      tempSelection.splice(idx, 1);
    else
      tempSelection.push(cluster);

    this.clusterSelection = tempSelection;
  }

  initClusterSelection(options: {notify?: boolean} = {}): void {
    options.notify ??= true;

    if (options.notify)
      this.clusterSelection = [];
    else
      this._clusterSelection = [];
  }

  setWebLogoInteraction(): void {
    const sourceView = this.analysisView.grid;
    const eventAction = (ev: MouseEvent): void => {
      const cell = sourceView.hitTest(ev.offsetX, ev.offsetY);
      if (cell?.isColHeader && cell.tableColumn?.semType === C.SEM_TYPES.MONOMER) {
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

  requestBarchartAction(ev: MouseEvent, barPart: {position: string, monomer: string} | null): void {
    if (!barPart)
      return;
    const monomer = barPart.monomer;
    const position = barPart.position;
    if (ev.type === 'click') {
      if (!ev.shiftKey)
        this.initInvariantMapSelection({cleanInit: true, notify: false});

      this.modifyMonomerPositionSelection(monomer, position, true);
    } else {
      const bar = `${position} = ${monomer}`;
      if (this.cachedWebLogoTooltip.bar === bar)
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
          const stats = this.monomerPositionStats[col.name];
          //TODO: precalc on stats creation
          const sortedStatsOrder = Object.keys(stats).sort((a, b) => {
            if (a === '' || a === '-')
              return -1;
            else if (b === '' || b === '-')
              return +1;
            return 0;
          }).filter((v) => v !== 'general');

          this.webLogoBounds[col.name] = CR.drawLogoInBounds(ctx, bounds, stats, col.name, sortedStatsOrder,
            this.df.rowCount, this.cp, this.headerSelectedMonomers[col.name]);
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
        if (layout.view.id === this.analysisView.id) {
          // this.analysisView.grid.onCellRender.subscribe((gcArgs) => headerRenderer(gcArgs));
          this.updateGrid();
        }
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

  showMonomerTooltip(aar: string, x: number, y: number): boolean {
    const tooltipElements: HTMLDivElement[] = [];
    const monomerName = aar.toLowerCase();

    const mw = getMonomerWorksInstance();
    const mol = mw?.getCappedRotatedMonomer('PEPTIDE', aar);

    if (mol) {
      tooltipElements.push(ui.div(monomerName));
      const options = {autoCrop: true, autoCropMargin: 0, suppressChiralText: true};
      tooltipElements.push(grok.chem.svgMol(mol, undefined, undefined, options));
    } else if (aar !== '')
      tooltipElements.push(ui.div(aar));
    else
      return true;


    ui.tooltip.show(ui.divV(tooltipElements), x, y);

    return true;
  }

  //TODO: move out to viewer code
  showTooltipAt(aar: string, position: string, x: number, y: number): HTMLDivElement | null {
    const stats = this.monomerPositionStats[position][aar];
    if (!stats?.count)
      return null;

    const activityCol = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
    const posCol = this.df.getCol(position);
    const posColCategories = posCol.categories;
    const aarCategoryIndex = posColCategories.indexOf(aar);
    const posColData = posCol.getRawData();
    const mask = DG.BitSet.create(activityCol.length, (i) => posColData[i] === aarCategoryIndex);

    const distributionTable = DG.DataFrame.fromColumns(
      [activityCol, DG.Column.fromBitSet(C.COLUMNS_NAMES.SPLIT_COL, mask)]);
    const labels = getDistributionLegend(`${position} : ${aar}`, 'Other');
    const hist = getActivityDistribution(prepareTableForHistogram(distributionTable), true);
    const tableMap = getStatsTableMap(stats);
    const aggregatedColMap = this.getAggregatedColumnValues({mask: mask});

    const resultMap = {...tableMap, ...aggregatedColMap};
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

  modifyMonomerPositionSelection(aar: string, position: string, isInvariantMap: boolean): void {
    if (this.df.getCol(position).categories.indexOf(aar) === -1)
      return;

    const tempSelection = isInvariantMap ? this.invariantMapSelection : this.mutationCliffsSelection;
    const tempSelectionAt = tempSelection[position];
    const aarIndex = tempSelectionAt.indexOf(aar);
    if (aarIndex === -1)
      tempSelectionAt.push(aar);
    else
      tempSelectionAt.splice(aarIndex, 1);

    if (isInvariantMap)
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

    const changeSelectionBitset = (currentBitset: DG.BitSet, clustColCat: string[],
      clustColData: type.RawData, customClust: {[key: string]: BitArray}): void => {
      const indexes = new Set<number>();
      for (const [position, monomerList] of Object.entries(this.mutationCliffsSelection)) {
        for (const monomer of monomerList) {
          const substitutions = this.mutationCliffs?.get(monomer)?.get(position) ?? null;
          if (substitutions === null)
            continue;
          for (const [key, value] of substitutions.entries()) {
            indexes.add(key);
            for (const v of value)
              indexes.add(v);
          }
        }
      }

      const positionList = Object.keys(this.invariantMapSelection);
      const rowCount = this.df.rowCount;
      for (const position of positionList) {
        const positionCol: DG.Column<string> = this.df.getCol(position);
        const positionColData = positionCol.getRawData();
        const positionColCat = positionCol.categories;
        const aarList = this.invariantMapSelection[position];
        for (const aar of aarList) {
          const aarIndex = positionColCat.indexOf(aar);
          if (aarIndex === -1)
            continue;
          for (let i = 0; i < rowCount; ++i) {
            if (positionColData[i] === aarIndex)
              indexes.add(i);
          }
        }
      }

      const getBitAt = (i: number): boolean => {
        if (indexes.has(i))
          return true;

        //TODO: preprocess
        const currentOrigClust = clustColCat[clustColData[i]];
        if (typeof currentOrigClust === undefined)
          return false;

        for (const clust of this.clusterSelection) {
          if (clust === currentOrigClust)
            return true;

          if (Object.hasOwn(customClust, clust) && customClust[clust].getBit(i))
            return true;
        }

        return false;
      };
      currentBitset.init((i) => getBitAt(i), false);
    };

    selection.onChanged.subscribe(() => {
      if (this.isUserChangedSelection)
        return;

      const clustColCat = clusterCol?.categories ?? [];
      const clustColData = clusterCol?.getRawData() ?? new Int32Array(0);
      const customClust: {[key: string]: BitArray} = {};
      const rowCount = this.df.rowCount;
      for (const clust of this.customClusters)
        customClust[clust.name] = BitArray.fromUint32Array(rowCount, clust.getRawData() as Uint32Array);

      changeSelectionBitset(selection, clustColCat, clustColData, customClust);
    });

    filter.onChanged.subscribe(() => {
      const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
      if (lstViewer !== null && typeof lstViewer.model !== 'undefined') {
        lstViewer.createLogoSummaryTableGrid();
        lstViewer.render();
      }
    });
    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(fireFilterChanged: boolean = false): void {
    this.isUserChangedSelection = false;
    this.df.selection.fireChanged();
    if (fireFilterChanged)
      this.df.filter.fireChanged();
    this.modifyOrCreateSplitCol();
    this.headerSelectedMonomers = calculateSelected(this.df);

    const acc = this.createAccordion();
    if (acc !== null) {
      grok.shell.o = acc.root;
      for (const pane of acc.panes)
        pane.expanded = true;
    }
    this.isUserChangedSelection = true;
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

    if (!this.isRibbonSet && this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1') {
      //TODO: don't pass model, pass parameters instead
      const settingsButton = ui.iconFA('wrench', () => getSettingsDialog(this), 'Peptides analysis settings');
      this.analysisView.setRibbonPanels([[settingsButton]], false);
      this.isRibbonSet = true;
      this.updateGrid();
    }

    if (!this.isToolboxSet && this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1') {
      let template: ISeqSplitted;
      const sequencesCol = this.df.getCol(this.settings.sequenceColumnName!);
      const minTemplateLength = this.splitSeqDf.columns.toList()
        .filter((col) => col.stats.missingValueCount === 0).length;
      const calculateIdentityBtn = ui.button('Identity', async () => {
        let identityScoresCol = calculateIdentity(template, this.splitSeqDf);
        identityScoresCol.name = this.df.columns.getUnusedName(identityScoresCol.name);
        identityScoresCol = this.df.columns.add(identityScoresCol);
        identityScoresCol.setTag(C.TAGS.IDENTITY_TEMPLATE, new Array(template).join(' '));
      }, 'Calculate identity');
      const calculateSimilarityBtn = ui.button('Similarity', async () => {
        let similarityScoresCol = await calculateSimilarity(template, this.splitSeqDf);
        similarityScoresCol.name = this.df.columns.getUnusedName(similarityScoresCol.name);
        similarityScoresCol = this.df.columns.add(similarityScoresCol);
        similarityScoresCol.setTag(C.TAGS.SIMILARITY_TEMPLATE, new Array(template).join(' '));
      }, 'Calculate similarity');
      const templateInput = ui.stringInput('Template', this.identityTemplate, async () => {
        this.identityTemplate = templateInput.value;
        if (isNaN(parseInt(templateInput.value))) {
          if (templateInput.value.length === 0) {
            calculateIdentityBtn.disabled = true;
            calculateSimilarityBtn.disabled = true;
            return;
          }
          try {
            template ??= await getTemplate(this.identityTemplate, sequencesCol);
            if (template.length < minTemplateLength) {
              grok.shell.warning(`Template length should be at least ${minTemplateLength} amino acids.`);
              calculateIdentityBtn.disabled = true;
              calculateSimilarityBtn.disabled = true;
              return;
            } else if (new Array(template).includes('') || new Array(template).includes('-')) {
              grok.shell.warning('Template shouldn\'t contain gaps or empty cells.');
              calculateIdentityBtn.disabled = true;
              calculateSimilarityBtn.disabled = true;
              return;
            }
          } catch (e) {
            grok.shell.warning(`Only ${sequencesCol.getTag(DG.TAGS.UNITS)} sequence format is supported.`);
            grok.log.warning(e as string);
            calculateIdentityBtn.disabled = true;
            return;
          }
        } else {
          const rowIndex = parseInt(templateInput.value) - 1;
          const selectedIndexes = this.df.filter.getSelectedIndexes();
          if (rowIndex < 0 || rowIndex >= selectedIndexes.length) {
            grok.shell.warning('Invalid row index');
            calculateIdentityBtn.disabled = true;
            calculateSimilarityBtn.disabled = true;
            return;
          }
          this.identityTemplate = sequencesCol.get(selectedIndexes[rowIndex]);
        }
        try {
          template = await getTemplate(this.identityTemplate, sequencesCol);
        } catch (e) {
          grok.shell.warning('Couldn\'t recognize sequence format.');
          grok.log.warning(e as string);
          calculateIdentityBtn.disabled = true;
          calculateSimilarityBtn.disabled = true;
          return;
        }
        calculateIdentityBtn.disabled = false;
        calculateSimilarityBtn.disabled = false;
      }, {placeholder: 'Sequence or row index...'});
      templateInput.setTooltip('Template sequence. Can be row index, peptide ID or sequence.');
      templateInput.fireChanged();
      const acc = this.analysisView.toolboxPage.accordion;
      acc.addPane('Sequence Identity and Similarity',
        () => ui.divV([ui.form([templateInput]), calculateIdentityBtn, calculateSimilarityBtn]), true, acc.panes[0]);
      this.isToolboxSet = true;
    }

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
}
