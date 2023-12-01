import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {monomerToShort, pickUpPalette, TAGS as bioTAGS, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {calculateScores, SCORE} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {BitArrayMetrics, StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

import wu from 'wu';
import * as rxjs from 'rxjs';
import * as uuid from 'uuid';
import $ from 'cash-dom';

import * as C from './utils/constants';
import * as type from './utils/types';
import {calculateSelected, extractColInfo, scaleActivity} from './utils/misc';
import {MONOMER_POSITION_PROPERTIES, MonomerPosition, MostPotentResidues} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getDistributionWidget} from './widgets/distribution';
import {ClusterTypeStats, MonomerPositionStats, PositionStats} from './utils/statistics';
import {LogoSummaryTable} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {_package, getTreeHelperInstance} from './package';
import {calculateClusterStatistics, calculateMonomerPositionStatistics, findMutations} from './utils/algorithms';
import {createDistanceMatrixWorker} from './utils/worker-creator';
import {getSelectionWidget} from './widgets/selection';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {DimReductionMethods, ITSNEOptions, IUMAPOptions} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {showMonomerTooltip} from './utils/tooltips';

export enum CLUSTER_TYPE {
  ORIGINAL = 'original',
  CUSTOM = 'custom',
}
export type ClusterType = `${CLUSTER_TYPE}`;
export enum VIEWER_TYPE {
  MONOMER_POSITION = 'Monomer-Position',
  MOST_POTENT_RESIDUES = 'Most Potent Residues',
  LOGO_SUMMARY_TABLE = 'Logo Summary Table',
  DENDROGRAM = 'Dendrogram',
}

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
  webLogoBounds: CR.WebLogoBounds = {};
  cachedWebLogoTooltip: {bar: string, tooltip: HTMLDivElement | null} = {bar: '', tooltip: null};
  _dm!: DistanceMatrix;
  _layoutEventInitialized = false;

  subs: rxjs.Subscription[] = [];
  isHighlighting: boolean = false;
  controlFire: boolean = false;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    if (dataFrame.columns.contains(C.COLUMNS_NAMES.ACTIVITY_SCALED) && !dataFrame.columns.contains(C.COLUMNS_NAMES.ACTIVITY))
      dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).name = C.COLUMNS_NAMES.ACTIVITY;

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
    this._monomerPositionStats ??= calculateMonomerPositionStatistics(this.df, this.positionColumns.toArray());
    return this._monomerPositionStats!;
  }

  set monomerPositionStats(mps: MonomerPositionStats) {
    this._monomerPositionStats = mps;
  }

  get positionColumns(): wu.WuIterable<DG.Column> {
    return wu(this.df.columns.byTags({[C.TAGS.POSITION_COL]: `${true}`}));
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
    this._clusterStats ??= calculateClusterStatistics(this.df, this.settings.clustersColumnName!,
      this.customClusters.toArray());
    return this._clusterStats!;
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
        const posCols = this.positionColumns.toArray().map((col) => col.name);

        for (let colIdx = 1; colIdx < this._analysisView.grid.columns.length; ++colIdx) {
          const gridCol = this._analysisView.grid.columns.byIndex(colIdx)!;
          gridCol.visible =
            posCols.includes(gridCol.column!.name) || (gridCol.column!.name === C.COLUMNS_NAMES.ACTIVITY);
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
    return splitByPosFlag === '1';
  }

  set splitByPos(flag: boolean) {
    const splitByMonomerFlag = (this.df.tags['distributionSplit'] || '00')[1];
    this.df.tags['distributionSplit'] = `${flag ? 1 : 0}${splitByMonomerFlag}`;
  }

  get splitByMonomer(): boolean {
    const splitByPosFlag = (this.df.tags['distributionSplit'] || '00')[1];
    return splitByPosFlag === '1';
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

  get customClusters(): wu.WuIterable<DG.Column<boolean>> {
    const query: { [key: string]: string } = {};
    query[C.TAGS.CUSTOM_CLUSTER] = '1';
    return wu(this.df.columns.byTags(query));
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
    let updateViewersData = false;
    for (const variable of updateVars) {
      switch (variable) {
      case 'activity':
        this.createScaledCol();
        updateViewersData = true;
        break;
      case 'mutationCliffs':
        this.updateMutationCliffs().then(() => {
          (this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition)?.viewerGrid.invalidate();
          (this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues)?.viewerGrid.invalidate();
        }).catch((e) => _package.logger.debug(e));
        break;
      case 'stats':
        this.monomerPositionStats = calculateMonomerPositionStatistics(this.df, this.positionColumns.toArray());
        this.clusterStats = calculateClusterStatistics(this.df, this.settings.clustersColumnName!,
          this.customClusters.toArray());
        updateViewersData = true;
        break;
      case 'grid':
        this.setGridProperties();
        const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
        lstViewer?.createLogoSummaryTableGrid();
        lstViewer?.render();
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
    if (updateViewersData)
      mpViewer?.createMonomerPositionGrid();
    mpViewer?.render();
    const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    if (updateViewersData)
      mprViewer?.createMostPotentResiduesGrid();
    mprViewer?.render();
    const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    if (updateViewersData)
      lstViewer?.createLogoSummaryTableGrid();
    lstViewer?.render();
  }

  async updateMutationCliffs(notify: boolean = true): Promise<void> {
    const scaledActivityCol: DG.Column<number> = this.df.getCol(C.COLUMNS_NAMES.ACTIVITY);
    //TODO: set categories ordering the same to share compare indexes instead of strings
    const monomerCols: type.RawColumn[] = this.df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER).map(extractColInfo);
    const targetCol = typeof this.settings.targetColumnName !== 'undefined' ?
      extractColInfo(this.df.getCol(this.settings.targetColumnName)) : null;
    let mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    const currentTarget = mpViewer?.getProperty(MONOMER_POSITION_PROPERTIES.TARGET)?.get(mpViewer) as string | undefined;
    const targetOptions = {targetCol: targetCol, currentTarget: currentTarget};
    const mutationCliffs = await findMutations(scaledActivityCol.getRawData(), monomerCols, this.settings, targetOptions);
    if (notify)
      this.mutationCliffs = mutationCliffs;
    else
      this._mutationCliffs = mutationCliffs;

    mpViewer ??= this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    mpViewer?.render(true);
    const mostPotentViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    mostPotentViewer?.render(true);
  }

  buildSplitSeqDf(): DG.DataFrame {
    const sequenceCol = this.df.getCol(this.settings.sequenceColumnName!);
    return splitAlignedSequences(sequenceCol);
  }

  getVisibleSelection(): DG.BitSet {
    return this.df.selection.clone().and(this.df.filter);
  }

  createAccordion(): DG.Accordion | null {
    const trueModel: PeptidesModel | undefined = grok.shell.t?.temp[PeptidesModel.modelName];
    if (!trueModel)
      return null;

    const acc = ui.accordion('Peptides analysis panel');
    acc.root.style.width = '100%';
    const filterAndSelectionBs = trueModel.getVisibleSelection();
    const filteredTitlePart = trueModel.df.filter.anyFalse ? ` among ${trueModel.df.filter.trueCount} filtered` : '';
    const getSelectionString = (selection: type.Selection): string => {
      const selectedMonomerPositions: string[] = [];
      for (const [pos, monomerList] of Object.entries(selection)) {
        for (const monomer of monomerList)
          selectedMonomerPositions.push(`${pos}:${monomer}`);
      }
      return selectedMonomerPositions.join(', ');
    };

    const selectionDescription = [];
    const selectedClusters = trueModel.clusterSelection[CLUSTER_TYPE.ORIGINAL]
        .concat(trueModel.clusterSelection[CLUSTER_TYPE.CUSTOM]).join(', ');
    if (selectedClusters.length !== 0)
      ui.divText(`Selected clusters: ${selectedClusters}`);
    const selectedMonomerPositions = getSelectionString(trueModel.invariantMapSelection);
    if (selectedMonomerPositions.length !== 0)
      selectionDescription.push(ui.divText(`Selected monomer-positions: ${selectedMonomerPositions}`));
    const selectedMutationCliffs = getSelectionString(trueModel.mutationCliffsSelection);
    if (selectedMutationCliffs.length !== 0)
      selectionDescription.push(ui.divText(`Selected mutation cliffs pairs: ${selectedMutationCliffs}`));

    const descritionsHost = ui.div(ui.divV(selectionDescription));
    acc.addTitle(ui.divV([
      ui.h1(`${filterAndSelectionBs.trueCount} selected rows${filteredTitlePart}`),
      descritionsHost,
    ], 'css-gap-small'));

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
      }, true);
    }
    const table = trueModel.df.filter.anyFalse ? trueModel.df.clone(trueModel.df.filter, null, true) : trueModel.df;
    acc.addPane('Mutation Cliffs pairs', () => mutationCliffsWidget(trueModel.df, trueModel).root, true);
    acc.addPane('Distribution', () => getDistributionWidget(table, trueModel).root, true);
    acc.addPane('Selection', () => getSelectionWidget(trueModel.df, trueModel), true);

    return acc;
  }

  updateGrid(): void {
    this.joinDataFrames();
    this.createScaledCol();
    // this.setWebLogoInteraction();
    this.webLogoBounds = {};

    CR.setWebLogoRenderer(this.analysisView.grid, this);
    if (!this._layoutEventInitialized) {
      grok.events.onViewLayoutApplied.subscribe((layout) => {
        if (layout.view.id === this.analysisView.id)
          this.updateGrid();
      });
      this._layoutEventInitialized = true;
    }

    this.setTooltips();
    this.setBitsetCallback();
    this.setGridProperties();
  }

  initInvariantMapSelection(options: {notify?: boolean} = {}): void {
    options.notify ??= true;

    const tempSelection: type.Selection = {};
    const positionColumns = this.positionColumns.toArray().map((col) => col.name);
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
    const positionColumns = this.positionColumns.toArray().map((col) => col.name);
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
    const splitSeqDf = this.buildSplitSeqDf();
    const positionColumns = splitSeqDf.columns.names();
    for (const colName of positionColumns) {
      let col = this.df.col(colName);
      const newCol = splitSeqDf.getCol(colName);
      if (col !== null)
        cols.remove(colName);

      const newColCat = newCol.categories;
      const newColData = newCol.getRawData();
      col = cols.addNew(newCol.name, newCol.type).init((i) => newColCat[newColData[i]]);
      col.setTag(C.TAGS.ANALYSIS_COL, `${true}`);
      col.setTag(C.TAGS.POSITION_COL, `${true}`);
      CR.setMonomerRenderer(col, this.alphabet);
    }
    this.df.name = name;
  }

  createScaledCol(): void {
    const sourceGrid = this.analysisView.grid;
    const scaledCol = scaleActivity(this.df.getCol(this.settings.activityColumnName!), this.settings.scaling);
    //TODO: make another func
    this.df.columns.replace(C.COLUMNS_NAMES.ACTIVITY, scaledCol);

    sourceGrid.columns.setOrder([scaledCol.name]);
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

  setTooltips(): void {
    this.analysisView.grid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader && cell.tableColumn!.semType === C.SEM_TYPES.MONOMER)
        return true;
      if (!(cell.isTableCell && cell.tableColumn!.semType === C.SEM_TYPES.MONOMER))
        return false;

      showMonomerTooltip(cell.cell.value, x, y);
      return true;
    });
  }

  getCombinedSelection(): DG.BitSet {
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
  }

  setBitsetCallback(): void {
    if (this.isBitsetChangedInitialized)
      return;
    const selection = this.df.selection;
    const filter = this.df.filter;

    const showAccordion = (): void => {
      const acc = this.createAccordion();
      if (acc === null)
        return;
      grok.shell.o = acc.root;
    };

    selection.onChanged.subscribe(() => {
      if (this.controlFire) {
        this.controlFire = false;
        return;
      }
      try {
        if (!this.isUserChangedSelection)
          selection.copyFrom(this.getCombinedSelection(), false);
      } catch (e) {
        _package.logger.debug('Peptides: Error on selection changed');
        _package.logger.debug(e as string);
      } finally {
        showAccordion();
      }
    });

    filter.onChanged.subscribe(() => {
      try {
        if (this.controlFire) {
          this.controlFire = false;
          return;
        }
        const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
        if (lstViewer !== null && typeof lstViewer.model !== 'undefined') {
          lstViewer.createLogoSummaryTableGrid();
          lstViewer.render();
        }
      } catch (e) {
        _package.logger.debug('Peptides: Error on filter changed');
        _package.logger.debug(e as string);
      } finally {
        showAccordion();
      }
    });

    this.isBitsetChangedInitialized = true;
  }

  fireBitsetChanged(fireFilterChanged: boolean = false): void {
    this.isUserChangedSelection = false;
    this.df.selection.fireChanged();
    if (fireFilterChanged)
      this.df.filter.fireChanged();

    // Fire bitset changed event again to update UI
    this.controlFire = true;
    this.df.selection.fireChanged();
    if (fireFilterChanged)
      this.df.filter.fireChanged();

    this.isUserChangedSelection = true;
    this.headerSelectedMonomers = calculateSelected(this.df);
  }

  setGridProperties(props?: DG.IGridLookSettings): void {
    const sourceGrid = this.analysisView.grid;
    const sourceGridProps = sourceGrid.props;
    sourceGridProps.allowColSelection = props?.allowColSelection ?? false;
    sourceGridProps.allowEdit = props?.allowEdit ?? false;
    sourceGridProps.showCurrentRowIndicator = props?.showCurrentRowIndicator ?? false;
    this.df.temp[C.EMBEDDING_STATUS] = false;
    const positionCols = this.positionColumns.toArray();
    let maxWidth = 10;
    const canvasContext = sourceGrid.canvas.getContext('2d');
    for (const positionCol of positionCols) {
      // Longest category
      const maxCategory = monomerToShort(positionCol.categories.reduce((a, b) => a.length > b.length ? a : b), 6);
      // Measure text width of longest category
      const width = Math.ceil(canvasContext!.measureText(maxCategory).width);
      maxWidth = Math.max(maxWidth, width);
    }
    setTimeout(() => {
      for (const positionCol of positionCols)
        sourceGrid.col(positionCol.name)!.width = maxWidth + 15;
    }, 100);
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
        calculateScores(this.df, seqCol, seqCol.get(this.df.currentRowIdx), SCORE.IDENTITY)
            .then((col: DG.Column<number>) => col.setTag(C.TAGS.IDENTITY_TEMPLATE, seqCol.get(this.df.currentRowIdx)))
            .catch((e) => _package.logger.debug(e));
      };
      actionsHost.append(ui.span([calculateIdentity], 'd4-markdown-row'));

      const calculateSimilarity = ui.label('Calculate similarity');
      calculateSimilarity.classList.add('d4-link-action');
      ui.tooltip.bind(calculateSimilarity, 'Adds a column with sequence similarity scores against sequence in the current row');
      calculateSimilarity.onclick = (): void => {
        const seqCol = this.df.getCol(this.settings.sequenceColumnName!);
        calculateScores(this.df, seqCol, seqCol.get(this.df.currentRowIdx), SCORE.SIMILARITY)
            .then((col: DG.Column<number>) => col.setTag(C.TAGS.SIMILARITY_TEMPLATE, seqCol.get(this.df.currentRowIdx)))
            .catch((e) => _package.logger.debug(e));
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
    if (typeof this.settings.targetColumnName === 'undefined') {
      this.updateMutationCliffs().then(() => {
        (this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition)?.viewerGrid.invalidate();
        (this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues)?.viewerGrid.invalidate();
      }).catch((e) => _package.logger.debug(e));
    }

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
    const newClusterCol = DG.Column.fromBitSet(clusterName, this.getVisibleSelection());
    newClusterCol.setTag(C.TAGS.CUSTOM_CLUSTER, '1');
    newClusterCol.setTag(C.TAGS.ANALYSIS_COL, `${true}`);
    this.df.columns.add(newClusterCol);
    this.analysisView.grid.col(newClusterCol.name)!.visible = false;
  }

  createNewView(): string {
    const rowMask = this.getVisibleSelection();
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

  async addSequenceSpace(): Promise<void> {
    let seqCol = this.df.getCol(this.settings.sequenceColumnName!);
    const uh = UnitsHandler.getOrCreate(seqCol);
    const isHelm = uh.isHelm();
    if (isHelm) {
      try {
        grok.shell.warning('Column is in HELM notation. Sequences space will linearize sequences from position 0 prior to analysis');
        const linearCol = uh.convert(NOTATION.SEPARATOR, '/');
        const newName = this.df.columns.getUnusedName(`Separator(${seqCol.name})`);
        linearCol.name = newName;
        this.df.columns.add(linearCol, true);
        this.analysisView.grid.col(newName)!.visible = false;
        seqCol = linearCol;
      } catch (e) {
        grok.shell.error('Error on converting HELM notation to linear notation');
        grok.shell.error(e as string);
        return;
      }
    }
    const seqSpaceParams: {table: DG.DataFrame, molecules: DG.Column, methodName: DimReductionMethods,
      similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames, plotEmbeddings: boolean,
      sparseMatrixThreshold?: number, options?: (IUMAPOptions | ITSNEOptions) & Options} =
      {table: this.df, molecules: seqCol,
        methodName: DimReductionMethods.UMAP, similarityMetric: MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH,
        plotEmbeddings: true, sparseMatrixThreshold: 0.3, options: {'bypassLargeDataWarning': true}};

    // Use counter to unsubscribe when 2 columns are hidden
    let counter = 0;
    const columnAddedSub = this.df.onColumnsAdded.subscribe((colArgs: DG.ColumnsArgs) => {
      for (const col of colArgs.columns) {
        if (col.name.startsWith('Embed_')) {
          this.analysisView.grid.col(col.name)!.visible = false;
          counter++;
        }
      }
      if (counter === 2)
        columnAddedSub.unsubscribe();
    });

    const seqSpaceViewer: DG.ScatterPlotViewer | undefined = await grok.functions.call('Bio:sequenceSpaceTopMenu', seqSpaceParams);
    if (!(seqSpaceViewer instanceof DG.ScatterPlotViewer))
      return;
    seqSpaceViewer.props.colorColumnName = C.COLUMNS_NAMES.ACTIVITY;
    seqSpaceViewer.props.showXSelector = false;
    seqSpaceViewer.props.showYSelector = false;
  }
}
