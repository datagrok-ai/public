import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {monomerToShort, pickUpPalette, TAGS as bioTAGS, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {calculateScores, SCORE} from '@datagrok-libraries/bio/src/utils/macromolecule/scoring';
import {Options} from '@datagrok-libraries/utils/src/type-declarations';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {BitArrayMetrics, StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {ITreeHelper} from '@datagrok-libraries/bio/src/trees/tree-helper';
import {TAGS as treeTAGS} from '@datagrok-libraries/bio/src/trees';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import wu from 'wu';
import * as rxjs from 'rxjs';
import $ from 'cash-dom';

import * as C from './utils/constants';
import {COLUMNS_NAMES} from './utils/constants';
import * as type from './utils/types';
import {PeptidesSettings} from './utils/types';
import {
  areParametersEqual,
  calculateSelected,
  getSelectionBitset,
  highlightMonomerPosition,
  initSelection,
  modifySelection,
  mutationCliffsToMaskInfo,
  scaleActivity,
} from './utils/misc';
import {ISARViewer, MonomerPosition, MostPotentResidues, SARViewer} from './viewers/sar-viewer';
import * as CR from './utils/cell-renderer';
import {mutationCliffsWidget} from './widgets/mutation-cliffs';
import {getDistributionWidget, PeptideViewer} from './widgets/distribution';
import {CLUSTER_TYPE, ILogoSummaryTable, LogoSummaryTable, LST_PROPERTIES} from './viewers/logo-summary';
import {getSettingsDialog} from './widgets/settings';
import {_package, getTreeHelperInstance} from './package';
import {calculateMonomerPositionStatistics} from './utils/algorithms';
import {createDistanceMatrixWorker} from './utils/worker-creator';
import {getSelectionWidget} from './widgets/selection';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {DimReductionMethods, ITSNEOptions, IUMAPOptions} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {showMonomerTooltip} from './utils/tooltips';
import {AggregationColumns, MonomerPositionStats} from './utils/statistics';
import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';

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
  _dm: DistanceMatrix | null = null;
  isInitialized = false;
  isRibbonSet = false;
  webLogoSelectedMonomers: type.SelectionStats = {};
  webLogoBounds: CR.WebLogoBounds = {};
  cachedWebLogoTooltip: {
    bar: string,
    tooltip: HTMLDivElement | null
  } = {bar: '', tooltip: null};
  _layoutEventInitialized = false;
  subs: rxjs.Subscription[] = [];
  isHighlighting: boolean = false;
  controlFire: boolean = false;
  accordionSource: VIEWER_TYPE | null = null;
  _scaledActivityColumn: DG.Column | null = null;

  private constructor(dataFrame: DG.DataFrame) {
    this.df = dataFrame;
  }

  _monomerPositionStats: MonomerPositionStats | null = null;

  get monomerPositionStats(): MonomerPositionStats | null {
    if (this._monomerPositionStats !== null)
      return this._monomerPositionStats;

    if (this.positionColumns === null)
      return null;

    this._monomerPositionStats ??= calculateMonomerPositionStatistics(this.df.getCol(this.settings!.activityColumnName),
      this.df.filter, this.positionColumns);
    return this._monomerPositionStats;
  }

  _analysisView?: DG.TableView;

  get analysisView(): DG.TableView {
    if (this._analysisView === undefined) {
      this._analysisView = wu(grok.shell.tableViews).find(({dataFrame}) => dataFrame?.getTag(DG.TAGS.ID) === this.id);
      if (typeof this._analysisView === 'undefined')
        this._analysisView = grok.shell.addTableView(this.df);

      const posCols = this.positionColumns?.map((col) => col.name);
      if (posCols != null) {
        for (let colIdx = 1; colIdx < this._analysisView.grid.columns.length; ++colIdx) {
          const gridCol = this._analysisView.grid.columns.byIndex(colIdx);
          if (gridCol === null)
            throw new Error(`PeptidesError: Could not get analysis view: grid column with index '${colIdx}' is null`);
          else if (gridCol.column === null) {
            throw new Error(`PeptidesError: Could not get analysis view: grid column with index '${colIdx}' has null ` +
              `column`);
          }

          gridCol.visible = posCols.includes(gridCol.column.name) || (gridCol.column.name === C.COLUMNS_NAMES.ACTIVITY);
        }
      }
    }

    if (this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1' && !this._layoutEventInitialized)
      grok.shell.v = this._analysisView;

    this._analysisView.grid.invalidate();
    return this._analysisView;
  }

  _settings: type.PeptidesSettings | null = null;

  get settings(): type.PeptidesSettings | null {
    const settingsStr = this.df.getTag(C.TAGS.SETTINGS);
    if (settingsStr == null)
      return null;
    this._settings ??= JSON.parse(settingsStr);
    return this._settings!;
  }

  set settings(s: type.PartialPeptidesSettings) {
    const newSettingsEntries = Object.entries(s) as ([keyof type.PeptidesSettings, never])[];
    const updateVars: Set<string> = new Set();
    for (const [key, value] of newSettingsEntries) {
      this.settings![key] = value;
      switch (key) {
      case 'activityColumnName':
      case 'activityScaling':
        updateVars.add('activity');
        updateVars.add('stats');
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
      case 'stats':
        this.webLogoSelection = {};
        this.webLogoBounds = {};
        this.cachedWebLogoTooltip = {bar: '', tooltip: null};
        this._monomerPositionStats = null;
        break;
      case 'dendrogram':
          this.settings!.showDendrogram ? this.addDendrogram() : this.closeViewer(VIEWER_TYPE.DENDROGRAM);
        break;
      case 'logoSummaryTable':
          this.settings!.showLogoSummaryTable ? this.addLogoSummaryTable() :
            this.closeViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE);
        break;
      case 'monomerPosition':
          this.settings!.showMonomerPosition ? this.addMonomerPosition() :
            this.closeViewer(VIEWER_TYPE.MONOMER_POSITION);
        break;
      case 'mostPotentResidues':
          this.settings!.showMostPotentResidues ? this.addMostPotentResidues() :
            this.closeViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES);
        break;
      }
    }
  }

  _webLogoSelection: type.Selection | null = null;

  get webLogoSelection(): type.Selection {
    const tagSelection = this.df.getTag(`${C.SUFFIXES.WL}${C.TAGS.INVARIANT_MAP_SELECTION}`);
    this._webLogoSelection ??= tagSelection === null && this.positionColumns !== null ?
      initSelection(this.positionColumns) : JSON.parse(tagSelection ?? `{}`);
    return this._webLogoSelection!;
  }

  set webLogoSelection(selection: type.Selection) {
    this._webLogoSelection = selection;
    this.df.setTag(`${C.SUFFIXES.WL}${C.TAGS.INVARIANT_MAP_SELECTION}`, JSON.stringify(selection));
    this.fireBitsetChanged(null);
    this.analysisView.grid.invalidate();
  }

  get positionColumns(): DG.Column<string>[] | null {
    const positionColumns = wu(this.df.columns.byTags({[C.TAGS.POSITION_COL]: `${true}`})).toArray();
    if (positionColumns.length === 0)
      return null;
    return positionColumns;
  }

  get id(): string {
    return this.df.getTag(DG.TAGS.ID)!;
  }

  get alphabet(): string {
    const col = this.df.getCol(this.settings!.sequenceColumnName);
    return col.getTag(bioTAGS.alphabet);
  }

  static getInstance(dataFrame: DG.DataFrame): PeptidesModel {
    if (dataFrame.columns.contains(C.COLUMNS_NAMES.ACTIVITY_SCALED) &&
      !dataFrame.columns.contains(C.COLUMNS_NAMES.ACTIVITY))
      dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED).name = C.COLUMNS_NAMES.ACTIVITY;

    dataFrame.temp[PeptidesModel.modelName] ??= new PeptidesModel(dataFrame);
    return dataFrame.temp[PeptidesModel.modelName] as PeptidesModel;
  }

  modifyWebLogoSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {
    shiftPressed: false,
    ctrlPressed: false,
  }, notify: boolean = true): void {
    if (notify)
      this.webLogoSelection = modifySelection(this.webLogoSelection, monomerPosition, options);
    else
      this._webLogoSelection = modifySelection(this.webLogoSelection, monomerPosition, options);
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
    const filteredTitlePart = trueModel.df.filter.anyFalse ? ` among ${trueModel.df.filter.trueCount} filtered` :
      '';
    const getSelectionString = (selection: type.Selection): string => {
      const selectedMonomerPositions: string[] = [];
      for (const [pos, monomerList] of Object.entries(selection)) {
        for (const monomer of monomerList)
          selectedMonomerPositions.push(`${pos}:${monomer}`);
      }
      return selectedMonomerPositions.join(', ');
    };

    // Logo Summary Table viewer selection overview
    const trueLSTViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    const selectionDescription: HTMLElement[] = [];
    const selectedClusters: string = (trueLSTViewer === null ? [] :
      trueLSTViewer.clusterSelection[CLUSTER_TYPE.ORIGINAL].concat(trueLSTViewer.clusterSelection[CLUSTER_TYPE.CUSTOM]))
      .join(', ');
    if (selectedClusters.length !== 0) {
      selectionDescription.push(ui.h1('Logo summary table selection'));
      selectionDescription.push(ui.divText(`Selected clusters: ${selectedClusters}`));
    }

    // Monomer-Position viewer selection overview
    const trueMPViewer = trueModel.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    const selectedMonomerPositions = getSelectionString(trueMPViewer?.invariantMapSelection ?? {});
    const selectedMutationCliffs = getSelectionString(trueMPViewer?.mutationCliffsSelection ?? {});
    if (selectedMonomerPositions.length !== 0 || selectedMutationCliffs.length !== 0)
      selectionDescription.push(ui.h1('Monomer-Position viewer selection'));
    if (selectedMonomerPositions.length !== 0)
      selectionDescription.push(ui.divText(`Selected monomer-positions: ${selectedMonomerPositions}`));
    if (selectedMutationCliffs.length !== 0)
      selectionDescription.push(ui.divText(`Selected mutation cliffs pairs: ${selectedMutationCliffs}`));

    // Most Potent Residues viewer selection overview
    const trueMPRViewer = trueModel.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    const selectedMPRMonomerPositions = getSelectionString(trueMPRViewer?.mutationCliffsSelection ?? {});
    if (selectedMPRMonomerPositions.length !== 0) {
      selectionDescription.push(ui.h1('Most Potent Residues viewer selection'));
      selectionDescription.push(ui.divText(`Selected monomer-positions: ${selectedMPRMonomerPositions}`));
    }

    // WebLogo selection overview
    const selectedMonomers = getSelectionString(trueModel.webLogoSelection);
    if (selectedMonomers.length !== 0) {
      selectionDescription.push(ui.h1('WebLogo selection'));
      selectionDescription.push(ui.divText(`Selected monomers: ${selectedMonomers}`));
    }

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
        newView.onmouseover = (ev): void =>
          ui.tooltip.show('Creates a new view from current selection', ev.clientX + 5, ev.clientY + 5);
        if (trueLSTViewer === null)
          return ui.divV([newView]);

        const newCluster = ui.label('New cluster');
        $(newCluster).addClass('d4-link-action');
        newCluster.onclick = (): void => {
          if (trueLSTViewer === null)
            throw new Error('Logo summary table viewer is not found');
          trueLSTViewer.clusterFromSelection();
        };
        newCluster.onmouseover = (ev): void =>
          ui.tooltip.show('Creates a new cluster from selection', ev.clientX + 5, ev.clientY + 5);

        const removeCluster = ui.label('Remove cluster');
        $(removeCluster).addClass('d4-link-action');
        removeCluster.onclick = (): void => {
          const lstViewer = trueModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
          if (lstViewer === null)
            throw new Error('Logo summary table viewer is not found');
          lstViewer.removeCluster();
        };
        removeCluster.onmouseover = (ev): void =>
          ui.tooltip.show('Removes currently selected custom cluster', ev.clientX + 5, ev.clientY + 5);
        removeCluster.style.visibility = trueLSTViewer.clusterSelection[CLUSTER_TYPE.CUSTOM].length === 0 ? 'hidden' :
          'visible';

        return ui.divV([newView, newCluster, removeCluster]);
      }, true);
    }
    const table = trueModel.df.filter.anyFalse ? trueModel.df.clone(trueModel.df.filter, null, true) : trueModel.df;

    // Get the source of the bitset change and find viewers that share the same parameters as source
    let requestSource: SARViewer | LogoSummaryTable | PeptidesSettings | null = trueModel.settings;
    const viewers: (PeptideViewer | PeptidesSettings | null)[] = [trueMPViewer, trueMPRViewer, trueLSTViewer]
      .filter((v) => {
        if (v === null)
          return false;
        if (v.type !== this.accordionSource)
          return true;
        requestSource = v;
        return false;
      });

    if (requestSource === null)
      throw new Error('PeptidesError: Model is the source of accordion but is not initialized');
    if (requestSource !== trueModel.settings)
      viewers.push(trueModel.settings);

    const notEmpty = (v: PeptideViewer | PeptidesSettings | null): v is PeptideViewer | PeptidesSettings =>
      v !== null && areParametersEqual(requestSource!, v) && (v !== trueModel.settings || trueModel.isInitialized);
    const panelDataSources = viewers.filter(notEmpty);
    panelDataSources.push(requestSource);
    const combinedBitset: DG.BitSet | null = DG.BitSet.create(table.rowCount);
    for (const panelDataSource of panelDataSources) {
      const bitset =
        (panelDataSource === this.settings) ? getSelectionBitset(this.webLogoSelection, this.monomerPositionStats!) :
          (panelDataSource instanceof LogoSummaryTable) ?
            getSelectionBitset(panelDataSource.clusterSelection, panelDataSource.clusterStats) :
            (panelDataSource instanceof SARViewer) ?
              getSelectionBitset(panelDataSource.mutationCliffsSelection,
                mutationCliffsToMaskInfo(panelDataSource.mutationCliffs ?? new Map(), table.rowCount)) :
              null;
      if (bitset !== null)
        combinedBitset.or(bitset);
      if (panelDataSource instanceof MonomerPosition) {
        const invariantMapSelectionBitset = getSelectionBitset(panelDataSource.invariantMapSelection,
          panelDataSource.monomerPositionStats);
        if (invariantMapSelectionBitset !== null)
          combinedBitset.or(invariantMapSelectionBitset);
      }
    }

    const sarViewer = requestSource as any as SARViewer | LogoSummaryTable;
    if (requestSource !== trueModel.settings && !(sarViewer instanceof LogoSummaryTable) &&
      sarViewer.mutationCliffs !== null) {
      // MC and Selection are left
      acc.addPane('Mutation Cliffs pairs', () => mutationCliffsWidget(trueModel.df, {
        mutationCliffs: sarViewer.mutationCliffs!, mutationCliffsSelection: sarViewer.mutationCliffsSelection,
        gridColumns: trueModel.analysisView.grid.columns, sequenceColumnName: sarViewer.sequenceColumnName,
        positionColumns: sarViewer.positionColumns, activityCol: sarViewer.getScaledActivityColumn(),
      }).root, true);
    }
    const isModelSource = requestSource === trueModel.settings;
    let filteredActivityCol = isModelSource ? this.getScaledActivityColumn()! :
      (requestSource as unknown as PeptideViewer).getScaledActivityColumn();
    if (trueModel.df.filter.anyFalse) {
      filteredActivityCol = DG.DataFrame.fromColumns([filteredActivityCol]).clone(trueModel.df.filter)
        .getCol(filteredActivityCol.name) as DG.Column<number>;
    }
    acc.addPane('Distribution', () => getDistributionWidget(table, {
      peptideSelection: combinedBitset,
      columns: isModelSource ? trueModel.settings!.columns ?? {} :
        (requestSource as SARViewer | LogoSummaryTable).getAggregationColumns(),
      activityCol: filteredActivityCol,
    }), true);
    const areObjectsEqual = (o1?: AggregationColumns | null, o2?: AggregationColumns | null): boolean => {
      if (o1 == null || o2 == null)
        return false;
      for (const [key, value] of Object.entries(o1)) {
        if (value !== o2[key])
          return false;
      }
      return true;
    };
    acc.addPane('Selection', () => getSelectionWidget(trueModel.df, {
      positionColumns: isModelSource ? this.positionColumns! :
        (requestSource as SARViewer | LogoSummaryTable).positionColumns,
      columns: isModelSource ? trueModel.settings!.columns ?? {} :
        (requestSource as SARViewer | LogoSummaryTable).getAggregationColumns(),
      activityColumn: isModelSource ? this.getScaledActivityColumn()! :
        (requestSource as SARViewer | LogoSummaryTable).getScaledActivityColumn(),
      gridColumns: trueModel.analysisView.grid.columns,
      colorPalette: pickUpPalette(this.df.getCol(isModelSource ? this.settings!.sequenceColumnName :
        (requestSource as SARViewer | LogoSummaryTable).sequenceColumnName)),
      tableSelection: trueModel.getCombinedSelection(),
      isAnalysis: this.settings !== null && (isModelSource ||
        areObjectsEqual(this.settings.columns, (requestSource as PeptideViewer).getAggregationColumns())),
    }), true);

    return acc;
  }

  getScaledActivityColumn(isFiltered: boolean = false): DG.Column<number> | null {
    this._scaledActivityColumn ??= this.df.col(C.COLUMNS_NAMES.ACTIVITY);
    if (isFiltered && this._scaledActivityColumn !== null) {
      return DG.DataFrame.fromColumns([this._scaledActivityColumn]).clone(this.df.filter)
        .getCol(this._scaledActivityColumn.name) as DG.Column<number>;
    }
    return this._scaledActivityColumn as DG.Column<number> | null;
  }

  updateGrid(): void {
    this.joinDataFrames();
    this.createScaledCol();
    this.webLogoBounds = {};

    const cellRendererOptions: CR.CellRendererOptions = {
      selectionCallback: (monomerPosition: type.SelectionItem, options: type.SelectionOptions): type.Selection =>
        modifySelection(this.webLogoSelection, monomerPosition, options),
      unhighlightCallback: (): void => this.unhighlight(),
      colorPalette: pickUpPalette(this.df.getCol(this.settings!.sequenceColumnName)),
      webLogoBounds: this.webLogoBounds,
      cachedWebLogoTooltip: this.cachedWebLogoTooltip,
      highlightCallback: (mp: type.SelectionItem, df: DG.DataFrame, mpStats: MonomerPositionStats): void =>
        highlightMonomerPosition(mp, df, mpStats),
      isSelectionTable: false,
      headerSelectedMonomers: this.webLogoSelectedMonomers,
    };
    if (this.monomerPositionStats === null || this.positionColumns === null)
      throw new Error('PeptidesError: Could not updage grid: monomerPositionStats or positionColumns are null');
    CR.setWebLogoRenderer(this.analysisView.grid, this.monomerPositionStats, this.positionColumns,
      this.getScaledActivityColumn()!, cellRendererOptions);
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

  joinDataFrames(): void {
    // append splitSeqDf columns to source table and make sure columns are not added more than once
    const name = this.df.name;
    const cols = this.df.columns;
    const splitSeqDf = splitAlignedSequences(this.df.getCol(this.settings!.sequenceColumnName));
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
    const scaledCol = scaleActivity(this.df.getCol(this.settings!.activityColumnName),
      this.settings!.activityScaling);
    //TODO: make another func
    this.df.columns.replace(COLUMNS_NAMES.ACTIVITY, scaledCol);

    sourceGrid.columns.setOrder([scaledCol.name]);
  }

  unhighlight(): void {
    if (!this.isHighlighting)
      return;
    this.df.rows.highlight(null);
    this.isHighlighting = false;
  }

  setTooltips(): void {
    this.analysisView.grid.onCellTooltip((cell, x, y) => {
      if (cell.isColHeader && cell.tableColumn?.semType === C.SEM_TYPES.MONOMER)
        return true;
      if (!(cell.isTableCell && cell.tableColumn?.semType === C.SEM_TYPES.MONOMER))
        return false;

      showMonomerTooltip(cell.cell.value, x, y);
      return true;
    });
  }

  getCombinedSelection(): DG.BitSet {
    const combinedSelection = new BitArray(this.df.rowCount, false);
    // Invariant map selection
    const addInvariantMapSelection = (selection: type.Selection, stats: MonomerPositionStats | null): void => {
      for (const [position, monomerList] of Object.entries(selection)) {
        for (const monomer of monomerList) {
          const positionStats = stats?.[position];
          if (typeof positionStats === 'undefined')
            continue;

          const monomerPositionStats = positionStats[monomer];
          if (typeof monomerPositionStats === 'undefined')
            continue;

          combinedSelection.or(monomerPositionStats.mask);
        }
      }
    };

    addInvariantMapSelection(this.webLogoSelection, this.monomerPositionStats);
    const mpViewer = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    addInvariantMapSelection(mpViewer?.invariantMapSelection ?? {}, mpViewer?.monomerPositionStats ?? null);

    // Mutation cliffs selection
    const addMutationCliffsSelection = (selection: type.Selection, mc: type.MutationCliffs | null): void => {
      for (const [position, monomerList] of Object.entries(selection)) {
        for (const monomer of monomerList) {
          const substitutions = mc?.get(monomer)?.get(position) ?? null;
          if (substitutions === null)
            continue;
          for (const [key, value] of substitutions.entries()) {
            combinedSelection.setTrue(key);
            for (const v of value)
              combinedSelection.setTrue(v);
          }
        }
      }
    };
    addMutationCliffsSelection(mpViewer?.mutationCliffsSelection ?? {}, mpViewer?.mutationCliffs ?? null);

    const mprViewer = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    addMutationCliffsSelection(mprViewer?.mutationCliffsSelection ?? {}, mprViewer?.mutationCliffs ?? null);

    // Cluster selection
    const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    for (const clustType of Object.keys(lstViewer?.clusterSelection ?? {})) {
      for (const clust of lstViewer!.clusterSelection[clustType] ?? []) {
        const clusterStats = lstViewer!.clusterStats[clustType as CLUSTER_TYPE][clust];
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

  fireBitsetChanged(source: VIEWER_TYPE | null, fireFilterChanged: boolean = false): void {
    this.accordionSource = source;
    if (!this.isBitsetChangedInitialized)
      this.setBitsetCallback();
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
    this.webLogoSelectedMonomers = calculateSelected(this.df);
  }

  setGridProperties(props?: DG.IGridLookSettings): void {
    const sourceGrid = this.analysisView.grid;
    const sourceGridProps = sourceGrid.props;
    sourceGridProps.allowColSelection = props?.allowColSelection ?? false;
    sourceGridProps.allowEdit = props?.allowEdit ?? false;
    sourceGridProps.showCurrentRowIndicator = props?.showCurrentRowIndicator ?? false;
    this.df.temp[C.EMBEDDING_STATUS] = false;
    const positionCols = this.positionColumns;
    if (positionCols === null)
      throw new Error('PeptidesError: Could not set grid properties: positionColumns are null');
    let maxWidth = 10;
    const canvasContext = sourceGrid.canvas.getContext('2d');
    if (canvasContext === null)
      throw new Error('PeptidesError: Could not set grid properties: canvas context is null');
    for (const positionCol of positionCols) {
      // Longest category
      const maxCategory = monomerToShort(positionCol.categories.reduce((a, b) => a.length > b.length ? a : b), 6);
      // Measure text width of longest category
      const width = Math.ceil(canvasContext.measureText(maxCategory).width);
      maxWidth = Math.max(maxWidth, width);
    }
    setTimeout(() => {
      for (const positionCol of positionCols) {
        const gridCol = sourceGrid.col(positionCol.name);
        if (gridCol === null)
          throw new Error(`PeptidesError: Could not set column width: grid column '${positionCol.name}' is null`);
        gridCol.width = maxWidth + 15;
      }
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
      const pepColValues: string[] = this.df.getCol(this.settings!.sequenceColumnName).toList();
      this._dm ??= new DistanceMatrix(await createDistanceMatrixWorker(pepColValues, StringMetricsNames.Levenshtein));
      const leafCol = this.df.col('~leaf-id') ?? this.df.columns.addNewString('~leaf-id').init((i) => i.toString());
      const treeHelper: ITreeHelper = getTreeHelperInstance();
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
  init(settings: type.PeptidesSettings): void {
    if (this.isInitialized)
      return;
    this.settings = settings;
    this.isInitialized = true;

    if (!this.isRibbonSet && this.df.getTag(C.TAGS.MULTIPLE_VIEWS) !== '1') {
      //TODO: don't pass model, pass parameters instead
      const settingsButton = ui.iconFA('wrench', () => getSettingsDialog(this), 'Peptides analysis settings');
      this.analysisView.setRibbonPanels([[settingsButton]], false);
      this.isRibbonSet = true;
      this.updateGrid();
    }

    this.subs.push(grok.events.onAccordionConstructed.subscribe((acc) => {
      if (!(grok.shell.o instanceof DG.SemanticValue || (grok.shell.o instanceof DG.Column &&
        this.df.columns.toList().includes(grok.shell.o))))
        return;

      const actionsPane = acc.getPane('Actions');

      const actionsHost = $(actionsPane.root).find('.d4-flex-col');
      const calculateIdentity = ui.label('Calculate identity');
      calculateIdentity.classList.add('d4-link-action');
      ui.tooltip.bind(calculateIdentity,
        'Adds a column with fractions of matching monomers against sequence in the current row');
      calculateIdentity.onclick = (): void => {
        const seqCol = this.df.getCol(this.settings!.sequenceColumnName);
        calculateScores(this.df, seqCol, seqCol.get(this.df.currentRowIdx), SCORE.IDENTITY)
          .then((col: DG.Column<number>) => col.setTag(C.TAGS.IDENTITY_TEMPLATE, seqCol.get(this.df.currentRowIdx)))
          .catch((e) => _package.logger.debug(e));
      };
      actionsHost.append(ui.span([calculateIdentity], 'd4-markdown-row'));

      const calculateSimilarity = ui.label('Calculate similarity');
      calculateSimilarity.classList.add('d4-link-action');
      ui.tooltip.bind(calculateSimilarity,
        'Adds a column with sequence similarity scores against sequence in the current row');
      calculateSimilarity.onclick = (): void => {
        const seqCol = this.df.getCol(this.settings!.sequenceColumnName);
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

    this.fireBitsetChanged(null, true);

    this.analysisView.grid.invalidate();
  }

  findViewer(viewerType: VIEWER_TYPE): DG.Viewer | null {
    return wu(this.analysisView.viewers).find((v) => v.type === viewerType) || null;
  }

  async addLogoSummaryTable(viewerProperties?: ILogoSummaryTable): Promise<void> {
    viewerProperties ??= {
      sequenceColumnName: this.settings!.sequenceColumnName,
      clustersColumnName: wu(this.df.columns.categorical).next().value,
      activityColumnName: this.settings!.activityColumnName, activityScaling: this.settings!.activityScaling,
    };
    const logoSummaryTable = await this.df.plot
      .fromType(VIEWER_TYPE.LOGO_SUMMARY_TABLE, viewerProperties) as LogoSummaryTable;
    this.analysisView.dockManager.dock(logoSummaryTable, DG.DOCK_TYPE.RIGHT, null, VIEWER_TYPE.LOGO_SUMMARY_TABLE);

    logoSummaryTable.viewerGrid.invalidate();
  }

  async addMonomerPosition(viewerProperties?: ISARViewer): Promise<void> {
    viewerProperties ??= {
      maxMutations: 1,
      activityScaling: this.settings!.activityScaling,
      activityColumnName: this.settings!.activityColumnName,
      sequenceColumnName: this.settings!.sequenceColumnName,
      minActivityDelta: 0,
    };
    const monomerPosition = await this.df.plot
      .fromType(VIEWER_TYPE.MONOMER_POSITION, viewerProperties) as MonomerPosition;
    const mostPotentResidues = this.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    const dm = this.analysisView.dockManager;
    const [dockType, refNode, ratio] = mostPotentResidues === null ? [DG.DOCK_TYPE.DOWN, null, undefined] :
      [DG.DOCK_TYPE.LEFT, this.findViewerNode(VIEWER_TYPE.MOST_POTENT_RESIDUES), 0.7];
    dm.dock(monomerPosition, dockType, refNode, VIEWER_TYPE.MONOMER_POSITION, ratio);
  }

  async addMostPotentResidues(viewerProperties?: ISARViewer): Promise<void> {
    viewerProperties ??= {
      activityScaling: this.settings!.activityScaling,
      activityColumnName: this.settings!.activityColumnName,
      sequenceColumnName: this.settings!.sequenceColumnName,
      minActivityDelta: 0,
      maxMutations: 1,
    };
    const mostPotentResidues =
      await this.df.plot.fromType(VIEWER_TYPE.MOST_POTENT_RESIDUES, viewerProperties) as MostPotentResidues;
    const monomerPosition = this.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    const dm = this.analysisView.dockManager;
    const [dockType, refNode, ratio] = monomerPosition === null ? [DG.DOCK_TYPE.DOWN, null, undefined] :
      [DG.DOCK_TYPE.RIGHT, this.findViewerNode(VIEWER_TYPE.MONOMER_POSITION), 0.3];
    dm.dock(mostPotentResidues, dockType, refNode, VIEWER_TYPE.MOST_POTENT_RESIDUES, ratio);
  }

  createNewView(): string {
    const rowMask = this.getVisibleSelection();
    const newDf = this.df.clone(rowMask);
    for (const [tag, value] of newDf.tags)
      newDf.setTag(tag, tag === C.TAGS.SETTINGS ? value : '');

    newDf.name = 'Peptides Multiple Views';
    newDf.setTag(C.TAGS.MULTIPLE_VIEWS, '1');

    const view = grok.shell.addTableView(newDf);
    const lstViewer = this.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    if (lstViewer != null) {
      view.addViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE, {
        [LST_PROPERTIES.SEQUENCE]: lstViewer.sequenceColumnName,
        [LST_PROPERTIES.CLUSTERS]: lstViewer?.clustersColumnName,
      });
    }

    return newDf.getTag(DG.TAGS.ID)!;
  }

  async addSequenceSpace(): Promise<void> {
    let seqCol = this.df.getCol(this.settings!.sequenceColumnName!);
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
    const seqSpaceParams: {
      table: DG.DataFrame,
      molecules: DG.Column,
      methodName: DimReductionMethods,
      similarityMetric: BitArrayMetrics | MmDistanceFunctionsNames,
      plotEmbeddings: boolean,
      sparseMatrixThreshold?: number,
      options?: (IUMAPOptions | ITSNEOptions) & Options
    } =
      {
        table: this.df, molecules: seqCol,
        methodName: DimReductionMethods.UMAP, similarityMetric: MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH,
        plotEmbeddings: true, sparseMatrixThreshold: 0.3, options: {'bypassLargeDataWarning': true},
      };

    // Use counter to unsubscribe when 2 columns are hidden
    let counter = 0;
    const columnAddedSub = this.df.onColumnsAdded.subscribe((colArgs: DG.ColumnsArgs) => {
      for (const col of colArgs.columns) {
        if (col.name.startsWith('Embed_')) {
          const gridCol = this.analysisView.grid.col(col.name);
          if (gridCol == null)
            continue;
          gridCol.visible = false;
          counter++;
        }
      }
      if (counter === 2)
        columnAddedSub.unsubscribe();
    });

    const seqSpaceViewer: DG.ScatterPlotViewer | undefined =
      await grok.functions.call('Bio:sequenceSpaceTopMenu', seqSpaceParams);
    if (!(seqSpaceViewer instanceof DG.ScatterPlotViewer))
      return;
    seqSpaceViewer.props.colorColumnName = this.getScaledActivityColumn()!.name;
    seqSpaceViewer.props.showXSelector = false;
    seqSpaceViewer.props.showYSelector = false;
  }
}
