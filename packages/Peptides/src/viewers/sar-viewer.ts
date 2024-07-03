/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import {COLUMN_NAME} from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import wu from 'wu';
import * as type from '../utils/types';
import {MutationCliffs, PeptidesSettings, SelectionItem} from '../utils/types';
import {
  AggregationColumns,
  getAggregatedColName,
  getAggregatedColumnValuesFromDf,
  getAggregatedValue,
  MonomerPositionStats,
  PositionStats,
  StatsItem,
} from '../utils/statistics';
import {_package} from '../package';
import {showTooltip} from '../utils/tooltips';
import {calculateCliffsStatistics, calculateMonomerPositionStatistics, findMutations, MutationCliffsOptions} from '../utils/algorithms';
import {
  debounce,
  extractColInfo,
  getTotalAggColumns,
  highlightMonomerPosition,
  initSelection,
  isApplicableDataframe,
  isSelectionEmpty,
  modifySelection,
  scaleActivity,
} from '../utils/misc';
import {splitAlignedSequences} from '@datagrok-libraries/bio/src/utils/splitter';
import {LogoSummaryTable} from './logo-summary';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';

export enum SELECTION_MODE {
  MUTATION_CLIFFS = 'Mutation Cliffs',
  INVARIANT_MAP = 'Invariant Map',
}

export enum SAR_PROPERTIES {
  SEQUENCE = 'sequence',
  ACTIVITY = 'activity',
  ACTIVITY_SCALING = 'activityScaling',
  TARGET = 'target',
  TARGET_CATEGORY = 'targetCategory',
  MIN_ACTIVITY_DELTA = 'minActivityDelta',
  MAX_MUTATIONS = 'maxMutations',
  COLUMNS = 'columns',
  AGGREGATION = 'aggregation',
  ACTIVITY_TARGET = 'activityTarget',
}

export enum MONOMER_POSITION_PROPERTIES {
  COLOR = 'color',
  COLOR_AGGREGATION = 'colorAggregation',
}

export enum PROPERTY_CATEGORIES {
  GENERAL = 'General',
  INVARIANT_MAP = 'Invariant Map',
  MUTATION_CLIFFS = 'Mutation Cliffs',
  AGGREGATION = 'Aggregation',
}

const MUTATION_CLIFFS_CELL_WIDTH = 40;
const AAR_CELL_WIDTH = 30;

export interface ISARViewer {
  sequenceColumnName: string;
  activityColumnName: string;
  activityScaling: C.SCALING_METHODS;
  minActivityDelta: number;
  maxMutations: number;
  activityTarget: C.ACTIVITY_TARGET;
}

/** Abstract class for MonomerPosition and MostPotentResidues viewers. */
export abstract class SARViewer extends DG.JsViewer implements ISARViewer {
  keyPressed: boolean = false;
  sequenceColumnName: string;
  activityColumnName: string;
  activityScaling: C.SCALING_METHODS;
  columns: string[];
  aggregation: string;
  targetColumnName: string;
  targetCategory: string;
  minActivityDelta: number;
  maxMutations: number;
  _scaledActivityColumn: DG.Column | null = null;
  doRender: boolean = true;
  activityTarget: C.ACTIVITY_TARGET;
  mutationCliffsDebouncer: (
    activityArray: type.RawData, monomerInfoArray: type.RawColumn[], options?: MutationCliffsOptions
    ) => Promise<type.MutationCliffs>;

  /** Sets common properties for inheritor viewers. */
  protected constructor() {
    super();
    // General properties
    this.sequenceColumnName = this.column(SAR_PROPERTIES.SEQUENCE,
      {category: PROPERTY_CATEGORIES.GENERAL, semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.activityColumnName = this.column(SAR_PROPERTIES.ACTIVITY,
      {category: PROPERTY_CATEGORIES.GENERAL, nullable: false});
    this.activityScaling = this.string(SAR_PROPERTIES.ACTIVITY_SCALING, C.SCALING_METHODS.NONE,
      {category: PROPERTY_CATEGORIES.GENERAL, choices: Object.values(C.SCALING_METHODS), nullable: false},
    ) as C.SCALING_METHODS;
    this.activityTarget = this.string(SAR_PROPERTIES.ACTIVITY_TARGET, C.ACTIVITY_TARGET.HIGH,
      {category: PROPERTY_CATEGORIES.GENERAL, choices: Object.values(C.ACTIVITY_TARGET), nullable: false},
    ) as C.ACTIVITY_TARGET;
    // Mutation Cliffs properties
    this.targetColumnName = this.column(SAR_PROPERTIES.TARGET, {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS});
    this.targetCategory = this.string(SAR_PROPERTIES.TARGET_CATEGORY, null,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, choices: []});
    this.minActivityDelta = this.float(SAR_PROPERTIES.MIN_ACTIVITY_DELTA, 0,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, min: 0, max: 100});
    this.maxMutations = this.int(SAR_PROPERTIES.MAX_MUTATIONS, 1,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, min: 1, max: 20});

    this.columns = this.columnList(SAR_PROPERTIES.COLUMNS, [], {category: PROPERTY_CATEGORIES.AGGREGATION});
    this.aggregation = this.string(SAR_PROPERTIES.AGGREGATION, DG.AGG.AVG,
      {category: PROPERTY_CATEGORIES.AGGREGATION, choices: C.AGGREGATION_TYPES});

    this.mutationCliffsDebouncer = debounce(
      async (activityArray: type.RawData, monomerInfoArray: type.RawColumn[], options?: MutationCliffsOptions) => {
        return await findMutations(activityArray, monomerInfoArray, options);
      });
  }

  _viewerGrid: DG.Grid | null = null;

  /**
   * Returns SARViewer grid. Creates a new one if it is null.
   * @return - SARViewer grid.
   */
  get viewerGrid(): DG.Grid {
    this._viewerGrid ??= this.createViewerGrid();
    return this._viewerGrid;
  }

  /**
   * Returns sequence column alphabet.
   * @return - sequence column alphabet.
   */
  get alphabet(): string {
    const col = this.dataFrame.getCol(this.sequenceColumnName);
    return col.getTag(bioTAGS.alphabet) ?? ALPHABET.UN;
  }

  /**
   * Returns PeptidesModel instance that belongs to the attached dataframe.
   * @return - PeptidesModel instance.
   */
  get model(): PeptidesModel {
    return PeptidesModel.getInstance(this.dataFrame);
  }

  _positionColumns: DG.Column<string>[] | null = null;

  /**
   * Returns position columns. If position columns are not attached to the viewer, it tries to get them from
   * other viewers or analysis (given that relevant parameters are the same), or creates own position columns.
   * @return - position columns.
   */
  get positionColumns(): DG.Column<string>[] {
    if (this._positionColumns != null)
      return this._positionColumns;


    const getSharedPositionColumns = (viewerType: VIEWER_TYPE): DG.Column<string>[] | null => {
      const viewer = this.model.findViewer(viewerType) as SARViewer | LogoSummaryTable | null;
      if (this.sequenceColumnName === viewer?.sequenceColumnName)
        return viewer._positionColumns;


      return null;
    };

    if (this.model.positionColumns != null && this.sequenceColumnName === this.model.settings?.sequenceColumnName)
      this._positionColumns = this.model.positionColumns;
    else if (this instanceof MonomerPosition)
      this._positionColumns = getSharedPositionColumns(VIEWER_TYPE.MOST_POTENT_RESIDUES);
    else if (this instanceof MostPotentResidues)
      this._positionColumns = getSharedPositionColumns(VIEWER_TYPE.MONOMER_POSITION);


    this._positionColumns ??= getSharedPositionColumns(VIEWER_TYPE.LOGO_SUMMARY_TABLE) ??
      splitAlignedSequences(this.dataFrame.getCol(this.sequenceColumnName)).columns.toList();
    return this._positionColumns!;
  }

  _monomerPositionStats: MonomerPositionStats | null = null;

  /**
   * Gets monomer-position statistics. If monomer-position statistics are not attached to the viewer, it tries to get
   * them from other viewers or analysis (given that relevant parameters are the same), or calculates own statistics.
   * @return - monomer-position statistics.
   */
  get monomerPositionStats(): MonomerPositionStats {
    if (this._monomerPositionStats != null)
      return this._monomerPositionStats;


    const isMonomerPositionStatsEqual = (other: SARViewer | PeptidesSettings | null): boolean =>
      this.sequenceColumnName === other?.sequenceColumnName &&
      this.activityColumnName === other?.activityColumnName &&
      this.activityScaling === other?.activityScaling;

    const getSharedStats = (viewerType: VIEWER_TYPE): MonomerPositionStats | null => {
      const viewer = this.model.findViewer(viewerType) as SARViewer | null;
      if (isMonomerPositionStatsEqual(viewer))
        return viewer!._monomerPositionStats;


      return null;
    };

    if (this.model.monomerPositionStats !== null && isMonomerPositionStatsEqual(this.model.settings))
      this._monomerPositionStats = this.model.monomerPositionStats;
    else if (this instanceof MonomerPosition)
      this._monomerPositionStats = getSharedStats(VIEWER_TYPE.MOST_POTENT_RESIDUES);
    else if (this instanceof MostPotentResidues)
      this._monomerPositionStats = getSharedStats(VIEWER_TYPE.MONOMER_POSITION);


    this._monomerPositionStats ??= calculateMonomerPositionStatistics(this.getScaledActivityColumn(),
      this.dataFrame.filter, this.positionColumns);
    return this._monomerPositionStats;
  }

  _mutationCliffs: type.MutationCliffs | null = null;
  _mutationCliffStats: type.MutationCliffStats | null = null;

  /**
   * Gets mutation cliffs. If mutation cliffs are not attached to the viewer, it tries to get them from other viewers,
   * or calculates its own.
   * @return - mutation cliffs.
   */
  get mutationCliffs(): type.MutationCliffs | null {
    if (this._mutationCliffs != null)
      return this._mutationCliffs;


    const isMutationCliffsEqual = (v1: SARViewer, v2: SARViewer | null): boolean =>
      v1.sequenceColumnName === v2?.sequenceColumnName &&
      v1.activityColumnName === v2.activityColumnName &&
      v1.activityScaling === v2.activityScaling &&
      v1.targetColumnName === v2?.targetColumnName &&
      v1.targetCategory === v2?.targetCategory &&
      v1.minActivityDelta === v2?.minActivityDelta &&
      v1.maxMutations === v2?.maxMutations;

    const getSharedMutationCliffs = (viewerType: VIEWER_TYPE): type.MutationCliffs | null => {
      const viewer = this.model.findViewer(viewerType) as SARViewer | null;
      if (isMutationCliffsEqual(this, viewer))
        return viewer!._mutationCliffs;


      return null;
    };

    if (this instanceof MonomerPosition)
      this._mutationCliffs = getSharedMutationCliffs(VIEWER_TYPE.MOST_POTENT_RESIDUES);
    else if (this instanceof MostPotentResidues)
      this._mutationCliffs = getSharedMutationCliffs(VIEWER_TYPE.MONOMER_POSITION);


    return this._mutationCliffs;
  }

  /**
   * Sets mutation cliffs.
   * @param mc - mutation cliffs to set.
   */
  set mutationCliffs(mc: type.MutationCliffs) {
    this._mutationCliffs = mc;
    this.viewerGrid.invalidate();
  }

  get cliffStats(): type.MutationCliffStats | null {
    if (!this._mutationCliffStats && this.mutationCliffs && this.activityColumnName &&
      this.sequenceColumnName && this.dataFrame.col(this.activityColumnName)) {
      const activityCol = this.dataFrame.col(this.activityColumnName)!;
      const activityData = activityCol.getRawData();
      this._mutationCliffStats = calculateCliffsStatistics(this.mutationCliffs, activityData);
    }
    return this._mutationCliffStats;
  }
  set cliffStats(stats: type.MutationCliffStats | null) {
    this._mutationCliffStats = stats;
    this.viewerGrid.invalidate;
  }

  _mutationCliffsSelection: type.Selection | null = null;

  /**
   * Gets mutation cliffs selection.
   * @return - mutation cliffs selection.
   */
  get mutationCliffsSelection(): type.Selection {
    const tagSuffix = this instanceof MonomerPosition ? C.SUFFIXES.MP : C.SUFFIXES.MPR;
    const tagSelection = this.dataFrame.getTag(`${tagSuffix}${C.TAGS.MUTATION_CLIFFS_SELECTION}`);
    this._mutationCliffsSelection ??= tagSelection === null ? initSelection(this.positionColumns) :
      JSON.parse(tagSelection);
    return this._mutationCliffsSelection!;
  }

  /**
   * Sets mutation cliffs selection.
   * @param selection - selection to set.
   */
  set mutationCliffsSelection(selection: type.Selection) {
    this._mutationCliffsSelection = selection;
    const tagSuffix = this instanceof MonomerPosition ? C.SUFFIXES.MP : C.SUFFIXES.MPR;
    this.dataFrame.setTag(`${tagSuffix}${C.TAGS.MUTATION_CLIFFS_SELECTION}`, JSON.stringify(selection));
    this.model.fireBitsetChanged(this instanceof MonomerPosition ? VIEWER_TYPE.MONOMER_POSITION :
      VIEWER_TYPE.MOST_POTENT_RESIDUES);

    const mpViewer = this.model.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition | null;
    mpViewer?.viewerGrid.invalidate();
    const mprViewer = this.model.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues | null;
    mprViewer?.viewerGrid.invalidate();

    this.model.analysisView.grid.invalidate();
  }

  /**
   * Gets scaled activity column.
   * @param isFiltered - flag indicating if only filtered rows should be taken into account.
   * @return - scaled activity column.
   */
  getScaledActivityColumn(isFiltered: boolean = false): DG.Column<number> {
    if (this.model.settings?.activityColumnName === this.activityColumnName &&
      this.model.settings?.activityScaling === this.activityScaling)
      this._scaledActivityColumn = this.model.getScaledActivityColumn(isFiltered);


    this._scaledActivityColumn ??= scaleActivity(this.dataFrame.getCol(this.activityColumnName),
      this.activityScaling);
    if (isFiltered) {
      return DG.DataFrame.fromColumns([this._scaledActivityColumn]).clone(this.dataFrame.filter)
        .getCol(this._scaledActivityColumn.name) as DG.Column<number>;
    }
    return this._scaledActivityColumn as DG.Column<number>;
  }

  /**
   * Modifies mutation cliffs selection. If shift and ctrl keys are both pressed, it removes mutation cliffs from
   * selection. If only shift key is pressed, it adds mutation cliffs to selection. If only ctrl key is pressed, it
   * changes mutation cliffs presence in selection. If none of the keys is pressed, it sets the mutation cliffs as the
   * only selected one.
   * @param monomerPosition - monomer-position to modify selection with.
   * @param options - selection options.
   * @param notify - flag indicating if bitset changed event should fire.
   */
  modifyMutationCliffsSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {
    shiftPressed: false,
    ctrlPressed: false,
  }, notify: boolean = true): void {
    if (notify)
      this.mutationCliffsSelection = modifySelection(this.mutationCliffsSelection, monomerPosition, options);
    else
      this._mutationCliffsSelection = modifySelection(this.mutationCliffsSelection, monomerPosition, options);
  }

  /**
   * Processes property changes.
   * @param property - changed property.
   */
  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.doRender = true;
    switch (property.name) {
    case `${SAR_PROPERTIES.SEQUENCE}${COLUMN_NAME}`:
      this._positionColumns = null;
      this._monomerPositionStats = null;
      this._mutationCliffs = null;
      this._mutationCliffStats = null;
      this._mutationCliffsSelection = null;
      this._viewerGrid = null;
      break;
    case `${SAR_PROPERTIES.ACTIVITY}${COLUMN_NAME}`:
    case SAR_PROPERTIES.ACTIVITY_SCALING:
      this._monomerPositionStats = null;
      this._mutationCliffs = null;
      this._mutationCliffStats = null;
      this._mutationCliffsSelection = null;
      this._viewerGrid = null;
      this._scaledActivityColumn = null;
      break;
    case `${SAR_PROPERTIES.TARGET}${COLUMN_NAME}`:
    case SAR_PROPERTIES.TARGET_CATEGORY:
    case SAR_PROPERTIES.MIN_ACTIVITY_DELTA:
    case SAR_PROPERTIES.MAX_MUTATIONS:
      this._mutationCliffs = null;
      this._mutationCliffStats = null;
      this._mutationCliffsSelection = null;
      this.doRender = false;
      break;
    case SAR_PROPERTIES.COLUMNS:
    case SAR_PROPERTIES.AGGREGATION:
      if (this instanceof MostPotentResidues)
        this._viewerGrid = null;
      break;
    case SAR_PROPERTIES.ACTIVITY_TARGET:
      if (this instanceof MostPotentResidues || this instanceof MonomerPosition)
        this._viewerGrid = null;
      break;
    }
    if (this._mutationCliffs === null && this.sequenceColumnName && this.activityColumnName)
      this.calculateMutationCliffs().then((mc) => {this.mutationCliffs = mc.cliffs; this.cliffStats = mc.cliffStats;});
  }

  /**
   * Gets a map of columns and aggregations.
   * @return - map of columns and aggregations.
   */
  getAggregationColumns(): AggregationColumns {
    return Object.fromEntries(
      this.columns.map((colName) => [colName, this.aggregation] as [string, DG.AGG]));
  }

  /**
   * Gets total viewer aggregated columns from viewer properties and analysis settings if applicable.
   * @return - total viewer aggregated columns.
   */
  getTotalViewerAggColumns(): [string, DG.AggregationType][] {
    const aggrCols = this.getAggregationColumns();
    return getTotalAggColumns(this.columns, aggrCols, this.model?.settings?.columns);
  }

  /** Creates viewer grid. */
  createViewerGrid(): DG.Grid {
    throw new Error('Not implemented');
  }

  /** Removes all the active subscriptions. */
  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  /** Processes attached table and sets viewer properties. */
  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    if (isApplicableDataframe(this.dataFrame)) {
      this.getProperty(`${SAR_PROPERTIES.SEQUENCE}${COLUMN_NAME}`)
        ?.set(this, this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)!.name);
      this.getProperty(`${SAR_PROPERTIES.ACTIVITY}${COLUMN_NAME}`)
        ?.set(this, wu(this.dataFrame.columns.numerical).next().value.name);
      if (this.mutationCliffs === null && this.sequenceColumnName && this.activityColumnName)
        this.calculateMutationCliffs().then((mc) => {this.mutationCliffs = mc.cliffs; this.cliffStats = mc.cliffStats;});
    } else {
      const msg = 'PeptidesError: dataframe is missing Macromolecule or numeric columns';
      grok.log.error(msg);
      grok.shell.warning(msg);
    }
  }

  /**
   * Calculates Mutation Cliffs.
   * @return - mutation cliffs.
   */
  async calculateMutationCliffs(): Promise<{cliffs: MutationCliffs, cliffStats: type.MutationCliffStats}> {
    const scaledActivityCol: DG.Column<number> = this.dataFrame.getCol(this.activityColumnName);
    //TODO: set categories ordering the same to share compare indexes instead of strings
    const monomerCols: type.RawColumn[] = this.positionColumns.map(extractColInfo);
    const targetCol = this.targetColumnName ? extractColInfo(this.dataFrame.getCol(this.targetColumnName)) : null;

    const options: MutationCliffsOptions = {
      maxMutations: this.maxMutations, minActivityDelta: this.minActivityDelta,
      targetCol, currentTarget: this.targetCategory,
    };
    const activityRawData = scaledActivityCol.getRawData();

    const mutRes = await this.mutationCliffsDebouncer(activityRawData, monomerCols, options);
    const mutStatistics = calculateCliffsStatistics(mutRes, activityRawData);
    return {cliffs: mutRes, cliffStats: mutStatistics};
  }
}

/** Structure-activity relationship viewer */
export class MonomerPosition extends SARViewer {
  colorColumnName: string;
  colorAggregation: string;
  currentGridCell: DG.GridCell | null = null;

  /** Sets MonomerPosition properties. */
  constructor() {
    super();

    const colorChoices = wu(grok.shell.t.columns.numerical).toArray().map((col) => col.name);
    this.colorColumnName = this.column(MONOMER_POSITION_PROPERTIES.COLOR,
      {category: PROPERTY_CATEGORIES.INVARIANT_MAP, choices: colorChoices, nullable: false});
    this.colorAggregation = this.string(MONOMER_POSITION_PROPERTIES.COLOR_AGGREGATION, DG.AGG.AVG,
      {category: PROPERTY_CATEGORIES.INVARIANT_MAP, choices: C.AGGREGATION_TYPES});
  }

  /**
   * Returns viewer name.
   * @return - viewer name.
   */
  get name(): string {
    return VIEWER_TYPE.MONOMER_POSITION;
  }

  /**
   * Returns current viewer operation mode.
   * @return - current viewer operation mode.
   */
  get mode(): SELECTION_MODE {
    return this.dataFrame.getTag(C.TAGS.MONOMER_POSITION_MODE) as SELECTION_MODE ??
      SELECTION_MODE.MUTATION_CLIFFS;
  }

  /**
   * Sets viewer operation mode.
   * @param mode - viewer operation mode to set.
   */
  set mode(mode: SELECTION_MODE) {
    this.dataFrame.setTag(C.TAGS.MONOMER_POSITION_MODE, mode);
    this.viewerGrid.invalidate();
    // setTimeout(() => this.viewerGrid.invalidate(), 300);
  }

  _invariantMapSelection: type.Selection | null = null;

  /**
   * Gets invariant map selection. Initializes it if it is null.
   * @return - invariant map selection.
   */
  get invariantMapSelection(): type.Selection {
    const tagSelection = this.dataFrame.getTag(`${C.SUFFIXES.MP}${C.TAGS.INVARIANT_MAP_SELECTION}`);
    this._invariantMapSelection ??= tagSelection === null ? initSelection(this.positionColumns) :
      JSON.parse(tagSelection);
    return this._invariantMapSelection!;
  }

  /**
   * Sets invariant map selection and notifies the model.
   * @param selection - selection to set.
   */
  set invariantMapSelection(selection: type.Selection) {
    this._invariantMapSelection = selection;
    this.dataFrame.setTag(`${C.SUFFIXES.MP}${C.TAGS.INVARIANT_MAP_SELECTION}`, JSON.stringify(selection));
    this.model.fireBitsetChanged(VIEWER_TYPE.MONOMER_POSITION);
    this.model.analysisView.grid.invalidate();
  }

  /** Processes attached table and sets viewer properties. */
  onTableAttached(): void {
    super.onTableAttached();
    if (isApplicableDataframe(this.dataFrame)) {
      this.getProperty(`${MONOMER_POSITION_PROPERTIES.COLOR}${COLUMN_NAME}`)
        ?.set(this, this.activityColumnName);
    } else {
      const msg = 'PeptidesError: dataframe is missing Macromolecule or numeric columns';
      grok.log.error(msg);
      grok.shell.warning(msg);
    }
    this.render();
  }

  /**
   * Modifies invariant map selection. If shift and ctrl keys are both pressed, it removes invariant map from
   * selection. If only shift key is pressed, it adds invariant map to selection. If only ctrl key is pressed, it
   * changes invariant map presence in selection. If none of the keys is pressed, it sets the invariant map as the
   * only selected one.
   * @param monomerPosition - monomer-position to modify selection with.
   * @param options - selection options.
   * @param notify - flag indicating if bitset changed event should fire.
   */
  modifyInvariantMapSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {
    shiftPressed: false,
    ctrlPressed: false,
  }, notify: boolean = true): void {
    if (notify)
      this.invariantMapSelection = modifySelection(this.invariantMapSelection, monomerPosition, options);
    else
      this._invariantMapSelection = modifySelection(this.invariantMapSelection, monomerPosition, options);
  }

  /**
   * Processes property changes.
   * @param property - changed property.
   */
  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    switch (property.name) {
    case MONOMER_POSITION_PROPERTIES.COLOR:
    case MONOMER_POSITION_PROPERTIES.COLOR_AGGREGATION:
      this.viewerGrid.invalidate();
      break;
    case SAR_PROPERTIES.SEQUENCE:
      this._invariantMapSelection = null;
      break;
    }

    // this will cause colors to recalculate
    this.model.df.columns.toList().forEach((col) => {
      col.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = null;
    });
    if (this.doRender)
      this.render();
  }

  /**
   * Creates monomer-position dataframe to be used in MonomerPosition grid.
   * @return - monomer-position dataframe.
   */
  createMonomerPositionDf(): DG.DataFrame {
    const uniqueMonomers = new Set<string>();
    const splitSeqCols = this.positionColumns;
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

  cacheInvariantMapColors(): void {
  // calculate and cache colors for invariant map
    if (this.colorColumnName && this.dataFrame.col(this.colorColumnName)) {
      const colorCol = this.dataFrame.getCol(this.colorColumnName);
      const colorColData = colorCol!.getRawData();
      let minColorVal = 9999999;
      let maxColorVal = -9999999;
      for (const pCol of this.positionColumns) {
        pCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = {};
        const colorCache = pCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE];
        const colName = pCol.name;
        const positionColData = pCol.getRawData();
        const positionColCategories = pCol.categories;
        if (!this.monomerPositionStats[colName])
          continue;
        const pStats = this.monomerPositionStats[colName]!;
        for (const pMonomer of Object.keys(pStats)) {
          if (pMonomer === 'general')
            continue;
          //const pStatItem = pStats[pMonomer]!;
          const colorValuesIndexes: number[] = [];
          for (let i = 0; i < pCol.length; ++i) {
            if (positionColCategories[positionColData[i]] === pMonomer)
              colorValuesIndexes.push(i);
          }
          const cellColorDataCol = DG.Column.float('color', colorValuesIndexes.length)
            .init((i) => colorColData[colorValuesIndexes[i]]);
          const aggColor = cellColorDataCol.aggregate(this.colorAggregation as DG.AGG);
          colorCache[pMonomer] = aggColor;
          minColorVal = Math.min(minColorVal, aggColor);
          maxColorVal = Math.max(maxColorVal, aggColor);
        }
        pCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = colorCache;
      }

      // do another swing to normalize colors
      for (const pCol of this.positionColumns) {
        const colorCache = pCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE];
        if (!colorCache)
          continue;
        for (const pMonomer of Object.keys(colorCache)) {
          if (this.activityTarget === C.ACTIVITY_TARGET.LOW)
            colorCache[pMonomer] = maxColorVal - colorCache[pMonomer] + minColorVal;
          colorCache[pMonomer] = DG.Color.scaleColor(colorCache[pMonomer], minColorVal, maxColorVal);
        }
        pCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = colorCache;
      }
    }
  }


  /**
   * Builds MonomerPosition interactive grid.
   * @return - MonomerPosition grid.
   */
  createViewerGrid(): DG.Grid {
    const monomerPositionDf = this.createMonomerPositionDf();
    const grid = monomerPositionDf.plot.grid();
    grid.sort([C.COLUMNS_NAMES.MONOMER]);
    const positionColumns = this.positionColumns.map((col) => col.name);
    grid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...positionColumns]);
    const monomerCol = monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setMonomerRenderer(monomerCol, this.alphabet);
    this.cacheInvariantMapColors();

    grid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this,
      this.mode === SELECTION_MODE.INVARIANT_MAP, this.dataFrame.getCol(this.colorColumnName),
      this.colorAggregation as DG.AGG));

    grid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      highlightMonomerPosition(monomerPosition, this.dataFrame, this.monomerPositionStats);
      this.model.isHighlighting = true;
      const columnEntries = this.getTotalViewerAggColumns();

      return showTooltip(this.model.df, this.getScaledActivityColumn(), columnEntries, {
        fromViewer: true,
        isMutationCliffs: this.mode === SELECTION_MODE.MUTATION_CLIFFS, monomerPosition, x, y,
        mpStats: this.monomerPositionStats, cliffStats: this.cliffStats?.stats ?? undefined,
      });
    });
    grid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    DG.debounce(grid.onCurrentCellChanged, 500).subscribe((gridCell: DG.GridCell) => {
      try {
        if (!this.keyPressed)
          return;


        if (this.currentGridCell !== null) {
          const previousMonomerPosition = this.getMonomerPosition(this.currentGridCell);
          if (this.mode === SELECTION_MODE.INVARIANT_MAP) {
            this.modifyInvariantMapSelection(previousMonomerPosition, {
              shiftPressed: true,
              ctrlPressed: true,
            }, false);
          } else {
            const hasMutationCliffs = this.mutationCliffs
              ?.get(previousMonomerPosition.monomerOrCluster)?.get(previousMonomerPosition.positionOrClusterType)?.size;
            if (hasMutationCliffs) {
              this.modifyMutationCliffsSelection(previousMonomerPosition, {
                shiftPressed: true,
                ctrlPressed: true,
              }, false);
            }
          }
        }
        const monomerPosition = this.getMonomerPosition(gridCell);
        if (this.mode === SELECTION_MODE.INVARIANT_MAP)
          this.modifyInvariantMapSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, true);
        else {
          const hasMutationCliffs = this.mutationCliffs
            ?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size;
          if (hasMutationCliffs) {
            this.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false},
              true);
          }
        }

        grid.invalidate();
        setTimeout(() => grid?.invalidate(), 300);
      } finally {
        this.keyPressed = false;
        this.currentGridCell = gridCell;
      }
    });
    grid.root.addEventListener('keydown', (ev) => {
      this.keyPressed = ev.key.startsWith('Arrow');
      if (this.keyPressed)
        return;


      if (ev.key === 'Escape' || (ev.code === 'KeyA' && ev.ctrlKey && ev.shiftKey)) {
        if (this.mode === SELECTION_MODE.INVARIANT_MAP)
          this._invariantMapSelection = initSelection(this.positionColumns);
        else
          this._mutationCliffsSelection = initSelection(this.positionColumns);
      } else if (ev.code === 'KeyA' && ev.ctrlKey) {
        const positions = Object.keys(this.monomerPositionStats).filter((pos) => pos !== 'general');
        for (const position of positions) {
          const monomers = Object.keys(this.monomerPositionStats[position]!)
            .filter((monomer) => monomer !== 'general');
          for (const monomer of monomers) {
            const monomerPosition = {monomerOrCluster: monomer, positionOrClusterType: position};
            if (this.mode === SELECTION_MODE.INVARIANT_MAP)
              this.modifyInvariantMapSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, false);
            else {
              this.modifyMutationCliffsSelection(monomerPosition, {
                shiftPressed: true,
                ctrlPressed: false,
              }, false);
            }
          }
        }
      } else
        return;


      this.model.fireBitsetChanged(VIEWER_TYPE.MONOMER_POSITION);
      grid.invalidate();
    });
    grid.root.addEventListener('click', (ev) => {
      const gridCell = grid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell?.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        return;


      const monomerPosition = this.getMonomerPosition(gridCell);
      if (this.mode === SELECTION_MODE.INVARIANT_MAP) {
        this.modifyInvariantMapSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
        if (isSelectionEmpty(this.invariantMapSelection))
          monomerPositionDf.currentRowIdx = -1;
      } else {
        const hasMutationCliffs = this.mutationCliffs?.get(monomerPosition.monomerOrCluster)
          ?.get(monomerPosition.positionOrClusterType)?.size;
        if (hasMutationCliffs) {
          this.modifyMutationCliffsSelection(monomerPosition,
            {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
        }
      }
      grid.invalidate();

      this.showHelp();
    });

    // const columnWidths = new Map<string, number>();
    // columnWidths.set(C.COLUMNS_NAMES.MONOMER, AAR_CELL_WIDTH);
    // for (const posCol of positionColumns)
    //   columnWidths.set(posCol, MUTATION_CLIFFS_CELL_WIDTH);

    // grid.onColumnResized.subscribe(() => {
    //   try {
    //     const resizedCol = Array.from(columnWidths.keys()).find((col) => grid.col(col)?.width !== columnWidths.get(col));
    //     if (!resizedCol)
    //       return;
    //     const resizedColWidth = Math.max(grid.col(resizedCol)?.width ?? MUTATION_CLIFFS_CELL_WIDTH / 2, MUTATION_CLIFFS_CELL_WIDTH / 2);
    //     Array.from(columnWidths.keys()).forEach((col) => {
    //       columnWidths.set(col, resizedColWidth);
    //       grid.col(col) && (grid.col(col)!.width = resizedColWidth);
    //     });
    //     grid.props.rowHeight = resizedColWidth / 1.5;
    //   } catch (e) {
    //     grok.log.error(e);
    //   }
    // });
    setViewerGridProps(grid);

    // Monomer cell renderer overrides width settings. This way I ensure is "initially" set.
    const afterDraw = grid.onAfterDrawContent.subscribe(() => {
      const monomerGCol = grid.col(C.COLUMNS_NAMES.MONOMER)!;
      if (monomerGCol.width === AAR_CELL_WIDTH) {
        afterDraw.unsubscribe();
        return;
      }
      monomerGCol.width = AAR_CELL_WIDTH;
      for (const posCol of positionColumns) {
        const gcCol = grid.col(posCol)!;
        gcCol.width = MUTATION_CLIFFS_CELL_WIDTH;
      }
    });

    return grid;
  }

  /** Shows viewer context help. */
  showHelp(): void {
    _package.files.readAsText('help/monomer-position.md').then((text) => {
      grok.shell.windows.help.showHelp(ui.markdown(text));
    }).catch((e) => grok.log.error(e));
  }

  /**
   * Gets monomer-position from MonomerPosition grid cell.
   * @param gridCell - MonomerPosition grid cell.
   * @return - monomer-position object.
   */
  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {
      monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!) as string,
      positionOrClusterType: gridCell!.tableColumn!.name,
    };
  }

  /** Renders the MonomerPosition viewer body. */
  render(): void {
    $(this.root).empty();
    if (!this.activityColumnName || !this.sequenceColumnName) {
      this.root.appendChild(ui.divText('Please, select a sequence and activity columns in the viewer properties'));
      return;
    }
    // Backward compatability with 1.16.0
    const columnProperty = this.getProperty(MONOMER_POSITION_PROPERTIES.COLOR);
    if (columnProperty) {
      columnProperty.choices = wu(grok.shell.t.columns.numerical).toArray().map((col) => col.name);
      if (columnProperty.get(this) === C.COLUMNS_NAMES.ACTIVITY_SCALED)
        columnProperty.set(this, C.COLUMNS_NAMES.ACTIVITY);
    }

    $(this.root).empty();
    let switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
    if (this.name === VIEWER_TYPE.MONOMER_POSITION) {
      const mutationCliffsMode = ui.input.bool(SELECTION_MODE.MUTATION_CLIFFS,
        {value: this.mode === SELECTION_MODE.MUTATION_CLIFFS});
      mutationCliffsMode.root.addEventListener('click', () => {
        invariantMapMode.value = false;
        mutationCliffsMode.value = true;
        this.mode = SELECTION_MODE.MUTATION_CLIFFS;
        this.showHelp();
      });
      mutationCliffsMode.setTooltip('Statistically significant changes in activity');
      const invariantMapMode = ui.input.bool(SELECTION_MODE.INVARIANT_MAP, {value: this.mode === SELECTION_MODE.INVARIANT_MAP});
      invariantMapMode.root.addEventListener('click', () => {
        mutationCliffsMode.value = false;
        invariantMapMode.value = true;
        this.mode = SELECTION_MODE.INVARIANT_MAP;
        this.showHelp();
      });
      invariantMapMode.setTooltip('Number of sequences having monomer-position');
      const setDefaultProperties = (input: DG.InputBase): void => {
        $(input.root).find('.ui-input-editor').css('margin', '0px').attr('type', 'radio');
        $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-right', '16px');
        $(input.captionLabel).addClass('ui-label-right');
      };
      setDefaultProperties(mutationCliffsMode);
      setDefaultProperties(invariantMapMode);

      switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root], {id: 'pep-viewer-title'});
      $(switchHost).css('width', 'auto').css('align-self', 'center');
    }
    const viewerRoot = this.viewerGrid.root;
    viewerRoot.style.width = 'auto';
    // expand button
    const expand = ui.iconFA('expand-alt', () => {
      const dialog = ui.dialog('Monomer Position');
      dialog.add(ui.divV([switchHost, viewerRoot], {style: {height: '100%'}}));
      dialog.onCancel(() => this.render());
      dialog.showModal(true);
      this.viewerGrid.invalidate();
    }, 'Show Monomer Position Table in full screen');
    $(expand).addClass('pep-help-icon');

    const header = ui.divH([expand, switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
    this.root.appendChild(ui.divV([header, viewerRoot]));
    this.viewerGrid?.invalidate();
  }
}

/** Vertical structure activity relationship viewer */
export class MostPotentResidues extends SARViewer {
  currentGridRowIdx: number | null = null;

  /** Sets MostPotentResidues properties. */
  constructor() {
    super();
  }

  /**
   * Returns viewer name.
   * @return - viewer name.
   */
  get name(): string {
    return VIEWER_TYPE.MOST_POTENT_RESIDUES;
  }

  /** Processes attached table and sets viewer properties. */
  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    this.render();
  }

  /**
   * Processes property changes.
   * @param property - changed property.
   */
  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (this.doRender)
      this.render();
  }

  /**
   * Creates most potent residues dataframe to be used in MostPotentResidues grid.
   * @return - most potent residues dataframe.
   */
  createMostPotentResiduesDf(): DG.DataFrame {
    const monomerPositionStatsEntries = Object.entries(this.monomerPositionStats) as [string, PositionStats][];
    const posData: number[] = new Array(monomerPositionStatsEntries.length - 1);
    const monomerData: string[] = new Array(posData.length);
    const mdData: number[] = new Array(posData.length);
    const pValData: (number | null)[] = new Array(posData.length);
    const countData: number[] = new Array(posData.length);
    const ratioData: number[] = new Array(posData.length);
    const meanData: number[] = new Array(posData.length);
    const aggrColumnEntries = this.getTotalViewerAggColumns();
    const aggColNames = aggrColumnEntries.map(([colName, aggFn]) => getAggregatedColName(aggFn, colName));
    const aggrColsData = new Array<Array<number>>(aggColNames.length);
    aggColNames.forEach(((_, idx) => aggrColsData[idx] = new Array<number>(posData.length)));
    let i = 0;
    for (const [position, positionStats] of monomerPositionStatsEntries) {
      const generalPositionStats = positionStats.general;
      if (!generalPositionStats)
        continue;


      if (Object.entries(positionStats).length === 1)
        continue;


      const filteredMonomerStats: [string, StatsItem][] = [];
      for (const [monomer, monomerStats] of Object.entries(positionStats)) {
        if (monomer === 'general')
          continue;


        if ((monomerStats as StatsItem).count > 1 && (monomerStats as StatsItem).pValue === null)
          filteredMonomerStats.push([monomer, monomerStats as StatsItem]);


        if ((monomerStats as StatsItem).pValue === generalPositionStats.minPValue)
          filteredMonomerStats.push([monomer, monomerStats as StatsItem]);
      }

      if (filteredMonomerStats.length === 0)
        continue;


      let maxEntry: [string, StatsItem] | null = null;
      // depending on the chosen target for activity, we might want to prioritize higher or lower activity.
      for (const [monomer, monomerStats] of filteredMonomerStats) {
        if (maxEntry === null ||
          (this.activityTarget === C.ACTIVITY_TARGET.HIGH && maxEntry[1].meanDifference < monomerStats.meanDifference) || // if we want to prioritize higher activity
          (this.activityTarget === C.ACTIVITY_TARGET.LOW && maxEntry[1].meanDifference > monomerStats.meanDifference) // if we want to prioritize lower activity
        )
          maxEntry = [monomer, monomerStats];
      }

      if (maxEntry === null)
        continue;


      posData[i] = parseInt(position);
      monomerData[i] = maxEntry![0];
      mdData[i] = maxEntry![1].meanDifference;
      pValData[i] = maxEntry![1].pValue;
      countData[i] = maxEntry![1].count;
      ratioData[i] = maxEntry![1].ratio;
      meanData[i] = maxEntry![1].mean;

      const stats = this.monomerPositionStats[position][maxEntry![0]];
      const mask = DG.BitSet.fromBytes(stats.mask.buffer.buffer, this.model.df.col(this.activityColumnName)!.length);
      for (let j = 0; j < aggColNames.length; j++) {
        const [colName, aggFn] = aggrColumnEntries[j];
        aggrColsData[j][i] = getAggregatedValue(this.model.df.getCol(colName), aggFn, mask);
      }
      ++i;
    }

    posData.length = i;
    monomerData.length = i;
    mdData.length = i;
    pValData.length = i;
    countData.length = i;
    ratioData.length = i;
    meanData.length = i;

    const mprDf = DG.DataFrame.create(i);
    const mprDfCols = mprDf.columns;
    mprDfCols.add(DG.Column.fromList(DG.TYPE.INT, C.COLUMNS_NAMES.POSITION, posData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.STRING, C.COLUMNS_NAMES.MONOMER, monomerData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.MEAN_DIFFERENCE, mdData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.MEAN, meanData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.P_VALUE, pValData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.INT, C.COLUMNS_NAMES.COUNT, countData));
    mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, C.COLUMNS_NAMES.RATIO, ratioData));
    aggColNames.forEach((it, idx) => mprDfCols.add(DG.Column.fromList(DG.TYPE.FLOAT, it, aggrColsData[idx])));

    return mprDf;
  }

  /**
   * Builds MostPotentResidues interactive grid.
   * @return - MostPotentResidues grid.
   */
  createViewerGrid(): DG.Grid {
    const mprDf = this.createMostPotentResiduesDf();
    const grid = mprDf.plot.grid();
    grid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = grid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = '#.000';
    pValGridCol.name = 'P-value';
    const monomerCol = mprDf.getCol(C.COLUMNS_NAMES.MONOMER);

    // Setting Monomer column renderer
    CR.setMonomerRenderer(monomerCol, this.alphabet);
    grid.onCellRender.subscribe(
      (args: DG.GridCellRenderArgs) => renderCell(args, this, false, undefined, undefined));

    grid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      highlightMonomerPosition(monomerPosition, this.dataFrame, this.monomerPositionStats);
      this.model.isHighlighting = true;

      if (gridCell.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        monomerPosition.positionOrClusterType = C.COLUMNS_NAMES.MONOMER;
      else if (gridCell.tableColumn?.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return false;


      const columnEntries = this.getTotalViewerAggColumns();
      const aggrValues = gridCell.tableRowIndex == null ? undefined :
        getAggregatedColumnValuesFromDf(gridCell.grid.dataFrame!, gridCell.tableRowIndex, columnEntries, {});
      return showTooltip(this.model.df, this.getScaledActivityColumn(), columnEntries,
        {
          fromViewer: true, isMostPotentResidues: true, monomerPosition, x, y, mpStats: this.monomerPositionStats,
          aggrColValues: aggrValues,
        });
    });
    DG.debounce(grid.onCurrentCellChanged, 500).subscribe((gridCell: DG.GridCell) => {
      try {
        if ((this.keyPressed && mprDf.currentCol.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE) || !this.keyPressed)
          return;


        const monomerPosition = this.getMonomerPosition(gridCell);
        if (this.currentGridRowIdx !== null) {
          const previousMonomerPosition = this.getMonomerPosition(grid.cell('Diff', this.currentGridRowIdx));
          this.modifyMutationCliffsSelection(previousMonomerPosition, {
            shiftPressed: true,
            ctrlPressed: true,
          }, false);
        }
        const hasMutationCliffs = this.mutationCliffs?.get(monomerPosition.monomerOrCluster)
          ?.get(monomerPosition.positionOrClusterType)?.size;
        if (hasMutationCliffs)
          this.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false});


        grid.invalidate();
      } finally {
        this.keyPressed = false;
        this.currentGridRowIdx = gridCell.gridRow;
      }
    });
    grid.root.addEventListener('keydown', (ev) => {
      this.keyPressed = ev.key.startsWith('Arrow');
      if (this.keyPressed)
        return;


      if (ev.key === 'Escape' || (ev.code === 'KeyA' && ev.ctrlKey && ev.shiftKey))
        this._mutationCliffsSelection = initSelection(this.positionColumns);
      else if (ev.code === 'KeyA' && ev.ctrlKey) {
        for (let rowIdx = 0; rowIdx < mprDf.rowCount; ++rowIdx) {
          const monomerPosition = this.getMonomerPosition(grid.cell('Diff', rowIdx));
          this.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: true, ctrlPressed: false}, false);
        }
      } else
        return;


      this.model.fireBitsetChanged(VIEWER_TYPE.MOST_POTENT_RESIDUES);
      grid.invalidate();
    });
    grid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    grid.root.addEventListener('click', (ev) => {
      const gridCell = grid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell!.tableColumn!.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;


      const monomerPosition = this.getMonomerPosition(gridCell);
      const hasMutationCliffs = this.mutationCliffs?.get(monomerPosition.monomerOrCluster)
        ?.get(monomerPosition.positionOrClusterType)?.size;
      if (!hasMutationCliffs)
        return;


      this.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      grid.invalidate();

      _package.files.readAsText('help/most-potent-residues.md').then((text) => {
        grok.shell.windows.help.showHelp(ui.markdown(text));
      }).catch((e) => grok.log.error(e));
    });

    setViewerGridProps(grid);
    const mdCol: DG.GridColumn = grid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';
    // Monomer cell renderer overrides width settings. This way I ensure is "initially" set.
    const afterDraw = grid.onAfterDrawContent.subscribe(() => {
      const monomerGCol = grid.col(C.COLUMNS_NAMES.MONOMER)!;
      if (monomerGCol.width === AAR_CELL_WIDTH) {
        afterDraw.unsubscribe();
        return;
      }
      monomerGCol.width = AAR_CELL_WIDTH;
      mdCol.width = MUTATION_CLIFFS_CELL_WIDTH;
    });

    return grid;
  }

  /**
   * Gets monomer-position from MostPotentResidues grid cell.
   * @param gridCell - MostPotentResidues grid cell.
   * @return - monomer-position object.
   */
  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {
      monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!),
      positionOrClusterType: `${gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.POSITION, gridCell!.tableRowIndex!)}`,
    };
  }

  /** Renders the MostPotentResidues viewer body. */
  render(): void {
    $(this.root).empty();
    if (!this.activityColumnName || !this.sequenceColumnName) {
      this.root.appendChild(ui.divText('Please, select a sequence and activity columns in the viewer properties'));
      return;
    }
    const switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
    const viewerRoot = this.viewerGrid.root;
    viewerRoot.style.width = 'auto';
    const header = ui.divH([switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
    this.root.appendChild(ui.divV([header, viewerRoot]));
    this.viewerGrid?.invalidate();
  }
}

/**
 * Renders SARViewer grid cells.
 * @param args - grid cell render arguments.
 * @param viewer - viewer instance.
 * @param isInvariantMap - flag indicating if viewr is operating in invariant map mode.
 * @param colorCol - column to use for color coding in invariant map mode.
 * @param colorAgg - aggregation function to use for color coding in invariant map mode.
 */
function renderCell(args: DG.GridCellRenderArgs, viewer: SARViewer, isInvariantMap?: boolean,
  _colorCol?: DG.Column<number>, _colorAgg?: DG.AGG): void {
  const renderColNames = [...viewer.positionColumns.map((col) => col.name), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
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
    canvasContext.restore();
    return;
  }

  const tableColName = cell.tableColumn?.name;
  const tableRowIndex = cell.tableRowIndex!;
  if (!cell.isTableCell || renderColNames.indexOf(tableColName!) === -1) {
    canvasContext.restore();
    return;
  }

  const gridTable = cell.grid.table;
  const currentMonomer: string = gridTable.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);
  const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
    gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex).toFixed();
  const currentPosStats = viewer.monomerPositionStats[currentPosition];

  if (!currentPosStats![currentMonomer]) {
    args.preventDefault();
    canvasContext.restore();
    return;
  }

  if (isInvariantMap) {
    const value = currentPosStats![currentMonomer]!.count;
    const positionCol = viewer.positionColumns.find((col) => col.name === currentPosition)!;
    const colorCache: { [_: string]: number } = positionCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] ?? {};
    let color: number = DG.Color.white;
    // const colorColStats = colorCol!.stats;
    if (colorCache[currentMonomer])
      color = colorCache[currentMonomer];
    else if (viewer instanceof MonomerPosition) {
      viewer.cacheInvariantMapColors();
      color = (positionCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] ?? {})[currentMonomer] ?? DG.Color.white;
    }
    // } else {
    //   const positionColData = positionCol.getRawData();
    //   const positionColCategories = positionCol.categories;

    //   const colorColData = colorCol!.getRawData();
    //   const colorValuesIndexes: number[] = [];
    //   for (let i = 0; i < positionCol.length; ++i) {
    //     if (positionColCategories[positionColData[i]] === currentMonomer)
    //       colorValuesIndexes.push(i);
    //   }
    //   const cellColorDataCol = DG.Column.float('color', colorValuesIndexes.length)
    //     .init((i) => colorColData[colorValuesIndexes[i]]);
    //   let aggColor = cellColorDataCol.aggregate(colorAgg!);
    //   // if the target activity is low, we might wat to reverse the color scale
    //   if (viewer.activityTarget === C.ACTIVITY_TARGET.LOW)
    //     aggColor = colorColStats.max - aggColor + colorColStats.min;
    //   color = DG.Color.scaleColor(aggColor, colorColStats.min, colorColStats.max);
    //   colorCache[currentMonomer] = aggColor;
    //   positionCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = colorCache;
    // }

    CR.renderInvariantMapCell(canvasContext, currentMonomer, currentPosition,
      (viewer as MonomerPosition).invariantMapSelection, value, bound, color);
  } else
    CR.renderMutationCliffCell(canvasContext, currentMonomer, currentPosition, viewer, bound);


  args.preventDefault();
  canvasContext.restore();
}

/**
 * Sets common grid properties.
 * @param grid - viewer grid.
 */
function setViewerGridProps(grid: DG.Grid): void {
  const gridProps = grid.props;
  gridProps.allowEdit = false;
  gridProps.allowRowSelection = false;
  gridProps.allowBlockSelection = false;
  gridProps.allowColSelection = false;
  gridProps.showRowHeader = false;
  gridProps.showCurrentRowIndicator = false;
  gridProps.showReadOnlyNotifications = false;

  gridProps.rowHeight = 20;
}
