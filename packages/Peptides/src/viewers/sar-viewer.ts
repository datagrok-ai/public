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
import {AggregationColumns, MonomerPositionStats, PositionStats, StatsItem, getAggregatedColName, getAggregatedColumnValuesFromDf, getAggregatedValue} from '../utils/statistics';
import {_package} from '../package';
import {showTooltip} from '../utils/tooltips';
import {calculateMonomerPositionStatistics, findMutations, MutationCliffsOptions} from '../utils/algorithms';
import {
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
}

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

  constructor() {
    super();
    // General properties
    this.sequenceColumnName = this.column(SAR_PROPERTIES.SEQUENCE,
      {category: PROPERTY_CATEGORIES.GENERAL, semType: DG.SEMTYPE.MACROMOLECULE, nullable: false});
    this.activityColumnName = this.column(SAR_PROPERTIES.ACTIVITY, {
      category: PROPERTY_CATEGORIES.GENERAL,
      nullable: false,
    });
    this.activityScaling = this.string(SAR_PROPERTIES.ACTIVITY_SCALING, C.SCALING_METHODS.NONE,
      {
        category: PROPERTY_CATEGORIES.GENERAL,
        choices: Object.values(C.SCALING_METHODS),
        nullable: false,
      }) as C.SCALING_METHODS;

    // Mutation Cliffs properties
    this.targetColumnName = this.column(SAR_PROPERTIES.TARGET, {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS});
    this.targetCategory = this.string(SAR_PROPERTIES.TARGET_CATEGORY, null,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, choices: []});
    this.minActivityDelta = this.float(SAR_PROPERTIES.MIN_ACTIVITY_DELTA, 0,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, min: 0, max: 100});
    this.maxMutations = this.int(SAR_PROPERTIES.MAX_MUTATIONS, 1,
      {category: PROPERTY_CATEGORIES.MUTATION_CLIFFS, min: 1, max: 50});

    this.columns = this.columnList(SAR_PROPERTIES.COLUMNS, [], {category: PROPERTY_CATEGORIES.AGGREGATION});
    const aggregationChoices = C.AGGREGATION_TYPES;
    this.aggregation = this.string(SAR_PROPERTIES.AGGREGATION, DG.AGG.AVG,
      {category: PROPERTY_CATEGORIES.AGGREGATION, choices: aggregationChoices});
  }

  _viewerGrid: DG.Grid | null = null;

  get viewerGrid(): DG.Grid {
    this._viewerGrid ??= this.createViewerGrid();
    return this._viewerGrid;
  }

  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  get alphabet(): string {
    const col = this.dataFrame.getCol(this.sequenceColumnName);
    return col.getTag(bioTAGS.alphabet) ?? ALPHABET.UN;
  }

  _model?: PeptidesModel;

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  _positionColumns: DG.Column<string>[] | null = null;

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

  set monomerPositionStats(mps: MonomerPositionStats) {
    this._monomerPositionStats = mps;
  }

  _mutationCliffs: type.MutationCliffs | null = null;

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

  set mutationCliffs(mc: type.MutationCliffs) {
    this._mutationCliffs = mc;
    this.viewerGrid.invalidate();
  }

  _mutationCliffsSelection: type.Selection | null = null;

  get mutationCliffsSelection(): type.Selection {
    const tagSuffix = this instanceof MonomerPosition ? C.SUFFIXES.MP : C.SUFFIXES.MPR;
    const tagSelection = this.dataFrame.getTag(`${tagSuffix}${C.TAGS.MUTATION_CLIFFS_SELECTION}`);
    this._mutationCliffsSelection ??= tagSelection === null ? initSelection(this.positionColumns) :
      JSON.parse(tagSelection);
    return this._mutationCliffsSelection!;
  }

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

  modifyMutationCliffsSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {
    shiftPressed: false,
    ctrlPressed: false,
  }, notify: boolean = true): void {
    if (notify)
      this.mutationCliffsSelection = modifySelection(this.mutationCliffsSelection, monomerPosition, options);
    else
      this._mutationCliffsSelection = modifySelection(this.mutationCliffsSelection, monomerPosition, options);
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.doRender = true;
    switch (property.name) {
    case `${SAR_PROPERTIES.SEQUENCE}${COLUMN_NAME}`:
      this._positionColumns = null;
      this._monomerPositionStats = null;
      this._mutationCliffs = null;
      this._mutationCliffsSelection = null;
      this._viewerGrid = null;
      break;
    case `${SAR_PROPERTIES.ACTIVITY}${COLUMN_NAME}`:
    case SAR_PROPERTIES.ACTIVITY_SCALING:
      this._monomerPositionStats = null;
      this._mutationCliffs = null;
      this._mutationCliffsSelection = null;
      this._viewerGrid = null;
      this._scaledActivityColumn = null;
      break;
    case `${SAR_PROPERTIES.TARGET}${COLUMN_NAME}`:
    case SAR_PROPERTIES.TARGET_CATEGORY:
    case SAR_PROPERTIES.MIN_ACTIVITY_DELTA:
    case SAR_PROPERTIES.MAX_MUTATIONS:
      this._mutationCliffs = null;
      this._mutationCliffsSelection = null;
      this.doRender = false;
      break;
    case SAR_PROPERTIES.COLUMNS:
    case SAR_PROPERTIES.AGGREGATION:
      if (this instanceof MostPotentResidues)
        this._viewerGrid = null;
      break;
    }
    if (this.mutationCliffs === null && this.sequenceColumnName && this.activityColumnName)
      this.calculateMutationCliffs().then((mc) => this.mutationCliffs = mc);
  }

  getAggregationColumns(): AggregationColumns {
    return Object.fromEntries(
      this.columns.map((colName) => [colName, this.aggregation] as [string, DG.AGG]));
  }

  getTotalViewerAggColumns(): [string, DG.AggregationType][] {
    const aggrCols = this.getAggregationColumns();
    return getTotalAggColumns(this.columns, aggrCols, this.model?.settings?.columns);
  }

  createViewerGrid(): DG.Grid {
    throw new Error('Not implemented');
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    if (isApplicableDataframe(this.dataFrame)) {
      this.getProperty(`${SAR_PROPERTIES.SEQUENCE}${COLUMN_NAME}`)
        ?.set(this, this.dataFrame.columns.bySemType(DG.SEMTYPE.MACROMOLECULE)!.name);
      this.getProperty(`${SAR_PROPERTIES.ACTIVITY}${COLUMN_NAME}`)
        ?.set(this, wu(this.dataFrame.columns.numerical).next().value.name);
      if (this.mutationCliffs === null && this.sequenceColumnName && this.activityColumnName)
        this.calculateMutationCliffs().then((mc) => this.mutationCliffs = mc);
    } else {
      const msg = 'PeptidesError: dataframe is missing Macromolecule or numeric columns';
      grok.log.error(msg);
      grok.shell.warning(msg);
    }
  }

  async calculateMutationCliffs(): Promise<MutationCliffs> {
    const scaledActivityCol: DG.Column<number> = this.dataFrame.getCol(this.sequenceColumnName);
    //TODO: set categories ordering the same to share compare indexes instead of strings
    const monomerCols: type.RawColumn[] = this.positionColumns.map(extractColInfo);
    const targetCol = this.targetColumnName ? extractColInfo(this.dataFrame.getCol(this.targetColumnName)) : null;

    const options: MutationCliffsOptions = {
      maxMutations: this.maxMutations, minActivityDelta: this.minActivityDelta,
      targetCol, currentTarget: this.targetCategory,
    };
    return await findMutations(scaledActivityCol.getRawData(), monomerCols, options);
  }
}

/** Structure-activity relationship viewer */
export class MonomerPosition extends SARViewer {
  colorColumnName: string;
  colorAggregation: string;
  currentGridCell: DG.GridCell | null = null;

  constructor() {
    super();

    const colorChoices = wu(grok.shell.t.columns.numerical).toArray().map((col) => col.name);
    this.colorColumnName = this.column(MONOMER_POSITION_PROPERTIES.COLOR,
      {category: PROPERTY_CATEGORIES.INVARIANT_MAP, choices: colorChoices, nullable: false});
    const aggregationChoices = C.AGGREGATION_TYPES;
    this.colorAggregation = this.string(MONOMER_POSITION_PROPERTIES.COLOR_AGGREGATION, DG.AGG.AVG,
      {category: PROPERTY_CATEGORIES.INVARIANT_MAP, choices: aggregationChoices});
  }

  get name(): string {
    return VIEWER_TYPE.MONOMER_POSITION;
  }

  get mode(): SELECTION_MODE {
    return this.dataFrame.getTag(C.TAGS.MONOMER_POSITION_MODE) as SELECTION_MODE ??
      SELECTION_MODE.MUTATION_CLIFFS;
  }

  set mode(mode: SELECTION_MODE) {
    this.dataFrame.setTag(C.TAGS.MONOMER_POSITION_MODE, mode);
    this.viewerGrid.invalidate();
  }

  _invariantMapSelection: type.Selection | null = null;

  get invariantMapSelection(): type.Selection {
    const tagSelection = this.dataFrame.getTag(`${C.SUFFIXES.MP}${C.TAGS.INVARIANT_MAP_SELECTION}`);
    this._invariantMapSelection ??= tagSelection === null ? initSelection(this.positionColumns) :
      JSON.parse(tagSelection);
    return this._invariantMapSelection!;
  }

  set invariantMapSelection(selection: type.Selection) {
    this._invariantMapSelection = selection;
    this.dataFrame.setTag(`${C.SUFFIXES.MP}${C.TAGS.INVARIANT_MAP_SELECTION}`, JSON.stringify(selection));
    this.model.fireBitsetChanged(VIEWER_TYPE.MONOMER_POSITION);
    this.model.analysisView.grid.invalidate();
  }

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

  modifyInvariantMapSelection(monomerPosition: type.SelectionItem, options: type.SelectionOptions = {
    shiftPressed: false,
    ctrlPressed: false,
  }, notify: boolean = true): void {
    if (notify)
      this.invariantMapSelection = modifySelection(this.invariantMapSelection, monomerPosition, options);
    else
      this._invariantMapSelection = modifySelection(this.invariantMapSelection, monomerPosition, options);
  }

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

  createViewerGrid(): DG.Grid {
    const monomerPositionDf = this.createMonomerPositionDf();
    const grid = monomerPositionDf.plot.grid();
    grid.sort([C.COLUMNS_NAMES.MONOMER]);
    const positionColumns = this.positionColumns.map((col) => col.name);
    grid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...positionColumns]);
    const monomerCol = monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setMonomerRenderer(monomerCol, this.alphabet);
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
        mpStats: this.monomerPositionStats,
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

  showHelp(): void {
    _package.files.readAsText('help/monomer-position.md').then((text) => {
      grok.shell.windows.help.showHelp(ui.markdown(text));
    }).catch((e) => grok.log.error(e));
  }

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {
      monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!) as string,
      positionOrClusterType: gridCell!.tableColumn!.name,
    };
  }

  render(refreshOnly = false): void {
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

    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      if (this.name === VIEWER_TYPE.MONOMER_POSITION) {
        const mutationCliffsMode = ui.boolInput(SELECTION_MODE.MUTATION_CLIFFS,
          this.mode === SELECTION_MODE.MUTATION_CLIFFS);
        mutationCliffsMode.root.addEventListener('click', () => {
          invariantMapMode.value = false;
          mutationCliffsMode.value = true;
          this.mode = SELECTION_MODE.MUTATION_CLIFFS;
          this.showHelp();
        });
        mutationCliffsMode.setTooltip('Statistically significant changes in activity');
        const invariantMapMode = ui.boolInput(SELECTION_MODE.INVARIANT_MAP, this.mode === SELECTION_MODE.INVARIANT_MAP);
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
      const header = ui.divH([switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

/** Vertical structure activity relationship viewer */
export class MostPotentResidues extends SARViewer {
  currentGridRowIdx: number | null = null;

  constructor() {
    super();
  }

  get name(): string {
    return VIEWER_TYPE.MOST_POTENT_RESIDUES;
  }

  detach(): void {
    this.subs.forEach((sub) => sub.unsubscribe());
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.helpUrl = 'https://datagrok.ai/help/datagrok/solutions/domains/bio/peptides-sar';
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (this.doRender)
      this.render();
  }

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
    aggColNames.forEach(((_, idx) => aggrColsData[idx] = new Array(posData.length)));
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
      for (const [monomer, monomerStats] of filteredMonomerStats) {
        if (maxEntry === null || maxEntry[1].meanDifference < monomerStats.meanDifference)
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
        const value = getAggregatedValue(this.model.df.getCol(colName), aggFn, mask);
        aggrColsData[j][i] = value;
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
      const aggrValues = gridCell.tableRowIndex ? 
        getAggregatedColumnValuesFromDf(gridCell.grid.dataFrame!, gridCell.tableRowIndex, columnEntries, {}) : undefined;
      return showTooltip(this.model.df, this.getScaledActivityColumn(), columnEntries,
        {fromViewer: true, isMutationCliffs: true, monomerPosition, x, y, mpStats: this.monomerPositionStats, 
          aggrColValues: aggrValues});
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

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {
      monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!),
      positionOrClusterType: `${gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.POSITION, gridCell!.tableRowIndex!)}`,
    };
  }

  render(refreshOnly = false): void {
    $(this.root).empty();
    if (!this.activityColumnName || !this.sequenceColumnName) {
      this.root.appendChild(ui.divText('Please, select a sequence and activity columns in the viewer properties'));
      return;
    }
    if (!refreshOnly) {
      const switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      const header = ui.divH([switchHost], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

function renderCell(args: DG.GridCellRenderArgs, viewer: SARViewer, isInvariantMap?: boolean,
  colorCol?: DG.Column<number>, colorAgg?: DG.AGG): void {
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
    let color: number | null = null;
    if (colorCache[currentMonomer])
      color = colorCache[currentMonomer];
    else {
      const positionColData = positionCol.getRawData();
      const positionColCategories = positionCol.categories;

      const colorColData = colorCol!.getRawData();
      const colorValuesIndexes: number[] = [];
      for (let i = 0; i < positionCol.length; ++i) {
        if (positionColCategories[positionColData[i]] === currentMonomer)
          colorValuesIndexes.push(i);
      }
      const cellColorDataCol = DG.Column.float('color', colorValuesIndexes.length)
        .init((i) => colorColData[colorValuesIndexes[i]]);
      const colorColStats = colorCol!.stats;
      color = DG.Color.scaleColor(cellColorDataCol.aggregate(colorAgg!), colorColStats.min, colorColStats.max);
      colorCache[currentMonomer] = color;
      positionCol.temp[C.TAGS.INVARIANT_MAP_COLOR_CACHE] = colorCache;
    }

    CR.renderInvariantMapCell(canvasContext, currentMonomer, currentPosition,
      (viewer as MonomerPosition).invariantMapSelection, value, bound, color);
  } else
    CR.renderMutationCliffCell(canvasContext, currentMonomer, currentPosition, viewer, bound);

  args.preventDefault();
  canvasContext.restore();
}

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
