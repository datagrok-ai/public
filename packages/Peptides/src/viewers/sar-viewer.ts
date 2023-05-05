import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {PeptidesModel, VIEWER_TYPE} from '../model';

export enum MONOMER_POSITION_MODE {
  MUTATION_CLIFFS = 'Mutation Cliffs',
  INVARIANT_MAP = 'Invariant Map',
}

/** Structure-activity relationship viewer */
export class MonomerPosition extends DG.JsViewer {
  _titleHost = ui.divText('Mutation Cliffs', {id: 'pep-viewer-title'});
  _viewerGrid!: DG.Grid;
  _model!: PeptidesModel;
  colorColumnName: string;
  aggregation: string;

  constructor() {
    super();
    this.colorColumnName = this.string('colorColumnName', C.COLUMNS_NAMES.ACTIVITY_SCALED, {category: 'Invariant Map'});
    this.aggregation = this.string('aggregation', 'avg', {category: 'Invariant Map', choices: Object.values(DG.AGG)});
  }

  get name(): string {return VIEWER_TYPE.MONOMER_POSITION;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMonomerPositionGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  get mode(): MONOMER_POSITION_MODE {
    return this.dataFrame.getTag(C.TAGS.MONOMER_POSITION_MODE) as MONOMER_POSITION_MODE ??
      MONOMER_POSITION_MODE.MUTATION_CLIFFS;
  }
  set mode(mode: MONOMER_POSITION_MODE) {
    this.dataFrame.setTag(C.TAGS.MONOMER_POSITION_MODE, mode);
    this.viewerGrid.invalidate();
  }

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  onTableAttached(): void {
    super.onTableAttached();
    this.subs.push(this.model.onMonomerPositionSelectionChanged.subscribe(() => this.viewerGrid.invalidate()));
    this.helpUrl = '/help/domains/bio/peptides.md';
    this.subs.push(this.model.onSettingsChanged.subscribe(() => {
      this.createMonomerPositionGrid();
      this.render();
    }));
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.render();
  }

  createMonomerPositionGrid(): void {
    this.viewerGrid = this.model.monomerPositionDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.MONOMER]);
    this.viewerGrid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...this.model.splitSeqDf.columns.names()]);
    const monomerCol = this.model.monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setAARRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model,
      this.mode === MONOMER_POSITION_MODE.INVARIANT_MAP, this.dataFrame.getCol(this.colorColumnName),
      this.aggregation as DG.AggregationType));
    this.viewerGrid.onCellTooltip((cell: DG.GridCell, x: number, y: number) => showTooltip(cell, x, y, this.model));
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell?.tableColumn?.name == C.COLUMNS_NAMES.MONOMER)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = monomerCol.get(gridCell!.tableRowIndex!);
      chooseAction(aar, position, ev.shiftKey, this.mode === MONOMER_POSITION_MODE.INVARIANT_MAP, this.model);
      this.viewerGrid.invalidate();
      this.model.fireBitsetChanged();
    });
    this.viewerGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(this.model.monomerPositionDf, this.model));

    setViewerGridProps(this.viewerGrid, false);
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      if (this.name == VIEWER_TYPE.MONOMER_POSITION) {
        const mutationCliffsMode = ui.boolInput('', this.mode === MONOMER_POSITION_MODE.MUTATION_CLIFFS);
        mutationCliffsMode.root.addEventListener('click', () => {
          invariantMapMode.value = false;
          mutationCliffsMode.value = true;
          this.mode = MONOMER_POSITION_MODE.MUTATION_CLIFFS;
        });
        mutationCliffsMode.addPostfix('Mutation Cliffs');
        const invariantMapMode = ui.boolInput('', this.mode === MONOMER_POSITION_MODE.INVARIANT_MAP);
        invariantMapMode.root.addEventListener('click', () => {
          mutationCliffsMode.value = false;
          invariantMapMode.value = true;
          this.mode = MONOMER_POSITION_MODE.INVARIANT_MAP;
        });
        invariantMapMode.addPostfix('Invariant Map');
        const setDefaultProperties = (input: DG.InputBase): void => {
          $(input.root).find('.ui-input-editor').css('margin', '0px').attr('type', 'radio');
          $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-left', '5px');
        };
        setDefaultProperties(mutationCliffsMode);
        setDefaultProperties(invariantMapMode);
        $(mutationCliffsMode.root).css('padding-right', '10px').css('padding-left', '5px');

        switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root], {id: 'pep-viewer-title'});
        $(switchHost).css('width', 'auto').css('align-self', 'center');
      }
      const tips: HTMLElement = ui.iconFA('question');
      ui.tooltip.bind(tips,
        () => ui.divV([ui.divText('Color intensity - p-value'), ui.divText('Circle size - Mean difference')]));

      $(tips).addClass('pep-help-icon');

      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      const header = ui.divH([switchHost, tips], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

/** Vertical structure activity relationship viewer */
export class MostPotentResiduesViewer extends DG.JsViewer {
  _titleHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
  _viewerGrid!: DG.Grid;
  _model!: PeptidesModel;

  constructor() {
    super();
  }

  get name(): string {return VIEWER_TYPE.MOST_POTENT_RESIDUES;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMostPotentResiduesGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  get model(): PeptidesModel {
    this._model ??= PeptidesModel.getInstance(this.dataFrame);
    return this._model;
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  onTableAttached(): void {
    super.onTableAttached();
    this.subs.push(this.model.onMonomerPositionSelectionChanged.subscribe(() => this.viewerGrid.invalidate()));
    this.helpUrl = '/help/domains/bio/peptides.md';
    this.subs.push(this.model.onSettingsChanged.subscribe(() => {
      this.createMostPotentResiduesGrid();
      this.render();
    }));
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    this.render();
  }

  createMostPotentResiduesGrid(): void {
    this.viewerGrid = this.model.mostPotentResiduesDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = this.viewerGrid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = '#.000';
    pValGridCol.name = 'P-value';
    const monomerCol = this.model.mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.MONOMER);
    const positionCol: DG.Column<number> = this.model.mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.POSITION);

    // Setting Monomer column renderer
    CR.setAARRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model));
    this.viewerGrid.onCellTooltip((cell: DG.GridCell, x: number, y: number) => showTooltip(cell, x, y, this.model));
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell!.tableColumn!.name != C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const tableRowIdx = gridCell!.tableRowIndex!;
      const position = positionCol.get(tableRowIdx);
      const aar = monomerCol.get(tableRowIdx);
      chooseAction(aar, position!.toFixed(), ev.shiftKey, false, this.model);
      this.viewerGrid.invalidate();
      this.model.fireBitsetChanged();
    });
    this.viewerGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(this.model.mostPotentResiduesDf, this.model));
    const mdCol: DG.GridColumn = this.viewerGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';
    setViewerGridProps(this.viewerGrid, true);
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      const switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      const tips: HTMLElement = ui.iconFA('question');
      ui.tooltip.bind(tips,
        () => ui.divV([ui.divText('Color intensity - p-value'), ui.divText('Circle size - Mean difference')]));

      $(tips).addClass('pep-help-icon');

      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      const header = ui.divH([switchHost, tips], {style: {alignSelf: 'center', lineHeight: 'normal'}});
      this.root.appendChild(ui.divV([header, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }
}

function renderCell(args: DG.GridCellRenderArgs, model: PeptidesModel, isInvariantMap?: boolean,
  colorCol?: DG.Column<number>, colorAgg?: DG.AggregationType): void {
  const renderColNames = [...model.splitSeqDf.columns.names(), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
  // const mdCol = model.monomerPositionStats.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
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
  if (!cell.isTableCell || renderColNames.indexOf(tableColName!) == -1) {
    canvasContext.restore();
    return;
  }

  const gridTable = cell.grid.table;
  const currentMonomer: string = gridTable.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);
  const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
    gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex).toFixed();
  const currentPosStats = model.monomerPositionStats[currentPosition];

  if (!currentPosStats[currentMonomer]) {
    args.preventDefault();
    canvasContext.restore();
    return;
  }

  if (isInvariantMap) {
    const value = currentPosStats[currentMonomer].count;
    const positionCol = model.df.getCol(currentPosition);
    const positionColData = positionCol.getRawData();
    const positionColCategories = positionCol.categories;

    const colorColData = colorCol!.getRawData();
    const colorDataList: number[] = [];
    for (let i = 0; i < positionColData.length; ++i) {
      if (positionColCategories[positionColData[i]] === currentMonomer)
        colorDataList.push(colorColData[i]);
    }
    const cellColorDataCol = DG.Column.fromList('double', '', colorDataList);
    const colorColStats = colorCol!.stats;

    const color = DG.Color.scaleColor(cellColorDataCol.aggregate(colorAgg!), colorColStats.min, colorColStats.max);
    CR.renderInvaraintMapCell(
      canvasContext, currentMonomer, currentPosition, model.monomerPositionFilter, value, bound, color);
  } else {
    CR.renderMutationCliffCell(canvasContext, currentMonomer, currentPosition, model.monomerPositionStats, bound,
      model.monomerPositionSelection, model.mutationCliffs, model.settings.isBidirectional);
  }
  args.preventDefault();
  canvasContext.restore();
}

function showTooltip(cell: DG.GridCell, x: number, y: number, model: PeptidesModel): boolean {
  const renderColNames = [...model.splitSeqDf.columns.names(), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
  const tableCol = cell.tableColumn;
  const tableColName = tableCol?.name;
  const tableRowIndex = cell.tableRowIndex;

  if (!cell.isRowHeader && !cell.isColHeader && tableCol && tableRowIndex != null) {
    const table = cell.grid.table;
    const currentAAR = table.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);

    if (tableCol.semType == C.SEM_TYPES.MONOMER)
      model.showMonomerTooltip(currentAAR, x, y);
    else if (renderColNames.includes(tableColName!)) {
      const currentPosition = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
        table.get(C.COLUMNS_NAMES.POSITION, tableRowIndex).toFixed();

      model.showTooltipAt(currentAAR, currentPosition, x, y);
    }
  }
  return true;
}

function chooseAction(aar: string, position: string, isShiftPressed: boolean, isFilter: boolean,
  model: PeptidesModel): void {
  if (!isShiftPressed)
    model.initMonomerPositionSelection(true);

  model.modifyMonomerPositionSelection(aar, position, isFilter);
}

function cellChanged(table: DG.DataFrame, model: PeptidesModel): void {
  if (model.isCellChanging)
    return;
  model.isCellChanging = true;
  table.currentRowIdx = -1;
  model.isCellChanging = false;
}

function setViewerGridProps(grid: DG.Grid, isMostPotentResiduesGrid: boolean): void {
  const gridProps = grid.props;
  gridProps.allowEdit = false;
  gridProps.allowRowSelection = false;
  gridProps.allowBlockSelection = false;
  gridProps.allowColSelection = false;

  gridProps.rowHeight = 20;
  const girdCols = grid.columns;
  const colNum = girdCols.length;
  for (let i = 0; i < colNum; ++i) {
    const col = girdCols.byIndex(i)!;
    const colName = col.name;
    col.width = isMostPotentResiduesGrid && colName !== 'Diff' && colName !== C.COLUMNS_NAMES.MONOMER ? 50 :
      gridProps.rowHeight + 10;
  }
}
