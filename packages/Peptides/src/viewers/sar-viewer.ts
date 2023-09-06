import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import wu from 'wu';
import {SelectionItem} from '../utils/types';

export enum SELECTION_MODE {
  MUTATION_CLIFFS = 'Mutation Cliffs',
  INVARIANT_MAP = 'Invariant Map',
}

export enum MONOMER_POSITION_PROPERTIES {
  COLOR_COLUMN_NAME = 'color',
  AGGREGATION = 'aggregation',
  TARGET = 'target',
};

/** Structure-activity relationship viewer */
export class MonomerPosition extends DG.JsViewer {
  _titleHost = ui.divText(SELECTION_MODE.MUTATION_CLIFFS, {id: 'pep-viewer-title'});
  _viewerGrid!: DG.Grid;
  _model!: PeptidesModel;
  colorCol: string;
  aggregation: string;
  target: string;

  constructor() {
    super();
    this.target = this.string(MONOMER_POSITION_PROPERTIES.TARGET, null,
      {category: SELECTION_MODE.MUTATION_CLIFFS, choices: []});
    this.colorCol = this.string(MONOMER_POSITION_PROPERTIES.COLOR_COLUMN_NAME, C.COLUMNS_NAMES.ACTIVITY_SCALED,
      {category: SELECTION_MODE.INVARIANT_MAP,
        choices: wu(grok.shell.t.columns.numerical).toArray().map((col) => col.name)});
    this.aggregation = this.string(MONOMER_POSITION_PROPERTIES.AGGREGATION, DG.AGG.AVG,
      {category: SELECTION_MODE.INVARIANT_MAP,
        choices: Object.values(DG.AGG)
          .filter((agg) => ![DG.AGG.KEY, DG.AGG.PIVOT, DG.AGG.SELECTED_ROWS_COUNT].includes(agg))});
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

  get mode(): SELECTION_MODE {
    return this.dataFrame.getTag(C.TAGS.MONOMER_POSITION_MODE) as SELECTION_MODE ??
      SELECTION_MODE.MUTATION_CLIFFS;
  }
  set mode(mode: SELECTION_MODE) {
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
    this.helpUrl = '/help/domains/bio/peptides.md';
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
    if (property.name === MONOMER_POSITION_PROPERTIES.TARGET)
      this.model.updateMutationCliffs();

    this.render();
  }

  createMonomerPositionGrid(): void {
    this.viewerGrid = this.model.monomerPositionDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.MONOMER]);
    this.viewerGrid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...this.model.splitSeqDf.columns.names()]);
    const monomerCol = this.model.monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setAARRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model,
      this.mode === SELECTION_MODE.INVARIANT_MAP, this.dataFrame.getCol(this.colorCol),
      this.aggregation as DG.AggregationType));

    this.viewerGrid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      this.model.highlightMonomerPosition(monomerPosition);
      return showTooltip(monomerPosition, x, y, this.model);
    });
    this.viewerGrid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell?.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        return;

      const monomerPosition = this.getMonomerPosition(gridCell);
      if (this.mode === SELECTION_MODE.INVARIANT_MAP)
        this.model.modifyInvariantMapSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      else if (this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
        this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      this.viewerGrid.invalidate();
    });

    setViewerGridProps(this.viewerGrid, false);
  }

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!) as string,
      positionOrClusterType: gridCell!.tableColumn!.name};
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      if (this.name === VIEWER_TYPE.MONOMER_POSITION) {
        const mutationCliffsMode = ui.boolInput('', this.mode === SELECTION_MODE.MUTATION_CLIFFS);
        mutationCliffsMode.root.addEventListener('click', () => {
          invariantMapMode.value = false;
          mutationCliffsMode.value = true;
          this.mode = SELECTION_MODE.MUTATION_CLIFFS;
        });
        mutationCliffsMode.setTooltip('Statistically significant changes in activity');
        mutationCliffsMode.addPostfix(SELECTION_MODE.MUTATION_CLIFFS);
        const invariantMapMode = ui.boolInput('', this.mode === SELECTION_MODE.INVARIANT_MAP);
        invariantMapMode.root.addEventListener('click', () => {
          mutationCliffsMode.value = false;
          invariantMapMode.value = true;
          this.mode = SELECTION_MODE.INVARIANT_MAP;
        });
        invariantMapMode.setTooltip('Number of sequences having monomer-position');
        invariantMapMode.addPostfix(SELECTION_MODE.INVARIANT_MAP);
        const setDefaultProperties = (input: DG.InputBase): void => {
          $(input.root).find('.ui-input-editor').css('margin', '0px').attr('type', 'radio');
          $(input.root).find('.ui-input-description').css('padding', '0px').css('padding-right', '16px');
        };
        setDefaultProperties(mutationCliffsMode);
        setDefaultProperties(invariantMapMode);

        switchHost = ui.divH([mutationCliffsMode.root, invariantMapMode.root], {id: 'pep-viewer-title'});
        $(switchHost).css('width', 'auto').css('align-self', 'center');
      }
      const tips: HTMLElement = ui.iconFA('question');
      ui.tooltip.bind(tips, () => ui.divV([
        ui.divText('Color intensity - p-value'),
        ui.divText('Circle size - Mean difference'),
        ui.divText('Number - # of unique sequences that form mutation cliffs pairs'),
      ]));

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
export class MostPotentResidues extends DG.JsViewer {
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
    this.helpUrl = '/help/domains/bio/peptides.md';
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

    // Setting Monomer column renderer
    CR.setAARRenderer(monomerCol, this.model.alphabet);
    this.viewerGrid.onCellRender.subscribe(
      (args: DG.GridCellRenderArgs) => renderCell(args, this.model, false, undefined, undefined));

    this.viewerGrid.onCellTooltip((gridCell: DG.GridCell, x: number, y: number) => {
      if (!gridCell.isTableCell) {
        this.model.unhighlight();
        return true;
      }
      const monomerPosition = this.getMonomerPosition(gridCell);
      this.model.highlightMonomerPosition(monomerPosition);
      if (gridCell.tableColumn?.name === C.COLUMNS_NAMES.MONOMER)
        monomerPosition.positionOrClusterType = C.COLUMNS_NAMES.MONOMER;
      else if (gridCell.tableColumn?.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return false;
      return showTooltip(monomerPosition, x, y, this.model);
    });
    this.viewerGrid.root.addEventListener('mouseleave', (_ev) => this.model.unhighlight());
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (!gridCell?.isTableCell || gridCell!.tableColumn!.name !== C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const monomerPosition = this.getMonomerPosition(gridCell);
      if (!this.model.mutationCliffs?.get(monomerPosition.monomerOrCluster)?.get(monomerPosition.positionOrClusterType)?.size)
        return;
      this.model.modifyMutationCliffsSelection(monomerPosition, {shiftPressed: ev.shiftKey, ctrlPressed: ev.ctrlKey});
      this.viewerGrid.invalidate();
    });
    const mdCol: DG.GridColumn = this.viewerGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';
    setViewerGridProps(this.viewerGrid, true);
  }

  getMonomerPosition(gridCell: DG.GridCell): SelectionItem {
    return {monomerOrCluster: gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.MONOMER, gridCell!.tableRowIndex!),
      positionOrClusterType: `${gridCell.cell.dataFrame.get(C.COLUMNS_NAMES.POSITION, gridCell!.tableRowIndex!)}`};
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      const switchHost = ui.divText(VIEWER_TYPE.MOST_POTENT_RESIDUES, {id: 'pep-viewer-title'});
      const tips: HTMLElement = ui.iconFA('question');
      ui.tooltip.bind(tips, () => ui.divV([
        ui.divText('Color intensity - p-value'),
        ui.divText('Circle size - Mean difference'),
        ui.divText('Number - # of unique sequences that form mutation cliffs pairs'),
      ]));

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
  colorCol?: DG.Column<number>, colorAgg?: DG.AggregationType, renderNums?: boolean): void {
  const renderColNames = [...model.splitSeqDf.columns.names(), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
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
  const currentPosStats = model.monomerPositionStats[currentPosition];

  if (!currentPosStats![currentMonomer]) {
    args.preventDefault();
    canvasContext.restore();
    return;
  }

  if (isInvariantMap) {
    const value = currentPosStats![currentMonomer]!.count;
    const positionCol = model.df.getCol(currentPosition);
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

    const color = DG.Color.scaleColor(cellColorDataCol.aggregate(colorAgg!), colorColStats.min, colorColStats.max);
    CR.renderInvaraintMapCell(
      canvasContext, currentMonomer, currentPosition, model.invariantMapSelection, value, bound, color);
  } else {
    CR.renderMutationCliffCell(canvasContext, currentMonomer, currentPosition, model.monomerPositionStats, bound,
      model.mutationCliffsSelection, model.mutationCliffs, model.settings.isBidirectional, renderNums);
  }
  args.preventDefault();
  canvasContext.restore();
}

export function showTooltip(monomerPosition: SelectionItem, x: number, y: number, model: PeptidesModel): boolean {
  if (monomerPosition.positionOrClusterType === C.COLUMNS_NAMES.MONOMER)
    model.showMonomerTooltip(monomerPosition.monomerOrCluster, x, y);
  else
    model.showTooltipAt(monomerPosition, x, y);
  return true;
}

function setViewerGridProps(grid: DG.Grid, isMostPotentResiduesGrid: boolean): void {
  const gridProps = grid.props;
  gridProps.allowEdit = false;
  gridProps.allowRowSelection = false;
  gridProps.allowBlockSelection = false;
  gridProps.allowColSelection = false;
  gridProps.showCurrentCellOutline = false;
  gridProps.showCurrentRowIndicator = false;

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
