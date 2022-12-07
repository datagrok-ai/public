import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as C from '../utils/constants';
import * as CR from '../utils/cell-renderer';
import {PeptidesModel} from '../model';
import {isGridCellInvalid} from '../utils/misc';

export class SARViewerBase extends DG.JsViewer {
  tempName!: string;
  _viewerGrid!: DG.Grid;
  sourceGrid!: DG.Grid;
  model!: PeptidesModel;
  isPropertyChanging: boolean = false;
  _isVertical = false;

  constructor() {
    super();
  }

  get name(): string {return '';}

  get viewerGrid(): DG.Grid {
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.sourceGrid = this.view?.grid ?? (grok.shell.v as DG.TableView).grid;
    this.model = PeptidesModel.getInstance(this.dataFrame);
    this.subs.push(this.model.onMutationCliffsSelectionChanged.subscribe(() => this.viewerGrid.invalidate()));
    this.helpUrl = '/help/domains/bio/peptides.md';
  }

  detach(): void {this.subs.forEach((sub) => sub.unsubscribe());}

  get isMutationCliffsMode(): string {
    return this.dataFrame.getTag(C.TAGS.SAR_MODE) ?? '1';
  }
  set isMutationCliffsMode(s: string) {
    this.dataFrame.setTag(C.TAGS.SAR_MODE, s);
  }

  render(refreshOnly = false): void {
    if (!refreshOnly) {
      $(this.root).empty();
      let switchHost = ui.divText('Most Potent Residues', {id: 'pep-viewer-title'});
      if (this.name == 'MC') {
        const mutationCliffsMode = ui.boolInput('', this.isMutationCliffsMode === '1');
        mutationCliffsMode.root.addEventListener('click', () => {
          invariantMapMode.value = false;
          mutationCliffsMode.value = true;
          this.isMutationCliffsMode = '1';
          this.model.isInvariantMap = false;
          this.viewerGrid.invalidate();
        });
        mutationCliffsMode.addPostfix('Mutation Cliffs');
        const invariantMapMode = ui.boolInput('', this.isMutationCliffsMode === '0');
        invariantMapMode.root.addEventListener('click', () => {
          mutationCliffsMode.value = false;
          invariantMapMode.value = true;
          this.isMutationCliffsMode = '0';
          this.model.isInvariantMap = true;
          this.viewerGrid.invalidate();
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
      const viewerRoot = this.viewerGrid.root;
      viewerRoot.style.width = 'auto';
      this.root.appendChild(ui.divV([switchHost, viewerRoot]));
    }
    this.viewerGrid?.invalidate();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);

    this.render(true);
  }
}

/** Structure-activity relationship viewer */
export class MonomerPosition extends SARViewerBase {
  _titleHost = ui.divText('Mutation Cliffs', {id: 'pep-viewer-title'});
  _name = 'MC';
  _isVertical = false;

  constructor() {super();}

  get name(): string {return this._name;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMonomerPositionGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
  }

  createMonomerPositionGrid(): void {
    this.viewerGrid = this.model.monomerPositionDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.MONOMER]);
    this.viewerGrid.columns.setOrder([C.COLUMNS_NAMES.MONOMER, ...this.model.splitSeqDf.columns.names()]);
    const monomerCol = this.model.monomerPositionDf.getCol(C.COLUMNS_NAMES.MONOMER);
    CR.setAARRenderer(monomerCol, this.model.alphabet, this.viewerGrid);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model));
    this.viewerGrid.onCellTooltip((cell: DG.GridCell, x: number, y: number) => showTooltip(cell, x, y, this.model));
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name == C.COLUMNS_NAMES.MONOMER)
        return;

      const position = gridCell!.tableColumn!.name;
      const aar = monomerCol.get(gridCell!.tableRowIndex!);
      chooseAction(aar, position, ev.shiftKey, this.model.isInvariantMap, this.model);
    });
    this.viewerGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(this.model.monomerPositionDf, this.model));

    setViewerGridProps(this.viewerGrid, false);
  }
}

/** Vertical structure activity relationship viewer */
export class MostPotentResiduesViewer extends SARViewerBase {
  _name = 'MPR';
  _titleHost = ui.divText('Most Potent Residues', {id: 'pep-viewer-title'});
  _isVertical = true;

  constructor() {
    super();
  }

  get name(): string {return this._name;}

  get viewerGrid(): DG.Grid {
    if (!this._viewerGrid)
      this.createMostPotentResiduesGrid();
    return this._viewerGrid;
  }
  set viewerGrid(grid: DG.Grid) {
    this._viewerGrid = grid;
  }

  onTableAttached(): void {
    super.onTableAttached();
    this.render();
  }

  onPropertyChanged(property: DG.Property): void {
    super.onPropertyChanged(property);
  }

  createMostPotentResiduesGrid(): void {
    this.viewerGrid = this.model.mostPotentResiduesDf.plot.grid();
    this.viewerGrid.sort([C.COLUMNS_NAMES.POSITION]);
    const pValGridCol = this.viewerGrid.col(C.COLUMNS_NAMES.P_VALUE)!;
    pValGridCol.format = '#.000';
    pValGridCol.name = 'P-value';
    const monomerCol = this.model.mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.MONOMER);
    const positionCol = this.model.mostPotentResiduesDf.getCol(C.COLUMNS_NAMES.POSITION);

    // Setting Monomer column renderer
    CR.setAARRenderer(monomerCol, this.model.alphabet, this.viewerGrid);
    this.viewerGrid.onCellRender.subscribe((args: DG.GridCellRenderArgs) => renderCell(args, this.model));
    this.viewerGrid.onCellTooltip((cell: DG.GridCell, x: number, y: number) => showTooltip(cell, x, y, this.model));
    this.viewerGrid.root.addEventListener('click', (ev) => {
      const gridCell = this.viewerGrid.hitTest(ev.offsetX, ev.offsetY);
      if (isGridCellInvalid(gridCell) || gridCell!.tableColumn!.name != C.COLUMNS_NAMES.MEAN_DIFFERENCE)
        return;

      const tableRowIdx = gridCell!.tableRowIndex!;
      const position = monomerCol.get(tableRowIdx);
      const aar = positionCol.get(tableRowIdx);
      chooseAction(aar, position, ev.shiftKey, false, this.model);
      this.viewerGrid.invalidate();
    });
    this.viewerGrid.onCurrentCellChanged.subscribe((_gc) => cellChanged(this.model.mostPotentResiduesDf, this.model));
    const mdCol: DG.GridColumn = this.viewerGrid.col(C.COLUMNS_NAMES.MEAN_DIFFERENCE)!;
    mdCol.name = 'Diff';
    setViewerGridProps(this.viewerGrid, true);
  }
}

function renderCell(args: DG.GridCellRenderArgs, model: PeptidesModel): void {
  const renderColNames = [...model.splitSeqDf.columns.names(), C.COLUMNS_NAMES.MEAN_DIFFERENCE];
  const mdCol = model.monomerPositionStatsDf.getCol(C.COLUMNS_NAMES.MEAN_DIFFERENCE);
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
    return;
  }

  const tableColName = cell.tableColumn?.name;
  const tableRowIndex = cell.tableRowIndex!;
  if (cell.isTableCell && tableColName && tableRowIndex !== null && renderColNames.indexOf(tableColName) !== -1) {
    const cellValue: number | null = cell.cell.value;

    if (cellValue && cellValue !== DG.INT_NULL && cellValue !== DG.FLOAT_NULL) {
      const gridTable = cell.grid.table;
      const currentPosition: string = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ?
        tableColName : gridTable.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);
      const currentAAR: string = gridTable.get(C.COLUMNS_NAMES.MONOMER, tableRowIndex);

      if (model.isInvariantMap) {
        const value: number = model.monomerPositionStatsDf
          .groupBy([C.COLUMNS_NAMES.POSITION, C.COLUMNS_NAMES.MONOMER, C.COLUMNS_NAMES.COUNT])
          .where(`${C.COLUMNS_NAMES.POSITION} = ${currentPosition} and ${C.COLUMNS_NAMES.MONOMER} = ${currentAAR}`)
          .aggregate().get(C.COLUMNS_NAMES.COUNT, 0);
        CR.renderInvaraintMapCell(
          canvasContext, currentAAR, currentPosition, model.invariantMapSelection, value, bound);
      } else {
        CR.renderMutationCliffCell(canvasContext, currentAAR, currentPosition, model.monomerPositionStatsDf,
          mdCol, bound, cellValue, model.mutationCliffsSelection, model.substitutionsInfo,
          model.settings.isBidirectional);
      }
    }
    args.preventDefault();
  }
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
    else if (cell.cell.value && renderColNames.includes(tableColName!)) {
      const currentPosition = tableColName !== C.COLUMNS_NAMES.MEAN_DIFFERENCE ? tableColName :
        table.get(C.COLUMNS_NAMES.POSITION, tableRowIndex);

      model.showTooltipAt(currentAAR, currentPosition, x, y);
    }
  }
  return true;
}

function chooseAction(aar: string, position: string, isShiftPressed: boolean, isInvariantMapSelection: boolean,
  model: PeptidesModel): void {
  if (isShiftPressed)
    model.modifyMonomerPositionSelection(aar, position, isInvariantMapSelection);
  else
    model.initMonomerPositionSelection(aar, position, isInvariantMapSelection);
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
