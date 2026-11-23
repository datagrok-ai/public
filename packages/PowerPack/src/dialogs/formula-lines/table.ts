import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {escapeColName, formatAreaFormula, parseAreaFormula} from './formulas';
import {ItemList} from './item-list';
import {AxisNames, BTN_CAPTION, EditorItem, isAnnotationRegionType, ITEM_TYPE} from './types';

const VISIBLE_COL_WIDTH = 40;
const REMOVE_COL_WIDTH = 35;

/**
 * Table Helper for displaying and navigating Formula Lines list.
 */
export class Table {
  private grid: DG.Grid;

  // Stops the onValuesChanged handler while the grid changes itself
  private notify: boolean = true;

  public get root(): HTMLElement {return this.grid.root;}

  public get currentItem(): EditorItem | null {return this.items.at(this.currentItemIdx);}

  private get dataFrame(): DG.DataFrame {return this.grid.dataFrame!;}

  public get currentItemIdx(): number {return this.dataFrame.currentRowIdx;}
  // Setting only `currentRowIdx` leaves `currentCell.isCell` false (no column), and the
  // grid's `_renderCurrentCell` paints the green stripe only when the current cell is a
  // real cell — so without a column the indicator stays invisible until the user clicks.
  // Using `currentCell` here sets both row and column atomically.
  public set currentItemIdx(rowIdx: number) {
    if (rowIdx >= 0 && rowIdx < this.dataFrame.rowCount) {
      const sameRow = this.dataFrame.currentRowIdx === rowIdx;
      this.dataFrame.currentCell = this.dataFrame.cell(rowIdx, 'title');
      // Defer scroll: on dialog open the grid viewport is still laying out.
      // Skip when already current (e.g. editing the current row) — would jolt the viewport.
      if (!sameRow)
        requestAnimationFrame(() => this.grid.scrollToCell('title', rowIdx));
    } else
      this.dataFrame.currentRowIdx = rowIdx;
  }

  constructor(
    public readonly items: ItemList,
    private onItemChangedAction: (rowIdx: number) => void,
    private onItemRemovedAction: (item: EditorItem) => void,
    srcAxes?: AxisNames,
    setCurrentByAxes?: boolean,
  ) {
    this.grid = this.createGrid();

    this.dataFrame.onCurrentRowChanged.subscribe((_) => this.onItemChangedAction(this.currentItemIdx));
    this.dataFrame.onValuesChanged.subscribe((ed) => this.onValuesChanged(ed.args?.indexes));

    if (setCurrentByAxes && srcAxes?.x && srcAxes?.y) {
      const rowIdx = this.findRowByAxes(srcAxes.x, srcAxes.y);
      if (rowIdx !== -1) {
        DG.delay(1).then((_) => {
          this.currentItemIdx = rowIdx;
        });
      }
    }
  }

  public setFirstItemAsCurrent() {
    if (this.items.length > 0) {
      this.currentItemIdx = 0;
      this.onItemChangedAction(0);
    } else
      this.onItemChangedAction(-1);
  }

  /** Shows the changed title and formula of the item in the [rowIdx] row and makes the row current. */
  public update(rowIdx: number) {
    const {title, formula} = this.cellValues(rowIdx);
    this.notify = false;
    this.dataFrame.set('title', rowIdx, title);
    this.dataFrame.set('formula', rowIdx, formula);
    this.notify = true;
    this.currentItemIdx = rowIdx;
  }

  /** Adds the item to the list and makes its row current; returns the row. */
  public add(item: EditorItem): number {
    const rowIdx = this.items.add(item);
    const {title, formula, visible} = this.cellValues(rowIdx);
    this.dataFrame.rows.insertAt(rowIdx);
    this.notify = false;
    this.dataFrame.set('title', rowIdx, title);
    this.dataFrame.set('formula', rowIdx, formula);
    this.dataFrame.set('visible', rowIdx, visible);
    this.dataFrame.set(BTN_CAPTION.REMOVE, rowIdx, '');
    this.notify = true;
    this.currentItemIdx = rowIdx;
    return rowIdx;
  }

  // What the row of an item shows
  private cellValues(rowIdx: number): {title?: string, formula?: string, visible?: boolean} {
    if (this.items.isFormulaLine(rowIdx)) {
      const line = this.items.formulaLines[rowIdx];
      return {title: line.title, formula: line.formula, visible: line.visible};
    }
    const region = this.items.at(rowIdx) as DG.AnnotationRegion;
    return {title: region.header, formula: formatAreaFormula(region), visible: !region.hidden};
  }

  private createGrid(): DG.Grid {
    const dataFrame = DG.DataFrame.fromColumns([
      DG.Column.fromList('string', 'title', []),
      DG.Column.fromList('string', 'formula', []),
      DG.Column.fromList('bool', 'visible', []),
    ]);

    for (let i = 0; i < this.items.length; i++) {
      const {title, formula, visible} = this.cellValues(i);
      dataFrame.rows.addNew([title ?? '', formula, visible]);
    }

    /** Column for "trash" buttons */
    dataFrame.columns.addNewString(BTN_CAPTION.REMOVE);

    const grid = DG.Grid.create(dataFrame);
    grid.setOptions({
      showCurrentRowIndicator: true,
      showSelectedRows: false,
      allowRowResizing: false,
      allowRowSelection: false,
      allowBlockSelection: false,
      allowColSelection: false,
      allowRowReordering: false,
      showRowHeader: false,
      showCellTooltip: false,
      showColumnTooltip: false,
      showCurrentCellOutline: false,
      showContextMenu: false,
      showEditRow: false,
    });

    grid.columns.setVisible(['title', 'formula', 'visible', BTN_CAPTION.REMOVE]);
    grid.columns.setOrder(['title', 'formula', 'visible', BTN_CAPTION.REMOVE]);

    grid.col('title')!.width = 120;
    grid.col('formula')!.width = 220;
    grid.col('visible')!.width = VISIBLE_COL_WIDTH;

    const deleteBtnCol = grid.col(BTN_CAPTION.REMOVE)!;
    deleteBtnCol.width = REMOVE_COL_WIDTH;
    deleteBtnCol.cellType = 'html';

    // Stretch title + formula to fill the grid width; visible + remove keep fixed widths.
    // Buffer accounts for the vertical scrollbar + cell padding so the row doesn't
    // overflow horizontally when the grid is at its initial size.
    const stretchColumns = () => {
      const total = grid.root.clientWidth;
      if (total <= 0)
        return;
      const remaining = total - VISIBLE_COL_WIDTH - REMOVE_COL_WIDTH - 24;
      if (remaining < 200)
        return;
      const titleW = Math.floor(remaining * 0.35);
      grid.col('title')!.width = titleW;
      grid.col('formula')!.width = remaining - titleW;
    };
    new ResizeObserver(stretchColumns).observe(grid.root);
    DG.delay(0).then(stretchColumns);

    grid.onCellPrepare((cell) => {
      if (cell.isColHeader) {
        switch (cell.gridColumn.name) {
        case 'title': cell.customText = 'Title'; break;
        case 'formula': cell.customText = 'Formula'; break;
        case 'visible': cell.customText = 'Show'; break;
        case BTN_CAPTION.REMOVE: cell.style.textColor = 0; break;
        }
      } else if (cell.isTableCell && cell.gridColumn.name === BTN_CAPTION.REMOVE)
        cell.style.element = this.deleteBtn(cell.tableRowIndex!); // the grid may be sorted
    });

    return grid;
  }

  private deleteBtn(rowIdx: number): HTMLElement {
    const btn = ui.button(ui.iconFA('trash-alt'), () => {
      const item = this.items.at(rowIdx)!;
      this.items.removeAt(rowIdx);
      this.dataFrame.rows.removeAt(rowIdx);
      if (this.currentItemIdx > rowIdx)
        this.currentItemIdx--;

      this.onItemRemovedAction(item);
      this.onItemChangedAction(this.currentItemIdx);
    }, 'Delete');
    btn.style.textAlign = 'center';
    btn.style.height = '20px';
    return btn;
  }

  // Writes the cells edited in the grid back to the item
  private onValuesChanged(indexes?: ArrayLike<number>): void {
    if (!this.notify)
      return;

    if (indexes && indexes.length === 1 && indexes[0] !== -1)
      this.currentItemIdx = indexes[0];

    const rowIdx = this.currentItemIdx;
    const item = this.currentItem;
    if (!item)
      return;

    if (isAnnotationRegionType(item.type)) {
      const region = item as DG.AnnotationRegion;
      region.header = this.dataFrame.get('title', rowIdx);
      region.hidden = !this.dataFrame.get('visible', rowIdx);
      parseAreaFormula(region, this.dataFrame.get('formula', rowIdx));
    } else {
      const line = item as DG.FormulaLine;
      line.title = this.dataFrame.get('title', rowIdx);
      line.formula = this.dataFrame.get('formula', rowIdx);
      line.visible = this.dataFrame.get('visible', rowIdx);
    }
    this.onItemChangedAction(rowIdx);
  }

  // The row of the first item on the [x] and [y] axes, or -1
  private findRowByAxes(x: string, y: string): number {
    // Escape axis names so nested braces (e.g. `${AGE} + 10`) match their escaped form in the formula.
    const escapedX = escapeColName(x);
    const escapedY = escapeColName(y);
    const checkAxesInFormula = (formula?: string) =>
      (formula ?? '').includes(escapedX) && (formula ?? '').includes(escapedY);

    const lineIdx = this.items.formulaLines.findIndex((item) => checkAxesInFormula(item.formula));
    if (lineIdx !== -1)
      return lineIdx;

    const regionIdx = this.items.annotationRegions.findIndex((item) => {
      if (item.type === ITEM_TYPE.AREA_REGION_ANNOTATION) {
        const region = item as DG.AreaAnnotationRegion;
        return region.x === x && region.y === y;
      }
      const region = item as DG.FormulaAnnotationRegion;
      return checkAxesInFormula(region.formula1) || checkAxesInFormula(region.formula2);
    });
    return regionIdx === -1 ? -1 : this.items.rowIdx(regionIdx, false);
  }
}
