import * as DG from 'datagrok-api/dg';

import * as rxjs from 'rxjs';

import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {MultipleSelection} from './SAR-multiple-selection';
import {FilteringStatistics, Stats as FilterStats} from './filtering-statistics';

/**
 * Implements multiple filtering callbacks to be used in events subscription.
 */
export class SARMultipleFilter {
  protected selection: MultipleSelection;
  protected resources: Resources;
  protected filterMode: boolean;
  protected stats: FilteringStatistics;

  /**
   * Creates an instance of SARMultipleFilter.
   * @param {boolean} filterMode Whether to filter or to select.
   */
  constructor(filterMode: boolean) {
    this.filterMode = filterMode;
    this.resources = new Resources();
    this.selection = new MultipleSelection();
    this.stats = new FilteringStatistics();
  }

  /** Makes label to display filtered residue-positions on accordion. */
  get filterLabel() {
    return this.selection.toString();
  }

  /** Makes query to display data frame statistics on accordion. */
  get query(): string {
    this.resources.assert(['residueColumnName']);
    return this.selection.toQuery(this.resources.residueColumnName);
  }

  /**
   * Updates filtering mode.
   * @param {boolean} value Filtering if true or selection else.
   */
  set filteringMode(value: boolean) {
    this.filterMode = value;
  }

  /** Returns bit set obtained from selection mask. */
  get mask() {
    this.resources.assert(['dataFrame', 'groupMapping']);

    const df = this.resources.dataFrame!;
    const mask = this.selection.eval(df, this.resources.groupMapping!);
    return DG.BitSet.create(mask.length, (i) => mask[i]);
  }

  /** Resets selection mask. */
  resetSelection() {
    this.selection = new MultipleSelection();
  }

  /**
   * Adds/updates one of objects required.
   * @param {ResourceKey} key Resource key.
   * @param {any} value Resource value.
   */
  addResource(key: ResourceKey, value: any) {
    this.resources.add(key, value);
  }

  /**
   * Evaluates selection into mask.
   * @return {boolean[]} List of booleans where a true corresponds to the row selected.
   */
  eval(): boolean[] {
    this.resources.assert(['dataFrame', 'groupMapping']);
    return this.selection.eval(this.resources.dataFrame!, this.resources.groupMapping!);
  }

  /**
   * Reacts on mouse event.
   * @param {MouseEventType} eventType Mouse event type.
   * @param {CellType} cellType Type of the grid cell.
   * @param {(string | null)} colName Name of column the cell from.
   * @param {(number | null)} rowIdx Index of the row the cell at.
   * @param {boolean} ctrlPressed Whether CTRL-button was pressed while catching event.
   */
  onSARGridMouseEvent(
    eventType: MouseEventType,
    cellType: CellType,
    colName: string | null,
    rowIdx: number | null,
    ctrlPressed: boolean,
  ) {
    switch (eventType) {
    case 'mouseout':
      this.onSARCellMouseOut();
      break;
    case 'click':
      if (colName)
        this.onSARCellClicked(colName, rowIdx, ctrlPressed);
      break;
    case 'mousemove':
      this.onSARCellHover(colName, rowIdx);
      break;
    }
  }

  /**
   * Reacts on click to the SAR viewer cell.
   * @param {string} colName Name of the column which the cell is from.
   * @param {(number | null)} rowIdx Index of the row at which the cell is.
   * @param {boolean} ctrlPressed Whether the CTRL button was pressed.
   */
  onSARCellClicked(colName: string, rowIdx: number | null, ctrlPressed: boolean) {
    this.resources.assert(['dataFrame', 'residueColumnName', 'grid']);

    const residueColumnName = this.resources.residueColumnName!;
    const grid = this.resources.grid!;
    const isMultiposSel = colName == residueColumnName;
    const isMultiresSel = rowIdx === null;

    if (isMultiposSel && isMultiresSel)
      return;

    const changeFilter = ctrlPressed ? this.selection.input : this.selection.set;
    const df = grid.dataFrame!;
    const columns = df.columns as DG.ColumnList;
    let pos: string | undefined;
    let res: string | undefined;

    if (!isMultiposSel)
      pos = colName;
    if (!isMultiresSel)
      res = df.get(residueColumnName, rowIdx);

    if (isMultiposSel)
      this.selection.setRes(res!, columns.names().slice(1)); // All except the first residues column.
    else if (isMultiresSel)
      this.selection.setPos(pos!, df.col(residueColumnName)!.toList());
    else
      changeFilter.bind(this.selection)(pos!, res!);

    console.warn([ctrlPressed ? 'ctrl+click' : 'click', this.selection.filter]);

    if (this.filterMode)
      this.resources.dataFrame!.rows.requestFilter();
    else
      this.maskRows();

    grid.invalidate();
  }

  /**
   * Reacts on mouse hover over SAR viewer cell.
   * @param {(string | null)} colName Name of the column which the cell is from.
   * @param {(number | null)} rowIdx Index of the row at which the cell is.
   */
  onSARCellHover(colName: string | null, rowIdx: number | null) {
    this.resources.assert(['residueColumnName', 'grid', 'dataFrame']);

    const residueColumnName = this.resources.residueColumnName!;

    if (!colName || !rowIdx || colName == residueColumnName)
      return;

    const df = this.resources.grid?.dataFrame!;

    if (!(df.columns as DG.ColumnList).names().includes(colName))
      return;

    const res = df.get(residueColumnName, rowIdx);

    this.resources.dataFrame?.rows.match({[colName]: res}).highlight();
  }

  /** Reacts on mouse out from the SAR viewer */
  onSARCellMouseOut() {
    this.resources.assert(['dataFrame']);
    this.resources.dataFrame?.rows.match({}).highlight();
  }

  /** Sets up SAR viewer grid visualization. */
  setupGridVizualization() {
    this.resources.assert(['grid']);
    this.resources.grid?.setOptions({
      'showCurrentRowIndicator': false,
      'allowBlockSelection': false,
      'allowColSelection': false,
      'allowRowSelection': false,
      'showMouseOverRowIndicator': false,
    });
  }

  /**
   * Custom renderrer to be applied to SAR viewer cell rendering.
   * @param {DG.GridCellRenderArgs} args Cell rendering arguments.
   */
  selectionCellRenderrer(args: DG.GridCellRenderArgs) {
    if (args.cell.isTableCell) {
      this.resources.assert(['residueColumnName', 'groupMapping']);

      const residueColumnName = this.resources.residueColumnName!;
      const cell = args.cell;
      const pos = cell.gridColumn.name;

      if (pos !== residueColumnName) {
        const rowIdx = cell.tableRowIndex!;
        const res = cell.cell.dataFrame.get(residueColumnName, rowIdx);

        if (this.selection.test(pos, res, this.resources.groupMapping!)) {
          const ctx = args.g;
          ctx.lineWidth = 1;
          ctx.strokeStyle = 'black';
          ctx.strokeRect(args.bounds.x, args.bounds.y, args.bounds.width, args.bounds.height);
        }
      }
    }
  }

  /** Reacts on SAR grid update/reloading. */
  onSARGridChanged() {
    this.resources.assert(['grid']);

    const grid = this.resources.grid!;

    addGridMouseHandler(grid, this.onSARGridMouseEvent.bind(this));
    this.setupGridVizualization();
    grid.onCellRender.subscribe(this.selectionCellRenderrer.bind(this));
  }

  /** Selects or filters rows depending from the filtering mode set. */
  maskRows() {
    this.resources.assert(['dataFrame']);

    const df = this.resources.dataFrame!;

    if (this.filterMode) {
      df.filter.and(this.mask);
      df.selection.setAll(false, false);

      console.warn(['onRowsFiltering', this.selection.filter, df.filter.trueCount]);
    } else {
      df.selection.copyFrom(this.mask);
      df.filter.fireChanged();

      console.warn(['onSelectionChanged', this.selection.filter, df.selection.trueCount]);
    }
  }

  /**
   * Calculates statistics on the activity column.
   * @return {FilterStats} Statistics.
   */
  getStatistics(): FilterStats {
    this.resources.assert(['activityColumnName', 'dataFrame']);

    const df = this.resources.dataFrame!;

    this.stats.setData(df.col(this.resources.activityColumnName!)?.getRawData() as Float32Array);
    this.stats.setMask(this.mask);
    return this.stats.result;
  }
}

/** Declares resources needed by the filtering. */
interface FilterResources {
  residueColumnName?: string;
  activityColumnName?: string;
  grid?: DG.Grid;
  dataFrame?: DG.DataFrame;
  groupMapping?: StringDictionary;
};

type ResourceKey = keyof FilterResources;

/** Resources controlling helper */
class Resources {
  protected data: FilterResources;

  /** Creates an instance of Resources. */
  constructor() {
    this.data = {};
  }

  /**
   * Adds/updates the resource given.
   * @param {ResourceKey} resource Resource key.
   * @param {any} value Resource value.
   */
  add(resource: ResourceKey, value: any) {
    this.data[resource] = value;
  }

  /**
   * Checks if all requested resources were added previously.
   * @param {ResourceKey[]} resources Requested resources.
   * @return {boolean} True if all those resources are initialized.
   */
  enough(resources: ResourceKey[]): boolean {
    return resources.every((k) => this.data[k] !== undefined);
  }

  /**
   * Throws an error if one of the resources given was not set.
   * @param {ResourceKey[]} [resources] Optional set of resources to check.
   */
  assert(resources?: ResourceKey[]) {
    if (resources) {
      if (!this.enough(resources))
        throw new Error(`Not enough one of ${resources} or more.`);
    } else {
      if (!this.ready())
        throw new Error(`Resources are not ready.`);
    }
  }

  /**
   * Checks if all resources were initialized.
   * @return {boolean} True if the test was successful.
   */
  ready(): boolean {
    return Object.values(this.data).every((v) => v !== undefined);
  }

  /** Data frame resource. */
  get dataFrame() {
    return this.data.dataFrame;
  }

  /** Grid resource. */
  get grid() {
    return this.data.grid;
  }

  /** Group mapping. */
  get groupMapping() {
    return this.data.groupMapping;
  }

  /** Residue column name. */
  get residueColumnName() {
    return this.data.residueColumnName;
  }

  /** Activity column name. */
  get activityColumnName() {
    return this.data.activityColumnName;
  }
}

const MouseEventsSource = {
  click: 0,
  mousemove: 1,
  mouseout: 2,
};
const MouseEvents = Object.keys(MouseEventsSource);
export type MouseEventType = keyof typeof MouseEventsSource;
export type CellType = 'isTableCell' | 'isColHeader' | 'unknown';

  type MouseEventHandler = (
    eventType: MouseEventType,
    cellType: CellType,
    colName: string | null,
    rowIdx: number | null,
    ctrlPressed: boolean
  ) => void;

/**
   * Adds mouse event handler to the click event bus.
   * @param {DG.Grid} grid Grid to add to.
   * @param {MouseEventHandler} handler Event handler.
   */
export function addGridMouseHandler(grid: DG.Grid, handler: MouseEventHandler) {
  const onMouseEvent = (mouseEvent: MouseEvent) => {
    if (!MouseEvents.includes(mouseEvent.type))
      return;

    const mouseEventType: MouseEventType = mouseEvent.type as MouseEventType;
    const keyPressed = mouseEvent.ctrlKey || mouseEvent.metaKey;
    const cell = grid.hitTest(mouseEvent.offsetX, mouseEvent.offsetY);
    let pos: string | null = null;
    let rowIdx: number | null = null;
    let cellType: CellType = 'unknown';

    if (cell) {
      pos = cell.gridColumn.name;

      if (pos.length == 0)
        pos = null;

      if (cell.isTableCell) {
        rowIdx = cell.tableRowIndex!;
        cellType = 'isTableCell';
      } else if (cell.isColHeader)
        cellType = 'isColHeader';
    }

    handler(mouseEventType, cellType, pos, rowIdx, keyPressed);
  };

  for (const e of MouseEvents)
    rxjs.fromEvent<MouseEvent>(grid.overlay, e).subscribe(onMouseEvent);
}
