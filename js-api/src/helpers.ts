import {toJs} from "./wrappers";

let api = <any>window;

/**
 * FormulaLine available parameters.
 */
export interface FormulaLine {
  id?: string;
  type?: string;
  title?: string;
  description?: string;
  color?: string;
  visible?: boolean;
  opacity?: number;
  zIndex?: number;
  min?: number;
  max?: number;
  formula?: string;
  // Specific to lines:
  width?: number;
  spline?: number;
  style?: string;
  // Specific to bands:
  column2?: string;
}

/**
 * FormulaLine meta information.
 */
export interface FormulaLineMeta {
  expression?: string;
  funcName?: string;
  argName?: string;
}

/**
 * FormulaLinesHelper provides methods for working with Formula Lines.
 * Primarily intended for use in Viewer and DataFrame metadata.
 * Documentation:  {@link https://datagrok.ai/help/develop/how-to/show-formula-lines}
 */
export abstract class FormulaLinesHelper {
  // The place where Formula Lines are saved must be defined in the descendant class.
  abstract get storage(): string;
  abstract set storage(value: string);

  private _isBuffering: boolean = false;
  private _bufferedItems: FormulaLine[] = [];

  get items(): FormulaLine[] {
    if (this._isBuffering)
      return this._bufferedItems;
    let json = this.storage;
    return json ? JSON.parse(json) : [];
  }

  set items(value: FormulaLine[]) {
    if (this._isBuffering)
      this._bufferedItems = value;
    else
      this.storage = JSON.stringify(value);
  }

  add(item: FormulaLine) {
    let newItem = this.setDefaults(item);
    if (this._isBuffering)
      this._bufferedItems.push(item);
    else
      this.items = this.items.concat(newItem);
  }

  addAll(items: FormulaLine[]) {
    let newItems = items.map((item) => { return this.setDefaults(item); });
    if (this._isBuffering)
      this._bufferedItems.push(...newItems);
    else
      this.items = this.items.concat(newItems);
  }

  addLine(item: FormulaLine) {
    item.type = 'line';
    this.add(item);
  }

  addBand(item: FormulaLine) {
    item.type = 'band';
    this.add(item);
  }

  updateAt(idx: number, value: FormulaLine) {
    if (this._isBuffering)
      this._bufferedItems[idx] = value;
    else {
      let newItems = this.items;
      newItems[idx] = value;
      this.items = newItems;
    }
  }

  removeAt(idx: number, count: number = 1) {
    if (this._isBuffering)
      this._bufferedItems = this._bufferedItems.slice(idx, idx + count - 1);
    else
      this.items = this.items.slice(idx, idx + count - 1);
  }

  removeWhere(predicate: (value: FormulaLine, index: number, array: FormulaLine[]) => boolean) {
    if (this._isBuffering)
      this._bufferedItems = this._bufferedItems.filter(predicate);
    else
      this.items = this.items.filter(predicate);
  }

  clear() {
    if (this._isBuffering)
      this._bufferedItems = [];
    else
      this.items = [];
  }

  setDefaults(item: FormulaLine): FormulaLine { return toJs(api.grok_FormulaLineHelper_SetDefaultParams(item)); }

  getMeta(item: FormulaLine): FormulaLineMeta { return toJs(api.grok_FormulaLineHelper_GetMeta(item)); }

  /**
   * When buffering is enabled, all changes made are saved in the buffer.
   * This increases the performance of the viewer, because it does not processing formulas immediately after changes.
   * It is recommended to use buffering when the number of lines depends on the number of dataframe rows.
   * Note: buffering is disabled by default. And any changes are immediately taken into processing by the viewer.
   */
  startBuffering() {
    if (this._isBuffering)
      return;
    this._bufferedItems = this.items;
    this._isBuffering = true;
  }

  /**
   * Copies all changes made in the buffer to the Formula Lines store.
   * From that moment, all changes are taken into processing by the viewer.
   */
  stopBuffering() {
    if (!this._isBuffering)
      return;
    this._isBuffering = false;
    this.items = this._bufferedItems;
  }

  cancelBuffering() {
    this._isBuffering = false;
  }
}
