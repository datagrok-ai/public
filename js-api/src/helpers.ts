import {toJs} from './wrappers';

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
  showOnPlot?:boolean;
  showOnTooltip?:boolean;
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

  get items(): FormulaLine[] {
    let json = this.storage;
    return json ? JSON.parse(json) : [];
  }

  set items(value: FormulaLine[]) {
    this.storage = JSON.stringify(value);
  }

  add(item: FormulaLine) {
    let newItem = FormulaLinesHelper.setDefaults(item);
    this.items = this.items.concat(newItem);
  }

  addAll(items: FormulaLine[]) {
    let newItems = items.map((item) => { return FormulaLinesHelper.setDefaults(item); });
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
    let newItems = this.items;
    newItems[idx] = value;
    this.items = newItems;
  }

  removeAt(idx: number, count: number = 1) {
    this.items = this.items.slice(idx, idx + count - 1);
  }

  removeWhere(predicate: (value: FormulaLine, index: number, array: FormulaLine[]) => boolean) {
    this.items = this.items.filter(predicate);
  }

  clear(): FormulaLinesHelper {
    this.items = [];
    return this;
  }

  static setDefaults(item: FormulaLine): FormulaLine { return toJs(api.grok_FormulaLineHelper_SetDefaultParams(item)); }

  static getMeta(item: FormulaLine): FormulaLineMeta { return toJs(api.grok_FormulaLineHelper_GetMeta(item)); }
}

export class StringUtils {
  public static hashCode(s: string): number {
    let hash: number = 0;
    if (s.length === 0)
      return hash;
    for (let i: number = 0; i < s.length; i++) {
      const chr: number = s.charCodeAt(i);
      hash = ((hash << 5) - hash) + chr;
      hash |= 0; // Convert to 32bit integer
    }
    return hash;
  }
}