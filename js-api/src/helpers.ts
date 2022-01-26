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

  get items(): FormulaLine[] {
    let json = this.storage;
    return json ? JSON.parse(json) : [];
  }

  set items(value: FormulaLine[]) { this.storage = JSON.stringify(value); }

  add(item: FormulaLine): void {
    let newItem = this.addDefaults(item);
    this.items = this.items.concat(newItem);
  }

  addAll(items: FormulaLine[]): void {
    let newItems = items.map((item) => { return this.addDefaults(item); });
    this.items = this.items.concat(newItems);
  }

  addLine(item: FormulaLine): void {
    item.type = 'line';
    this.add(item);
  }

  addBand(item: FormulaLine): void {
    item.type = 'band';
    this.add(item);
  }

  clear(): void { this.items = []; }

  addDefaults(item: FormulaLine): FormulaLine { return toJs(api.grok_FormulaLineHelper_AddDefaultParams(item)); }

  getMeta(item: FormulaLine): FormulaLineMeta { return toJs(api.grok_FormulaLineHelper_GetMeta(item)); }
}
