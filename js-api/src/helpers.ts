import {toJs} from './wrappers';
import {IDartApi} from "./api/grok_api.g";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

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
  orientation?: string;
  width?: number;
  spline?: number;
  style?: string;
  // Specific to bands:
  column2?: string;
  xMap?: string;
  yMap?: string;
}

/**
 * Represents a core annotation region parameters.
 */
export interface AnnotationRegion {
  type?: string;
  header?: string;
  headerColor?: string;
  fillColor?: string;
  outlineWidth?: number;
  outlineColor?: number;
  xMap?: string;
  yMap?: string;
  hidden?: boolean;
  description?: string;
  opacity?: number;
  isDataFrameRegion?: boolean;
}

/**
 * Represents an area annotation region (defined by polygon points) parameters.
 */
export interface AreaAnnotationRegion extends AnnotationRegion {
  x?: string;
  y?: string;
  area?: [number, number][];
}

/**
 * Represents a formula annotation region (defined by two formula lines) parameters.
 */
export interface FormulaAnnotationRegion extends AnnotationRegion {
  formula1?: string;
  formula2?: string;
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

  static getMetaByFormula(formula: string, type: string): FormulaLineMeta { return toJs(api.grok_FormulaLineHelper_GetMetaByFormula(formula, type)); }
}

export abstract class AnnotationRegionsHelper {
  // The place where Annotation Regions are saved must be defined in the descendant class.
  abstract get storage(): string;
  abstract set storage(value: string);

  get items(): AnnotationRegion[] {
    let json = this.storage;
    return json ? JSON.parse(json) : [];
  }

  set items(value: AnnotationRegion[]) {
    this.storage = JSON.stringify(value);
  }

  add(item: AnnotationRegion) {
    this.items = this.items.concat(item);
  }

  addAll(items: AnnotationRegion[]) {
    this.items = this.items.concat(items);
  }

  updateAt(idx: number, value: AnnotationRegion) {
    let newItems = this.items;
    newItems[idx] = value;
    this.items = newItems;
  }

  removeAt(idx: number, count: number = 1) {
    this.items = this.items.slice(idx, idx + count - 1);
  }

  removeWhere(predicate: (value: AnnotationRegion, index: number, array: AnnotationRegion[]) => boolean) {
    this.items = this.items.filter(predicate);
  }

  clear(): AnnotationRegionsHelper {
    this.items = [];
    return this;
  }
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

  public static camelCaseToSentence(s: string, o?: {capitalizeFirst?: boolean, capitalizeNext?: boolean, capitalizeConjunctions?: boolean}) {
    return api.grok_StringUtils_CamelCaseToSentence(s, o?.capitalizeFirst ?? true, o?.capitalizeNext ?? false, o?.capitalizeConjunctions ?? false);
  }
}