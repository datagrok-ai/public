/**
 * Helper classes for column metadata, colors, markers, and dialogs.
 * @module dataframe/column-helpers
 */

import {TAGS, ColorCodingType, MarkerCodingType, LINK_CLICK_BEHAVIOR} from "../const";
import {Color} from '../color';
import {filter} from "rxjs/operators";
import {Tags} from "../api/ddt.api.g";
import {IDartApi} from "../api/grok_api.g";
import type {Column} from "./column";
import type {MatchType} from "./types";

declare let grok: any;
declare let DG: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export class ColumnDialogHelper {
  private readonly column: Column;
  constructor(column: Column) {
    this.column = column;
  }

  /** Opens an editor dialog with preview for a calculated column. */
  editFormula(): void {
    let formula = this.column.meta.formula;
    // let df = this.column.dataFrame;
    if (formula == null)
      formula = '';
    if (!(this.column.name && this.column.dataFrame?.columns.contains(this.column.name)))
      return;
    let params = { table: this.column.dataFrame, expression: formula, name: this.column.name, type: this.column.type };
    let call = DG.Func.byName('AddNewColumn').prepare(params);
    call.aux['addColumn'] = false;
    call.edit();
    let sub = grok.functions.onAfterRunAction.pipe(filter((c: any) => c == call)).subscribe(() => {
      let newCol = call.getOutputParamValue();
      for (let [key, value] of this.column.tags) if (key !== DG.TAGS.FORMULA) newCol.setTag(key, value);
      let name = newCol.name;
      newCol = this.column.dataFrame.columns.replace(this.column, newCol);
      newCol.name = name;
      api.grok_Column_UpdateDependentColumns(newCol.dart);
      sub.unsubscribe();
    });
  }
}

export class ColumnColorHelper {
  private readonly column: Column;

  constructor(column: Column) {
    this.column = column;
  }

  getType(): ColorCodingType {
    if (this.column.tags.has(DG.TAGS.COLOR_CODING_TYPE))
      return this.column.tags[DG.TAGS.COLOR_CODING_TYPE];
    else if (this.column.tags.has(DG.TAGS.COLOR_CODING_CATEGORICAL))
      return DG.COLOR_CODING_TYPE.CATEGORICAL;
    return DG.COLOR_CODING_TYPE.OFF;
  }

  private _setOutOfRangeLinearColors(belowMinColor?: string, aboveMaxColor?: string): void {
    if (belowMinColor != null)
      this.column.tags[DG.TAGS.COLOR_CODING_LINEAR_BELOW_MIN_COLOR] = belowMinColor;
    if (aboveMaxColor != null)
      this.column.tags[DG.TAGS.COLOR_CODING_LINEAR_ABOVE_MAX_COLOR] = aboveMaxColor;
  }

  /** Enables linear color-coding on a column.
   * @param range - list of palette colors (ARGB integers; see {@link Color}).
   * @param options - list of additional parameters, such as the minimum/maximum value to be used for scaling and the colors for values below the minimum and above the maximum.
   * Use the same numeric representation as [Column.min] and [Column.max].
   */
  setLinear(range: number[] | null = null, options: {min?: number, belowMinColor?: string, max?: number, aboveMaxColor?: string} | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.LINEAR;
    if (range != null)
      this.column.tags[DG.TAGS.COLOR_CODING_LINEAR] = JSON.stringify(range);
    if (options?.min != null)
      this.column.tags[DG.TAGS.COLOR_CODING_SCHEME_MIN] = `${options.min}`;
    if (options?.max != null)
      this.column.tags[DG.TAGS.COLOR_CODING_SCHEME_MAX] = `${options.max}`;
    this._setOutOfRangeLinearColors(options?.belowMinColor, options?.aboveMaxColor);
  }

  /** Enables linear color-coding on a column with absolute values.
   * @param valueColors - dictionary of numerical values and hex-colors.
   * @param options - list of additional parameters, such as the colors for values below the minimum and above the maximum.
   *
   * See samples: {@link https://public.datagrok.ai/js/samples/grid/color-coding/color-coding}}
   */
  setLinearAbsolute(valueColors: {[value: number]: string}, options: {belowMinColor?: string, aboveMaxColor?: string} | null = null): void {
    this.column.tags[TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.LINEAR;
    const orderedEntries = Object.entries(valueColors).sort(([a], [b]) => +a - +b);
    const colors = orderedEntries.map(([_, value]) => Color.fromHtml(value));
    const stringifiedObj = '{' + orderedEntries.map(([key, value]) => `"${key}":"${value}"`).join(',') + '}';
    this.column.tags[TAGS.COLOR_CODING_LINEAR_IS_ABSOLUTE] = 'true';
    this.column.tags[TAGS.COLOR_CODING_LINEAR_ABSOLUTE] = stringifiedObj;
    this.column.tags[TAGS.COLOR_CODING_LINEAR] = JSON.stringify(colors);
    this._setOutOfRangeLinearColors(options?.belowMinColor, options?.aboveMaxColor);
  }

  setCategorical(colorMap: {} | null = null, options: {fallbackColor: string | number, matchType?: MatchType} | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.CATEGORICAL;
    if (colorMap != null)
      this.column.tags[DG.TAGS.COLOR_CODING_CATEGORICAL] = JSON.stringify(colorMap);
    if (options?.fallbackColor != null)
      this.column.tags[DG.TAGS.COLOR_CODING_FALLBACK_COLOR] = JSON.stringify(options.fallbackColor);
    if (options?.matchType != null)
      this.column.tags[DG.TAGS.COLOR_CODING_MATCH_TYPE] = options.matchType;
  }

  setConditional(rules: {[index: string]: number | string}  | null = null): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.CONDITIONAL;
    if (rules != null) {
      for (let [rule, color] of Object.entries(rules)) {
        if (typeof color === 'number')
          rules[rule] = DG.Color.toHtml(color);
      }
      this.column.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = JSON.stringify(rules);
    }
  }

  setDisabled(): void {
    this.column.tags[DG.TAGS.COLOR_CODING_TYPE] = DG.COLOR_CODING_TYPE.OFF;
  }

  getColor(i: number): number {
    return api.grok_Column_GetColor(this.column.dart, i);
  }

  getColors(): Uint32Array {
    return api.grok_Column_GetColors(this.column.dart);
  }
}

export class ColumnMarkerHelper {
  private readonly column: Column;

  constructor(column: Column) {
    this.column = column;
  }

  assign(category: string, marker: MarkerCodingType): ColumnMarkerHelper {
    let jsonTxt: string | null = this.column.getTag(TAGS.MARKER_CODING);
    const jsonMap: {[key: string]: string} = jsonTxt ? JSON.parse(jsonTxt) : {};
    jsonMap[category] = marker;
    jsonTxt = JSON.stringify(jsonMap);
    this.column.setTag(TAGS.MARKER_CODING, jsonTxt);
    return this;
  }

  default(marker: MarkerCodingType): ColumnMarkerHelper {
    return this.assign('~DEFAULT', marker);
  }

  // Obsolete. Recommended method is "assign".
  setMarkerCoding(category: string, marker: MarkerCodingType): void {
    let jsonTxt: string | null = this.column.getTag(TAGS.MARKER_CODING);
    const jsonMap: {[key: string]: string} = jsonTxt ? JSON.parse(jsonTxt) : {};
    jsonMap[category] = marker;
    jsonTxt = JSON.stringify(jsonMap);
    this.column.setTag(TAGS.MARKER_CODING, jsonTxt);
  }

  setAll(categoryMarkerMap: {[key: string]: MarkerCodingType}): ColumnMarkerHelper {
    for (const [category, marker] of Object.entries(categoryMarkerMap))
      this.assign(category, marker);
    return this;
  }

  reset(): ColumnMarkerHelper {
    this.column.tags[TAGS.MARKER_CODING] = '{}';
    return this;
  }
}

export class ColumnMetaHelper {
  private readonly column: Column;
  private _dialogs: ColumnDialogHelper | undefined;
  private _colors: ColumnColorHelper | undefined;
  private _markers: ColumnMarkerHelper | undefined;

  constructor(column: Column) {
    this.column = column;
  }

  get dialogs(): ColumnDialogHelper {
    if (this._dialogs == undefined)
      this._dialogs = new ColumnDialogHelper(this.column);
    return this._dialogs;
  }

  get colors(): ColumnColorHelper {
    if (this._colors == undefined)
      this._colors = new ColumnColorHelper(this.column);
    return this._colors;
  }

  get markers(): ColumnMarkerHelper {
    if (this._markers == undefined)
      this._markers = new ColumnMarkerHelper(this.column);
    return this._markers;
  }

  private setNonNullTag(key: string, value: string | null) {
    if (value === null)
      api.grok_Column_Remove_Tag(this.column.dart, key);
    else
      this.column.setTag(key, value);
  }

  /** Specifies the name to be shown in the UI */
  get friendlyName(): string | null { return this.column.getTag(TAGS.FRIENDLY_NAME); }
  set friendlyName(x: string | null) { this.setNonNullTag(TAGS.FRIENDLY_NAME, x); }

  /** Column description (usually shown in tooltips) */
  get description(): string | null { return this.column.getTag(TAGS.DESCRIPTION); }
  set description(x: string | null) { this.setNonNullTag(TAGS.DESCRIPTION, x); }

  /** Specifies the data format of the dataframe column. See also [GridColumn.format] */
  get format(): string | null { return this.column.getTag(TAGS.FORMAT) ?? api.grok_Column_GetAutoFormat(this.column.dart); }
  set format(x: string | null) { this.setNonNullTag(TAGS.FORMAT, x); }

  /** Returns the maximum amount of significant digits detected in the column. */
  get sourcePrecision(): number | null { return this.column.getTag(TAGS.SOURCE_PRECISION) != null ? +this.column.getTag(TAGS.SOURCE_PRECISION) : null; }

  /** When set, uses the formula to calculate the column values. */
  get formula(): string | null { return this.column.getTag(TAGS.FORMULA); }
  set formula(x: string | null) { this.setNonNullTag(TAGS.FORMULA, x); }

  /** Specifies the units of the dataframe column. */
  get units(): string | null { return this.column.getTag(TAGS.UNITS); }
  set units(x: string | null) { this.setNonNullTag(TAGS.UNITS, x); }

  /** Specifies the units of the dataframe column. */
  get cellRenderer(): string | null { return this.column.getTag(TAGS.CELL_RENDERER); }
  set cellRenderer(x: string | null) { this.setNonNullTag(TAGS.CELL_RENDERER, x); }

  /** When set, switches the cell editor to a combo box that only allows to choose specified values.
   * Applicable for string columns only.
   * See also {@link autoChoices}. */
  get choices(): string[] | null { return JSON.parse(this.column.getTag(TAGS.CHOICES)); }
  set choices(x: string[] | null) { this.setNonNullTag(TAGS.CHOICES, x != null ? JSON.stringify(x) : null); }

  /** When set to 'true', switches the cell editor to a combo box that only allows to choose values
   * from a list of already existing values in the column.
   * Applicable for string columns only.
   * See also {@link choices}. */
  get autoChoices(): boolean { return this.column.getTag(TAGS.AUTO_CHOICES) == null ? false :
    this.column.getTag(TAGS.AUTO_CHOICES).toLowerCase() == 'true'; }
  set autoChoices(x: boolean) { this.column.setTag(TAGS.AUTO_CHOICES, x.toString()); }

  /** Specifies the behavior of link click (open in new tab, open in context panel, custom). Open in new tab is used by default. */
  get linkClickBehavior() : LINK_CLICK_BEHAVIOR {
    return this.column.getTag(TAGS.LINK_CLICK_BEHAVIOR) as LINK_CLICK_BEHAVIOR ?? LINK_CLICK_BEHAVIOR.OPEN_IN_NEW_TAB;
  }
  set linkClickBehavior(x: LINK_CLICK_BEHAVIOR) { this.column.setTag(TAGS.LINK_CLICK_BEHAVIOR, x); }

  /** Specifies whether the column is exported as part of the CSV file. Defaults to true. */
  get includeInCsvExport(): boolean { return this.column.getTag(Tags.IncludeInCsvExport) != 'false'; }
  set includeInCsvExport(x) { this.column.setTag(Tags.IncludeInCsvExport, x.toString()); }

  /** Specifies whether the column is exported as part of the binary file. Defaults to true. */
  get includeInBinaryExport(): boolean { return this.column.getTag(Tags.IncludeInBinaryExport) != 'false'; }
  set includeInBinaryExport(x) { this.column.setTag(Tags.IncludeInBinaryExport, x.toString()); }

  /** When specified, column filter treats the split strings as separate values.
   * See also: https://datagrok.ai/help/visualize/viewers/filters#column-tags */
  get multiValueSeparator(): string { return this.column.getTag(Tags.MultiValueSeparator); }
  set multiValueSeparator(x) { this.setNonNullTag(Tags.MultiValueSeparator, x); }

  get allowColorPicking(): boolean { return this.column.getTag(Tags.AllowColorPicking) != 'false'; }
  set allowColorPicking(x) { this.column.setTag(Tags.AllowColorPicking, x.toString()); }
}
