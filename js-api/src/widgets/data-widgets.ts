/**
 * Data-related widgets: RangeSlider, HtmlTable, ColumnComboBox, Legend, PropertyGrid.
 * @module widgets/data-widgets
 */

import {toDart, toJs} from "../wrappers";
import * as rxjs from "rxjs";
import {Observable} from "rxjs";
import {__obs, observeStream} from "../events";
import {Property} from "../entities";
import {Column, DataFrame} from "../dataframe";
import {LegendPosition} from "../const";
import {IDartApi} from "../api/grok_api.g";
import {DartWidget} from "./base";
import {RangeSliderStyle, SliderOptions} from "./types";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** A slider that lets user control both min and max values. */
export class RangeSlider extends DartWidget {

  static create(vertical: boolean = false, optionsStyle: RangeSliderStyle | SliderOptions = 'barbell'): RangeSlider {
    let options: SliderOptions = {};
    if (typeof optionsStyle !== 'string') {
      options.style = optionsStyle.style;
    } else {
      options.style = optionsStyle;
    }
    return toJs(api.grok_RangeSlider(options.style, vertical));
  }

  /** Minimum range value. */
  get minRange(): number { return api.grok_RangeSlider_Get_MinRange(this.dart); };

  /** Gets maximum range value. */
  get maxRange(): number { return api.grok_RangeSlider_Get_MaxRange(this.dart); };

  /** Gets minimum value. */
  get min(): number { return api.grok_RangeSlider_Get_Min(this.dart); };

  /** Gets maximum value. */
  get max(): number { return api.grok_RangeSlider_Get_Max(this.dart); };

  /** Sets values to range slider.
   * @param {number} minRange
   * @param {number} maxRange
   * @param {number} min
   * @param {number} max */
  setValues(minRange: number, maxRange: number, min: number, max: number): void {
    api.grok_RangeSlider_SetValues(this.dart, minRange, maxRange, min, max);
  }

  /** Sets showHandles in range slider.
   * @param {boolean} value */
  setShowHandles(value: boolean): void {
    api.grok_RangeSlider_SetShowHandles(this.dart, value);
  }

  /** Sets the specified min value, preserving the range (scrolls). */
  scrollTo(newMinValue: number): void { return api.grok_RangeSlider_ScrollTo(this.dart, newMinValue); }

  /** Shifts min and max values by the specified delta. */
  scrollBy(delta: number): void { return api.grok_RangeSlider_ScrollTo(this.dart, delta); }

  /** @returns {Observable} */
  get onValuesChanged(): Observable<any> {
    return observeStream(api.grok_RangeSlider_Get_OnValuesChanged(this.dart));
  }
}


export class HtmlTable extends DartWidget {

  /** @constructs {HtmlTable} */
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a visual table based on [items], [renderer], and [columnNames] */
  static create(items: any, renderer: any, columnNames: any = null): HtmlTable {
    return toJs(api.grok_HtmlTable(items, renderer !== null ? (object: any, ind: number) => renderer(toJs(object), ind) : null, columnNames));
  }

  /** Removes item */
  remove(item: any): void {
    api.grok_HtmlTable_Remove(this.dart, item);
  }
}


/** A combo box with columns as item.
 * Supports sorting, searching, custom tooltips, column-specific rendering, drag-and-drop, etc */
export class ColumnComboBox extends DartWidget {
  /** @constructs {ColumnComboBox} */
  constructor(dart: any) {
    super(dart);
  }

  /** Creates a column combo box with specified [dataframe] and [predicate]. */
  static create(dataframe: DataFrame, predicate: Function): ColumnComboBox {
    return toJs(api.grok_ColumnComboBox(dataframe.dart, (x: any) => predicate(toJs(x))));
  }

  /** @type {boolean} */
  get vertical(): boolean {
    return api.grok_ColumnComboBox_Get_Vertical(this.dart);
  }

  /** Text to be shown before the combo box */
  get caption(): string { return api.grok_ColumnComboBox_Get_Caption(this.dart); }
  set caption(c: string) { api.grok_ColumnComboBox_Set_Caption(this.dart, c); }

  /** @type {Property} */
  get property(): Property {
    return toJs(api.grok_ColumnComboBox_Get_Property(this.dart));
  }

  onEvent(eventId: string): rxjs.Observable<any> {
    return __obs(eventId, this.dart);
  }

  /** Occurs when the value is changed. */
  get onChanged(): rxjs.Observable<String> { return this.onEvent('d4-column-box-column-changed'); }
}


/** Column legend for viewers */
export class Legend extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

  static create(column: Column): Legend {
    return toJs(api.grok_Legend(column.dart));
  }

  /** Column for the legend */
  get column(): Column { return toJs(api.grok_Legend_Get_Column(this.dart)); }
  set column(column: Column) { api.grok_Legend_Set_Column(this.dart, column.dart); }

  /** Whether to show empty categories */
  get showNulls(): boolean { return api.grok_Legend_Get_ShowNulls(this.dart); }
  set showNulls(show: boolean) { api.grok_Legend_Set_ShowNulls(this.dart, show); }

  /** Position (left / right / top / bottom) */
  get position(): LegendPosition { return api.grok_Legend_Get_Position(this.dart); }
  set position(pos: LegendPosition) { api.grok_Legend_Set_Position(this.dart, pos); }

  set onViewerLegendChanged(handler: Function) {
    api.grok_Legend_Set_OnViewerLegendChanged(this.dart, handler);
  }

  /** Mapped indices of filtered (selected) categories. */
  get selectedCategories(): number[] | null {
    return api.grok_Legend_Get_SelectedCategories(this.dart);
  }

  /** Mapped indices of extra (selected) categories. */
  get selectedExtraCategories(): number[] | null {
    return api.grok_Legend_Get_SelectedExtraCategories(this.dart);
  }

  /** Whether the legend is in tooltip (collapsed) mode. */
  get isTooltipMode(): boolean { return api.grok_Legend_Get_IsToolTipMode(this.dart); }
}


export class PropertyGrid extends DartWidget {

  constructor() {
    super(api.grok_PropertyGrid());
  }

  update(src: any, props: Property[]) {
    api.grok_PropertyGrid_Update(this.dart, src, props.map((x) => toDart(x)));
  }
}

export class GridFilterBase extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

  saveState(): Record<string, any> {
    return toJs(api.grok_GridFilterBase_SaveState(this.dart));
  }

  applyState(state: Record<string, any>): void {
    api.grok_GridFilterBase_ApplyState(this.dart, toDart(state));
  }

  get filterColumnName(): string {
    return api.grok_GridFilterBase_Get_FilterColumnName(this.dart);
  }

  get filterType(): string {
    return api.grok_GridFilterBase_Get_FilterType(this.dart);
  }

  get isActive(): boolean {
    return api.grok_GridFilterBase_Get_IsActive(this.dart);
  }

  set isActive(value: boolean) {
    api.grok_GridFilterBase_Set_IsActive(this.dart, value);
  }

  get isFiltering(): boolean {
    return api.grok_GridFilterBase_Get_IsFiltering(this.dart);
  }

  setSorting(columnName: 'name' | 'count', ascending: boolean): void {
    api.grok_GridFilterBase_Set_Sorting(this.dart, columnName, ascending);
  }

  get caption(): string {
    return api.grok_GridFilterBase_Get_Caption(this.dart);
  }

  set caption(value: string) {
    api.grok_GridFilterBase_Set_Caption(this.dart, value);
  }


}
