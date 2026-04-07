/**
 * Filter abstract base class.
 * @module widgets/filter
 */

import {filter} from 'rxjs/operators';
import {Column, DataFrame} from "../dataframe";
import {Widget} from "./base";

declare let ui: any;

/** Base class for DataFrame-bound filtering controls.
 * Supports collaborative filtering by efficiently working together with
 * other filters. */
export abstract class Filter extends Widget {

  /** An indicator icon on the left of the filter header, before the name.
   * A filter is responsible for hiding or showing it, depending on its state. */
  indicator: HTMLDivElement;

  /** Group of control icons on the right of the filter header, after the name.
   * FilterGroup takes care of its visibility. */
  controls: HTMLDivElement;

  /** A DataFrame this filter is associated with. */
  dataFrame: DataFrame | null;

  /** A column this filter is associated with. */
  column: Column | null = null;

  columnName?: string;

  /** Caption to be shown in the filter group. */
  get caption(): string { return this.columnName ?? ''}

  constructor() {
    super(ui.div());

    this.dataFrame = null;
    this.indicator = ui.div([], 'd4-filter-indicator');
    this.controls = ui.div([], 'd4-flex-row');

    this.indicator.style.display = 'none';
  }

  // static create(type: FILTER_TYPE | null, column: Column): Filter {
  //   return api.grok_Filter_ForColumn(column, type);
  // }
  //
  // static forColumn(column: Column) { return Filter.create(null, column); }
  //
  // static histogram(column: Column) { return Filter.create(FILTER_TYPE.HISTOGRAM, column); }
  // static categorical(column: Column) { return Filter.create(FILTER_TYPE.CATEGORICAL, column); }
  // static multiValue(column: Column) { return Filter.create(FILTER_TYPE.MULTI_VALUE, column); }
  // static freeText(column: Column) { return Filter.create(FILTER_TYPE.FREE_TEXT, column); }

  /** Override to indicate whether the filter actually filters something (most don't in the initial state).
   * This is used to minimize the number of unnecessary computations.
   * Make sure to call super.isFiltering to check whether the filter has been disabled by the user */
  get isFiltering(): boolean {
    return !(this.root.parentElement?.classList?.contains('d4-filter-disabled') == true) &&
      !(this.root.parentElement?.parentElement?.parentElement?.classList?.contains('d4-filters-disabled') == true);
  }

  /** Whether a filter is ready to apply the filtering mask synchronously. */
  get isReadyToApplyFilter(): boolean { return true; }

  /** Override to provide short filter summary that might be shown on viewers or in the context panel. */
  abstract get filterSummary(): string;

  /** Override to filter the dataframe.
   * The method should work with `this.dataFrame.filter`, should disregard
   * false values (these are filtered out already by other filters), and should
   * filter out corresponding indexes. */
  abstract applyFilter(): void;

  /** Override to save filter state. */
  saveState(): any {
    return {
      column: this.columnName,
      columnName: this.columnName
    };
  }

  /** Override to load filter state. */
  applyState(state: any): void {
    this.columnName = state.columnName;
    this.column = this.columnName && this.dataFrame ? this.dataFrame.col(this.columnName) : null;
  }

  /** Gets called when a data frame is attached.
   * Make sure to call super.attach(dataFrame) when overriding.
   * @param {DataFrame} dataFrame*/
  attach(dataFrame: DataFrame): void {
    this.dataFrame = dataFrame;
    this.subs.push(this.dataFrame.onRowsFiltering
      .pipe(filter((_) => this.isFiltering && this.isReadyToApplyFilter))
      .subscribe((_) => {
        this.applyFilter();
        this.dataFrame!.rows.filters.push(`${this.columnName}: ${this.filterSummary}`);
      })
    );
  }

  /** Gets called when a previously used filter gets moved in the DOM tree.
   * Normally, you don't need to do anything, but this might be handy for
   * the iframe-based filters. */
  refresh() {}

  detach(): void {
    super.detach();
    if (this.isFiltering)
      this.dataFrame?.rows?.requestFilter();
  }
}
