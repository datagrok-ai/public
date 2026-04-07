/**
 * Formula and annotation helpers for DataFrames.
 * @module dataframe/formula-helpers
 */

import {FormulaLinesHelper, AnnotationRegionsHelper} from "../helpers";
import type {DataFrame} from "./data-frame";

declare let DG: any;

export class DataFrameFormulaLinesHelper extends FormulaLinesHelper {
  readonly df: DataFrame;

  get storage(): string { return this.df.getTag(DG.TAGS.FORMULA_LINES) ?? ''; }
  set storage(value: string) { this.df.setTag(DG.TAGS.FORMULA_LINES, value); }

  constructor(df: DataFrame) {
    super();
    this.df = df;
  }
}

export class DataFrameAnnotationRegionsHelper extends AnnotationRegionsHelper {
  readonly df: DataFrame;

  get storage(): string { return this.df.getTag(DG.TAGS.ANNOTATION_REGIONS) ?? ''; }
  set storage(value: string) { this.df.setTag(DG.TAGS.ANNOTATION_REGIONS, value); }

  constructor(df: DataFrame) {
    super();
    this.df = df;
  }
}
