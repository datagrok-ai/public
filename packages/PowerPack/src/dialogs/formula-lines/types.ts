import * as DG from 'datagrok-api/dg';

/** Formula Line types */
export const enum ITEM_TYPE {
  LINE = 'line',
  BAND = 'band',
  AREA_REGION_ANNOTATION = 'area',
  FORMULA_REGION_ANNOTATION = 'formula',
}

/** Formula Line captions */
export const enum ITEM_CAPTION {
  LINE = 'Line',
  BAND = 'Band',
  CONST_LINE = 'Constant line',
  VERT_LINE = 'Line - Vertical',
  HORZ_LINE = 'Line - Horizontal',
  HORZ_BAND = 'Band - Horizontal',
  VERT_BAND = 'Band - Vertical',
  FORMULA_REGION = 'Region - Formula Lines',
  RECT_REGION = 'Region - Draw Rectangle',
  POLYGON_REGION = 'Region - Draw Lasso',
  // Axis-bounds formula region — mirrors the axis context-menu's `Annotations > Region`.
  // VERT bounds the X axis (vertical band-shape); HORZ bounds the Y axis.
  VERT_RANGE_REGION = 'Region - Vertical Range',
  HORZ_RANGE_REGION = 'Region - Horizontal Range',
}

export const enum ITEM_ORIENTATION {
  VERTICAL = 'Vertical',
  HORIZONTAL = 'Horizontal',
}

/** Formula Line sources */
export const enum ITEM_SOURCE {
  VIEWER = 'Viewer',
  DATAFRAME = 'DataFrame',
}

export const enum BTN_CAPTION {
  ADD_NEW = 'Add new',
  CLONE = 'Clone',
  REMOVE = 'Remove',
  FORMULA_LINES_HISTORY = 'Formula Lines History',
  ANNOTATION_REGIONS_HISTORY = 'Annotation Regions History',
}

export type EditorItem = DG.FormulaLine | DG.AnnotationRegion;

export interface AxisNames {
  y?: string,
  x?: string,
  yMap?: string,
  xMap?: string,
}

export interface AxisColumns {
  y: DG.Column | null,
  x: DG.Column | null,
  yMap?: string,
  xMap?: string,
  /** Optional stats override. When set, default formula values (q1/q2/q3) come
   *  from here instead of `y.stats` / `x.stats`. Aggregated hosts (bar chart)
   *  use this so defaults reflect the displayed value-axis range, while `y` / `x`
   *  still carry the raw column whose name is referenced in the formula. */
  yStats?: DG.Stats,
  xStats?: DG.Stats,
}

export interface EditorOptions {
  allowEditDFLines: boolean,
  [propertyName: string]: any,
}

export const DEFAULT_OPTIONS: EditorOptions = {
  allowEditDFLines: true,
};

export const isAnnotationRegionType = (type: string | undefined): boolean =>
  type === ITEM_TYPE.AREA_REGION_ANNOTATION || type === ITEM_TYPE.FORMULA_REGION_ANNOTATION;
