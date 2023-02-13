import * as DG from 'datagrok-api/dg';


export enum VISIBILITY_MODE {
  ALWAYS = 'Always',
  AUTO = 'Auto',
  NEVER = 'Never',
};

export type markerType = 'circle' | 'rect' | 'ring' | 'diamond';
export type markerPosition = 'main line' | 'above main line' | 'scatter';
export type visibilityMode = `${VISIBILITY_MODE}`;
export type timePoint = Date | number | null;

export interface Indexable {
  [key: string]: any,
}

export interface ColumnData {
  column: DG.Column,
  data: Int32Array | Float32Array | Float64Array | Uint32Array,
  categories: string[] | null,
}

export interface ColumnsData extends Indexable {
  splitByColumnName?: ColumnData | null,
  startColumnName?: ColumnData | null,
  endColumnName?: ColumnData | null,
  colorByColumnName?: ColumnData | null,
  eventColumnName?: ColumnData | null,
  eventsColumnNames?: {[key: string]: ColumnData},
}
