import * as DG from 'datagrok-api/dg';

export const LARGE_DF_THRESHOLD = 50_000;
export const MAX_INDEXES_IN_DESCRIBE = 100;
export const MAX_SAMPLE_COLUMNS = 5;

export type ColumnDescription = {
  name: string;
  type: DG.ColumnType;
  semType: string | null;
  length: number;
  missing: number;
  unique: number;
  units?: string;
  format?: string;
  friendlyName?: string;
  description?: string;
  numerical?: {min: number; max: number; avg: number; stdev: number; med: number};
  categorical?: {topCategories: Array<{value: string; count: number}>};
};

export type CloneDfOpts = {
  rows?: DG.BitSet | 'filtered' | 'selected';
  cols?: string[];
  withSelection?: boolean;
  withTags?: boolean;
};

export type ColumnMetaSpec = {
  semType?: string | null;
  format?: string | null;
  units?: string | null;
  friendlyName?: string | null;
  description?: string | null;
  choices?: string[] | null;
  colorCoding?: ColorCodeSpec;
  tags?: Record<string, string | null>;
};

export type AddColumnSpec = {
  name: string;
  type?: DG.ColumnType;
  formula?: string;
  values?: any[];
  init?: (i: number) => any;
  virtual?: (i: number) => any;
  insertAt?: number;
  ensureUniqueName?: boolean;
  meta?: ColumnMetaSpec;
};

export type FilterCriteria = {
  min?: number;
  max?: number;
  values?: string[];
  contains?: string;
  regex?: string;
  substructure?: string;
  searchType?: 'substructure' | 'similar' | 'exact';
};

export type FilterRowsOpts = {
  ui?: boolean;
  molBlockFailover?: string;
};

export type FilteredDfOpts = {
  cols?: string[];
  withSelection?: boolean;
};

export type FilterSubstructureOpts = {
  ui?: boolean | 'auto';
  molBlockFailover?: string;
  searchType?: 'substructure' | 'similar' | 'exact';
};

export type SelectMode = 'replace' | 'add' | 'remove' | 'intersect';

export type SelectInput =
  | DG.BitSet
  | ((i: number) => boolean)
  | ArrayLike<number>;

export type SelectRowsOpts = {
  mode?: SelectMode;
};

export type SelectedDfOpts = {
  cols?: string[];
  withSelection?: boolean;
};

export type SelectionDescription = {
  count: number;
  total: number;
  indexes: number[];
  currentRowIdx: number;
  sample?: Record<string, unknown>;
};

export type SortOrder = 'asc' | 'desc' | boolean;

export type WidthPolicy = 'Minimal' | 'Compact' | 'Optimal' | 'Maximal';

export type ConfigureGridOpts = {
  hide?: string[];
  show?: string[];
  order?: string[];
  widths?: Record<string, number>;
  widthPolicy?: WidthPolicy;
  pin?: string[];
  unpin?: string[];
  formats?: Record<string, string>;
  rowHeight?: number;
  frozenColumns?: number;
};

export type ColorCodeSpec =
  | {
      kind: 'linear';
      range?: (string | number)[];
      min?: number;
      max?: number;
      belowMinColor?: string;
      aboveMaxColor?: string;
    }
  | {
      kind: 'categorical';
      colors?: Record<string, string | number>;
      fallbackColor?: string | number;
      matchType?: 'exact' | 'regex';
    }
  | {
      kind: 'conditional';
      rules: Record<string, string | number>;
    }
  | {kind: 'off'};

export type ResetGridOpts = {
  visibility?: boolean;
  widths?: boolean;
  sort?: boolean;
  colors?: boolean;
};

export interface JoinOptions {
  leftKey: string;
  rightKey: string;
  joinType?: DG.JoinType;
  normalizeKey?: 'canonicalSmiles' | 'lowercase' | 'trim';
}
