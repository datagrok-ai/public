export type MwxColumnType = 'numeric' | 'text' | 'datetime';

export interface MwxColumn {
  name: string;
  description: string;
  type: MwxColumnType;
  values: (number | string | null)[];
  categories?: {[key: string]: number};
}

export interface MwxWorksheet {
  name: string;
  version: number;
  columns: MwxColumn[];
  rowCount: number;
}
