export type MwxColumnType = 'numeric' | 'text' | 'datetime';

export interface MwxColumn {
  name: string;
  description: string;
  type: MwxColumnType;
  values: (number | string | null)[];
  categories?: {[key: string]: number};
  format?: MwxColumnFormat;
}

export interface MwxColumnFormat {
  formatKey: number;
  autoFormat?: boolean;
  numDecPlaces?: number;
  minValue?: number;
  maxValue?: number;
  charCount?: number;
}

export interface MwxWorksheet {
  name: string;
  version: number;
  worksheetId?: string;
  columns: MwxColumn[];
  rowCount: number;
}
