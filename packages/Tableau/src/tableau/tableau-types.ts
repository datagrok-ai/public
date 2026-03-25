export interface TwbFile {
  version: string;
  datasources: TwbDatasource[];
  worksheets: TwbWorksheet[];
}

export interface TwbDatasource {
  name: string;
  caption: string;
  sourceFile: string;
  sourceDirectory: string;
  columns: TwbColumn[];
}

export interface TwbColumn {
  name: string;
  caption: string;
  datatype: string;
  ordinal: number;
  role: string;
  type: string;
  aggregation: string;
  containsNull: boolean;
}

export interface TwbWorksheet {
  name: string;
  datasourceName: string;
  rows: string;
  cols: string;
  markClass: string;
  usedColumns: string[];
}
