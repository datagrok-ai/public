import * as DG from 'datagrok-api/dg';

export type ReferenceObject = {
  [tableName: string]: { [columnName: string]: { refTable: string; refColumn: string } };
};
export type ReferencedByObject = {
  [tableName: string]: { [columnName: string]: { refTable: string; refColumn: string }[] };
};
export type SchemaInfo = { references: ReferenceObject; referencedBy: ReferencedByObject };

export type EntryPointOptions = {
  valueConverter: (value: string | number) => string | number;
  joinOptions: QueryJoinOptions[];
};
export type QueryJoinOptions = {
  fromTable: string;
  tableName: string;
  columnName: string;
  onColumn: string;
  select: string[];
};

export type DBExplorerConfig = {
    connectionName: string;
    schemaName: string;
    nqName?: string;
    dataSourceName?: string;
    entryPoints: {[semtype: string]: {
        table: string;
        column: string;
    }}
    joinOptions?: QueryJoinOptions[];
    headerNames?: {[tableName: string]: string};
    uniqueColumns?: {[tableName: string]: string};
    customSelectedColumns?: {[tableName: string]: string[]};
}

export class DBValueObject {
  constructor(
    public table: string,
    public column: string,
    public value: string | number,
    public semValue?: DG.SemanticValue
  ) {}
}
