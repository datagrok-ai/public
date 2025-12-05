import * as DG from 'datagrok-api/dg';

export type ReferenceObject = {
  [tableName: string]: { [columnName: string]: { refTable: string; refColumn: string } };
};
export type ReferencedByObject = {
  [tableName: string]: { [columnName: string]: { refTable: string; refColumn: string }[] };
};
export type SchemaInfo = { references: ReferenceObject; referencedBy: ReferencedByObject };

export type SchemaAndConnection = {
  schema: SchemaInfo;
  connection: DG.DataConnection;
}

export type EntryPointOptions = {
  valueConverter: (value: string | number) => string | number;
  joinOptions: QueryJoinOptions[];
  regexpExample?: (typeof DG.ObjectHandler.prototype.regexpExample)
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
        regexpExample?: EntryPointOptions['regexpExample'];
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

type COLUMN_TYPE = string;
type TABLE_COLUMN_SEPARATED_BY_DOT = `${string}.${string}`;
export type DBConnectionDescriptor = {
  connectionID: string;
  name: string;
  provider?: string;
  description: string;
  schemas: {
    [schemaName: string]: {
      description: string;
      tables: {[tableName: string]: {
        columns: {[columnName: string]: COLUMN_TYPE};
        }
      }
      references: {
        from: TABLE_COLUMN_SEPARATED_BY_DOT;
        to: TABLE_COLUMN_SEPARATED_BY_DOT;
      }[];
    }
  }
}
