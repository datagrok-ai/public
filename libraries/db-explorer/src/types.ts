import * as DG from 'datagrok-api/dg';
import {DB_EXPLORER_OBJ_HANDLER_TYPE} from './object-handlers';

export type ReferenceObject = {
  [schemaName: string]: {
    [tableName: string]: {
      [columnName: string]: { refSchema: string, refTable: string; refColumn: string }
    }
  }
};
export type ReferencedByObject = {
  [schemaName: string]: {
    [tableName: string]: {
      [columnName: string]: { refSchema: string, refTable: string; refColumn: string }[]
    }
  }
};
export type SchemaInfo = { references: ReferenceObject; referencedBy: ReferencedByObject };

export type SchemaAndConnection = {
  schema: SchemaInfo;
  connection: DG.DataConnection;
}

export type SupportedRenderer = 'rawImage' | 'imageURL' | 'molecule' | 'helm';

export type EntryPointOptions = {
  valueConverter: (value: string | number) => string | number;
  joinOptions: QueryJoinOptions[];
  regexpExample?: (typeof DG.ObjectHandler.prototype.regexpExample)
  matchRegexp?: string;
};
export type QueryJoinOptions = {
  fromSchema?: string;
  onSchema?: string;
  fromTable: string;
  tableName: string;
  columnName: string;
  onColumn: string;
  select: string[];
};

export type ExplicitReference = {
  schema?: string;
  table: string;
  column: string;
  refSchema?: string;
  refTable: string;
  refColumn: string;
}

export type DBExplorerConfig = {
    connectionName: string;
    /** Primary schema, that can be overriden by schemas in entry points */
    schemaName: string;
    nqName?: string;
    dataSourceName?: string;
    entryPoints: {[semtype: string]: {
        schama?: string;
        table: string;
        column: string;
        regexpExample?: EntryPointOptions['regexpExample'];
        matchRegexp?: string;
    }}
    joinOptions?: QueryJoinOptions[];
    headerNames?: {[tableName: string]: string};
    uniqueColumns?: {[tableName: string]: string};
    customSelectedColumns?: {[tableName: string]: string[]};
    customRenderers?: {
        table: string,
        column: string,
        renderer: SupportedRenderer;
    }[];
    explicitReferences?: ExplicitReference[];
}

export class DBValueObject {
  name: string = DB_EXPLORER_OBJ_HANDLER_TYPE;
  constructor(
    public connectionNqName: string,
    public schemaName: string,
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
