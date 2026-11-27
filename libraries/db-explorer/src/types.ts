import * as DG from 'datagrok-api/dg';
import {DB_EXPLORER_OBJ_HANDLER_TYPE} from './object-handlers';

export type ReferenceObject = {
  [schemaName: string]: {
    [tableName: string]: {
      /** [label] names the row the column points at, as seen from the column's row (ExplicitReference.label) */
      [columnName: string]: { refSchema: string, refTable: string; refColumn: string, label?: string }
    }
  }
};
export type ReferencedByObject = {
  [schemaName: string]: {
    [tableName: string]: {
      /** [reverseLabel] names the rows that point back, in the Links section (ExplicitReference.reverseLabel) */
      [columnName: string]: { refSchema: string, refTable: string; refColumn: string, reverseLabel?: string }[]
    }
  }
};
export type SchemaInfo = {
  references: ReferenceObject;
  referencedBy: ReferencedByObject;
  /** Table comments, from the database or from Datagrok's annotations of it, by schema and table */
  comments?: {[schemaName: string]: {[tableName: string]: string}};
};

/** A column that points at a key, or the key itself */
export type LinkRef = { schema: string; table: string; column: string; reverseLabel?: string };

/** A group of the Links section: the rows that point at one key value of the card */
export type LinkGroup = {
  /** The key is the card's own (This drug), not one a column of the card points at (Same target) */
  own: boolean;
  /** The key the rows point at */
  target: LinkRef;
  value: string | number;
  /** The card's column the value is read from */
  field: string;
  /** ExplicitReference.label of the reference the value is read through */
  label?: string;
  links: LinkRef[];
};

/** How the Links section names and shows the rows of a table. Every field is optional: the rest is inferred */
export type EntityOptions = {
  /** One record, e.g. 'Compound' (default: the identifier's semantic type, or the table's name) */
  label?: string;
  /** Code shown after the name, e.g. 'CPD' (default: the fixed part of the identifier's pattern) */
  short?: string;
  /** More than one record (default: the label with an s) */
  plural?: string;
  /** Shown in the tooltip of the table's links (default: the table's comment) */
  description?: string;
  /** Fold the table's links under More, or keep them in the list (default: folded when the configuration names
   * card columns for several tables and none, nor an identifier, a key or a header column, for this one) */
  more?: boolean;
  /** Leave the table's links out */
  hidden?: boolean;
};

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
  /** The row the column points at, from the column's row: 'Parent compound' (default: read from the column's name) */
  label?: string;
  /** The rows that point back, in the Links section of the row pointed at: 'Child compounds' */
  reverseLabel?: string;
}

export type DBExplorerConfig = {
    connectionName: string;
    /** Catalog the schemas live in, for a database that has catalogs (Databricks, Snowflake, MS SQL).
     * When omitted, tables are addressed as schema.table and resolve in the connection's own catalog */
    catalogName?: string;
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
    /** Names and display of tables in the Links section, by table (or schema.table) */
    entities?: {[tableName: string]: EntityOptions};
}

export class DBValueObject {
  name: string = DB_EXPLORER_OBJ_HANDLER_TYPE;
  constructor(
    public connectionNqName: string,
    public schemaName: string,
    public table: string,
    public column: string,
    public value: string | number,
    public semValue?: DG.SemanticValue,
    public catalogName?: string
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
