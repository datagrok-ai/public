/**
 * DataConnection, DataQuery, TableQuery, TableQueryBuilder, DataJob, and Credentials classes.
 * @module entities/data-connection
 */

import {AGG, AggregationType, JOIN_TYPE, JoinType, Type, TYPE} from "../const";
import {toDart, toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import {DataFrame} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import {Func} from "./func";
import {TableInfo} from "./table-info";
import {DataConnectionProperties, FieldOrder, FieldPredicate, GroupAggregation, TableJoin} from "./types";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** Represents a data query
 * @extends Func
 * {@link https://datagrok.ai/help/access/access#data-query}
 * */
export class DataQuery extends Func {
  /** @constructs DataQuery*/
  constructor(dart: any) {
    super(dart);
  }

  /** @deprecated Use FuncCall.adHoc instead **/
  get adHoc(): boolean { return api.grok_Query_Get_AdHoc(this.dart); }
  /** @deprecated Use FuncCall.adHoc instead **/
  set adHoc(a: boolean) { api.grok_Query_Set_AdHoc(this.dart, a); }

  /** Query text */
  get query(): string { return api.grok_Query_Query(this.dart); }
  set query(q: string) { api.grok_Query_Set_Query(this.dart, q); }

  get connection(): DataConnection { return toJs(api.grok_Query_Get_Connection(this.dart)); }
  set connection(c: DataConnection) { api.grok_Query_Set_Connection(this.dart, toDart(c)); }

  get postProcessScript(): string {return api.grok_Query_Get_PostProcessScript(this.dart);}
  set postProcessScript(script: string) { api.grok_Query_Set_PostProcessScript(this.dart, script);}
  /** Executes query
   * @returns {Promise<DataFrame>} */
  async executeTable(): Promise<DataFrame> { return toJs(await api.grok_Query_ExecuteTable(this.dart)); }
}

/** Represents a table query
 * @extends DataQuery */
export class TableQuery extends DataQuery {
  /** @constructs TableQuery */
  constructor(dart: any) { super(dart); }

  /** Creates a TableQuery
   * @param {DataConnection} connection - DataConnection to query table from
   * @returns {TableQuery} */
  static create(connection: DataConnection): TableQuery {
    return toJs(api.grok_TableQuery_Create(connection.dart));
  }

  /** Table name
   * @type {string} */
  get table(): string { return api.grok_TableQuery_GetTable(this.dart); }
  set table(tableName: string) { api.grok_TableQuery_SetTable(this.dart, tableName); }

  /** Fields array
   * @type {string[]} */
  get fields(): string[] { return api.grok_TableQuery_GetFields(this.dart); }
  set fields(fields: string[]) { api.grok_TableQuery_SetFields(this.dart, fields); }

  /** Executes query */
  async executeTable(): Promise<DataFrame> { return toJs(await api.grok_TableQuery_ExecuteTable(this.dart)); }

  /** Where clauses */
  get where(): FieldPredicate[] { return toJs(api.grok_TableQuery_Get_WhereClausesDB(this.dart)); }
  set where(wl: FieldPredicate[]) { api.grok_TableQuery_Set_WhereClausesDB(this.dart, wl.map(param => toDart(param))); }

  /** Aggregation clauses {queryPartParams} */
  get aggregations(): GroupAggregation[] { return toJs(api.grok_TableQuery_Get_AggregationsDB(this.dart)); }
  set aggregations(wl: GroupAggregation[]) { api.grok_TableQuery_Set_AggregationsDB(this.dart, wl.map(param => toDart(param))); }

  get joins(): TableJoin[] { return toJs(api.grok_TableQuery_Get_Joins(this.dart)); }
  set joins(j: TableJoin[]) { api.grok_TableQuery_Set_Joins(this.dart, j.map(join => toDart(join))); }

  /** Having clauses */
  get having(): FieldPredicate[] { return toJs(api.grok_TableQuery_Get_HavingDB(this.dart)); }
  set having(wl: FieldPredicate[]) { api.grok_TableQuery_Set_HavingDB(this.dart, wl.map(param => toDart(param))); }

  /** Order By clauses */
  get orderBy(): FieldOrder[] { return toJs(api.grok_TableQuery_Get_OrderByDB(this.dart)); }
  set orderBy(wl: FieldOrder[]) { api.grok_TableQuery_Set_OrderByDB(this.dart, wl.map(param => toDart(param))); }

  set limit(rows: number | undefined) { api.grok_TableQuery_Set_Limit(this.dart, rows); }
  get limit(): number | undefined { return api.grok_TableQuery_Get_Limit(this.dart); }

  /** Creates {@link TableQueryBuilder} from table name
   * @param {string} table - Table name
   * @returns {TableQueryBuilder} */
  static from(table: string): TableQueryBuilder {return toJs(api.grok_TableQuery_From(table)); }

  /** Creates {@link TableQueryBuilder} from {@link TableInfo}
   * @param {TableInfo} table - TableInfo object
   * @returns {TableQueryBuilder} */
  static fromTable(table: TableInfo): TableQueryBuilder {return toJs(api.grok_TableQuery_FromTable(table.dart)); }
}

/** Table query builder that works with database tables */
export class TableQueryBuilder {
  dart: any;

  /** @constructs TableQueryBuilder */
  constructor(dart: any) { this.dart = dart; }

  /** Creates {@link TableQueryBuilder} from table name
   * @param {string} table - Table name
   * @param connection - {@link DataConnection} that {@link TableQuery} will use after the build. Can be passed lately directly to {@link TableQuery}.
   * @returns {TableQueryBuilder} */
  static from(table: string, connection?: DataConnection | string): TableQueryBuilder {
    let conn: any = connection === undefined ? connection : connection instanceof DataConnection ? connection.dart : connection;
    return toJs(api.grok_DbTableQueryBuilder_From(table, conn));
  }

  /** Creates {@link TableQueryBuilder} from {@link TableInfo}
   * @param {TableInfo} table - TableInfo object
   * @returns {TableQueryBuilder} */
  static fromTable(table: TableInfo): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_FromTable(table.dart));
  }

  /** Selects all fields of the table
   * @returns {TableQueryBuilder} */
  selectAll(): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_SelectAll(this.dart)); }

  /** Selects specified fields of the table. All fields will be selected if not provided.
   * @param {string[]} fields - Array of fields to select
   * @returns {TableQueryBuilder} */
  select(fields: string[]): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_Select(this.dart, fields)); }

  /**
   * Performs aggregation
   * @param {AggregationType} type - Aggregation type.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  selectAggr(type: AggregationType, field: string | null, fieldAlias: string | null): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_SelectAggr(this.dart, type, field, fieldAlias));
  }

  /**
   * Adds an aggregation that calculates average value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  avg(field: string, fieldAlias: string | null = 'avg'): TableQueryBuilder {
    return this.selectAggr(AGG.AVG, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates minimum value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  min(field: string, fieldAlias: string | null = 'min'): TableQueryBuilder {
    return this.selectAggr(AGG.MIN, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates maximum value for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  max(field: string, fieldAlias: string | null = 'max'): TableQueryBuilder {
    return this.selectAggr(AGG.MAX, field, fieldAlias);
  }

  /**
   * Adds an aggregation that calculates sum of the values for the specified column.
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  sum(field: string, fieldAlias: string | null = 'sum'): TableQueryBuilder {
    return this.selectAggr(AGG.SUM, field, fieldAlias);
  }

  /**
   * Adds an aggregation that counts rows. Equivalent to count(*).
   * @param {string} fieldAlias - Name of the resulting column. Default value is count.
   */
  count(fieldAlias: string = 'count'): TableQueryBuilder {
    return this.selectAggr(AGG.TOTAL_COUNT, null, fieldAlias);
  }

  /** Adds an aggregation that counts rows for the specified column. Equivalent to count(field).
   * @param field - Column name.
   * @param fieldAlias Name of the resulting column. Default value is count.
   * @returns {GroupByBuilder} */
  valueCount(field: string, fieldAlias: string = 'count'): TableQueryBuilder {
    return this.selectAggr(AGG.VALUE_COUNT, field, fieldAlias);
  }

  /**
   * Adds an aggregation that counts rows with null values
   * @param {string} field - Column name.
   * @param {string} fieldAlias - Name of the resulting column. Default value is agg(colName).
   */
  nulls(field: string, fieldAlias: string | null = 'nulls'): TableQueryBuilder {
    return this.selectAggr(AGG.MISSING_VALUE_COUNT, field, fieldAlias);
  }

  /**
   * Performs join operation of main table or table specified in {@link leftTable} to table specified in {@link rightTable}.
   * Specify joining fields of the main table (or table specified in {@link leftTable}) in {@link leftTableKeys} and joining fields of {@link rightTable}
   * in {@link rightTableKeys}.
   * @param rightTable {string}
   * @param joinType {JoinType}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param rightTableAlias {string} - use to specify desired alias for the joining table and apply this alias to fields specified in {@link select}.
   * @param leftTable {string}
   */
  join(rightTable: string, joinType: JoinType, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Join(this.dart, rightTable, joinType, leftTableKeys, rightTableKeys, rightTableAlias ?? null, leftTable ?? null));
  }

  /**
   * Performs left join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  leftJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.LEFT, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs right join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  rightJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.RIGHT, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs inner join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  innerJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.INNER, leftTableKeys, rightTableKeys, leftTable)
  }

  /**
   * Performs outer join operation.
   * @param rightTable {string}
   * @param rightTableAlias {string}
   * @param leftTableKeys {string[]}
   * @param rightTableKeys {string[]}
   * @param leftTable {string}
   */
  outerJoin(rightTable: string, leftTableKeys: string[], rightTableKeys: string[], rightTableAlias?: string, leftTable?: string): TableQueryBuilder {
    return this.join(rightTable, JOIN_TYPE.OUTER, leftTableKeys, rightTableKeys, leftTable)
  }

  /** Groups rows that have the same values into summary values
   * @param {string[]} fields - Array of fields to group by
   * @returns {TableQueryBuilder} */
  groupBy(fields: string[]): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_GroupBy(this.dart, fields));
  }

  /** Rotates a table-valued expression by turning the unique values from one column in the expression into multiple
   * columns in the output
   * @param {string[]} fields - Array of fields to pivot on
   * @returns {TableQueryBuilder} */
  pivotOn(fields: string[]): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_PivotOn(this.dart, fields));
  }

  /** Adds a where clause to the query.
   * @param {string} field - Field name. If you join to other tables use correct alias or table name as prefix, e.g. <table alias>.<field>
   * @param {string} pattern - Pattern to test field values against
   * @param columnType Should be provided, if this builder wasn't created from {@link TableInfo} or alias is used. Defaults to {@link TYPE.STRING}
   * @returns {TableQueryBuilder} */
  where(field: string, pattern: string, columnType?: Type): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Where(this.dart, field, pattern, columnType ?? null));
  }

  /** Adds a having clause to the query. Use only with {@link groupBy}
   * @param {string} field - Field name. If you join to other tables use correct alias or table name as prefix, e.g. <table alias>.<field>
   * @param {string} pattern - Pattern to test field values against
   * @param columnType Should be provided, if this builder wasn't created from {@link TableInfo} or alias is used. Defaults to {@link TYPE.STRING}
   * @returns {TableQueryBuilder} */
  having(field: string, pattern: string, columnType?: Type): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_Having(this.dart, field, pattern, columnType ?? null));
  }

  /** Sorts results in ascending or descending order
   * @param {string} field - Field to sort based on
   * @param {boolean} asc - Sort in ascending order
   * @returns {TableQueryBuilder} */
  sortBy(field: string, asc: boolean = true): TableQueryBuilder {
    return toJs(api.grok_DbTableQueryBuilder_SortBy(this.dart, field, asc));
  }

  /** Selects limited number of records
   * @param {TableQueryBuilder} n - Number of records to select
   * @returns {TableQueryBuilder} */
  limit(n: number): TableQueryBuilder { return toJs(api.grok_DbTableQueryBuilder_Limit(this.dart, n)); }

  /** Builds a query
   * @returns {TableQuery} */
  build(): TableQuery { return toJs(api.grok_DbTableQueryBuilder_Build(this.dart)); }
}

/** Represents a data job
 * @extends Func
 * {@link https://datagrok.ai/help/access}
 * */
export class DataJob extends Func {
  /** @constructs DataJob */
  constructor(dart: any) {
    super(dart);
  }
}

/** Represents a data connection
 * @extends Entity
 * {@link https://datagrok.ai/help/access/access#data-connection}
 * */
export class DataConnection extends Entity {
  parameters: any;

  /** @constructs DataConnection */
  constructor(dart: any) {
    super(dart);
    this.parameters = new MapProxy(api.grok_DataConnection_Get_Parameters(this.dart), 'parameters');
  }

  get credentials(): Credentials {
    return toJs(api.grok_DataConnection_Get_Credentials(this.dart));
  }

  get dataSource(): string {
    return api.grok_DataConnection_Get_DataSource(this.dart);
  }

  /** Collection of parameters: server, database, endpoint, etc. */
  // get parameters(): DataConnectionParams { return api.grok_DataConnection_Parameters(this.dart); }

  /** Tests the connection, returns "ok" on success or an error message on error
   * @returns {Promise<string>}*/
  test(): Promise<string> {
    return api.grok_DataConnection_Test(this.dart);
  }

  /**
   * Creates {@link DataQuery} using this connection. Can be used only with database connections.
   * @param name - name of the query
   * @param sql - text of the query
   */
  query(name: string, sql: string): DataQuery {
    return toJs(api.grok_DataConnection_Query(this.dart, name, sql));
  }

  buildQuery(table: string): TableQueryBuilder {
    return TableQueryBuilder.from(table, this);
  }

  tableQuery(): TableQuery {
    return toJs(api.grok_TableQuery_Create(this.dart));
  }

  /** Creates a data connection. Note that in order to be used, it has to be saved first using {@link DataConnectionsDataSource}
   * @param {string} name - Connection name
   * @param {DataConnectionProperties} parameters - Connection properties
   * */
   static create(name: string, parameters: DataConnectionProperties): DataConnection {
    return toJs(api.grok_DataConnection_Create(name, parameters.dataSource, parameters));
  }
}

/** Represents connection credentials
 *  Usually it is a login and a password pair
 *  Passwords are stored in the secured credentials storage
 *  See also: {@link https://datagrok.ai/help/datagrok/solutions/enterprise/security#credentials}
 *  */
export class Credentials extends Entity {
  /**
   * Represents credentials parameters that are hidden from other users and when they are used real values
   * of them will be replaced by values from {@link openParameters}.
   * Parameters can be filled with key-value entries that are suitable for the entity that owns that Credentials.
   */
  public parameters: any;

  constructor(dart: any) {
    super(dart);
    this.parameters = new MapProxy(api.grok_Credentials_Parameters(this.dart), 'parameters');
  }

  /**
   * Represents opened credential parameters. They will be showed, for example, in the connection edit form.
   * They are updated when {@link parameters} are changed and instance is saved.
   * See also {@link CredentialsDataSource}.
   */
  get openParameters(): Record<string, string> { return api.grok_Credentials_OpenParameters(this.dart); }
}
