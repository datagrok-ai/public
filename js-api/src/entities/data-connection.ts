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

/** One failed row of a {@link MutationResult}. */
export interface RowError {
  /** 0-based index of the failing row within the payload (-1 for statement-level errors). */
  index: number;
  /** Column the error is attributed to, when the driver reports it. */
  column?: string;
  /** Machine code: `unique`, `fk`, `notnull`, `check`, `value`, `conflict`, `missing`, `db-error`. */
  code: string;
  message: string;
}

/** The structured result of a database write ({@link DbTable} insert/upsert/update/delete).
 * Mutations never return a DataFrame. `inserted`/`updated`/`skipped` are best-effort — some
 * dialects cannot tell inserts from updates on an upsert and leave the split `null`. */
export interface MutationResult {
  /** Total rows affected across the operation. */
  affectedRows: number;
  /** Per-row failures (populated on partial mode / batch errors). */
  errors?: RowError[];
  inserted?: number | null;
  updated?: number | null;
  skipped?: number | null;
  errorCount?: number | null;
  /** Top-level SQL error message, when the whole operation failed. */
  errorMessage?: string;
  perStatement?: {statementIndex: number, affectedRows: number}[];
  generatedKeys?: Record<string, any>[];
  /** The executed (or, under `dryRun`, planned) DDL plan; DDL and dryRun results only. */
  plan?: DdlPlan;
}

/** Fluent builder for structured `UPDATE`/`DELETE` mutations. A thin wrapper over
 * {@link DbTable} that accumulates `SET` assignments and `WHERE` conditions. In phase A a
 * `WHERE` condition is **equality by value**; a string value may use string patterns
 * (e.g. `contains foo`). Numeric/date range grammar (`> 5`, `10-20`) in `where` is phase B.
 * Build via `grok.data.db.buildMutation(...)`. */
export class TableMutationBuilder {
  private readonly _set: Record<string, any> = {};
  private readonly _where: Record<string, any> = {};
  private _allowFullTable = false;

  constructor(public readonly connectionId: string, public readonly tableName: string) {}

  /** Creates a builder for [tableName] against the [connectionId] connection (nqName). */
  static from(tableName: string, connectionId: string): TableMutationBuilder {
    return new TableMutationBuilder(connectionId, tableName);
  }

  /** Assigns `column = value` in the `SET` clause (update only). */
  set(column: string, value: any): TableMutationBuilder { this._set[column] = value; return this; }

  /** Adds a `WHERE` condition: equality by value (a string value may use string patterns). */
  where(column: string, pattern: any): TableMutationBuilder { this._where[column] = pattern; return this; }

  /** Allows a `DELETE` with no `WHERE` conditions (a full-table wipe — off by default). */
  allowFullTable(): TableMutationBuilder { this._allowFullTable = true; return this; }

  /** Runs the accumulated `UPDATE`. */
  async update(): Promise<MutationResult> {
    return JSON.parse(await api.grok_DbTable_Update(this.connectionId, this.tableName,
      JSON.stringify(this._set), JSON.stringify(this._where)));
  }

  /** Runs the accumulated `DELETE`. */
  async delete(): Promise<MutationResult> {
    return JSON.parse(await api.grok_DbTable_Delete(this.connectionId, this.tableName,
      JSON.stringify(this._where), this._allowFullTable));
  }
}

/** One destructive (or refused) action of a {@link DdlPlan}. */
export interface DestructiveAction {
  /** `table` or `columns.<name>`. */
  path: string;
  /** Machine code: `drop-table`, `truncate-table`, `drop-column`, `type-change`,
   * `not-null-violation`, `required-without-default`. */
  code: string;
  message: string;
  /** Live-data count from the pre-check (rows for table-level actions, non-null
   * values for column-level ones). */
  count?: number;
}

/** The execution plan of a DDL operation ({@link DdlCommand.dryRun}, or riding a
 * {@link DdlConfirmationRequiredError}): the exact per-provider SQL plus live-data
 * pre-checks. The apply path always recomputes the plan server-side — an earlier
 * dry run is informational, never trusted. */
export interface DdlPlan {
  /** The exact SQL statements the operation will run, in order. */
  statements: string[];
  /** Destructive actions that require `confirmDestructive: true` to proceed. */
  destructive: DestructiveAction[];
  /** Whether the provider can roll DDL back on failure (Postgres/MSSQL true;
   * MySQL/Oracle auto-commit DDL, so a failed multi-statement op cannot roll back). */
  transactionalDdl?: boolean;
}

/** Thrown by {@link DdlCommand.execute} when the recomputed plan has destructive
 * actions and `confirmDestructive` was not set. Carries the {@link plan} so the
 * caller can render a confirmation without a second round trip, then re-run
 * `execute({confirmDestructive: true})`. */
export class DdlConfirmationRequiredError extends Error {
  constructor(message: string, public readonly plan: DdlPlan, public readonly provider?: string) {
    super(message);
    this.name = 'DdlConfirmationRequiredError';
  }
}

/** JSON-parse-probes a thrown mutation error: a `destructive-confirmation-required`
 * refusal (whose 400 body rides the thrown message verbatim) becomes a typed
 * {@link DdlConfirmationRequiredError}; anything else is returned untouched. */
function toDdlError(e: any): any {
  const msg: string = e?.message ?? `${e}`;
  const start = msg.indexOf('{'), end = msg.lastIndexOf('}');
  if (start >= 0 && end > start) {
    let body: any = null;
    try {
      body = JSON.parse(msg.substring(start, end + 1));
    } catch (_) {}
    if (body?.error === 'destructive-confirmation-required')
      return new DdlConfirmationRequiredError(body.message ?? msg, body.plan, body.provider);
  }
  return e;
}

/** Column definition for {@link DdlBuilder.createTable} / `alterTable(...).addColumn`.
 * `type` is a dg type (`string`, `int`, `bigint`, `float`, `bool`, `datetime`) mapped
 * to the provider's native type server-side. */
export interface DdlColumn {
  name: string;
  /** Dg type: `string`, `int`, `bigint`, `float`, `bool`, `datetime`. */
  type: string;
  /** Defaults to true. */
  nullable?: boolean;
  /** Always a string; typed by `type` and embedded as an escaped literal
   * (DDL cannot bind parameters). */
  defaultValue?: string;
}

/** Index definition for {@link DdlBuilder.createTable}. */
export interface DdlIndex {
  columns: string[];
  unique?: boolean;
  /** Auto-generated (`ix_<table>_<col>`) when absent. */
  name?: string;
}

/** A single addressed DDL operation, built by {@link DdlBuilder}. `dryRun()` returns
 * the plan without executing; `execute()` recomputes the plan server-side and refuses
 * destructive actions with {@link DdlConfirmationRequiredError} unless confirmed. */
export class DdlCommand {
  constructor(public readonly connectionId: string, private readonly _mutation: Record<string, any>) {}

  /** Returns the execution plan — exact per-provider SQL and live-data destructive
   * pre-checks — without executing anything. Passes the same privilege and capability
   * gates as a real run. */
  async dryRun(): Promise<DdlPlan> {
    return (await this._run({...this._mutation, dryRun: true})).plan!;
  }

  /** Executes the operation. When the recomputed plan has destructive actions
   * (drop/truncate table, drop column, type change) and `confirmDestructive` is not
   * set, rejects with {@link DdlConfirmationRequiredError} carrying the plan. */
  async execute(options: {confirmDestructive?: boolean} = {}): Promise<MutationResult> {
    return this._run(options?.confirmDestructive ? {...this._mutation, confirmDestructive: true} : this._mutation);
  }

  private async _run(mutation: Record<string, any>): Promise<MutationResult> {
    try {
      return JSON.parse(await api.grok_Db_DdlExecute(this.connectionId, JSON.stringify(mutation)));
    } catch (e: any) {
      throw toDdlError(e);
    }
  }
}

/** Single-action `ALTER TABLE` operations, obtained via {@link DdlBuilder.alterTable}.
 * Each method returns a {@link DdlCommand} for exactly one alteration (many databases
 * auto-commit DDL, so batching alterations would fake an atomicity that isn't there) —
 * run several alterations as sequential commands. */
export class AlterTableBuilder {
  constructor(public readonly connectionId: string, private readonly _tableName: string,
    private readonly _schema?: string) {}

  /** Adds a column. A non-nullable column without a `defaultValue` is refused on a
   * populated table (`required-without-default`). */
  addColumn(column: DdlColumn): DdlCommand { return this._op({action: 'addColumn', column}); }

  /** Drops a column — destructive: requires `confirmDestructive` when the column has
   * non-null values. */
  dropColumn(columnName: string): DdlCommand { return this._op({action: 'dropColumn', columnName}); }

  /** Renames a column (native `RENAME COLUMN` — non-destructive). */
  renameColumn(columnName: string, newName: string): DdlCommand {
    return this._op({action: 'renameColumn', columnName, newName});
  }

  /** Changes the column's type to a dg type — destructive: requires `confirmDestructive`
   * when the column has non-null values (castability is not pre-verified; a failing cast
   * surfaces as the database's own error). MySQL and MSSQL restate the whole column on a
   * type change, so `options.nullable` is required there — the server refuses a MySQL/MSSQL
   * `changeType` without an explicit nullability. */
  changeType(columnName: string, newType: string, options?: {nullable?: boolean}): DdlCommand {
    return this._op(options?.nullable !== undefined
      ? {action: 'changeType', columnName, newType, nullable: options.nullable}
      : {action: 'changeType', columnName, newType});
  }

  /** Sets or drops NOT NULL. `setNullable(col, false)` with NULLs present is a hard
   * refusal (`not-null-violation`) — confirmation cannot help. */
  setNullable(columnName: string, nullable: boolean): DdlCommand {
    return this._op({action: 'setNullable', columnName, nullable});
  }

  private _op(fields: Record<string, any>): DdlCommand {
    const m: Record<string, any> = {'#type': 'AlterTable', tableName: this._tableName, ...fields};
    if (this._schema != null)
      m.schema = this._schema;
    return new DdlCommand(this.connectionId, m);
  }
}

/** Fluent builder for structured DDL against one connection, obtained via
 * `grok.data.db.ddl(connectionId)`. Mirrors {@link TableMutationBuilder}; every
 * terminal returns a {@link DdlCommand} carrying exactly one operation.
 *
 * Requires a provider with DDL support (Postgres, MySQL/MariaDB, MSSQL, Oracle —
 * there is no generic-JDBC DDL fallback; others reject with a structured capability
 * error) and the matching fine-grained saved-connection privilege
 * (`DataConnection.CreateTable` for createTable, `DropTable` for dropTable,
 * `TruncateTable` for truncateTable, `AlterSchema` for the rest), enforced
 * server-side. */
export class DdlBuilder {
  private _schema?: string;

  constructor(public readonly connectionId: string) {}

  /** Sets the database schema for subsequently built operations. */
  schema(schema: string): DdlBuilder { this._schema = schema; return this; }

  /** `CREATE TABLE` with optional primary key and indexes. Foreign keys are not part
   * of the create — create first, then {@link addKey}. */
  createTable(tableName: string, columns: DdlColumn[],
    options?: {primaryKey?: string[], indexes?: DdlIndex[], ifNotExists?: boolean}): DdlCommand {
    const fields: Record<string, any> = {columns};
    if (options?.primaryKey !== undefined) fields.primaryKey = options.primaryKey;
    if (options?.indexes !== undefined) fields.indexes = options.indexes;
    if (options?.ifNotExists !== undefined) fields.ifNotExists = options.ifNotExists;
    return this._op('CreateTable', tableName, fields);
  }

  /** Single-action `ALTER TABLE` operations on [tableName]. */
  alterTable(tableName: string): AlterTableBuilder {
    return new AlterTableBuilder(this.connectionId, tableName, this._schema);
  }

  /** `CREATE [UNIQUE] INDEX` on [columns]; the name is auto-generated
   * (`ix_<table>_<col>`) when absent. */
  createIndex(tableName: string, columns: string[], options?: {unique?: boolean, name?: string}): DdlCommand {
    const fields: Record<string, any> = {columns};
    if (options?.unique !== undefined) fields.unique = options.unique;
    if (options?.name !== undefined) fields.indexName = options.name;
    return this._op('CreateIndex', tableName, fields);
  }

  /** `DROP INDEX` by name. */
  dropIndex(tableName: string, indexName: string): DdlCommand {
    return this._op('DropIndex', tableName, {indexName});
  }

  /** Adds a primary- or foreign-key constraint; `refTable`/`refColumns` are foreign-only.
   * The constraint name is auto-generated (`pk_<table>` / `fk_<table>_<col>`) when absent. */
  addKey(tableName: string, spec: {keyType: 'primary' | 'foreign', columns: string[],
    refTable?: string, refColumns?: string[], name?: string}): DdlCommand {
    const fields: Record<string, any> = {keyType: spec.keyType, columns: spec.columns};
    if (spec.refTable !== undefined) fields.refTable = spec.refTable;
    if (spec.refColumns !== undefined) fields.refColumns = spec.refColumns;
    if (spec.name !== undefined) fields.constraintName = spec.name;
    return this._op('AddKey', tableName, fields);
  }

  /** `DROP TABLE` — destructive: requires `execute({confirmDestructive: true})` (the
   * refusal carries the plan with the live row count). */
  dropTable(tableName: string): DdlCommand { return this._op('DropTable', tableName, {}); }

  /** `TRUNCATE TABLE` — destructive, same confirmation contract as {@link dropTable}. */
  truncateTable(tableName: string): DdlCommand { return this._op('TruncateTable', tableName, {}); }

  private _op(type: string, tableName: string, fields: Record<string, any>): DdlCommand {
    const m: Record<string, any> = {'#type': type, tableName, ...fields};
    if (this._schema != null)
      m.schema = this._schema;
    return new DdlCommand(this.connectionId, m);
  }
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
