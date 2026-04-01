/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BuiltinDBInfoMeta, getDBColumnMetaData, getDBTableMetaData} from './query-meta-utils';

const suspiciousSQlPatterns = ['DROP ', 'DELETE ', 'UPDATE ', 'INSERT ', 'ALTER ', 'CREATE ', 'TRUNCATE ', 'EXEC ', 'MERGE '];

// Tool definitions are registered in the MCP server (dockerfiles/mcp-server/src/index.ts).
// The Claude runtime intercepts calls to mcp__datagrok__db_* tools and forwards them
// to the browser, where handleToolCall() dispatches to the methods below.

/**
 * Context class that manages tool execution for SQL generation.
 * Supports multiple catalogs (e.g. for Databricks) with a default catalog.
 */
export class SQLGenerationContext {
  private dbInfoMap: Map<string, BuiltinDBInfoMeta> = new Map();

  constructor(
    private connectionID: string,
    private defaultCatalog: string,
    dbInfos: BuiltinDBInfoMeta[],
    private connection: DG.DataConnection
  ) {
    for (const dbInfo of dbInfos)
      this.dbInfoMap.set(dbInfo.name, dbInfo);
  }

  static async create(connectionID: string, defaultCatalog: string): Promise<SQLGenerationContext> {
    const connection = await grok.dapi.connections.find(connectionID);
    if (!connection)
      throw new Error(`Connection with ID ${connectionID} not found.`);
    const dbInfos = await BuiltinDBInfoMeta.allFromConnection(connection);
    return new SQLGenerationContext(connectionID, defaultCatalog, dbInfos, connection);
  }

  private getDbInfo(catalogName?: string): BuiltinDBInfoMeta {
    const name = catalogName ?? this.defaultCatalog;
    const dbInfo = this.dbInfoMap.get(name);
    if (!dbInfo)
      throw new Error(`Catalog '${name}' not found. Available catalogs: ${Array.from(this.dbInfoMap.keys()).join(', ')}`);
    return dbInfo;
  }

  /**
   * Parses a table reference that may include catalog, schema, and table name.
   * Supports: catalog.schema.table, schema.table, or just table
   */
  private parseTableRef(tableRef: string): {catalog: string, schema: string, table: string} {
    const parts = tableRef.trim().split('.');
    if (parts.length >= 3)
      return {catalog: parts[0], schema: parts[1], table: parts.slice(2).join('.')};
    if (parts.length === 2)
      return {catalog: this.defaultCatalog, schema: parts[0], table: parts[1]};
    return {catalog: this.defaultCatalog, schema: '', table: parts[0]};
  }

  async listCatalogs(): Promise<string[]> {
    return Array.from(this.dbInfoMap.keys());
  }

  async listAllSchemas(catalogName?: string): Promise<string[]> {
    return (await this.getDbInfo(catalogName).getSchemas()).map((s) => s.name);
  }

  /**
   * Tool: list_tables
   * Returns all tables with their descriptions
   */
  async listTables(includeAllDescriptions: boolean, schemaName: string, catalogName?: string): Promise<{description: string, tableCount: number}> {
    const dbInfo = this.getDbInfo(catalogName);

    const schema = (await dbInfo.getSchemas()).find((s) => s.name === schemaName);
    if (!schema)
      return {description: `Schema ${schemaName} not found in catalog ${catalogName ?? this.defaultCatalog}.`, tableCount: 0};

    const tables = (await schema.getTables());
    const tableDescs = tables.map((table) => {
      const parts = [`${table.friendlyName ?? table.name}`];
      const tableMeta = getDBTableMetaData(table);
      if (includeAllDescriptions) {
        if (tableMeta?.llmComment)
          parts.push(`  LLM Comment: ${tableMeta.llmComment.substring(0, 400)}...`);
        else if (tableMeta.comment)
          parts.push(`  Description: ${tableMeta.comment.substring(0, 400)}...`);
      }
      parts.push(`  Columns: ${table.columns.length}, ${tableMeta?.rowCount ? `Rows: ${tableMeta.rowCount}` : ''}`);
      return parts.join('\n');
    });
    return {description: `Tables in ${schemaName}:\n\n${tableDescs.join('\n')}`, tableCount: tableDescs.length};
  }

  /**
   * Tool: describe_tables
   * Returns detailed column information for specified tables.
   * Supports catalog.schema.table and schema.table formats.
   */
  async describeTables(tableNames: string[]): Promise<string> {
    const refs = tableNames.map((t) => this.parseTableRef(t));
    const descriptions: string[] = [];

    // Group by catalog+schema for efficient fetching
    const groupKey = (r: {catalog: string, schema: string}) => `${r.catalog}\0${r.schema}`;
    const schemaTableCache: Map<string, DG.TableInfo[]> = new Map();

    for (const key of new Set(refs.map(groupKey))) {
      const [catalog, schema] = key.split('\0');
      try {
        const dbInfo = this.getDbInfo(catalog);
        const schemaInfo = (await dbInfo.getSchemas()).find((s) => s.name === schema);
        if (schemaInfo)
          schemaTableCache.set(key, await schemaInfo.getTables());
      } catch { /* catalog not found */ }
    }

    for (const ref of refs) {
      const tables = schemaTableCache.get(groupKey(ref));
      if (!tables) {
        descriptions.push(`Schema ${ref.schema} not found in catalog ${ref.catalog}`);
        continue;
      }
      const table = tables.find((t) => t.friendlyName === ref.table || t.name === ref.table);
      if (!table) {
        descriptions.push(`Table ${ref.table} not found in ${ref.schema}`);
        continue;
      }
      descriptions.push(this.buildDetailedTableDescription(table));
    }

    return descriptions.join('\n\n---\n\n');
  }

  /**
   * Tool: list_joins
   * Returns all foreign key relationships involving the specified tables.
   * Supports catalog.schema.table and schema.table formats.
   */
  async listJoins(tableNames: string[]): Promise<string> {
    const refs = tableNames.map((t) => this.parseTableRef(t));

    const schemaQualifiedTableSet = new Set<string>();
    for (const ref of refs)
      schemaQualifiedTableSet.add(`${ref.schema}.${ref.table}`);

    // Collect unique catalog+schema pairs for efficient relation fetching
    const catalogSchemaMap = new Map<string, Set<string>>();
    for (const ref of refs) {
      if (!catalogSchemaMap.has(ref.catalog))
        catalogSchemaMap.set(ref.catalog, new Set());
      catalogSchemaMap.get(ref.catalog)!.add(ref.schema);
    }

    const allRelations: DG.DbRelationInfo[] = [];
    for (const [catalog, schemas] of catalogSchemaMap) {
      const dbInfo = this.getDbInfo(catalog);
      for (const schema of schemas) {
        const relations = await dbInfo.getRelationsForSchema(schema);
        allRelations.push(...relations);
      }
    }

    const filteredRelations = allRelations
      .filter((rel) => schemaQualifiedTableSet.has(`${rel.fromSchema}.${rel.fromTable}`) || schemaQualifiedTableSet.has(`${rel.toSchema}.${rel.toTable}`));
    const joins = filteredRelations.map((rel) => this.buildJoinDescription(rel));

    if (joins.length === 0)
      return `No foreign key relationships found for tables: ${tableNames.join(', ')}`;

    return `Foreign key relationships involving ${tableNames.join(', ')}:\n\n${joins.join('\n')}`;
  }

  static _lastFc: DG.FuncCall | null = null;
  /**
   * Tool: try_sql
   * Executes SQL and returns row count and column names
   */
  async trySql(sql: string, description: string): Promise<string> {
    console.log('trySql called with description:', description);
    try {
      // Add LIMIT if not present
      let testSql = sql.trim();
      if (testSql.endsWith(';'))
        testSql = testSql.substring(0, testSql.length - 1);
      if (!testSql.toUpperCase().includes('LIMIT'))
        testSql += ' LIMIT 10';

      // Basic safety check to prevent destructive queries
      const upperSql = sql.toUpperCase();
      if (suspiciousSQlPatterns.some((pattern) => upperSql.includes(pattern))) {
        // prompt user confirmation

        const p = new Promise<boolean>((resolve) => {
          ui.dialog('Confirm SQL Execution')
            .add(ui.divText('Our AI is trying to execute a potentially destructive SQL command. For your safety, please confirm if you want to proceed.'))
            .add(ui.markdown(`\`\`\`sql\n${sql}\n\`\`\``))
            .onOK(() => resolve(true))
            .onCancel(() => resolve(false))
            .show();
        });
        const userConfirmed = await p;
        if (!userConfirmed)
          return 'SQL Execution Not allowed by User: Destructive commands are not allowed. Please revise the query.';
      }
      const queryName = `test-query-${Date.now()}`;
      const wrappedSql = `--name: ${queryName}\n${testSql}`;

      const fc = this.connection.query(queryName, wrappedSql).prepare({});
      SQLGenerationContext._lastFc = fc;
      await fc.call(false, undefined, {processed: true, report: false});
      const df = await fc.getOutputParamValue();
      SQLGenerationContext._lastFc = null;
      if (!df)
        throw new Error('No data returned');
      const result = [
        `Query executed successfully!`,
        `Row count: ${df.rowCount}`,
        `Columns (${df.columns.length}): ${df.columns.toList().map((c: DG.Column) => `${c.name} (${c.type})`).join(', ')}`,
      ];

      // If there are rows, show a sample
      if (df.rowCount > 0) {
        const sample = df.columns.toList().map((col: DG.Column) => {
          const val = col.get(0);
          return `${col.name}=${this.formatColValue(val)}`; // save some tokens by formatting (molfiles, images, etc.)
        }).join(', ');
        result.push(`Sample row: ${sample}`);
      }

      return result.join('\n');
    } catch (error) {
      SQLGenerationContext._lastFc = null;
      const message = error instanceof Error ? error.message : String(error);
      return `SQL Error: ${message}\n\nThis query failed to execute. Please revise the SQL based on the schema information.`;
    }
  }

  /**
   * Dispatches a client-side tool call to the appropriate handler.
   * Returns the tool result as a string.
   */
  async handleToolCall(toolName: string, input: any): Promise<string> {
    switch (toolName) {
    case 'db_list_catalogs':
      return (await this.listCatalogs()).join(', ');
    case 'db_list_schemas':
      return (await this.listAllSchemas(input.catalogName)).join(', ');
    case 'db_list_tables': {
      const res = await this.listTables(false, input.schemaName, input.catalogName);
      return res.description;
    }
    case 'db_describe_tables':
      return this.describeTables(input.tables);
    case 'db_list_joins':
      return this.listJoins(input.tables);
    case 'db_try_sql':
      return this.trySql(input.sql, input.description);
    default:
      return `Unknown tool: ${toolName}`;
    }
  }

  private formatColValue(value: string | number | boolean | null | undefined | object): string {
    if (value === null || value === undefined || value === '' || value === DG.FLOAT_NULL || value === DG.INT_NULL)
      return 'NULL';
    try {
      const valueString = value?.toString();
      return valueString.length > 50 ? `'${valueString.substring(0, 47)}...'` : `'${valueString}'`;
    } catch {
      return 'NON-STRINGIFIABLE-VALUE';
    }
  }

  /**
   * Helper: Build detailed table description
   */
  private buildDetailedTableDescription(table: DG.TableInfo): string {
    const lines: string[] = [];
    const tableMeta = getDBTableMetaData(table);
    lines.push(`TABLE: ${table.friendlyName ?? table.name} (${tableMeta.rowCount ?? 'Unknown ammount of'} rows)`);
    if (tableMeta.comment)
      lines.push(`Comment: ${tableMeta.comment}`);
    if (tableMeta.llmComment)
      lines.push(`LLM Comment: ${tableMeta.llmComment}`);

    lines.push('\nColumns:');

    for (const col of table.columns) {
      const colParts: string[] = [`  ${col.name} (${col.type})`];
      const colMeta = getDBColumnMetaData(col);


      // not working currently, add back later
      // if (col.semanticType)
      //   colParts.push(`[Semantic: ${col.semanticType}]`);
      if (colMeta.isUnique)
        colParts.push('[UNIQUE]');
      if (colMeta.llmComment)
        colParts.push(`- ${colMeta.llmComment}`);
      else if (colMeta.comment)
        colParts.push(`- ${colMeta.comment}`);
      if (colMeta.min !== undefined && colMeta.max !== undefined)
        colParts.push(`Range: ${colMeta.min} to ${colMeta.max}`);
      if ((colMeta.values && colMeta.values.length > 0) || (colMeta.sampleValues && colMeta.sampleValues.length > 0)) {
        const values = (Array.isArray(colMeta.values) ? colMeta.values : undefined) ?? colMeta.sampleValues!;
        colParts.push(`Values: ${values.slice(0, 15).join(', ')}${values.length > 15 ? ', ...' : ''}`);
      }
      if (colMeta.uniqueCount)
        colParts.push(`Unique Values Count: ${colMeta.uniqueCount}`);

      lines.push(colParts.join(' '));
    }

    return lines.join('\n');
  }

  /**
   * Helper: Build join description
   */
  private buildJoinDescription(relation: DG.DbRelationInfo): string {
    const parts: string[] = [];

    parts.push(`${relation.fromSchema}.${relation.fromTable}(${(relation.fromColumns ?? ['Unknown Column']).join(', ')}) -> ${relation.toSchema}.${relation.toTable}(${(relation.toColumns ?? ['Unknown Column']).join(', ')})`);

    if (relation.cardinality)
      parts.push(`[${relation.cardinality}]`);
    if (relation.isPrimaryPath === false)
      parts.push('[LEGACY - prefer other paths if available]');
    if (relation.llmComment)
      parts.push(`- ${relation.llmComment}`);
    else if (relation.comment)
      parts.push(`- ${relation.comment}`);

    return parts.join(' ');
  }
}
