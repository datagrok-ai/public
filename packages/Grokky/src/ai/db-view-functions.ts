import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {SQLGenerationContext} from '../db/sql-tools';

// Implementations of the `meta.viewType: DataQueryView` package functions (package.ts).
// They complement the Dart-native functions the query view returns from getFunctions()
// (get_query_info / set_query_and_run) with SQL schema exploration and test execution.
// Each takes the generic `view` (the current DataQueryView) as its first argument.

const sqlContexts = new Map<string, Promise<SQLGenerationContext>>();

async function getSqlContext(connectionId: string): Promise<SQLGenerationContext> {
  if (!sqlContexts.has(connectionId)) {
    sqlContexts.set(connectionId, (async () => {
      const connection = await grok.dapi.connections.find(connectionId);
      const catalog = connection?.parameters?.['catalog'] ?? connection?.parameters?.['db'] ?? '';
      return SQLGenerationContext.create(connectionId, catalog);
    })());
  }
  return sqlContexts.get(connectionId)!;
}

/** The connection is discovered through the view's own Dart-native get_query_info function. */
async function queryViewContext(view: any): Promise<SQLGenerationContext> {
  const info = (view?.getFunctions?.() as DG.Func[] | undefined)?.find((f) => f.name === 'get_query_info');
  const connectionId = info ? ((await info.apply({})) ?? {})['connectionId'] : null;
  if (!connectionId)
    throw new Error('This view has no database connection');
  return getSqlContext(connectionId);
}

const splitTables = (tables: string) => (tables ?? '').split(',').map((t) => t.trim()).filter((t) => t.length > 0);

export async function listDbCatalogs(view: any): Promise<string> {
  return (await (await queryViewContext(view)).listCatalogs()).join(', ');
}

export async function listDbSchemas(view: any, catalogName?: string): Promise<string> {
  return (await (await queryViewContext(view)).listAllSchemas(catalogName)).join(', ');
}

export async function listDbTables(view: any, schemaName: string, catalogName?: string): Promise<string> {
  return (await (await queryViewContext(view)).listTables(false, schemaName, catalogName)).description;
}

export async function getDbTableDetails(view: any, tables: string): Promise<string> {
  return (await queryViewContext(view)).describeTables(splitTables(tables));
}

export async function listDbJoins(view: any, tables: string): Promise<string> {
  return (await queryViewContext(view)).listJoins(splitTables(tables));
}

export async function getSqlTestResult(view: any, sql: string, description: string): Promise<string> {
  return (await queryViewContext(view)).trySql(sql, description ?? '');
}
