import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {SQLGenerationContext} from '../db/sql-tools';

// viewAIToolsProvider implementations registered in package.ts. These complement the
// Dart-native tools views declare via getAITools() (e.g. DataQueryView's get_query_info /
// set_query_and_run) with tools that are easier to build in JS.

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

const stringArraySchema = (desc: string) => ({
  type: 'object',
  properties: {tables: {type: 'array', items: {type: 'string'}, description: desc}},
  required: ['tables'],
});

/** Schema-exploration and SQL-testing tools for the database query editor.
 * The connection is discovered through the view's Dart-native get_query_info tool. */
export async function buildQueryViewAITools(view: DG.ViewBase): Promise<DG.AIViewTool[]> {
  const native: DG.AIViewTool[] = (view as any).getAITools?.() ?? [];
  const info = native.find((t) => t.name === 'get_query_info');
  if (!info)
    return [];
  const {connectionId} = (await info.run({})) ?? {};
  if (!connectionId)
    return [];
  const ctx = await getSqlContext(connectionId);
  return [
    {
      name: 'list_db_catalogs',
      description: 'List the catalogs available on this connection.',
      run: async () => (await ctx.listCatalogs()).join(', '),
    },
    {
      name: 'list_db_schemas',
      description: 'List schemas of a catalog (defaults to the connection default catalog).',
      inputSchema: {type: 'object', properties: {catalogName: {type: 'string'}}},
      run: async (args: any) => (await ctx.listAllSchemas(args?.catalogName)).join(', '),
    },
    {
      name: 'list_db_tables',
      description: 'List tables of a schema with row counts.',
      inputSchema: {
        type: 'object',
        properties: {schemaName: {type: 'string'}, catalogName: {type: 'string'}},
        required: ['schemaName'],
      },
      run: async (args: any) => (await ctx.listTables(false, args.schemaName, args?.catalogName)).description,
    },
    {
      name: 'get_db_table_details',
      description: 'Detailed column info (types, comments, ranges, sample values) for the given tables. ' +
        'Table refs: catalog.schema.table, schema.table, or table.',
      inputSchema: stringArraySchema('Table references to describe.'),
      run: async (args: any) => ctx.describeTables(args.tables ?? []),
    },
    {
      name: 'list_db_joins',
      description: 'Foreign-key relationships involving the given tables — use to build correct JOINs.',
      inputSchema: stringArraySchema('Table references to find relationships for.'),
      run: async (args: any) => ctx.listJoins(args.tables ?? []),
    },
    {
      name: 'get_sql_test_result',
      description: 'Test-execute a SELECT (auto-LIMITed) and report row count, columns, and a sample row. ' +
        'Use to validate SQL before setting it on the editor with set_query_and_run.',
      inputSchema: {
        type: 'object',
        properties: {
          sql: {type: 'string', description: 'The SQL to test.'},
          description: {type: 'string', description: 'One line describing what the query does.'},
        },
        required: ['sql', 'description'],
      },
      run: async (args: any) => ctx.trySql(args.sql, args.description ?? ''),
    },
  ];
}

/** Read/write tools for the script editor. The language is part of the script header
 * (`#language:` / `//language:` comment) inside the code itself. */
export function buildScriptViewAITools(view: DG.ViewBase): DG.AIViewTool[] {
  const sv = new DG.ScriptView(DG.toDart(view));
  return [
    {
      name: 'get_script_code',
      description: 'Read the current content of this script editor, including the metadata header ' +
        '(name, language, inputs, outputs).',
      run: () => ({code: sv.code ?? ''}),
    },
    {
      name: 'set_script_code',
      description: 'Replace the content of this script editor. Include the full metadata header ' +
        '(#name, #language, #input/#output annotations) appropriate for the chosen language.',
      inputSchema: {
        type: 'object',
        properties: {code: {type: 'string', description: 'The full script text to set.'}},
        required: ['code'],
      },
      run: (args: any) => {
        if (!args?.code)
          return {success: false, error: 'code is empty'};
        sv.code = args.code;
        return {success: true};
      },
    },
  ];
}
