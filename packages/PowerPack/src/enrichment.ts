/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

/** Cache for full-schema calls, keyed by "nqName/schema/db". Evicts on rejection; capped at 50 entries. */
const schemaCache = new Map<string, Promise<DG.TableInfo[]>>();
const SCHEMA_CACHE_MAX = 50;

function getCachedSchema(conn: DG.DataConnection, schema: string | null, db: string | null): Promise<DG.TableInfo[]> {
  const key = `${conn.nqName}/${schema ?? ''}/${db ?? ''}`;
  if (!schemaCache.has(key)) {
    if (schemaCache.size >= SCHEMA_CACHE_MAX)
      schemaCache.delete(schemaCache.keys().next().value!);
    const p = grok.dapi.connections.getSchema(conn, schema, null, db);
    p.catch(() => schemaCache.delete(key));
    schemaCache.set(key, p);
  }
  return schemaCache.get(key)!;
}

export interface FkRef { schema: string; table: string; column: string; }

/** The foreign key chain of the column: the column it refers to, the one that one refers to, and so on. */
export async function getFkRefs(conn: DG.DataConnection, schemaName: string | null, table: string, column: string, db: string | null): Promise<FkRef[]> {
  const fkRefs: FkRef[] = [];
  try {
    const allTables: DG.TableInfo[] = await getCachedSchema(conn, schemaName, db);
    const visited = new Set<string>();
    let curTable = table;
    let curColumn = column;
    while (true) {
      const tableInfo = allTables.find((t) => t.friendlyName === curTable);
      const colInfo = tableInfo?.columns.find((c) => c.name === curColumn);
      const refSchema = colInfo?.tags[DG.Tags.ReferencesSchema];
      const refTable = colInfo?.tags[DG.Tags.ReferencesTable];
      const refColumn = colInfo?.tags[DG.Tags.ReferencesColumn];
      if (!refSchema || !refTable || !refColumn)
        break;
      const fkKey = `${refSchema}.${refTable}.${refColumn}`;
      if (visited.has(fkKey))
        break;
      visited.add(fkKey);
      fkRefs.push({schema: refSchema, table: refTable, column: refColumn});
      curTable = refTable;
      curColumn = refColumn;
    }
  }
  catch (_) { }
  return fkRefs;
}

export async function getEnrichConfigsNames(connectionNqName: string, schema: string, table: string, column: string, db: string | null = null): Promise<string[]> {
  try {
    const files = await grok.dapi.files
      .list(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}${db ? `/${db}` : ''}/${schema}/${table}/${column}`);
    return files.map((f) => f.name.substring(0, f.name.lastIndexOf('.')));
  } catch (_: any) {
    return [];
  }
}

export async function deleteEnrichment(connectionNqName: string, schema: string, table: string, column: string, name: string, db: string | null = null): Promise<boolean> {
  try {
    await grok.dapi.files
      .delete(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}${db ? `/${db}` : ''}/${schema}/${table}/${column}/${name}.json`);
    return true;
  } catch (_: any) {
    new DG.Balloon().error('Something went wrong while deleting enrichment. Please, try again.');
    return false;
  }
}

export async function readEnrichConfig(connectionNqName: string, schema: string, table: string, column: string, name: string, db: string | null = null): Promise<Enrichment | null> {
  try {
    const data = await grok.dapi.files
      .readAsText(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}${db ? `/${db}` : ''}/${schema}/${table}/${column}/${name}.json`);
    return JSON.parse(data);
  } catch (_: any) {
    return null;
  }
}

export async function saveEnrichment(e: Enrichment, connectionNqName: string): Promise<void> {
  try {
    await grok.dapi.files
      .writeAsText(`System:AppData/PowerPack/enrichments/${nqNameToPath(connectionNqName)}${e.keyDb ? `/${e.keyDb}` : ''}/${e.keySchema}/${e.keyTable}/${e.keyColumn}/${e.name}.json`,
        JSON.stringify(e));
  } catch (_: any) {
    new DG.Balloon().error('Could not save enrichment config. Please try again.');
  }
}

function nqNameToPath(connectionNqName: string) {
  return connectionNqName.replace(':', '_');
}

export function convertTableQueryToEnrichment(query: DG.TableQuery, name: string, keySchema: string, keyTable: string, keyColumn: string, keyDb: string | null = null): Enrichment {
  return {name: name, keySchema: keySchema, keyTable: keyTable, keyColumn: keyColumn, fields: query.fields, joins: query.joins, keyDb: keyDb};
}

function convertEnrichmentToQuery(e: Enrichment, conn: DG.DataConnection): DG.TableQuery {
  const query = conn.tableQuery();
  query.table = `${e.keySchema}.${e.keyTable}`;
  query.fields = e.fields;
  query.joins = e.joins;
  return query;
}

export async function runEnrichmentFromConfig(conn: DG.DataConnection, schema: string, table: string, column: string, name: string, df: DG.DataFrame, db: string | null = null, localColumn: string | null = null): Promise<void> {
  const config = await readEnrichConfig(conn.nqName, schema, table, column, name, db);
  if (!config) {
    new DG.Balloon().error('Could not find enrichment. Please try again or create a new one.');
    return;
  }
  const query = convertEnrichmentToQuery(config, conn);
  await executeEnrichQuery(query, df, column, localColumn ?? undefined);
}

export async function executeEnrichQuery(query: DG.TableQuery, df: DG.DataFrame, keyCol: string, localKeyCol?: string): Promise<void> {
  try {
    query.limit = undefined;
    const dfKeyCol = df.columns.firstWhere(col => col.tags.get(DG.Tags.DbColumn) === (localKeyCol ?? keyCol));
    if (!dfKeyCol) {
      grok.shell.error('Could not resolve the column.');
      return;
    }
    const keyColValues = dfKeyCol.toList();
    const dataType = dfKeyCol.type;
    let res: DG.DataFrame | undefined;
    for (let i = 0; i < keyColValues.length; i += 200) {
      const inValues = keyColValues.slice(i, i + 200);
      query.where = [{field: `${query.table}.${keyCol}`, pattern: `in (${inValues.join(',')})`, dataType}];
      const queryCall = query.prepare();
      const run = await queryCall.call(false, undefined,
        {processed: true, report: false});
      const queryResult = run.getOutputParamValue();
      if (!queryResult)
        break;
      if (!res)
        res = queryResult;
      else
        res.append(queryResult, true);
    }
    if (!res) {
      grok.shell.error('Something went wrong and query result returned null.');
      return;
    }

    await grok.data.detectSemanticTypes(res);
    const joinTable = DG.Func.byName('JoinTables');
    const previousName = df.name;
    joinTable.applySync({
      'table1': df,
      'table2': res,
      'keys1': [dfKeyCol.name],
      'keys2': [keyCol],
      'values1': df.columns.names(),
      'values2': res.columns.names().filter((n) => n !== keyCol),
      'joinType': 'left',
      'inPlace': true
    });
    df.name = previousName;
  } catch (e: any) {
    grok.shell.error(`Failed to enrich:\n ${e}`);
  }
}

export interface Enrichment {
  name: string;
  keyDb: string | null;
  keySchema: string;
  keyTable: string,
  keyColumn: string;
  fields: string[];
  joins: DG.TableJoin[];
}
