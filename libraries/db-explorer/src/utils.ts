/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBConnectionDescriptor, ReferencedByObject, ReferenceObject, SchemaInfo} from './types';

/** Providers where one query reads any catalog as catalog.schema.table. Postgres reports catalogs too, but each
 * is a separate database that a query cannot name */
const CATALOG_SOURCES: string[] = [DG.DataSourceType.Databricks, DG.DataSourceType.Snowflake, DG.DataSourceType.MsSql];

export function hasCatalogs(connection: DG.DataConnection): boolean {
  return CATALOG_SOURCES.includes(connection.dataSource) && !!DG.ConnectionDataSource.byType(connection.dataSource)?.supportCatalogs;
}

export function ownCatalog(connection: DG.DataConnection): string | null {
  return connection.parameters['db'] || null;
}

/** The catalog a configuration has to name: none for the connection's own, where schema.table already resolves
 * and whose name the connection may spell in another case than the database does */
export function explicitCatalog(connection: DG.DataConnection, catalogName?: string | null): string | undefined {
  return !catalogName || catalogName.toLowerCase() === ownCatalog(connection)?.toLowerCase() ? undefined : catalogName;
}

export async function listCatalogs(connection: DG.DataConnection): Promise<string[]> {
  return (await grok.data.db.getInfo(connection)).map((info) => info.name)
    .sort((a, b) => a.localeCompare(b, undefined, {sensitivity: 'base'}));
}

/** Foreign keys of a schema, as what each column refers to and what refers to it */
export async function getSchemaReferences(connection: DG.DataConnection, schemaName: string, catalogName?: string): Promise<SchemaInfo> {
  const catalog = explicitCatalog(connection, catalogName);
  const tables = await grok.dapi.connections.getSchema(connection, schemaName, null, catalog ?? null);
  if (!tables)
    throw new Error(`Schema ${schemaName} not found`);
  const references: ReferenceObject = {[schemaName]: {}};
  const referencedBy: ReferencedByObject = {[schemaName]: {}};
  const schemaRefs = references[schemaName];
  const schemaRefBy = referencedBy[schemaName];
  // the platform reads foreign keys from the connection's own catalog, whichever catalog the tables come from
  const keysKnown = !catalog;
  for (const table of tables) {
    const tableName = table.friendlyName ?? table.name;
    schemaRefs[tableName] = {};
    if (!keysKnown)
      continue;
    for (const column of table.columns) {
      const ref = column.referenceInfo;
      if (!ref || !ref.table || !ref.column)
        continue;
      schemaRefs[tableName][column.name] = {refTable: ref.table, refColumn: ref.column, refSchema: schemaName};
      schemaRefBy[ref.table] ??= {};
      schemaRefBy[ref.table][ref.column] ??= [];
      schemaRefBy[ref.table][ref.column].push({refTable: tableName, refColumn: column.name, refSchema: schemaName});
    }
  }
  return {references, referencedBy};
}

export async function getSchemaDescriptor(connectionID: string, schemaName: string, description?: string): Promise<DBConnectionDescriptor['schemas'][string]> {
  const out: DBConnectionDescriptor['schemas'][string] = {
    tables: {},
    description: description ?? '',
    references: [],
  };

  try {
    const connection = await grok.dapi.connections.find(connectionID);
    if (!connection)
      throw new Error('Connection not found');
    const schema = await grok.dapi.connections.getSchema(connection, schemaName);
    if (!schema)
      throw new Error(`Schema not found in connection ${connection.name}`);

    schema.forEach((table) => {
      const tableName = table.friendlyName ?? table.name;
      out.tables[tableName] = {
        columns: {}
      };
      table.columns.forEach((column) => {
        out.tables[tableName].columns[column.name] = column.type;
      });
      table.columns.forEach((column) => {
        const ref = column.referenceInfo;
        if (ref && ref.table && ref.column)
          out.references.push({from: `${tableName}.${column.name}`, to: `${ref.table}.${ref.column}`});
      });
    });
  } catch (e) {
    console.error(`Failed to get schema descriptor for ${schemaName}`, e);
  }
  return out;
}
