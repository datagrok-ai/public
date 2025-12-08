/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DBConnectionDescriptor} from './types';


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
