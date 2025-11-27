/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenAIHelpClient} from './openAI-client';
import {modelName} from '../package';

type DBID = string;
type SchemaInfo = {
    tableInfos: string[];
    referenceInfos: string[];
}
class DBSchemaInfo {
  private static _cache: Map<DBID, SchemaInfo> = new Map();

  static async getSchemaInfo(connectionID: string, schemaName: string): Promise<SchemaInfo> {
    if (this._cache.has(connectionID))
      return this._cache.get(connectionID)!;
    try {
      const connection = await grok.dapi.connections.find(connectionID);
      if (!connection)
        throw new Error('Connection not found');
      const schema = await grok.dapi.connections.getSchema(connection, schemaName);
      if (!schema)
        throw new Error(`Schema not found in connection ${connection.name}`);
      const tableInfos: string[] = [];
      const referenceInfos: string[] = [];
      schema.forEach((table) => {
        const tableName = table.friendlyName ?? table.name;
        // this.references[tableName] = {};
        // const t = this.references[tableName];
        tableInfos.push(`${tableName}: (${table.columns.map((c) => `${c.name} (${c.type})`).join(', ')})`);
        table.columns.forEach((column) => {
          const ref = column.referenceInfo;
          if (ref && ref.table && ref.column)
            referenceInfos.push(`${tableName}.${column.name} -> ${ref.table}.${ref.column}`);
        });
      });
      const schemaInfo = {tableInfos, referenceInfos};
      this._cache.set(connectionID, schemaInfo);
      return schemaInfo;
    } catch (e) {
      console.error('Failed to get schema info', e);
      this._cache.set(connectionID, {tableInfos: [], referenceInfos: []});
      return {tableInfos: [], referenceInfos: []};
    }
  }
}

export async function generateAISqlQuery(prompt: string, connectionID: string, schemaName: string): Promise<string> {
  const client = OpenAIHelpClient.getInstance();
  // get the schema info
  const schema = await DBSchemaInfo.getSchemaInfo(connectionID, schemaName);
  const tableContents: string[] = schema.tableInfos;
  const referenceInfos: string[] = schema.referenceInfos;
  if (tableContents.length === 0)
    throw new Error('No schema information available to generate SQL query.');


  const systemPrompt = `You are an expert SQL generator. Given the database schema below, with tables, their columns, types and reference info (foreign keys), write a SQL query for the user question.
    Only return one SQL statement, with no explanation, with no semicolon in the end, no markdown or any extra text, NO \`\`\`sql prefix. ONLY THE SQL AND NOTHING ELSE!!!.
    DO NOT USE 'to' alias (or any prohibited alieas) for any table, otherwise, try to use aliases as much as possible to shorten the sql. when referencing tables, always use the schema name as prefix. Schema name is: ${schemaName}.
    Here is the database schema:
    Tables:
    ${tableContents.join('\n')}
    References:
    ${referenceInfos.join('\n')}


    User prompt:
    ${prompt}

    DO NOT MAKE UP TABLES OR COLUMNS THAT DO NOT EXIST IN THE SCHEMA PROVIDED ABOVE! IF THE QUERY IS NOT POSSIBLE WITH THE GIVEN SCHEMA, RETURN A QUERY THAT RETURNS NO RESULTS, E.G., "SELECT * FROM table WHERE 1=0"
    DO NOT USE 'to' alias (or any prohibited alieas) for any table.
    MAKE SURE THE SQL QUERY IS NICELY FORMATTED FOR BETTER READABILITY.
    `;
  const sqlQuery = await client.generalPromptCached('gpt-4.1-nano', '', systemPrompt);
  if (sqlQuery.startsWith('```sql') && sqlQuery.endsWith('```'))
    return sqlQuery.substring(6, sqlQuery.length - 3).trim();

  return sqlQuery;
}
