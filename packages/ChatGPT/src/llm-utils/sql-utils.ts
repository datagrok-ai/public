/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenAIHelpClient} from './openAI-client';
import {getSchemaDescriptor} from '@datagrok-libraries/db-explorer/src/utils';
import {modelName} from '../package';

type DBIDWithSchema = `${string}.${string}`;
type SchemaInfo = {
    tableInfos: string[];
    referenceInfos: string[];
}
class DBSchemaInfo {
  private static _cache: Map<DBIDWithSchema, SchemaInfo> = new Map();

  static async getSchemaInfo(connectionID: string, schemaName: string): Promise<SchemaInfo> {
    const cacheKey: DBIDWithSchema = `${connectionID}.${schemaName}`;
    if (this._cache.has(cacheKey))
      return this._cache.get(cacheKey)!;
    const schemaDescriptor = await getSchemaDescriptor(connectionID, schemaName);
    const tableInfos: string[] = [];
    const referenceInfos: string[] = [];
    for (const [tableName, tableInfo] of Object.entries(schemaDescriptor.tables))
      tableInfos.push(`${tableName}: (${Object.entries(tableInfo.columns).map(([colName, colType]) => `${colName} (${colType})`).join(', ')})`);
    for (const ref of schemaDescriptor.references)
      referenceInfos.push(`${ref.from} -> ${ref.to}`);
    const schemaInfo = {tableInfos, referenceInfos};
    this._cache.set(cacheKey, schemaInfo);
    return schemaInfo;
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

CRITICAL JOIN RULES:
- ALWAYS follow the foreign key relationships provided in the References section
- NEVER create direct joins between tables that don't have a direct foreign key relationship
- If table A needs data from table C, but only A->B and B->C relationships exist, you MUST join through B (the intermediate table)
- Before writing any JOIN, verify the relationship exists in the References section
- DO NOT assume relationships exist between tables based on column name similarity (e.g., user_id in different tables)
- ONLY use relationships explicitly listed in the References section

QUERY PLANNING PROCESS (follow these steps mentally before writing SQL):
1. Identify all tables mentioned in the user question
2. Trace the foreign key path between these tables using ONLY the References section
3. List all intermediate tables needed for joins
4. Verify each join relationship exists in the References section
5. Write the query following this foreign key path

EXAMPLE - Multi-hop join:
Question: "Get events for a specific user"
CORRECT Path: events -> users_sessions -> users (follow the FK chain through intermediate table)
WRONG Path: events -> users (no direct relationship exists - this will fail)

OUTPUT REQUIREMENTS:
- Only return one SQL statement
- No explanation, no semicolon at the end
- No markdown or any extra text
- NO \`\`\`sql prefix
- ONLY THE SQL AND NOTHING ELSE!!!
- DO NOT USE 'to' alias (or any prohibited aliases) for any table
- Use aliases as much as possible to shorten the SQL (except for prohibited ones)
- When referencing tables, always use the schema name as prefix
- Schema name is: ${schemaName}
- Make sure the SQL query is nicely formatted for better readability

CRITICAL VALIDATION:
- DO NOT MAKE UP TABLES OR COLUMNS THAT DO NOT EXIST IN THE SCHEMA PROVIDED BELOW
- DO NOT CREATE JOINS THAT ARE NOT LISTED IN THE REFERENCES SECTION
- For each JOIN you write, verify the exact relationship exists in References
- If the query is not possible with the given schema and available foreign keys, RETURN: SELECT * FROM table WHERE 1=0
    
    Here is the database schema:
    Tables:
    ${tableContents.join('\n')}


    References (Foreign Key Relationships - USE ONLY THESE FOR JOINS):
    ${referenceInfos.join('\n')}


    User prompt:
    ${prompt}

    DO NOT MAKE UP TABLES OR COLUMNS THAT DO NOT EXIST IN THE SCHEMA PROVIDED ABOVE! IF THE QUERY IS NOT POSSIBLE WITH THE GIVEN SCHEMA, RETURN A QUERY THAT RETURNS NO RESULTS, E.G., "SELECT * FROM table WHERE 1=0"
    DO NOT USE 'to' alias (or any prohibited alieas) for any table.
    MAKE SURE THE SQL QUERY IS NICELY FORMATTED FOR BETTER READABILITY.
    `;
  const sqlQuery = await client.generalPromptCached(modelName, '', systemPrompt);
  if (sqlQuery.startsWith('```sql') && sqlQuery.endsWith('```'))
    return sqlQuery.substring(6, sqlQuery.length - 3).trim();

  return sqlQuery;
}
