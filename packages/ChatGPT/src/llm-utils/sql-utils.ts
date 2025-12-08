/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {OpenAIHelpClient} from './openAI-client';
import {modelName} from '../package';
import {DBSchemaInfo, DBConnectionMeta, DBTableMeta, DBRelationMeta} from './db-index-tools';
import {JsonSchema} from '../prompt-engine/interfaces';
import {chemblIndex} from './indexes/chembl-index';
import {biologicsIndex} from './indexes/biologics-index';

// Schema for table selection
const TableSelectionSchema: JsonSchema = {
  type: 'object',
  properties: {
    tables: {
      type: 'array',
      items: {type: 'string'},
      description: 'List of table names that are needed to answer the user query',
    },
    reasoning: {
      type: 'string',
      description: 'Brief explanation of why these tables were selected',
    },
  },
  required: ['tables', 'reasoning'],
  additionalProperties: false,
};

interface TableSelection {
  tables: string[];
  reasoning: string;
}

/**
 * Generates an SQL query using a two-step approach:
 * 1. First, identify relevant tables needed for the query
 * 2. Then, generate SQL with focused context on those tables
 *
 * @param prompt - User's natural language query
 * @param connectionID - Database connection ID
 * @param schemaName - Schema name
 * @param dbMeta - Optional DBConnectionMeta object with enriched metadata
 */
export async function generateAISqlQuery(
  prompt: string,
  connectionID: string,
  schemaName: string,
  dbMeta?: DBConnectionMeta
): Promise<string> {
  // Step 1: Identify relevant tables
  const connection = await grok.dapi.connections.find(connectionID);
  if (connection.name.toLowerCase() === 'chembl')
    dbMeta = chemblIndex;
  else if (connection.name.toLowerCase() === 'biologics')
    dbMeta = biologicsIndex;
  console.log('Step 1: Identifying relevant tables for query...');
  const tableSelection = await identifyRelevantTables(connectionID, prompt, schemaName, dbMeta);
  console.log(`Selected tables: ${tableSelection.tables.join(', ')}`);
  console.log(`Reasoning: ${tableSelection.reasoning}`);

  // Step 2: Get focused schema context
  let focusedContext: FocusedSchemaContext;
  if (dbMeta)
    focusedContext = extractFocusedContext(dbMeta, schemaName, tableSelection.tables);
  else {
    // Fallback to schema descriptor approach
    const schema = await DBSchemaInfo.getSchemaInfo(connectionID, schemaName);
    focusedContext = buildFocusedContextFromDescriptor(schema, tableSelection.tables);
  }

  if (focusedContext.tableInfos.length === 0)
    throw new Error('No schema information available for selected tables.');

  // Step 3: Generate SQL query with focused context
  console.log('Step 2: Generating SQL query with focused context...');
  const sqlQuery = await generateSqlWithContext(prompt, schemaName, focusedContext);

  if (sqlQuery.startsWith('```sql') && sqlQuery.endsWith('```'))
    return sqlQuery.substring(6, sqlQuery.length - 3).trim();

  return sqlQuery;
}

interface FocusedSchemaContext {
  tableInfos: string[];
  relationInfos: string[];
  columnDetails: string[];
}

/**
 * First LLM call: Identify which tables are needed for the query
 */
async function identifyRelevantTables(
  connectionID: string,
  prompt: string,
  schemaName: string,
  dbMeta?: DBConnectionMeta
): Promise<TableSelection> {
  const client = OpenAIHelpClient.getInstance();

  // Build table overview for selection
  let tableOverview: string;
  if (dbMeta) {
    const schema = dbMeta.schemas.find((s) => s.name === schemaName);
    if (!schema)
      throw new Error(`Schema ${schemaName} not found in DBConnectionMeta`);

    tableOverview = schema.tables.map((table) => {
      const parts = [`${table.name}`];
      if (table.LLMComment) parts.push(`  ${table.LLMComment}`);
      else if (table.comment) parts.push(`- ${table.comment}`);
      parts.push(`  Columns: ${table.columns.map((c) => c.name).join(', ')}`);
      return parts.join('\n');
    }).join('\n\n');
  } else {
    const schema = await DBSchemaInfo.getSchemaInfo(connectionID, schemaName);
    tableOverview = schema.tableInfos.join('\n');
  }

  const systemPrompt = `You are a database expert. Given a list of available tables in a database schema, identify which tables are needed to answer the user's query.

Consider:
- Direct tables mentioned or implied in the query
- Related tables needed for joins (check foreign key relationships)
- Intermediate tables required for multi-hop joins
- Tables containing relevant attributes even if not explicitly mentioned

Be inclusive but focused - include tables that might be needed, it is much better to have something extra rather than lack it, but avoid including the entire schema unless necessary.`;

  const userPrompt = `Schema: ${schemaName}

Available tables:
${tableOverview}

User query: "${prompt}"

Identify which tables are needed to answer this query. Return ONLY the table names (without schema prefix).`;

  const response = await client.generalPromptCached(
    modelName,
    systemPrompt,
    userPrompt,
    TableSelectionSchema
  );

  return JSON.parse(response) as TableSelection;
}

/**
 * Extract focused context from DBConnectionMeta for selected tables
 */
function extractFocusedContext(
  dbMeta: DBConnectionMeta,
  schemaName: string,
  selectedTables: string[]
): FocusedSchemaContext {
  const schema = dbMeta.schemas.find((s) => s.name === schemaName);
  if (!schema)
    throw new Error(`Schema ${schemaName} not found in DBConnectionMeta`);

  const tableInfos: string[] = [];
  const columnDetails: string[] = [];
  const selectedTableSet = new Set(selectedTables);

  // Extract table and column information for selected tables
  for (const table of schema.tables) {
    if (!selectedTableSet.has(table.name)) continue;

    // Build table info with comments
    const tableInfo = buildTableInfo(table);
    tableInfos.push(tableInfo);

    // Build detailed column info
    const colDetail = buildColumnDetails(table);
    columnDetails.push(colDetail);
  }

  // Extract relevant relations (only between selected tables or involving selected tables)
  const relationInfos: string[] = [];
  for (const relation of dbMeta.relations) {
    if (relation.fromSchema !== schemaName && relation.toSchema !== schemaName)
      continue;

    // Include relation if either table is selected, or if it's a bridge between selected tables
    const fromSelected = selectedTableSet.has(relation.fromTable);
    const toSelected = selectedTableSet.has(relation.toTable);

    if (fromSelected || toSelected) {
      // Add the intermediate table if this relation connects to selected tables
      if (fromSelected && !toSelected)
        selectedTableSet.add(relation.toTable);
      else if (!fromSelected && toSelected)
        selectedTableSet.add(relation.fromTable);

      const relInfo = buildRelationInfo(relation);
      relationInfos.push(relInfo);
    }
  }

  // Second pass: add any missing tables that were identified as intermediates
  for (const table of schema.tables) {
    if (selectedTableSet.has(table.name) && !selectedTables.includes(table.name)) {
      const tableInfo = buildTableInfo(table);
      tableInfos.push(tableInfo);

      const colDetail = buildColumnDetails(table);
      columnDetails.push(colDetail);
    }
  }

  return {tableInfos, relationInfos, columnDetails};
}

/**
 * Build table info string with comments
 */
function buildTableInfo(table: DBTableMeta): string {
  const columns = table.columns.map((c) => {
    const parts = [c.name, c.type];
    if (c.semanticType) parts.push(`sem:${c.semanticType}`);
    if (c.isUnique) parts.push('UNIQUE');
    return parts.join(' ');
  }).join(', ');

  const parts = [`${table.name}: (${columns})`];
  if (table.comment) parts.push(`  // ${table.comment}`);
  if (table.LLMComment) parts.push(`  // AI Context: ${table.LLMComment}`);

  return parts.join('\n');
}

/**
 * Build detailed column information
 */
function buildColumnDetails(table: DBTableMeta): string {
  const details = table.columns.map((col) => {
    const parts = [`  ${col.name} (${col.type})`];
    if (col.LLMComment) parts.push(`- ${col.LLMComment}`);
    else if (col.comment) parts.push(`- ${col.comment}`);
    if (col.semanticType) parts.push(`- Semantic type: ${col.semanticType}`);
    if (col.isUnique) parts.push('- UNIQUE');
    if (col.min !== undefined && col.max !== undefined)
      parts.push(`- Range: ${col.min} to ${col.max}`);
    if (col.categoryValues && col.categoryValues.length > 0) {
      const values = col.categoryValues.slice(0, 10);
      parts.push(`- Values: ${values.join(', ')}${col.categoryValues.length > 10 ? '...' : ''}`);
    }
    return parts.join(' ');
  }).join('\n');

  return `Table: ${table.name} (${table.rowCount} rows)\n${details}`;
}

/**
 * Build relation info string with cardinality and comments
 */
function buildRelationInfo(relation: DBRelationMeta): string {
  const relStr = `${relation.fromTable}(${relation.fromColumns.join(', ')}) -> ${relation.toTable}(${relation.toColumns.join(', ')})`;
  const parts = [relStr];

  if (relation.cardinality) parts.push(`[${relation.cardinality}]`);
  if (relation.IsPrimaryPath === false) parts.push('[LEGACY]');
  if (relation.LLMComment) parts.push(`// ${relation.LLMComment}`);
  else if (relation.comment) parts.push(`// ${relation.comment}`);

  return parts.join(' ');
}

/**
 * Fallback: Build focused context from schema descriptor (when DBConnectionMeta is not available)
 */
function buildFocusedContextFromDescriptor(
  schema: {tableInfos: string[]; referenceInfos: string[]},
  selectedTables: string[]
): FocusedSchemaContext {
  const selectedTableSet = new Set(selectedTables);

  // Filter table infos to only selected tables
  const tableInfos = schema.tableInfos.filter((info) => {
    const tableName = info.split(':')[0].trim();
    return selectedTableSet.has(tableName);
  });

  // Filter references to only those involving selected tables
  const relationInfos = schema.referenceInfos.filter((ref) => {
    const [from, to] = ref.split('->').map((s) => s.trim());
    const fromTable = from.split('.')[0];
    const toTable = to.split('.')[0];
    return selectedTableSet.has(fromTable) || selectedTableSet.has(toTable);
  });

  return {
    tableInfos,
    relationInfos,
    columnDetails: [], // Not available in descriptor mode
  };
}

/**
 * Second LLM call: Generate SQL query with focused context
 */
async function generateSqlWithContext(
  prompt: string,
  schemaName: string,
  context: FocusedSchemaContext
): Promise<string> {
  const client = OpenAIHelpClient.getInstance();

  const systemPrompt = `You are an expert SQL generator. Given the database schema below, with tables, their columns, types and reference info (foreign keys), write a SQL query for the user question.

CRITICAL JOIN RULES:
- ALWAYS follow the foreign key relationships provided in the References section
- NEVER create direct joins between tables that don't have a direct foreign key relationship
- If table A needs data from table C, but only A->B and B->C relationships exist, you MUST join through B (the intermediate table)
- Before writing any JOIN, verify the relationship exists in the References section
- DO NOT assume relationships exist between tables based on column name similarity (e.g., user_id in different tables)
- ONLY use relationships explicitly listed in the References section
- Pay attention to [LEGACY] markers - prefer non-legacy relationships when multiple paths exist
- Pay attention to relationship cardinality (one-to-many, many-to-one, etc.) when writing JOINs

QUERY PLANNING PROCESS (follow these steps mentally before writing SQL):
1. Identify all needed tables mentioned in the user question
2. Trace the foreign key path between these tables using ONLY the References section
3. List all intermediate tables needed for joins
4. Verify each join relationship exists in the References section
5. Write the query following this foreign key path

EXAMPLE - Multi-hop join:
Question: "Get events for a specific user"
CORRECT Path: events -> users_sessions -> users (follow the FK chain through intermediate table)
WRONG Path: events -> users (no direct relationship exists - this will fail)

COLUMN SELECTION:
- Use column comments and AI context to understand what columns mean
- Pay attention to semantic types (e.g., Molecule, Macromolecule) when filtering or selecting
- Use category values when filtering on categorical columns
- Consider value ranges when writing WHERE clauses

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
- If the query is not possible with the given schema and available foreign keys, RETURN: SELECT * FROM ${schemaName}.${context.tableInfos[0]?.split(':')[0] || 'table'} WHERE 1=0

Here is the focused database schema for your query:

Tables:
${context.tableInfos.join('\n')}

${context.columnDetails.length > 0 ? `\nDetailed Column Information:\n${context.columnDetails.join('\n\n')}` : ''}

References (Foreign Key Relationships - USE ONLY THESE FOR JOINS):
${context.relationInfos.join('\n')}

User prompt:
${prompt}

DO NOT MAKE UP TABLES OR COLUMNS THAT DO NOT EXIST IN THE SCHEMA PROVIDED ABOVE! 
DO NOT USE 'to' alias (or any prohibited aliases) for any table.
MAKE SURE THE SQL QUERY IS NICELY FORMATTED FOR BETTER READABILITY.`;

  return await client.generalPromptCached(modelName, '', systemPrompt);
}
