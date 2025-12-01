/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import OpenAI from 'openai';
import {LLMCredsManager} from './creds';
import {modelName} from '../package';
import {DBSchemaInfo, DBConnectionMeta, DBTableMeta, DBRelationMeta, DBSchemaMeta} from './db-index-tools';
import {chemblIndex} from './indexes/chembl-index';
import {biologicsIndex} from './indexes/biologics-index';
import {AbortPointer} from '../utils';
import * as rxjs from 'rxjs';
import {getTopKSimilarQueries, getVectorEmbedding} from './embeddings';

export const AI_SQL_QUERY_ABORT_EVENT = 'd4-ai-sql-query-abort';
/**
 * Generates SQL query using function calling approach where LLM can explore schema interactively
 * @param prompt - User's natural language query
 * @param connectionID - Database connection ID
 * @param schemaName - Schema name
 * @param dbMeta - Optional DBConnectionMeta object with enriched metadata
 */
export async function generateAISqlQueryWithTools(
  prompt: string,
  connectionID: string,
  schemaName: string,
  abortSignal: AbortPointer,
  dbMeta?: DBConnectionMeta
): Promise<string> {
  // Load pre-indexed metadata if available
  const connection = await grok.dapi.connections.find(connectionID);
  let dbQueryEmbeddings: {query: string, embedding: number[]}[] = [];
  if (connection.name.toLowerCase() === 'chembl') {
    dbMeta = chemblIndex;
    // kind of lazy load to save memory
    dbQueryEmbeddings = (await import('./indexes/chembl-query-embeddings')).chemblQueryEmbeddings;
  } else if (connection.name.toLowerCase() === 'biologics') {
    dbMeta = biologicsIndex;
    dbQueryEmbeddings = (await import('./indexes/biologics-query-embeddings')).biologicsQueryEmbeddings;
  }
  if (connection.name.toLowerCase() === 'datagrok')
    dbQueryEmbeddings = (await import('./indexes/datagrok-query-embeddings')).datagrokQueryEmbeddings;

  // Create tool execution context
  const context = new SQLGenerationContext(connectionID, schemaName, dbMeta, connection);

  // Initialize OpenAI client
  const openai = new OpenAI({
    apiKey: LLMCredsManager.getApiKey(),
    dangerouslyAllowBrowser: true,
  });

  // Get initial table list
  const initialTableList = await context.listTables();

  let similarQueriesInfo = '';
  if ((dbQueryEmbeddings?.length ?? 0) > 0) {
    const userQueryEmbedding = await getVectorEmbedding(openai, prompt);
    const topSimilarQueries = getTopKSimilarQueries(userQueryEmbedding, dbQueryEmbeddings, 3);
    console.log('Top similar queries:', topSimilarQueries);
    similarQueriesInfo = `Here are some similar previously executed queries that might help you. 
    Note that these queries might have annotations/comments like name descriptions and so on that can help you understand the schema better.\n\n
    Here they are:
    ${topSimilarQueries.map((q) => `- Query: ${q.query} (similarity score: ${q.score.toFixed(4)})`).join('\n')}
    `;
  }
  // Initial system message
  const systemMessage = `You are an expert SQL query generator. You have access to tools to explore a database schema and generate SQL queries.

WORKFLOW:
1. You are provided with a list of available tables and their descriptions
2. Use describe_tables to get detailed column information for relevant tables
3. Use list_joins to understand relationships between tables
4. Use try_sql to test your generated queries (it returns row count and column names)
5. Iterate and refine your query based on test results
6. Once satisfied, provide the final SQL query

CRITICAL RULES:
- ALWAYS use list_joins to verify relationships before writing JOIN clauses
- NEVER assume joins exist based on column name similarity (unless no other option, or explicitly instructed by table/column comments)
- ONLY use relationships explicitly listed by list_joins
- Use try_sql to validate your query before finalizing
- Pay attention to semantic types, value ranges, and category values from describe_tables
- Schema name is: ${schemaName}
- Always prefix table names with schema name in SQL
- DO NOT USE 'to' as a table alias

When you have the final SQL query ready, respond with ONLY the SQL query text (no markdown, no explanation, no semicolon at the end).`;

  const messages: OpenAI.Chat.ChatCompletionMessageParam[] = [
    {role: 'system', content: systemMessage},
    {
      role: 'user',
      content: `Available tables in schema ${schemaName}:\n\n${initialTableList}\n\nUser Query: ${prompt}\n\n
      Explore the schema using the available tools and generate an SQL query to answer this question.
        
      ${similarQueriesInfo}
      `,
    },
  ];

  // Define available tools
  const tools: OpenAI.Chat.ChatCompletionTool[] = [
    {
      type: 'function',
      function: {
        name: 'list_tables',
        description: 'Returns the list of all table names along with their descriptions/comments. Use this to understand what data is available.',
        parameters: {
          type: 'object',
          properties: {},
          required: [],
        },
      },
    },
    {
      type: 'function',
      function: {
        name: 'describe_tables',
        description: 'Get detailed information about specific table(s) including all columns, their types, semantic types, comments, value ranges, and category values. Essential for understanding what data is in each table.',
        parameters: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names (without schema prefix) to describe',
            },
          },
          required: ['tables'],
        },
      },
    },
    {
      type: 'function',
      function: {
        name: 'list_joins',
        description: 'Lists all foreign key relationships (joins) that involve the specified table(s). CRITICAL: Use this before writing any JOIN clause to verify the relationship exists.',
        parameters: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names (without schema prefix) to find joins for',
            },
          },
          required: ['tables'],
        },
      },
    },
    {
      type: 'function',
      function: {
        name: 'try_sql',
        description: 'Execute an SQL query to test it. Returns row count and column names (limited to 10 rows). Use this to validate your query before providing the final answer. If row count is 0, the query might be wrong or the data might genuinely be absent.',
        parameters: {
          type: 'object',
          properties: {
            sql: {
              type: 'string',
              description: 'The SQL query to test (will be automatically limited to 10 rows)',
            },
          },
          required: ['sql'],
        },
      },
    },
  ];

  // Function calling loop
  let iterations = 0;
  const maxIterations = 15; // Prevent infinite loops

  while (iterations < maxIterations) {
    if (abortSignal.aborted) {
      grok.shell.info('SQL generation aborted by user.');
      try {
        SQLGenerationContext._lastFc?.cancel();
        SQLGenerationContext._lastFc = null;
      } catch (_) {
      }

      return '';
    }
    iterations++;
    console.log(`\n=== Iteration ${iterations} ===`);

    const response = await openai.chat.completions.create({
      model: modelName,
      messages,
      tools,
      temperature: 0,
    });

    const choice = response.choices[0];
    if (!choice) break;

    const message = choice.message;
    messages.push(message);

    // Check if the model wants to call functions
    if (message.tool_calls && message.tool_calls.length > 0) {
      const functionToolCalls = message.tool_calls.filter((tc) => tc.type === 'function');
      console.log(`Tool calls requested: ${functionToolCalls.map((tc) => tc.function.name).join(', ')}`);

      // Execute each tool call
      for (const toolCall of functionToolCalls) {
        if (toolCall.type !== 'function') continue;
        const functionName = toolCall.function.name;
        const functionArgs = JSON.parse(toolCall.function.arguments);

        console.log(`Executing ${functionName} with args:`, functionArgs);

        let result: string;
        try {
          switch (functionName) {
          case 'list_tables':
            result = await context.listTables();
            break;
          case 'describe_tables':
            result = await context.describeTables(functionArgs.tables);
            break;
          case 'list_joins':
            result = await context.listJoins(functionArgs.tables);
            break;
          case 'try_sql':
            result = await context.trySql(functionArgs.sql);
            break;
          default:
            result = `Error: Unknown function ${functionName}`;
          }
        } catch (error: any) {
          result = `Error executing ${functionName}: ${error.message}`;
        }

        console.log(`Result: ${result.substring(0, 200)}${result.length > 200 ? '...' : ''}`);

        // Add function result to messages
        messages.push({
          role: 'tool',
          tool_call_id: toolCall.id,
          content: result,
        });
      }
    } else if (message.content) {
      // Model provided a text response (hopefully the final SQL)
      console.log('Model response:', message.content);

      // Check if this looks like SQL (basic heuristic)
      const content = message.content.trim();
      if (content.toUpperCase().includes('SELECT') || content.toUpperCase().includes('WITH')) {
        // Clean up the response
        let sql = content;
        if (sql.startsWith('```sql'))
          sql = sql.substring(6);
        if (sql.startsWith('```'))
          sql = sql.substring(3);
        if (sql.endsWith('```'))
          sql = sql.substring(0, sql.length - 3);

        return sql.trim().replace(/;+$/, ''); // Remove trailing semicolons
      }

      // Model is explaining something, let it continue
      if (choice.finish_reason === 'stop')
        break;
    } else {
      // No content and no tool calls
      break;
    }

    // Check finish reason
    if (choice.finish_reason === 'stop')
      break;
  }

  // If we exhausted iterations or didn't get a proper response
  throw new Error('Failed to generate SQL query: Maximum iterations reached or no valid SQL returned');
}

/**
 * Context class that manages tool execution for SQL generation
 */
class SQLGenerationContext {
  private schema?: DBSchemaMeta;

  constructor(
    private connectionID: string,
    private schemaName: string,
    private dbMeta: DBConnectionMeta | undefined,
    private connection: DG.DataConnection
  ) {
    if (dbMeta)
      this.schema = dbMeta.schemas.find((s) => s.name === schemaName);
  }

  /**
   * Tool: list_tables
   * Returns all tables with their descriptions
   */
  async listTables(): Promise<string> {
    if (this.schema && this.dbMeta) {
      const tables = this.schema.tables.map((table) => {
        const parts = [`${table.name}`];
        if (table.LLMComment)
          parts.push(`  Description: ${table.LLMComment}`);
        else if (table.comment)
          parts.push(`  Description: ${table.comment}`);
        parts.push(`  Columns: ${table.columns.length}, Rows: ${table.rowCount}`);
        return parts.join('\n');
      });
      return `Tables in ${this.schemaName}:\n\n${tables.join('\n\n')}`;
    }

    // Fallback to schema descriptor
    const schemaInfo = await DBSchemaInfo.getSchemaInfo(this.connectionID, this.schemaName);
    return `Tables in ${this.schemaName}:\n\n${schemaInfo.tableInfos.join('\n')}`;
  }

  /**
   * Tool: describe_tables
   * Returns detailed column information for specified tables
   */
  async describeTables(tableNames: string[]): Promise<string> {
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    const descriptions: string[] = [];

    for (const tableName of tableNames) {
      if (this.schema) {
        const table = this.schema.tables.find((t) => t.name === tableName);
        if (!table) {
          descriptions.push(`Table ${tableName} not found`);
          continue;
        }

        const tableDesc = this.buildDetailedTableDescription(table);
        descriptions.push(tableDesc);
      } else {
        // Fallback: basic info from schema descriptor
        const schemaInfo = await DBSchemaInfo.getSchemaInfo(this.connectionID, this.schemaName);
        const tableInfo = schemaInfo.tableInfos.find((info) => info.startsWith(`${tableName}:`));
        if (tableInfo)
          descriptions.push(tableInfo);
        else
          descriptions.push(`Table ${tableName} not found`);
      }
    }

    return descriptions.join('\n\n---\n\n');
  }

  /**
   * Tool: list_joins
   * Returns all foreign key relationships involving the specified tables
   */
  async listJoins(tableNames: string[]): Promise<string> {
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    const tableSet = new Set(tableNames);
    const joins: string[] = [];

    if (this.dbMeta) {
      for (const relation of this.dbMeta.relations) {
        if (relation.fromSchema !== this.schemaName && relation.toSchema !== this.schemaName)
          continue;

        // Include if any of the specified tables is involved
        if (tableSet.has(relation.fromTable) || tableSet.has(relation.toTable)) {
          const joinDesc = this.buildJoinDescription(relation);
          joins.push(joinDesc);
        }
      }
    } else {
      // Fallback to schema descriptor
      const schemaInfo = await DBSchemaInfo.getSchemaInfo(this.connectionID, this.schemaName);
      for (const ref of schemaInfo.referenceInfos) {
        const [from, to] = ref.split('->').map((s) => s.trim());
        const fromTable = from.split('.')[0];
        const toTable = to.split('.')[0];

        if (tableSet.has(fromTable) || tableSet.has(toTable))
          joins.push(ref);
      }
    }

    if (joins.length === 0)
      return `No foreign key relationships found for tables: ${tableNames.join(', ')}`;

    return `Foreign key relationships involving ${tableNames.join(', ')}:\n\n${joins.join('\n')}`;
  }

  static _lastFc: DG.FuncCall | null = null;
  /**
   * Tool: try_sql
   * Executes SQL and returns row count and column names
   */
  async trySql(sql: string): Promise<string> {
    let sub: rxjs.Subscription | null = null;
    try {
      // Add LIMIT if not present
      let testSql = sql.trim();
      if (testSql.endsWith(';'))
        testSql = testSql.substring(0, testSql.length - 1);
      if (!testSql.toUpperCase().includes('LIMIT'))
        testSql += ' LIMIT 10';

      const queryName = `test-query-${Date.now()}`;
      const wrappedSql = `--name: ${queryName}\n${testSql}`;

      const fc = this.connection.query(queryName, wrappedSql).prepare({});
      SQLGenerationContext._lastFc = fc;
      sub = grok.events.onCustomEvent(AI_SQL_QUERY_ABORT_EVENT).subscribe(() => {
        try {
          fc.cancel();
        } catch (_) {
        } finally {
          SQLGenerationContext._lastFc = null;
          sub?.unsubscribe();
        }
      });
      await fc.call(false, undefined, {processed: true, report: false});
      const df = await fc.getOutputParamValue();
      SQLGenerationContext._lastFc = null;
      sub.unsubscribe();
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
    } catch (error: any) {
      return `SQL Error: ${error?.message}\n\nThis query failed to execute. Please revise the SQL based on the schema information.`;
    }
  }

  private formatColValue(value: any): string {
    if (value === null || value === undefined || value === '' || value === DG.FLOAT_NULL || value === DG.INT_NULL)
      return 'NULL';
    try {
      const valueString = value?.toString();
      return valueString.length > 50 ? `'${valueString.substring(0, 47)}...'` : `'${valueString}'`;
    } catch (_e) {
      return 'NON-STRINGIFIABLE-VALUE';
    }
  }

  /**
   * Helper: Build detailed table description
   */
  private buildDetailedTableDescription(table: DBTableMeta): string {
    const lines: string[] = [];

    lines.push(`TABLE: ${table.name} (${table.rowCount} rows)`);
    if (table.comment)
      lines.push(`Description: ${table.comment}`);
    if (table.LLMComment)
      lines.push(`AI Context: ${table.LLMComment}`);

    lines.push('\nColumns:');

    for (const col of table.columns) {
      const colParts: string[] = [`  ${col.name} (${col.type})`];

      if (col.semanticType)
        colParts.push(`[Semantic: ${col.semanticType}]`);
      if (col.isUnique)
        colParts.push('[UNIQUE]');
      if (col.LLMComment)
        colParts.push(`- ${col.LLMComment}`);
      else if (col.comment)
        colParts.push(`- ${col.comment}`);
      if (col.min !== undefined && col.max !== undefined)
        colParts.push(`Range: ${col.min} to ${col.max}`);
      if (col.categoryValues && col.categoryValues.length > 0) {
        const values = col.categoryValues.slice(0, 15);
        colParts.push(`Values: ${values.join(', ')}${col.categoryValues.length > 15 ? ', ...' : ''}`);
      }

      lines.push(colParts.join(' '));
    }

    return lines.join('\n');
  }

  /**
   * Helper: Build join description
   */
  private buildJoinDescription(relation: DBRelationMeta): string {
    const parts: string[] = [];

    parts.push(`${relation.fromTable}(${relation.fromColumns.join(', ')}) -> ${relation.toTable}(${relation.toColumns.join(', ')})`);

    if (relation.cardinality)
      parts.push(`[${relation.cardinality}]`);
    if (relation.IsPrimaryPath === false)
      parts.push('[LEGACY - prefer other paths if available]');
    if (relation.LLMComment)
      parts.push(`- ${relation.LLMComment}`);
    else if (relation.comment)
      parts.push(`- ${relation.comment}`);

    return parts.join(' ');
  }
}
