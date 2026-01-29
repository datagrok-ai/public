/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getAIAbortSubscription} from '../utils';
import * as _rxjs from 'rxjs';
import {SemValueObjectHandler} from '@datagrok-libraries/db-explorer/src/object-handlers';
import {ChatModel} from 'openai/resources/index';
import {AIPanelFuncs, MessageType} from './panel';
import {BuiltinDBInfoMeta, getDBColumnMetaData, getDBTableMetaData} from './query-meta-utils';
import {ModelType, OpenAIClient} from './openAI-client';
import {OpenAIResponsesProvider} from './AI-API-providers/openai-responses-provider';
import {OpenAIChatCompletionsProvider} from './AI-API-providers/openai-chat-completions-provider';
import {AIProvider, FunctionToolSpec} from './AI-API-providers/types';


const suspiciousSQlPatterns = ['DROP ', 'DELETE ', 'UPDATE ', 'INSERT ', 'ALTER ', 'CREATE ', 'TRUNCATE ', 'EXEC ', 'MERGE '];

type JsonPrimitive = string | number | boolean | null;
type JsonValue = JsonPrimitive | JsonObject | JsonValue[];
type JsonObject = { [key: string]: JsonValue };

function isJsonObject(value: JsonValue | object | null | undefined): value is JsonObject {
  return !!value && typeof value === 'object' && !Array.isArray(value);
}

function getStringProp(obj: JsonObject | null | undefined, key: string): string | null {
  const v = obj?.[key];
  return typeof v === 'string' ? v : null;
}

function getStringArrayProp(obj: JsonObject | null | undefined, key: string): string[] {
  const v = obj?.[key];
  if (!Array.isArray(v))
    return [];
  return v.filter((item): item is string => typeof item === 'string');
}

function parseJsonObject(text: string): JsonObject | null {
  try {
    const parsed = JSON.parse(text) as JsonValue;
    return isJsonObject(parsed) ? parsed : null;
  } catch {
    return null;
  }
}

/**
 * Generates SQL query using function calling approach where LLM can explore schema interactively
 * @throws Error if SQL generation fails
 * @param prompt - User's natural language query
 * @param connectionID - Database connection ID
 * @param schemaName - Schema name
 * @param dbMeta - Optional DBConnectionMeta object with enriched metadata
 */
export async function generateAISqlQueryWithTools(
  prompt: string,
  connectionID: string,
  schemaName: string,
  options: {
    oldMessages?: MessageType[],
    aiPanel?: AIPanelFuncs<MessageType>,
    modelName?: ChatModel,
    disableVerbose?: boolean,
  } = {}
): Promise<string> {
  let aborted = false;

  options.aiPanel?.addUiMessage(prompt, true);

  // Load pre-indexed metadata if available
  const connection = await grok.dapi.connections.find(connectionID);
  const dbInfo = await BuiltinDBInfoMeta.fromConnection(connection);
  // Create tool execution context
  const context = new SQLGenerationContext(connectionID, schemaName, dbInfo, connection);
  const abortSub = getAIAbortSubscription().subscribe(() => {
    aborted = true;
    console.log('Aborting SQL generation as per user request');
    try {
      SQLGenerationContext._lastFc?.cancel();
      SQLGenerationContext._lastFc = null;
    // eslint-disable-next-line @typescript-eslint/no-unused-vars
    } catch (_) {
    }
    abortSub.unsubscribe();
  });
  try {
  // Initialize OpenAI client provider
    const provider = AIProvider.getProvider();

    // Initial system message
    const systemMessage = `You are an expert SQL query generator. You have access to tools to explore a database schema and generate SQL queries.

    Your name is Datagrok SQL Expert.

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
    - ALWAYS Use try_sql to validate your query before finalizing and submitting it
    - Pay attention to semantic types, value ranges, and category values from describe_tables
    - Current Schema name is: ${schemaName}
    - ALWAYS prefix table names with schema name in SQL or other functions!!!
    - DO NOT USE 'to' as a table alias
    - If some categorical value is supplied to match (e.g. status value or measurement units), use the information from corresponding column's category values (if any). otherwise, try to use multiple options using OR (for example for micromolar units, try 'uM', 'μM', ets).
    - Some queries might require cartridge use (like RDKIT functions for similarity / substructure search). Use provided query examples and your knowledge of such functions as needed.
    - When working with large tables, try not to use order unless explicitly requested.
    - VERY IMPORTANT: IN some cases, user might want to just get a plain reply without sql query. In such cases, Always begin your response with "REPLY ONLY:"!!!
    - If you fail to generate a valid working query, admit it politely instead of making something up and ask for more information from the user.
    - If the user prompt is ambiguous, ask for clarifications instead of guessing.
    - If user is plainly trying to communicate or chat or ask general questions, use REPLY ONLY mode to answer them appropriately.
    - Use the reply_to_user tool to communicate with the user, explain your reasoning during chain of though and execution and so on.
    When you have the final SQL query ready, respond with ONLY the SQL query text (no markdown, no explanation, no semicolon at the end).`;

    const input: MessageType[] = options.oldMessages ? [...options.oldMessages] : [];
    const addEngineMessage = (message: MessageType) => {
      options.aiPanel?.addEngineMessage(message);
    };
    if (input.length === 0) {
      // Get initial table list
      const initialTableList = (await context.listTables(false, schemaName)).description;
      let semTypeWithinPromptInfo = '';
      // try to also match some identifiers within the user prompt to provide even more context
      const parsedPrompt = DG.SemanticValue.parse(prompt);
      if (parsedPrompt?.semType && parsedPrompt.value) {
        // @ts-ignore
        const objHandler = DG.ObjectHandler.list().find((oh) => (oh as SemValueObjectHandler)?.columnName && (oh as SemValueObjectHandler)?.tableName && oh.type === parsedPrompt.semType) as SemValueObjectHandler | undefined;
        if (objHandler) {
          // TODO: remove ts ignores after updating db-explorer library in plugins
          // @ts-ignore
          const [entryPointTable, entryPointColumn] = [objHandler.tableName, objHandler.columnName];
          semTypeWithinPromptInfo = `\n\n Note: The ${parsedPrompt.value} mentioned in the prompt is identified as a ${parsedPrompt.semType} type, typically found in the ${entryPointTable}.${entryPointColumn} column. Use this information to guide your SQL generation.`;
          console.log(semTypeWithinPromptInfo);
        }
      }

      if (grok.ai.config.indexEntities) {
        const similarQueries = await findSimilarQueriesToPrompt(prompt);
        console.log(`found ${similarQueries.n} similar queries for the prompt`);
        if (similarQueries.n > 0 && similarQueries.text) {
          semTypeWithinPromptInfo += `\n\n Additionally, here are some similar previously executed queries that might help you:
          \n${similarQueries.text}
          `;
        }
      }

      const systemMsg = provider.createSystemMessage(systemMessage);
      input.push(systemMsg);
      addEngineMessage(systemMsg);
      const initialUserMsg = provider.createUserMessage(`
        Available schemas: ${(await context.listAllSchemas()).join(', ')}\n\n
        Available tables in schema ${schemaName}:\n\n${initialTableList}\n\nUser Query: ${prompt}\n\n
        Explore the schema using the available tools and generate an SQL query to answer this question.
        ${semTypeWithinPromptInfo}`);
      addEngineMessage(initialUserMsg); // prompt will be shown in UI
      input.push(initialUserMsg);
    } else {
      const followUpUserMsg = provider.createUserMessage(`User Follow up (do modifications as/if needed based on the following): ${prompt}`);
      input.push(followUpUserMsg);
      addEngineMessage(followUpUserMsg);
    }

    // Define available tools
    const functionTools: FunctionToolSpec[] = [
      {
        type: 'function',
        name: 'list_tables_in_schema',
        description: 'Returns the list of all table names along with their descriptions/comments. Use this to understand what data is available. You must provide schemaName parameter to specify which schema to list tables from.',
        parameters: {
          type: 'object',
          properties: {
            schemaName: {
              type: 'string',
              description: 'Name of the schema to list tables from',
            },
          },
          required: ['schemaName'],
          additionalProperties: false,
        },
        strict: true,
      },
      {
        type: 'function',
        name: 'list_all_schemas',
        description: 'Retyurns the names all all schemas in the connection.',
        parameters: {
          type: 'object',
          properties: {},
          required: [],
          additionalProperties: false,
        },
        strict: true,
      },

      {
        type: 'function',
        name: 'describe_tables',
        description: 'Get detailed information about specific table(s) including all columns, their types, semantic types, comments, value ranges, and category values. Essential for understanding what data is in each table.',
        parameters: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names (MANDATORY WITH schema prefix) to describe. ALWAYS provide schema name as prefix, for example public.tableName...',
            },
          },
          required: ['tables'],
          additionalProperties: false,
        },
        strict: true,
      },
      {
        type: 'function',
        name: 'list_joins',
        description: 'Lists all foreign key relationships (joins) that involve the specified table(s). CRITICAL: Use this before writing any JOIN clause to verify the relationship exists.',
        parameters: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names (WITH schema prefix) to find joins for',
            },
          },
          required: ['tables'],
          additionalProperties: false,
        },
        strict: true,
      },
      {
        type: 'function',
        name: 'try_sql',
        description: 'Execute an SQL query to test it. Returns row count and column names (limited to 10 rows). Use this to validate your query before providing the final answer. If row count is 0, the query might be wrong or the data might genuinely be absent. TRY NOT TO USE ORDER BY in YOUR TEST QUERIES on large tables unless absolutely necessary.',
        parameters: {
          type: 'object',
          properties: {
            sql: {
              type: 'string',
              description: 'The SQL query to test (will be automatically limited to 10 rows)',
            },
            description: {
              type: 'string',
              description: 'Short description of what this SQL is trying to achieve in markdown format. This will be put in ui for user context.',
            }
          },
          required: ['sql', 'description'],
          additionalProperties: false,
        },
        strict: true,
      },
      // plain function tool for general purpose question answering
      {
        type: 'function',
        name: 'reply_to_user',
        description: 'Use this tool to esentially communicate with the user directly, give feedback, add context to what you are doing, explain reasoning(!!!) and so on. This is NOT for providing the final SQL query, but rather for intermediate communication. When using this tool, always provide a clear, informative and well formatted (markdown) message to the user in the reply parameter.',
        parameters: {
          type: 'object',
          properties: {
            reply: {
              type: 'string',
              description: 'Your text reply to the user question/prompt. this message will be shown in the UI, so make it clear, informative and formatted in markdown as needed.',
            },
          },
          required: ['reply'],
          additionalProperties: false,
        },
        strict: true,
      }
    ];

    // if we have some embeddings, then also add a tool so that LLM can use it to find similar queries
    if (grok.ai.config.indexEntities) {
      functionTools.push({
        type: 'function',
        name: 'find_similar_queries',
        description: 'Finds existing similar queries based on the prompt you provide (which should be based on user question). Use this to get inspiration from previously executed queries that might be similar to the user question. Return up to 3 similar queries with their descriptions/comments. Similarity is based on LLM embeddings and cosine similarity.',
        parameters: {
          type: 'object',
          properties: {
            prompt: {
              type: 'string',
              description: 'prompt derived from user question to find similar queries for',
            },
          },
          required: ['prompt'],
          additionalProperties: false,
        },
        strict: true,
      });
    }


    async function findSimilarQueriesToPrompt(aiPrompt: string) {
      try {
        const resultEntities = (await grok.ai.searchEntities(aiPrompt, 0.3, 10, ['DataQuery']))
          .filter((e) => !!e && e instanceof DG.DataQuery && e.connection?.id === connectionID)
          .filter((_, i) => i < 4) as DG.DataQuery[]; // limit to top 4

        if (resultEntities.length === 0)
          return {text: 'No similar queries found for this connection', n: 0};
        return {text: `Here are some similar previously executed queries that might help you.
          Note that these queries might have annotations/comments like name descriptions and so on that can help you understand the schema better.\n\n
          Here they are:
          ${resultEntities.map((q) => `- Query:\n ${q.query} \n`).join('\n\n')}`, n: resultEntities.length};
      } catch (error) {
        console.error('Error generating embedding for similar query search:', error);
        return {text: 'Error generating embedding for similar query search', n: 0};
      }
    };

    // Function calling loop
    let iterations = 0;
    let maxIterations = 15; // Prevent infinite loops

    while (iterations < maxIterations) {
      if (aborted) {
        grok.shell.info('SQL generation aborted by user.');
        return '';
      }
      iterations++;
      console.log(`\n=== Iteration ${iterations} ===`);
      const modelName = options?.modelName ?? ModelType.Fast;
      const response = await provider.create({
        model: modelName,
        messages: input,
        tools: functionTools,
        ...(modelName.startsWith('gpt-5') ? {reasoning: {effort: 'high'}} : {}),
      });

      const outputs = response.outputMessages;
      input.push(...outputs);
      outputs.forEach((item) => addEngineMessage(item));

      // Execute function calls (if any)
      let hadToolCalls = false;
      for (const call of response.toolCalls) {
        hadToolCalls = true;

        const functionName = call.name;
        const callId = call.id;
        const rawArgs = call.arguments ?? '{}';
        const args = parseJsonObject(rawArgs);

        if (!functionName || !callId) {
          const errorText = 'Error: malformed function call (missing name/call_id).';
          const outputItem = provider.createToolOutputMessage(callId ?? '', errorText);
          if (callId)
            input.push(outputItem);
          if (callId)
            addEngineMessage(outputItem);
          continue;
        }

        let result = '';
        try {
          switch (functionName) {
          case 'list_all_schemas':
            options.aiPanel?.addUiMessage(`📂 Listing all schemas in the database.`, false);
            result = (await context.listAllSchemas()).join(', ');
            !options.disableVerbose && options.aiPanel?.addUiMessage(`✅ Schemas found: *${result}*`, false);
            break;
          case 'reply_to_user':
            options.aiPanel?.addUiMessage(`💬 ${getStringProp(args, 'reply') ?? ''}`, false);
            result = 'Reply sent to user.';
            break;
          case 'list_tables_in_schema':
            options.aiPanel?.addUiMessage(`📑 Listing all tables in schema *${getStringProp(args, 'schemaName') ?? ''}*.`, false);
            const res = await context.listTables(false, getStringProp(args, 'schemaName') ?? schemaName);
            result = res.description;
            !options.disableVerbose && options.aiPanel?.addUiMessage(res.tableCount ? `✅ Found ${res.tableCount} tables` : '⚠️ No tables found', false);
            break;
          case 'describe_tables': {
            const tables = getStringArrayProp(args, 'tables');
            options.aiPanel?.addUiMessage(`📋 Getting content of following tables: *${tables.join(', ')}*.`, false);
            result = await context.describeTables(tables);
            break;
          }
          case 'list_joins': {
            const tables = getStringArrayProp(args, 'tables');
            options.aiPanel?.addUiMessage(`🔗 Listing relations for tables: *${tables.join(', ')}*.`, false);
            result = await context.listJoins(tables);
            break;
          }
          case 'try_sql': {
            const sql = getStringProp(args, 'sql') ?? '';
            const description = getStringProp(args, 'description') ?? '';
            options.aiPanel?.addUiMessage(`🧪 Testing SQL query:\n${description}\n\`\`\`sql\n${sql}\n\`\`\``, false);
            result = await context.trySql(sql, description);
            !options.disableVerbose && options.aiPanel?.addUiMessage(`📊 SQL Test Result: \n\n${result}`, false);
            break;
          }
          case 'find_similar_queries': {
            options.aiPanel?.addUiMessage(`🔍 Searching for similar queries to help with SQL generation.`, false);
            const similarQueries = await findSimilarQueriesToPrompt(getStringProp(args, 'prompt') ?? '');
            !options.disableVerbose && options.aiPanel?.addUiMessage(similarQueries.n > 0 ? `✅ Found ${similarQueries.n} similar queries.` : `⚠️ ${similarQueries.text}`, false);
            result = similarQueries.text;
            break;
          }
          default:
            result = `Error: Unknown function ${functionName}`;
          }
        } catch (error) {
          const message = error instanceof Error ? error.message : String(error);
          result = `Error executing ${functionName}: ${message}`;
        }

        console.log(`Result: ${result.substring(0, 200)}${result.length > 200 ? '...' : ''}`);

        const outputItem = provider.createToolOutputMessage(callId, result);
        input.push(outputItem);
        addEngineMessage(outputItem);
      }

      // If the model called tools, we must continue so it can consume tool outputs.
      if (hadToolCalls) {
        if (iterations >= maxIterations && options.aiPanel) {
          const cont = await options.aiPanel.addConfirmMessage('Maximum iterations reached. Do you want to continue generating the SQL query?');
          if (cont)
            maxIterations += 10; // reset max iterations to continue
        }
        continue;
      }

      const content = response.outputText.trim();
      if (content.length === 0)
        throw new Error('Model returned no text and no tool calls');

      console.log('Model response:', content);

      // Check if this looks like SQL (basic heuristic)
      let sql = content;
      if (sql.startsWith('```sql'))
        sql = sql.substring(6);
      if (sql.startsWith('```'))
        sql = sql.substring(3);
      if (sql.endsWith('```'))
        sql = sql.substring(0, sql.length - 3);
      const contentUpperCase = sql.toUpperCase().trim();
      if (contentUpperCase.length > 2 && (contentUpperCase.startsWith('SELECT') || contentUpperCase.startsWith('WITH'))) {
        // Final SQL detected
        const res = sql.trim().replace(/;+$/, ''); // Remove trailing semicolons
        const resUpperCase = res.toUpperCase();
        if (suspiciousSQlPatterns.some((pattern) => resUpperCase.includes(pattern))) {
          const out = await new Promise<string>((resolve) => {
            ui.dialog('Potentially Destructive SQL Detected')
              .add(ui.divText('The generated SQL query contains potentially destructive commands. For safety, please confirm if you want to proceed with this query.'))
              .add(ui.markdown(`\`\`\`sql\n${res}\n\`\`\``))
              .onOK(() => resolve(res))
              .onCancel(() => {
                console.log('User cancelled execution of potentially destructive SQL');
                resolve('');
              })
              .show();
          });
          options.aiPanel?.addUiMessage(!out ? 'User cancelled execution of potentially destructive SQL.' : 'User confirmed execution of potentially destructive SQL.', false);
          if (out)
            !options.disableVerbose && options.aiPanel?.addUiMessage(`Final SQL Query:\n\`\`\`sql\n${out}\n\`\`\``, false, {finalResult: out});
          return out;
        }
        if (res)
          !options.disableVerbose && options.aiPanel?.addUiMessage(`Final SQL Query:\n\`\`\`sql\n${res}\n\`\`\``, false, {finalResult: res});
        return res;
      } else {
        const replyText = contentUpperCase.startsWith('REPLY ONLY:') ? content.substring('REPLY ONLY:'.length).trim() : content;
        options.aiPanel?.addUiMessage(replyText, false);
        return '';
      }
    }
    // If we exhausted iterations or didn't get a proper response
    throw new Error('Failed to generate SQL query: Maximum iterations reached or no valid SQL returned');
  } finally {
    abortSub.unsubscribe();
  }
}

/**
 * Context class that manages tool execution for SQL generation
 */
class SQLGenerationContext {
  constructor(
    private connectionID: string,
    private schemaName: string,
    private dbMeta: BuiltinDBInfoMeta,
    private connection: DG.DataConnection
  ) {}

  async listAllSchemas(): Promise<string[]> {
    return (await this.dbMeta.getSchemas()).map((s) => s.name);
  }

  /**
   * Tool: list_tables
   * Returns all tables with their descriptions
   */
  async listTables(includeAllDescriptions: boolean, schemaName: string): Promise<{description: string, tableCount: number}> {
    schemaName ??= this.schemaName;

    const schema = (await this.dbMeta.getSchemas()).find((s) => s.name === schemaName);
    if (!schema)
      return {description: `Schema ${schemaName} not found in metadata.`, tableCount: 0};

    const tables = (await schema.getTables());
    //const parts: string[] = [];
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
   * Returns detailed column information for specified tables
   */
  async describeTables(tableNames: string[]): Promise<string> {
    const schemaNames = tableNames.map((t) => t.trim().indexOf('.') >= 0 ? t.split('.')[0] : this.schemaName);
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    const descriptions: string[] = [];

    const schemas = await this.dbMeta.getSchemas();
    const schemaTableMap: Map<string, DG.TableInfo[]> = new Map();
    for (const sn of Array.from(new Set(schemaNames))) {
      const schema = schemas.find((s) => s.name === sn);
      if (schema) {
        const tables = await schema.getTables();
        schemaTableMap.set(sn, tables);
      }
    }
    for (let i = 0; i < tableNames.length; i++) {
      const tables = schemaTableMap.get(schemaNames[i]);
      if (!tables) {
        descriptions.push(`Schema ${schemaNames[i]} not found`);
        continue;
      }
      const table = tables.find((t) => t.friendlyName === tableNames[i] || t.name === tableNames[i]);
      if (!table) {
        descriptions.push(`Table ${tableNames[i]} not found`);
        continue;
      }
      const tableDesc = this.buildDetailedTableDescription(table);
      descriptions.push(tableDesc);
    }

    return descriptions.join('\n\n---\n\n');
  }

  /**
   * Tool: list_joins
   * Returns all foreign key relationships involving the specified tables
   */
  async listJoins(tableNames: string[]): Promise<string> {
    const schemaNames = tableNames.map((t) => t.trim().indexOf('.') >= 0 ? t.split('.')[0] : this.schemaName);
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    // in case of cross-schema requests,
    const schemaQualifiedTableSet = new Set<string>();
    for (let i = 0; i < tableNames.length; i++)
      schemaQualifiedTableSet.add(`${schemaNames[i]}.${tableNames[i]}`);

    // const joins: string[] = [];
    const uniqueSchemas = Array.from(new Set(schemaNames));
    const relations = (await Promise.all(uniqueSchemas.map(async (sn) => await this.dbMeta.getRelationsForSchema(sn))))
      .flat()
      .filter((rel) => schemaQualifiedTableSet.has(`${rel.fromSchema}.${rel.fromTable}`) || schemaQualifiedTableSet.has(`${rel.toSchema}.${rel.toTable}`));
    const joins = relations.map((rel) => this.buildJoinDescription(rel));

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
