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
import {ModelOption, ModelType, LLMClient} from './LLM-client';
import {LanguageModelV3FunctionTool, LanguageModelV3Message} from '@ai-sdk/provider';
import {findLast} from '../utils';


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
 * @param options - Options including catalogName, old messages, AI panel, etc.
 */
export async function generateAISqlQueryWithTools(
  prompt: string,
  connectionID: string,
  options: {
    catalogName?: string,
    oldMessages?: MessageType[],
    aiPanel?: AIPanelFuncs<MessageType>,
    modelType?: ModelOption,
    disableVerbose?: boolean,
  } = {}
): Promise<string> {
  let aborted = false;

  options.aiPanel?.addUiMessage(prompt, true);

  // Load pre-indexed metadata if available
  const connection = await grok.dapi.connections.find(connectionID);
  const allDbInfos = await BuiltinDBInfoMeta.allFromConnection(connection);
  const catalogNames = allDbInfos.map((d) => d.name);
  const defaultCatalog = options.catalogName ??
    connection.parameters?.['catalog'] ??
    connection.parameters?.['db'] ??
    catalogNames[0] ??
    '';
  // Create tool execution context
  const context = new SQLGenerationContext(connectionID, defaultCatalog, allDbInfos, connection);
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
    const langTool = LLMClient.getInstance();
    const modelType = options.modelType ?? 'Fast';
    const client = langTool.aiModels[modelType];

    // Initial system message
    const systemMessage = `You are an expert SQL query generator. You have access to tools to explore a database schema and generate SQL queries.

    Your name is Datagrok SQL Expert.

    WORKFLOW:
    1. You are provided with the default catalog and available schemas
    2. Use list_tables_in_schema to see what tables are available in a schema
    3. Use describe_tables to get detailed column information for relevant tables
    4. Use list_joins to understand relationships between tables
    5. Use try_sql to test your generated queries (it returns row count and column names)
    6. Iterate and refine your query based on test results
    7. Once satisfied, provide the final SQL query

    CRITICAL RULES:
    - ALWAYS use list_joins to verify relationships before writing JOIN clauses
    - NEVER assume joins exist based on column name similarity (unless no other option, or explicitly instructed by table/column comments)
    - ONLY use relationships explicitly listed by list_joins
    - ALWAYS Use try_sql to validate your query before finalizing and submitting it
    - Pay attention to semantic types, value ranges, and category values from describe_tables
    - The default database catalog is: ${defaultCatalog}. Use list_catalogs if you need to explore other catalogs.
    - Use list_all_schemas and list_tables_in_schema to discover the database structure before writing queries.
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
      // Get initial schema list for the default catalog
      const defaultSchemas = await context.listAllSchemas(defaultCatalog);
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

      const systemMsg = langTool.createSystemMessage(systemMessage);
      input.push(systemMsg);
      addEngineMessage(systemMsg);
      const initialUserMsg = langTool.createUserMessage(`
        Default catalog: ${defaultCatalog}\n
        Available catalogs: ${catalogNames.join(', ')}\n\n
        Available schemas in catalog '${defaultCatalog}': ${defaultSchemas.join(', ')}\n\nUser Query: ${prompt}\n\n
        Explore the database using the available tools and generate an SQL query to answer this question.
        ${semTypeWithinPromptInfo}`);
      addEngineMessage(initialUserMsg); // prompt will be shown in UI
      input.push(initialUserMsg);
    } else {
      const followUpUserMsg = langTool.createUserMessage(`User Follow up (do modifications as/if needed based on the following): ${prompt}`);
      input.push(followUpUserMsg);
      addEngineMessage(followUpUserMsg);
    }

    // Define available tools
    const functionTools: LanguageModelV3FunctionTool[] = [
      {
        type: 'function',
        name: 'list_tables_in_schema',
        description: 'Returns the list of all table names along with their descriptions/comments. Use this to understand what data is available.',
        inputSchema: {
          type: 'object',
          properties: {
            schemaName: {
              type: 'string',
              description: 'Name of the schema to list tables from',
            },
            catalogName: {
              type: 'string',
              description: 'Name of the catalog containing the schema',
            },
          },
          required: ['schemaName', 'catalogName'],
          additionalProperties: false,
        },
        strict: true,
      },
      {
        type: 'function',
        name: 'list_all_schemas',
        description: 'Returns the names of all schemas in the specified catalog.',
        inputSchema: {
          type: 'object',
          properties: {
            catalogName: {
              type: 'string',
              description: 'Name of the catalog to list schemas from',
            },
          },
          required: ['catalogName'],
          additionalProperties: false,
        },
        strict: true,
      },
      {
        type: 'function',
        name: 'list_catalogs',
        description: 'Returns the names of all database catalogs available on this connection. Use this if you need to explore data across different catalogs.',
        inputSchema: {
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
        inputSchema: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names to describe. Use schema.table format (e.g. public.tableName) or catalog.schema.table format for cross-catalog queries.',
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
        inputSchema: {
          type: 'object',
          properties: {
            tables: {
              type: 'array',
              items: {type: 'string'},
              description: 'List of table names to find joins for. Use schema.table or catalog.schema.table format.',
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
        inputSchema: {
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
        inputSchema: {
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
        inputSchema: {
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
      const response = await client.doGenerate({
        prompt: input,
        tools: functionTools,
        providerOptions: {
          openai: {
            ...(ModelType.Coding.startsWith('gpt-5') ? {reasoning: {effort: 'medium'}} : {}),
          }
        }
      });

      // Fix tool_use input fields - Anthropic expects objects, not strings
      // The AI SDK's doGenerate puts item IDs in providerMetadata, but the
      // Responses API input converter only reads reasoning items from providerOptions
      // (tool-call has a fallback to providerMetadata, but reasoning does not).
      // Copy providerMetadata → providerOptions so reasoning item_references are created.
      const fixedContent = response.content.map((item) => {
        const patched = {
          ...item,
          providerOptions: (item as any).providerOptions ?? (item as any).providerMetadata,
        };
        if (item.type === 'tool-call') {
          return {
            ...patched,
            input: typeof item.input === 'string' ? parseJsonObject(item.input) ?? item.input : item.input
          };
        }
        return patched;
      });

      const formattedOutput: LanguageModelV3Message = {
        role: 'assistant',
        // @ts-ignore
        content: fixedContent
      };
      //const outputs: MessageType[] = response.content.map((item) => ({role: 'assistant', content: [item]} as MessageType));
      input.push(formattedOutput);
      options.aiPanel?.addEngineMessage(formattedOutput);

      // Execute function calls (if any)
      let hadToolCalls = false;
      for (const call of response.content.filter((c) => c.type === 'tool-call')) {
        hadToolCalls = true;

        const functionName = call.toolName;
        const callId = call.toolCallId;
        const rawArgs = call.input ?? '{}';
        const args = parseJsonObject(rawArgs);

        if (!functionName || !callId) {
          const errorText = 'Error: malformed function call (missing name/call_id).';
          const outputItem = langTool.createToolOutputMessage(callId ?? '', functionName ?? '', errorText);
          if (callId)
            input.push(outputItem);
          if (callId)
            addEngineMessage(outputItem);
          continue;
        }

        let result = '';
        try {
          switch (functionName) {
          case 'list_catalogs':
            options.aiPanel?.addUiMessage(`📂 Listing all catalogs in the database.`, false);
            result = (await context.listCatalogs()).join(', ');
            !options.disableVerbose && options.aiPanel?.addUiMessage(`✅ Catalogs found: *${result}*`, false);
            break;
          case 'list_all_schemas':
            options.aiPanel?.addUiMessage(`📂 Listing all schemas in catalog *${getStringProp(args, 'catalogName') ?? ''}*.`, false);
            result = (await context.listAllSchemas(getStringProp(args, 'catalogName') ?? defaultCatalog)).join(', ');
            !options.disableVerbose && options.aiPanel?.addUiMessage(`✅ Schemas found: *${result}*`, false);
            break;
          case 'reply_to_user':
            options.aiPanel?.addUiMessage(`💬 ${getStringProp(args, 'reply') ?? ''}`, false);
            result = 'Reply sent to user.';
            break;
          case 'list_tables_in_schema':
            options.aiPanel?.addUiMessage(`📑 Listing all tables in schema *${getStringProp(args, 'schemaName') ?? ''}* (catalog: *${getStringProp(args, 'catalogName') ?? ''}*).`, false);
            const res = await context.listTables(false, getStringProp(args, 'schemaName') ?? '', getStringProp(args, 'catalogName') ?? undefined);
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

        const outputItem = langTool.createToolOutputMessage(callId, functionName, result);
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

      const content = findLast(response.content, (c) => c.type === 'text')?.text!;
      if ((content?.length ?? 0) === 0)
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
      } else if (contentUpperCase.startsWith('REPLY ONLY:')) {
        const replyText = contentUpperCase.startsWith('REPLY ONLY:') ? content.substring('REPLY ONLY:'.length).trim() : content;
        options.aiPanel?.addUiMessage(replyText, false);
      } else {
        // Not SQL - add a message to AI to remind it that it needs to provide ONLY SQL Beggining with SELECT or WITH
        const reminderMsg = langTool.createUserMessage('The last response was not a valid SQL query. Please provide ONLY the SQL query text beginning with SELECT or WITH statement, without any explanations or additional text. Remember to prefix table names with schema name.');
        input.push(reminderMsg);
        addEngineMessage(reminderMsg);
      }
    }
    // If we exhausted iterations or didn't get a proper response
    throw new Error('Failed to generate SQL query: Maximum iterations reached or no valid SQL returned');
  } finally {
    abortSub.unsubscribe();
  }
}

/**
 * Context class that manages tool execution for SQL generation.
 * Supports multiple catalogs (e.g. for Databricks) with a default catalog.
 */
class SQLGenerationContext {
  private dbInfoMap: Map<string, BuiltinDBInfoMeta> = new Map();

  constructor(
    private connectionID: string,
    private defaultCatalog: string,
    dbInfos: BuiltinDBInfoMeta[],
    private connection: DG.DataConnection
  ) {
    for (const dbInfo of dbInfos)
      this.dbInfoMap.set(dbInfo.name, dbInfo);
  }

  private getDbInfo(catalogName?: string): BuiltinDBInfoMeta {
    const name = catalogName ?? this.defaultCatalog;
    const dbInfo = this.dbInfoMap.get(name);
    if (!dbInfo)
      throw new Error(`Catalog '${name}' not found. Available catalogs: ${Array.from(this.dbInfoMap.keys()).join(', ')}`);
    return dbInfo;
  }

  /**
   * Parses a table reference that may include catalog, schema, and table name.
   * Supports: catalog.schema.table, schema.table, or just table
   */
  private parseTableRef(tableRef: string): {catalog: string, schema: string, table: string} {
    const parts = tableRef.trim().split('.');
    if (parts.length >= 3)
      return {catalog: parts[0], schema: parts[1], table: parts.slice(2).join('.')};
    if (parts.length === 2)
      return {catalog: this.defaultCatalog, schema: parts[0], table: parts[1]};
    return {catalog: this.defaultCatalog, schema: '', table: parts[0]};
  }

  async listCatalogs(): Promise<string[]> {
    return Array.from(this.dbInfoMap.keys());
  }

  async listAllSchemas(catalogName?: string): Promise<string[]> {
    return (await this.getDbInfo(catalogName).getSchemas()).map((s) => s.name);
  }

  /**
   * Tool: list_tables
   * Returns all tables with their descriptions
   */
  async listTables(includeAllDescriptions: boolean, schemaName: string, catalogName?: string): Promise<{description: string, tableCount: number}> {
    const dbInfo = this.getDbInfo(catalogName);

    const schema = (await dbInfo.getSchemas()).find((s) => s.name === schemaName);
    if (!schema)
      return {description: `Schema ${schemaName} not found in catalog ${catalogName ?? this.defaultCatalog}.`, tableCount: 0};

    const tables = (await schema.getTables());
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
   * Returns detailed column information for specified tables.
   * Supports catalog.schema.table and schema.table formats.
   */
  async describeTables(tableNames: string[]): Promise<string> {
    const refs = tableNames.map((t) => this.parseTableRef(t));
    const descriptions: string[] = [];

    // Group by catalog+schema for efficient fetching
    const groupKey = (r: {catalog: string, schema: string}) => `${r.catalog}\0${r.schema}`;
    const schemaTableCache: Map<string, DG.TableInfo[]> = new Map();

    for (const key of new Set(refs.map(groupKey))) {
      const [catalog, schema] = key.split('\0');
      try {
        const dbInfo = this.getDbInfo(catalog);
        const schemaInfo = (await dbInfo.getSchemas()).find((s) => s.name === schema);
        if (schemaInfo)
          schemaTableCache.set(key, await schemaInfo.getTables());
      } catch { /* catalog not found */ }
    }

    for (const ref of refs) {
      const tables = schemaTableCache.get(groupKey(ref));
      if (!tables) {
        descriptions.push(`Schema ${ref.schema} not found in catalog ${ref.catalog}`);
        continue;
      }
      const table = tables.find((t) => t.friendlyName === ref.table || t.name === ref.table);
      if (!table) {
        descriptions.push(`Table ${ref.table} not found in ${ref.schema}`);
        continue;
      }
      descriptions.push(this.buildDetailedTableDescription(table));
    }

    return descriptions.join('\n\n---\n\n');
  }

  /**
   * Tool: list_joins
   * Returns all foreign key relationships involving the specified tables.
   * Supports catalog.schema.table and schema.table formats.
   */
  async listJoins(tableNames: string[]): Promise<string> {
    const refs = tableNames.map((t) => this.parseTableRef(t));

    const schemaQualifiedTableSet = new Set<string>();
    for (const ref of refs)
      schemaQualifiedTableSet.add(`${ref.schema}.${ref.table}`);

    // Collect unique catalog+schema pairs for efficient relation fetching
    const catalogSchemaMap = new Map<string, Set<string>>();
    for (const ref of refs) {
      if (!catalogSchemaMap.has(ref.catalog))
        catalogSchemaMap.set(ref.catalog, new Set());
      catalogSchemaMap.get(ref.catalog)!.add(ref.schema);
    }

    const allRelations: DG.DbRelationInfo[] = [];
    for (const [catalog, schemas] of catalogSchemaMap) {
      const dbInfo = this.getDbInfo(catalog);
      for (const schema of schemas) {
        const relations = await dbInfo.getRelationsForSchema(schema);
        allRelations.push(...relations);
      }
    }

    const filteredRelations = allRelations
      .filter((rel) => schemaQualifiedTableSet.has(`${rel.fromSchema}.${rel.fromTable}`) || schemaQualifiedTableSet.has(`${rel.toSchema}.${rel.toTable}`));
    const joins = filteredRelations.map((rel) => this.buildJoinDescription(rel));

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
