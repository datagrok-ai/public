/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import OpenAI from 'openai';
import {LLMCredsManager} from './creds';
import {DBSchemaInfo, DBConnectionMeta, DBTableMeta, DBRelationMeta, DBSchemaMeta} from './db-index-tools';
import {chemblIndex} from './indexes/chembl-index';
import {biologicsIndex} from './indexes/biologics-index';
import {getAIAbortSubscription} from '../utils';
import * as rxjs from 'rxjs';
import {getTopKSimilarQueries, getVectorEmbedding} from './embeddings';
import {SemValueObjectHandler} from '@datagrok-libraries/db-explorer/src/object-handlers';
import {ChatModel} from 'openai/resources/index';
import {UIMessageOptions} from './panel';

type AIPanelFuncs = {
  addUserMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam, msg: string) => void,
  addAIMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam, title: string, msg: string) => void,
  addEngineMessage: (aiMsg: OpenAI.Chat.ChatCompletionMessageParam) => void, // one that is not shown in the UI
  addUiMessage: (msg: string, fromUser: boolean, messageOptions?: UIMessageOptions) => void
}

const suspiciousSQlPatterns = ['DROP ', 'DELETE ', 'UPDATE ', 'INSERT ', 'ALTER ', 'CREATE ', 'TRUNCATE ', 'EXEC ', 'MERGE '];
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
    oldMessages?: OpenAI.Chat.ChatCompletionMessageParam[]
    aiPanel?: AIPanelFuncs,
    modelName?: ChatModel
  } = {}
): Promise<string> {
  let aborted = false;
  let dbMeta: DBConnectionMeta | undefined = undefined;
  options.aiPanel?.addUiMessage(prompt, true);

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
  const abortSub = getAIAbortSubscription().subscribe(() => {
    aborted = true;
    console.log('Aborting SQL generation as per user request');
    try {
      SQLGenerationContext._lastFc?.cancel();
      SQLGenerationContext._lastFc = null;
    } catch (_) {
    }
    abortSub.unsubscribe();
  });
  try {
  // Initialize OpenAI client
    const openai = new OpenAI({
      apiKey: LLMCredsManager.getApiKey(),
      dangerouslyAllowBrowser: true,
    });


    // Initial system message
    const systemMessage = `You are an expert SQL query generator. You have access to tools to explore a database schema and generate SQL queries.

    Your name is GrokBerg Databazawsky.

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
    - Always prefix table names with schema name in SQL
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


    const messages: OpenAI.Chat.ChatCompletionMessageParam[] = options.oldMessages ? options.oldMessages.slice() : [];
    if (messages.length === 0) {
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

      if ((dbQueryEmbeddings?.length ?? 0) > 0) {
        const similarQueries = await findSimilarQueriesToPrompt(prompt);
        if (similarQueries.n > 0 && similarQueries.text) {
          semTypeWithinPromptInfo += `\n\n Additionally, here are some similar previously executed queries that might help you:
          \n${similarQueries.text}
          `;
        }
      }

      const systemMsg: OpenAI.Chat.ChatCompletionMessageParam = {role: 'system', content: systemMessage};
      messages.push(systemMsg);
      options.aiPanel?.addEngineMessage(systemMsg);
      const initialUserMsg: OpenAI.Chat.ChatCompletionMessageParam = {
        role: 'user',
        content: `
        Available schemas: ${(await context.listAllSchemas()).join(', ')}\n\n
        Available tables in schema ${schemaName}:\n\n${initialTableList}\n\nUser Query: ${prompt}\n\n
        Explore the schema using the available tools and generate an SQL query to answer this question.
        ${semTypeWithinPromptInfo}`,
      };
      options.aiPanel?.addEngineMessage(initialUserMsg); // prompt will be shown in UI
      messages.push(initialUserMsg);
    } else {
      const followUpUserMsg: OpenAI.Chat.ChatCompletionMessageParam = {role: 'user', content: `User Follow up (do modifications as/if needed based on the following): ${prompt}`};
      messages.push(followUpUserMsg);
      options.aiPanel?.addEngineMessage(followUpUserMsg);
    }

    // Define available tools
    const tools: OpenAI.Chat.ChatCompletionTool[] = [
      {
        type: 'function',
        function: {
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
          },
        },
      },
      {
        type: 'function',
        function: {
          name: 'list_all_schemas',
          description: 'Retyurns the names all all schemas in the connection.',
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
                description: 'List of table names (WITH schema prefix) to describe',
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
                description: 'List of table names (WITH schema prefix) to find joins for',
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
            required: ['sql'],
          },
        },
      },
      // plain function tool for general purpose question answering
      {
        type: 'function',
        function: {
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
          }
        }
      }
    ];

    // if we have some embeddings, then also add a tool so that LLM can use it to find similar queries
    if ((dbQueryEmbeddings?.length ?? 0) > 0) {
      tools.push({
        type: 'function',
        function: {
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
          },
        }
      });
    }

    async function findSimilarQueriesToPrompt(aiPrompt: string) {
      if ((dbQueryEmbeddings?.length ?? 0) === 0)
        return {text: 'No queries found', n: 0};
      if ((dbQueryEmbeddings!.length < 4))
        return {text: `Only ${dbQueryEmbeddings.length} queries found:\n\n${dbQueryEmbeddings!.map((q) => `- ${q.query}`).join('\n\n')}`, n: dbQueryEmbeddings.length};
      try {
        const userQueryEmbedding = await getVectorEmbedding(openai, aiPrompt);
        const topSimilarQueries = getTopKSimilarQueries(userQueryEmbedding, dbQueryEmbeddings!, 3);
        console.log('Top similar queries:', topSimilarQueries);
        if (topSimilarQueries.length === 0)
          return {text: 'No similar queries found', n: 0};
        return {text: `Here are some similar previously executed queries that might help you.
          Note that these queries might have annotations/comments like name descriptions and so on that can help you understand the schema better.\n\n
          Here they are:
          ${topSimilarQueries.map((q) => `- Query:\n ${q.query} \n (similarity score: ${q.score.toFixed(4)})`).join('\n\n')}`, n: topSimilarQueries.length};
      } catch (error) {
        console.error('Error generating embedding for similar query search:', error);
        return {text: 'Error generating embedding for similar query search', n: 0};
      }
    };

    // Function calling loop
    let iterations = 0;
    const maxIterations = 15; // Prevent infinite loops

    while (iterations < maxIterations) {
      if (aborted) {
        grok.shell.info('SQL generation aborted by user.');
        return '';
      }
      iterations++;
      console.log(`\n=== Iteration ${iterations} ===`);

      const response = await openai.chat.completions.create({
        //   model: modelName,
        //   temperature: 0,
        model: options?.modelName ?? 'gpt-4o-mini',
        temperature: 1,
        reasoning_effort: options?.modelName?.startsWith('gpt-5') ? 'high' : undefined,
        messages,
        tools,
      });

      const choice = response.choices[0];
      if (!choice) break;

      const message = choice.message;
      messages.push(message);
      options.aiPanel?.addEngineMessage(message);
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
            case 'list_all_schemas':
              options.aiPanel?.addUiMessage(`Listing all schemas in the database.`, false);
              result = (await context.listAllSchemas()).join(', ');
              options.aiPanel?.addUiMessage(`Schemas found: *${result}*`, false);
              break;
            case 'reply_to_user':
              options.aiPanel?.addUiMessage(functionArgs.reply ?? '', false);
              result = 'Reply sent to user.';
              break;
            case 'list_tables_in_schema':
              options.aiPanel?.addUiMessage(`Listing all tables in schema *${functionArgs.schemaName}*.`, false);
              const res = await context.listTables(false, functionArgs.schemaName);
              result = res.description;
              options.aiPanel?.addUiMessage(res.tableCount ? `Found ${res.tableCount} tables` : 'No tables found', false);
              break;
            case 'describe_tables':
              options.aiPanel?.addUiMessage(`Getting content of following tables: *${functionArgs.tables.join(', ')}*.`, false);
              result = await context.describeTables(functionArgs.tables);
              break;
            case 'list_joins':
              options.aiPanel?.addUiMessage(`Listing relations for tables: *${functionArgs.tables.join(', ')}*.`, false);
              result = await context.listJoins(functionArgs.tables);
              break;
            case 'try_sql':
              options.aiPanel?.addUiMessage(`Testing SQL query:\n${functionArgs.description}\n\`\`\`sql\n${functionArgs.sql}\n\`\`\``, false);
              result = await context.trySql(functionArgs.sql, functionArgs.description);
              options.aiPanel?.addUiMessage(`SQL Test Result: \n\n${result}`, false);
              break;
            case 'find_similar_queries':
              options.aiPanel?.addUiMessage(`Searching for similar queries to help with SQL generation.`, false);
              const similarQueries = await findSimilarQueriesToPrompt(functionArgs.prompt);
              options.aiPanel?.addUiMessage(similarQueries.n > 0 ? `Found ${similarQueries.n} similar queries.` : similarQueries.text, false);
              result = similarQueries.text;
              break;
            default:
              result = `Error: Unknown function ${functionName}`;
            }
          } catch (error: any) {
            result = `Error executing ${functionName}: ${error.message}`;
          }

          console.log(`Result: ${result.substring(0, 200)}${result.length > 200 ? '...' : ''}`);

          // Add function result to messages
          const msg = {
            role: 'tool',
            tool_call_id: toolCall.id,
            content: result,
          } as const;
          messages.push(msg);
          options.aiPanel?.addEngineMessage({...msg});
        }
      } else if (message.content) {
      // Model provided a text response (hopefully the final SQL)
        console.log('Model response:', message.content);

        // Check if this looks like SQL (basic heuristic)
        const content = message.content.trim();
        // Clean up the response
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
              options.aiPanel?.addUiMessage(`Final SQL Query:\n\`\`\`sql\n${out}\n\`\`\``, false, {result: {finalResult: out}});
            return out;
          }
          if (res)
            options.aiPanel?.addUiMessage(`Final SQL Query:\n\`\`\`sql\n${res}\n\`\`\``, false, {result: {finalResult: res}});
          return res;
        } else {
          const replyText = contentUpperCase.startsWith('REPLY ONLY:') ? content.substring('REPLY ONLY:'.length).trim() : content;
          options.aiPanel?.addUiMessage(replyText, false);
          return '';
        }

        // Model is explaining something, let it continue
        // if (choice.finish_reason === 'stop')
        //   break;
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
  } finally {
    abortSub.unsubscribe();
  }
}

/**
 * Context class that manages tool execution for SQL generation
 */
class SQLGenerationContext {
  private schemas?: DBSchemaMeta[];

  constructor(
    private connectionID: string,
    private schemaName: string,
    private dbMeta: DBConnectionMeta | undefined,
    private connection: DG.DataConnection
  ) {
    if (dbMeta && (dbMeta.schemas?.length ?? 0) > 0)
      this.schemas = dbMeta.schemas;
  }

  private static schemasCache: {[connectionId: string]: string[]} = {};

  async listAllSchemas(): Promise<string[]> {
    if (this.schemas)
      return this.schemas.map((s) => s.name);
    if (SQLGenerationContext.schemasCache[this.connectionID])
      return SQLGenerationContext.schemasCache[this.connectionID];
    let out: string[] = [];
    try {
      const connection = await grok.dapi.connections.find(this.connectionID);
      const schemas = await grok.dapi.connections.getSchemas(connection);
      out = schemas;
    } catch (e) {
      console.log('Error fetching schemas:', e);
    }
    SQLGenerationContext.schemasCache[this.connectionID] = out;
    return out;
  }

  /**
   * Tool: list_tables
   * Returns all tables with their descriptions
   */
  async listTables(includeAllDescriptions: boolean, schemaName: string): Promise<{description: string, tableCount: number}> {
    schemaName ??= this.schemaName;
    if (this.schemas && this.dbMeta) {
      const schema = this.schemas.find((s) => s.name === schemaName);
      if (!schema)
        return {description: `Schema ${schemaName} not found in metadata.`, tableCount: 0};
      const tables = schema.tables.map((table) => {
        const parts = [`${table.name}`];
        if (includeAllDescriptions) {
          if (table.LLMComment)
            parts.push(`  Description: ${table.LLMComment.substring(0, 200)}...`);
          else if (table.comment)
            parts.push(`  Description: ${table.comment.substring(0, 200)}...`);
        }
        parts.push(`  Columns: ${table.columns.length}, Rows: ${table.rowCount}`);
        return parts.join('\n');
      });
      return {description: `Tables in ${this.schemaName}:\n\n${tables.join('\n')}`, tableCount: tables.length};
    }

    // Fallback to schema descriptor
    const schemaInfo = await DBSchemaInfo.getSchemaInfo(this.connectionID, schemaName);
    return {description: `Tables in ${this.schemaName}:\n\n${schemaInfo.tableInfos.join('\n')}`, tableCount: schemaInfo.tableInfos.length};
  }

  /**
   * Tool: describe_tables
   * Returns detailed column information for specified tables
   */
  async describeTables(tableNames: string[]): Promise<string> {
    const schemaNames = tableNames.map((t) => t.trim().indexOf('.') >= 0 ? t.split('.')[0] : this.schemaName);
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    const descriptions: string[] = [];

    for (let i = 0; i < tableNames.length; i++) {
      const tableName = tableNames[i];
      const schemaName = schemaNames[i];
      if (this.schemas) {
        const schema = this.schemas.find((s) => s.name === schemaName);
        if (!schema) {
          descriptions.push(`Schema ${schemaName} not found`);
          continue;
        }
        const table = schema.tables.find((t) => t.name === tableName);
        if (!table) {
          descriptions.push(`Table ${tableName} not found`);
          continue;
        }

        const tableDesc = this.buildDetailedTableDescription(table);
        descriptions.push(tableDesc);
      } else {
        // Fallback: basic info from schema descriptor
        const schemaInfo = await DBSchemaInfo.getSchemaInfo(this.connectionID, schemaNames[i]);
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
    const schemaNames = tableNames.map((t) => t.trim().indexOf('.') >= 0 ? t.split('.')[0] : this.schemaName);
    tableNames = tableNames.map((t) => t.trim()).map((t) => t.indexOf('.') >= 0 ? t.split('.')[1] : t); // Remove schema prefix if present
    // in case of cross-schema requests,
    const schemaQualifiedTableSet = new Set<string>();
    for (let i = 0; i < tableNames.length; i++)
      schemaQualifiedTableSet.add(`${schemaNames[i]}.${tableNames[i]}`);

    const joins: string[] = [];
    if (this.dbMeta) {
      for (let i = 0; i < this.dbMeta.relations.length; i++) {
        const relation = this.dbMeta.relations[i];
        const fromFullName = `${relation.fromSchema}.${relation.fromTable}`;
        const toFullName = `${relation.toSchema}.${relation.toTable}`;
        if (!schemaQualifiedTableSet.has(fromFullName) && !schemaQualifiedTableSet.has(toFullName))
          continue;
        const joinDesc = this.buildJoinDescription(relation);
        joins.push(joinDesc);
      }
    } else {
      // Fallback to schema descriptor
      const applicableSchemaNames = Array.from(new Set(schemaNames));
      const applicableSchemas = await Promise.all(applicableSchemaNames.map((sn) => DBSchemaInfo.getSchemaInfo(this.connectionID, sn)));
      for (let i = 0; i < applicableSchemas.length; i++) {
        const schemaInfo = applicableSchemas[i];
        const schemaName = applicableSchemaNames[i];
        for (const ref of schemaInfo.referenceInfos) {
          const [from, to] = ref.split('->').map((s) => s.trim());
          const fromTable = from.split('.')[0];
          const toTable = to.split('.')[0];
          const qualifiedFrom = `${schemaName}.${fromTable}`;
          const qualifiedTo = `${schemaName}.${toTable}`;

          if (schemaQualifiedTableSet.has(qualifiedFrom) || schemaQualifiedTableSet.has(qualifiedTo))
            joins.push(ref);
        }
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
    } catch (error: any) {
      SQLGenerationContext._lastFc = null;
      return `SQL Error: ${typeof error == 'string' ? error : error?.message}\n\nThis query failed to execute. Please revise the SQL based on the schema information.`;
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
        const values = col.categoryValues;
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
