/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getSchemaDescriptor} from '@datagrok-libraries/db-explorer/src/utils';
import {LLMClient} from './LLM-client';
import {JsonSchema} from '../prompt-engine/interfaces';

type DBIDWithSchema = `${string}.${string}`;
type DBReference = `${string}.${string} -> ${string}.${string}`;
type SchemaInfo = {
    tableInfos: string[];
    referenceInfos: DBReference[];
    descriptor: Awaited<ReturnType<typeof getSchemaDescriptor>>;
}
export class DBSchemaInfo {
  private static _cache: Map<DBIDWithSchema, SchemaInfo> = new Map();

  static async getSchemaInfo(connectionID: string, schemaName: string): Promise<SchemaInfo> {
    const cacheKey: DBIDWithSchema = `${connectionID}.${schemaName}`;
    if (this._cache.has(cacheKey))
      return this._cache.get(cacheKey)!;
    const schemaDescriptor = await getSchemaDescriptor(connectionID, schemaName);
    const tableInfos: string[] = [];
    const referenceInfos: DBReference[] = [];
    for (const [tableName, tableInfo] of Object.entries(schemaDescriptor.tables))
      tableInfos.push(`${tableName}: (${Object.entries(tableInfo.columns).map(([colName, colType]) => `${colName} (${colType})`).join(', ')})`);
    for (const ref of schemaDescriptor.references)
      referenceInfos.push(`${ref.from} -> ${ref.to}`);
    const schemaInfo: SchemaInfo = {tableInfos, referenceInfos, descriptor: schemaDescriptor};
    this._cache.set(cacheKey, schemaInfo);
    return schemaInfo;
  }
}

export type CommentedEntity = {
    comment?: string;
    LLMComment?: string;
}
export type DBColumnMeta = CommentedEntity & {
    name: string;
    table: string;
    schema: string;
    type: string;
    semanticType?: string; // for example molecule, macromolecule, curve, etc.
    isUnique?: boolean;
    min?: number;
    max?: number;
    categoryValues?: string[]; // if less than 100, list them all
    uniqueValueCount?: number;
};

export type DBTableMeta = CommentedEntity & {
    name: string;
    schema: string;
    columns: DBColumnMeta[];
    domains?: string[]; // list of domain names that are related to this table
    rowCount: number;
}
export type DBSchemaMeta = CommentedEntity & {
    name: string;
    tables: DBTableMeta[];
};

export type DBConnectionMeta = CommentedEntity & {
    name: string;
    schemas: DBSchemaMeta[];
    relations: DBRelationMeta[];
};


export type DBRelationMeta = CommentedEntity & {
    fromSchema: string;
    fromTable: string;
    fromColumns: string[];
    toSchema: string;
    toTable: string;
    toColumns: string[];
    cardinality?: typeof DG.DbRelationInfo.prototype.cardinality;
   // comment         // from DB FK comment if exists
 // LLMComment      // "Each activity belongs to one assay; one assay has many activities"
    IsPrimaryPath?: boolean // recommended join vs weird legacy join
}

async function getTableLength(connection: DG.DataConnection, schemaName: string, tableName: string): Promise<number> {
  const sql = `--name: count-query \n
    SELECT COUNT(*) AS cnt FROM ${schemaName}.${tableName}`;
  const df = await connection.query('count-query', sql).apply({});
  return df.get('cnt', 0) as number;
}

async function get100TableSample(connection: DG.DataConnection, schemaName: string, tableName: string): Promise<DG.DataFrame> {
  const sql = `--name: sample-query \n
    SELECT * FROM ${schemaName}.${tableName} LIMIT 100`;
  return await connection.query('sample-query', sql).apply({});
}

async function getColumnMinMax(connection: DG.DataConnection, schemaName: string, tableName: string, columnName: string): Promise<{min: number; max: number}> {
  const sql = `--name: min-max-query \n
    SELECT MIN(${columnName}) AS min_val, MAX(${columnName}) AS max_val FROM ${schemaName}.${tableName}`;
  const df = await connection.query('min-max-query', sql).apply({});
  return {min: df.get('min_val', 0) as number, max: df.get('max_val', 0) as number};
}

async function getMax100CategoryValues(connection: DG.DataConnection, schemaName: string, tableName: string, columnName: string): Promise<string[] | undefined> {
  const sql = `--name: category-values-query \n
    SELECT DISTINCT ${columnName} AS catVals FROM ${schemaName}.${tableName} LIMIT 102`;
  const df: DG.DataFrame = await connection.query('category-values-query', sql).apply({});
  if (df.rowCount > 100) // too many unique values
    return undefined;
  const values = df.col('catVals')?.toList().map((a) => a?.toString() || 'null') as string[];
  if (values.length === 0 || values.some((v) => v.length > 100)) // there might be some big text values, skip those
    return undefined;
  return values;
}

async function getIsColumnUnique(connection: DG.DataConnection, schemaName: string, tableName: string, columnName: string, totalRowCount: number): Promise<boolean> {
  const sql = `--name: unique-count-query \n
    SELECT COUNT(DISTINCT ${columnName}) AS unique_count FROM ${schemaName}.${tableName}`;
  const df = await connection.query('unique-count-query', sql).apply({});
  const uniqueCount = df.get('unique_count', 0) as number;
  return uniqueCount === totalRowCount;
}


/**
 * Indexes the whole connection for LLM usage
 * Pass empty schemas array to index all schemas
 * @param connection
 * @param schemas
 */
export async function genDBConnectionMeta(connection: DG.DataConnection, schemas: string[]) {
  const out: DBConnectionMeta = {
    name: connection.name,
    schemas: [],
    relations: [],
  };

  const allSchemaNames = await grok.dapi.connections.getSchemas(connection);
  const schemaNames = schemas.length > 0 ? schemas : allSchemaNames;
  console.log(`Indexing schemas: ${schemaNames.join(', ')}`);
  for (const schemaName of schemaNames) {
    console.log(`Indexing schema: ${schemaName}`);
    console.log('Fetching schema descriptor...');
    const dbSchemaDescriptor = (await DBSchemaInfo.getSchemaInfo(connection.id, schemaName)).descriptor;
    // first, fill in the relations data that is already there
    const relations = out.relations;
    for (const ref of dbSchemaDescriptor.references) {
      const [fromTable, fromCol] = ref.from.split('.');
      const [toTable, toCol] = ref.to.split('.');
      relations.push({
        fromSchema: schemaName,
        fromTable,
        fromColumns: [fromCol],
        toSchema: schemaName,
        toTable,
        toColumns: [toCol],
      }); // we will handle cardinality later using LLM and querying the data
    }
    console.log(`Indexed ${relations.length} relations.`);


    const tables = Object.keys(dbSchemaDescriptor.tables);
    const schemaMeta: DBSchemaMeta = {
      name: schemaName,
      tables: [],
    };
    console.log(`Indexing ${tables.length} tables...`);
    for (const tableName of tables) {
      console.log(`Indexing table: ${tableName}`);
      // query for row count
      const rowCount = await getTableLength(connection, schemaName, tableName);
      const tableMeta: DBTableMeta = {
        name: tableName,
        schema: schemaName,
        columns: [],
        rowCount,
      };
      // get table sample, detect semtypes
      const tableSample = await get100TableSample(connection, schemaName, tableName);
      await tableSample.meta.detectSemanticTypes();
      console.log(`Indexing ${tableSample.columns.length} columns...`);
      for (const col of tableSample.columns) {
        const colMeta: DBColumnMeta = {
          name: col.name,
          table: tableName,
          schema: schemaName,
          type: dbSchemaDescriptor.tables[tableName]?.columns[col.name] || col.type, // should be always present
        };
        if (rowCount <= 1_000_000) {
          // get min/max for numeric columns
          if (col.isNumerical && col.type !== DG.COLUMN_TYPE.DATE_TIME) {
            const {min, max} = await getColumnMinMax(connection, schemaName, tableName, col.name);
            colMeta.min = min;
            colMeta.max = max;
          } else if (col.isCategorical && !dbSchemaDescriptor.tables[tableName]?.columns[col.name]?.toLowerCase().startsWith('timestamp')) { // we have a bug with timestamps
            const catValues = await getMax100CategoryValues(connection, schemaName, tableName, col.name);
            if (catValues)
              colMeta.categoryValues = catValues;
          }
          if (rowCount > 0) {
          // check if unique
            const isUnique = await getIsColumnUnique(connection, schemaName, tableName, col.name, rowCount);
            colMeta.isUnique = isUnique;
          }
        } else
          console.log(`Skipping detailed indexing for column ${col.name} in table ${tableName} due to large row count (${rowCount}).`);


        // get semantic type
        if (col.semType)
          colMeta.semanticType = col.semType;
        tableMeta.columns.push(colMeta);
      }
      schemaMeta.tables.push(tableMeta);
    }
    console.log(`Finished indexing schema: ${schemaName}`);
    out.schemas.push(schemaMeta);
  }

  //NB: Uncomment the following lines to enable LLM-based post-processing of DB metadata
  //const enrichedMeta = await postProcessDBConnectionMeta(out, 'gpt-4o-mini');
  //return enrichedMeta;
  return out;
}

// Helper schemas for structured LLM output
const ConnectionCommentSchema: JsonSchema = {
  type: 'object',
  properties: {
    comment: {type: 'string', description: 'Brief overview of the database purpose and domain'},
    LLMComment: {type: 'string', description: 'AI-friendly description of database structure and usage patterns'},
  },
  required: ['comment', 'LLMComment'],
  additionalProperties: false,
};

const SchemaCommentSchema: JsonSchema = {
  type: 'object',
  properties: {
    comment: {type: 'string', description: 'Brief description of schema purpose'},
    LLMComment: {type: 'string', description: 'AI-friendly description for query generation'},
  },
  required: ['comment', 'LLMComment'],
  additionalProperties: false,
};

const TableCommentSchema: JsonSchema = {
  type: 'object',
  properties: {
    comment: {type: 'string', description: 'Brief table description'},
    LLMComment: {type: 'string', description: 'AI-friendly description with domain context'},
  },
  required: ['comment', 'LLMComment'],
  additionalProperties: false,
};

const ColumnCommentsSchema: JsonSchema = {
  type: 'object',
  properties: {
    columns: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: {type: 'string'},
          comment: {type: 'string', description: 'Brief column description or empty if self-evident'},
          LLMComment: {type: 'string', description: 'AI-friendly description or empty if self-evident'},
        },
        required: ['name', 'comment', 'LLMComment'],
        additionalProperties: false,
      },
    },
  },
  required: ['columns'],
  additionalProperties: false,
};

const RelationCommentSchema: JsonSchema = {
  type: 'object',
  properties: {
    comment: {type: 'string', description: 'Brief relation description'},
    LLMComment: {type: 'string', description: 'Natural language explanation (e.g., "Each activity belongs to one assay")'},
    cardinality: {type: 'string', enum: ['one-to-one', 'one-to-many', 'many-to-one', 'many-to-many']},
    IsPrimaryPath: {type: 'boolean', description: 'Is this a recommended join vs legacy join'},
  },
  required: ['comment', 'LLMComment', 'cardinality', 'IsPrimaryPath'],
  additionalProperties: false,
};

const HiddenRelationsSchema: JsonSchema = {
  type: 'object',
  properties: {
    hiddenRelations: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          fromSchema: {type: 'string'},
          fromTable: {type: 'string'},
          fromColumns: {type: 'array', items: {type: 'string'}},
          toSchema: {type: 'string'},
          toTable: {type: 'string'},
          toColumns: {type: 'array', items: {type: 'string'}},
          comment: {type: 'string'},
          LLMComment: {type: 'string'},
          cardinality: {type: 'string', enum: ['one-to-one', 'one-to-many', 'many-to-one', 'many-to-many']},
          IsPrimaryPath: {type: 'boolean'},
        },
        required: ['fromSchema', 'fromTable', 'fromColumns', 'toSchema', 'toTable', 'toColumns', 'comment', 'LLMComment', 'cardinality'],
        additionalProperties: false,
      },
    },
  },
  required: ['hiddenRelations'],
  additionalProperties: false,
};

/**
 * Post-processes DB metadata by enriching it with LLM-generated comments and insights
 * @param connectionMeta - The raw connection metadata
 * @param openAIModel - The OpenAI model to use (e.g., 'gpt-4o')
 */
export async function postProcessDBConnectionMeta(connectionMeta: DBConnectionMeta, openAIModel: string): Promise<DBConnectionMeta> {
  const llmClient = LLMClient.getInstance();
  const result = {...connectionMeta};

  // Step 1: Generate connection-level comment
  console.log('Generating connection-level comment...');
  const connectionOverview = generateConnectionOverview(result);
  const connectionCommentResponse = await llmClient.generalPromptCached(
    openAIModel,
    'You are a database documentation expert. Generate concise, informative descriptions for database connections.',
    `Database: ${result.name}\n\n${connectionOverview}\n\nGenerate a brief overview (comment) and an AI-friendly description (LLMComment) that helps AI models understand this database's purpose and structure.`,
    ConnectionCommentSchema
  );
  const connectionComment = JSON.parse(connectionCommentResponse);
  result.comment = connectionComment.comment;
  result.LLMComment = connectionComment.LLMComment;

  // Step 2: Process each schema
  for (const schema of result.schemas) {
    console.log(`Processing schema: ${schema.name}`);

    // Generate schema comment
    const schemaOverview = generateSchemaOverview(schema);
    const schemaCommentResponse = await llmClient.generalPromptCached(
      openAIModel,
      'You are a database documentation expert. Generate concise schema descriptions.',
      `Schema: ${schema.name}\n\n${schemaOverview}\n\nGenerate a brief description and AI-friendly context.`,
      SchemaCommentSchema
    );
    const schemaComment = JSON.parse(schemaCommentResponse);
    schema.comment = schemaComment.comment;
    schema.LLMComment = schemaComment.LLMComment;

    // Step 3: Process tables in batches (to provide context)
    for (const table of schema.tables) {
      console.log(`Processing table: ${table.name}`);

      // Generate table comment with full context
      const tableContext = generateTableContext(table, schema);
      const tableCommentResponse = await llmClient.generalPromptCached(
        openAIModel,
        'You are a database documentation expert. Generate table descriptions that help AI understand the domain.',
        `Table: ${schema.name}.${table.name}\n\n${tableContext}\n\nGenerate a brief description and AI-friendly context explaining this table's role.`,
        TableCommentSchema
      );
      const tableComment = JSON.parse(tableCommentResponse);
      table.comment = tableComment.comment;
      table.LLMComment = tableComment.LLMComment;

      // Step 4: Generate column comments (batch processing for efficiency)
      const columnContext = generateColumnContext(table);
      const columnCommentsResponse = await llmClient.generalPromptCached(
        openAIModel,
        `You are a database documentation expert. For each column, generate a brief description. 
If a column is self-evident (like 'id' in a table named 'users' -> 'User ID'), you can use a very short description or leave empty strings.
Focus on non-obvious columns, semantic types, and domain-specific meanings.`,
        `Table: ${schema.name}.${table.name}\n\n${columnContext}\n\nFor each column, provide comment and LLMComment. Leave both empty if the column is completely self-evident.`,
        ColumnCommentsSchema
      );
      const columnComments = JSON.parse(columnCommentsResponse);

      // Apply column comments
      for (const colComment of columnComments.columns) {
        const col = table.columns.find((c) => c.name === colComment.name);
        if (col) {
          col.comment = colComment.comment || undefined;
          col.LLMComment = colComment.LLMComment || undefined;
        }
      }
    }
  }

  // Step 5: Process relations
  console.log('Processing relations...');
  const relationsContext = generateRelationsContext(result);
  for (let i = 0; i < result.relations.length; i++) {
    const relation = result.relations[i];
    const relationContext = `${relation.fromSchema}.${relation.fromTable}(${relation.fromColumns.join(', ')}) -> ${relation.toSchema}.${relation.toTable}(${relation.toColumns.join(', ')})`;

    const relationCommentResponse = await llmClient.generalPromptCached(
      openAIModel,
      `You are a database expert. Analyze foreign key relationships and provide:
1. Brief comment
2. Natural language explanation (e.g., "Each assay_result belongs to one assay; one assay has many assay_results")
3. Cardinality (one-to-one, one-to-many, many-to-one, many-to-many)
4. Whether this is a primary/recommended join path (true) or a legacy/unusual join (false)`,
      `Relation: ${relationContext}\n\nDatabase context:\n${relationsContext}\n\nAnalyze this relationship.`,
      RelationCommentSchema
    );
    const relationComment = JSON.parse(relationCommentResponse);
    relation.comment = relationComment.comment;
    relation.LLMComment = relationComment.LLMComment;
    relation.cardinality = relationComment.cardinality;
    relation.IsPrimaryPath = relationComment.IsPrimaryPath;
  }

  // Step 6: Discover hidden relations
  console.log('Discovering hidden relations...');
  const hiddenRelationsResponse = await llmClient.generalPromptCached(
    openAIModel,
    `You are a database expert. Analyze the database schema and identify potential "hidden" relationships that are not defined as foreign keys but could be logically joined.
Examples:
- Tables with similar column names (e.g., user_id in different tables)
- Semantic relationships (e.g., name fields that could match)
- Domain-specific relationships

Only suggest high-confidence relationships that would be useful for queries. Return empty array if no hidden relations are found.`,
    `Database: ${result.name}\n\nSchemas and tables:\n${generateHiddenRelationsContext(result)}\n\nIdentify any hidden relationships.`,
    HiddenRelationsSchema
  );
  const hiddenRelations = JSON.parse(hiddenRelationsResponse);

  // Add hidden relations to the result
  if (hiddenRelations.hiddenRelations && hiddenRelations.hiddenRelations.length > 0) {
    result.relations.push(...hiddenRelations.hiddenRelations);
    console.log(`Discovered ${hiddenRelations.hiddenRelations.length} hidden relations`);
  }

  return result;
}

// Helper functions to generate context strings

function generateConnectionOverview(conn: DBConnectionMeta): string {
  const schemaCount = conn.schemas.length;
  const tableCount = conn.schemas.reduce((sum, s) => sum + s.tables.length, 0);
  const relationCount = conn.relations.length;

  return `Schemas: ${schemaCount}, Tables: ${tableCount}, Relations: ${relationCount}
Schema names: ${conn.schemas.map((s) => s.name).join(', ')}`;
}

function generateSchemaOverview(schema: DBSchemaMeta): string {
  const tableNames = schema.tables.map((t) => t.name).join(', ');
  return `Tables (${schema.tables.length}): ${tableNames}`;
}

function generateTableContext(table: DBTableMeta, schema: DBSchemaMeta): string {
  const colNames = table.columns.map((c) => `${c.name} (${c.type})`).join(', ');
  const relatedTables = schema.tables
    .filter((t) => t.name !== table.name)
    .map((t) => t.name)
    .join(', ');

  return `Row count: ${table.rowCount}
Columns: ${colNames}
Other tables in schema: ${relatedTables}`;
}

function generateColumnContext(table: DBTableMeta): string {
  const columnDetails = table.columns.map((col) => {
    const parts = [`${col.name}: ${col.type}`];
    if (col.semanticType) parts.push(`semantic: ${col.semanticType}`);
    if (col.isUnique) parts.push('unique');
    if (col.min !== undefined) parts.push(`range: ${col.min}-${col.max}`);
    if (col.categoryValues) parts.push(`categories: ${col.categoryValues.slice(0, 5).join(', ')}${col.categoryValues.length > 5 ? '...' : ''}`);
    return parts.join(', ');
  }).join('\n');

  return `Table: ${table.name} (${table.rowCount} rows)\nColumns:\n${columnDetails}`;
}

function generateRelationsContext(conn: DBConnectionMeta): string {
  // Provide overview of all tables and their columns for context
  const tablesOverview = conn.schemas.flatMap((schema) =>
    schema.tables.map((table) =>
      `${schema.name}.${table.name}: ${table.columns.map((c) => `${c.name}(${c.type})`).join(', ')}`
    )
  ).join('\n');

  return `Tables in database:\n${tablesOverview}`;
}

function generateHiddenRelationsContext(conn: DBConnectionMeta): string {
  // Generate a compact representation of all tables and columns
  const context = conn.schemas.flatMap((schema) =>
    schema.tables.map((table) => {
      const cols = table.columns.map((c) => {
        const parts = [c.name, c.type];
        if (c.semanticType) parts.push(`sem:${c.semanticType}`);
        if (c.isUnique) parts.push('unique');
        return parts.join(':');
      }).join(', ');
      return `${schema.name}.${table.name} [${cols}]`;
    })
  ).join('\n');

  // Also include existing relations for reference
  const existingRels = conn.relations.map((r) =>
    `${r.fromSchema}.${r.fromTable}(${r.fromColumns.join(',')}) -> ${r.toSchema}.${r.toTable}(${r.toColumns.join(',')})`
  ).join('\n');

  return `${context}\n\nExisting relations:\n${existingRels}`;
}


/**
 * ############################################################################################################################################
 * Up until now is the code for meta generation and post-processing using LLMs
 * ###########################################################################################################################################
 * Bellow is tool to upload these annotations to the sticky meta :3
 */

export async function moveDBMetaToStickyMetaOhCoolItEvenRhymes(connectionMeta: DBConnectionMeta) {
  const connection = await grok.dapi.connections.filter(`name="${connectionMeta.name}"`).first();
  if (!connection)
    throw new Error(`Connection ${connectionMeta.name} not found`);
  const dbInfo = await grok.data.db.getInfo(connection);

  await dbInfo.clearProperties();
  console.log(`Cleared existing metadata for connection ${connectionMeta.name}`);
  // connection-level comments
  if (connectionMeta.comment) {
    await failable('Setting connection comment', async () => {
      await dbInfo.setComment(connectionMeta.comment!);
    });
  }
  if (connectionMeta.LLMComment) {
    await failable('Setting connection LLMComment', async () => {
      await dbInfo.setLlmComment(connectionMeta.LLMComment!);
    });
  }
  // relation-level annotations
  const relations = connectionMeta.relations ?? [];

  for (const rel of relations) {
    await failable(`adding relation ${rel.fromSchema}.${rel.fromTable} -> ${rel.toSchema}.${rel.toTable}`, async () => {
      await dbInfo.addRelation(rel.fromTable, rel.fromColumns, rel.toTable, rel.toColumns, {
        fromSchema: rel.fromSchema,
        toSchema: rel.toSchema,
        comment: rel.comment,
        llmComment: rel.LLMComment,
        cardinality: rel.cardinality,
        isPrimaryPath: rel.IsPrimaryPath,
      });
    });
  }
  const schemas = await dbInfo.getSchemas();
  for (const schemaMeta of connectionMeta.schemas) {
    const dbSchema = schemas.find((s) => s.name === schemaMeta.name);
    if (!dbSchema) {
      console.warn(`Schema ${schemaMeta.name} not found in connection ${connectionMeta.name}, skipping`);
      continue;
    }
    // schema-level comments
    if (schemaMeta.comment) {
      await failable(`Setting schema ${schemaMeta.name} comment`, async () => {
        await dbSchema.setComment(schemaMeta.comment!);
      });
    }
    if (schemaMeta.LLMComment) {
      await failable(`Setting schema ${schemaMeta.name} LLMComment`, async () => {
        await dbSchema.setLlmComment(schemaMeta.LLMComment!);
      });
    }
    // moving to tables
    const tables = await dbSchema.getTables();
    for (const tableMeta of schemaMeta.tables) {
      const dbTable = tables.find((t) => t.friendlyName === tableMeta.name || t.name === tableMeta.name);
      if (!dbTable) {
        console.warn(`Table ${tableMeta.name} not found in schema ${schemaMeta.name}, skipping`);
        continue;
      }
      // table-level comments

      await failable(`Setting table ${schemaMeta.name}.${tableMeta.name} annotations`, async () => {
        await dbSchema.annotateTable(dbTable, removeUndefinesFromObject({
          comment: tableMeta.comment,
          llmComment: tableMeta.LLMComment,
          domains: tableMeta.domains?.join(',') ?? undefined,
          rowCount: tableMeta.rowCount ?? undefined,
        }));
      });
      // move to columns
      const columns = dbTable.columns;
      for (const columnMeta of tableMeta.columns) {
        const dbColumn = columns.find((c) => c.name === columnMeta.name);
        if (!dbColumn) {
          console.warn(`Column ${columnMeta.name} not found in table ${tableMeta.name}, skipping`);
          continue;
        }
        await failable(`Setting column ${schemaMeta.name}.${tableMeta.name}.${columnMeta.name} annotations`, async () => {
          await dbSchema.annotateColumn(dbTable, dbColumn, removeUndefinesFromObject({
            comment: columnMeta.comment ?? undefined,
            llmComment: columnMeta.LLMComment ?? undefined,
            isUnique: columnMeta.isUnique ?? undefined,
            sampleValues: columnMeta.categoryValues?.slice(0, 10) ?? undefined,
            values: columnMeta.categoryValues ?? undefined,
            max: columnMeta.max ?? undefined,
            min: columnMeta.min ?? undefined,
            uniqueCount: columnMeta.uniqueValueCount ?? undefined,
          }));
        });
      }
    }
  }
  console.log(`Finished moving metadata to sticky meta for connection ${connectionMeta.name}`);
}

async function failable(desc: string, f: () => Promise<void>) {
  try {
    console.log(desc);
    await f();
    console.log(`${desc} - done`);
    return;
  } catch (e) {
    console.error('Failable function error:', e);
    return null;
  }
}

function removeUndefinesFromObject<T>(obj: T): T {
  const keys = Object.keys(obj as Object) as (keyof T)[];
  for (const key of keys) {
    if (obj[key] == undefined)
      delete obj[key];
  }
  return obj;
}
