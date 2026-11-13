/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {QueryJoinOptions} from './types';

export type ValueQuery = {
  sql: string;
  /** Joins that apply to the queried table */
  joins: QueryJoinOptions[];
  /** Table and schema of every joined column, by its name in the result */
  columnTables: {[colName: string]: string};
  columnSchemas: {[colName: string]: string};
};

/**
 * SQL that reads the rows of a table whose [match] column equals [matchValue], with the joins configured for it.
 * Pass the @match and @matchValue as empty strings to get 1 random row (for example queries)
 */
export function buildValueQuery(
  dataSource: string | undefined,
  tableName: string,
  match: string,
  matchValue: string | number,
  schemaName: string,
  joinOptions: QueryJoinOptions[] = [],
  max1: boolean = false,
  catalogName?: string
): ValueQuery {
  const matchValueStr = typeof matchValue === 'string' ? `'${matchValue.replace(/'/g, '\'\'')}'` : matchValue;
  const startCharCode = 98; // ascii b
  const [nbs, nbe] = getNameBrackets(dataSource);
  const tableRef = (schema: string, table: string) =>
    `${catalogName ? `${nbs}${catalogName}${nbe}.` : ''}${nbs}${schema}${nbe}.${nbs}${table}${nbe}`;

  const joins = joinOptions.filter((opt) => opt.fromTable.toLowerCase() === tableName.toLowerCase() && (!opt.fromSchema || opt.fromSchema.toLowerCase() === schemaName.toLowerCase()));

  const tableAliases = joins.map((_, i) => String.fromCharCode(startCharCode + i));
  const columnTables: {[colName: string]: string} = {};
  const columnSchemas: {[colName: string]: string} = {};
  const otherCols =
    joins
      .map((opt, i) => {
        const alias = tableAliases[i];
        return opt.select.map((col) => {
          // if the col is confugured with alias, make sure it has quotes
          // split cace insensitive way
          const parts = col.split(/ as /i);
          if (parts.length > 1) {
            columnTables[parts[1].trim()] = opt.tableName;
            columnSchemas[parts[1].trim()] = opt.onSchema ?? schemaName;
            return `${alias}.${nbs}${parts[0].trim()}${nbe} as ${nbs}${parts[1].trim()}${nbe}`;
          }
          columnTables[col] = opt.tableName;
          columnSchemas[col] = opt.onSchema ?? schemaName;
          return `${alias}.${nbs}${col}${nbe}`;
        }).join(', ');
      })
      .join(', ');
  const joinStr = joins
    .map((opt, i) => {
      const alias = tableAliases[i];
      return `left join ${tableRef(opt.onSchema ?? schemaName, opt.tableName)} ${alias} on a.${nbs}${opt.columnName}${nbe} = ${alias}.${nbs}${opt.onColumn}${nbe}`;
    })
    .join(' ');
  const isExampleQuery = match === '' && matchValue === ''; // used for setup editor, will not happen otherwise if not intentionally
  const oneRow = max1 || isExampleQuery;
  const topClause = dataSource === DG.DataSourceType.MsSql;
  const selected = `a.* ${(joins.length > 0 ? ',' : '')} ${otherCols} `;

  const sql = `
        --name: getDBValueInfo
        --output: dataframe result
        select ${oneRow && topClause ? 'top 1 ' : ''}${selected}
        from ${tableRef(schemaName, tableName)} a
        ${joinStr}
        ${isExampleQuery ? '' : `where a.${match} = ${matchValueStr}`}
        ${oneRow && !topClause ? 'limit 1' : ''}
    `;
  return {sql, joins, columnTables, columnSchemas};
}

/**
 *
 * @param connection
 * @param tableName
 * @param match - Pass the @match and @matchValue as empty strings to get 1 random row (for example queries)
 * @param matchValue
 * @param schemaName
 * @param joinOptions
 * @param max1
 * @param catalogName - catalog of the schema, for a database that has catalogs
 * @returns
 */
export async function queryDB(
  connection: DG.DataConnection | null,
  tableName: string,
  match: string,
  matchValue: string | number,
  schemaName: string,
  joinOptions: QueryJoinOptions[] = [],
  max1: boolean = false,
  catalogName?: string
): Promise<DG.DataFrame> {
  if (connection == null)
    return DG.DataFrame.create(0);
  const {sql, joins, columnTables, columnSchemas} = buildValueQuery(connection.dataSource, tableName, match, matchValue, schemaName, joinOptions, max1, catalogName);
  const ciMapGet = (map: {[k: string]: string}, key: string): string | undefined => {
    if (map[key] !== undefined) return map[key];
    const lk = key.toLowerCase();
    for (const k of Object.keys(map)) {
      if (k.toLowerCase() === lk)
        return map[k];
    }
    return undefined;
  };

  const q = connection.query('getDBValueInfo', sql);

  const res: DG.DataFrame = await q.apply({});
  if (!res || !res.rowCount)
    return DG.DataFrame.create(0);
  // add correct tags for dbschema and dbTables.
  for (const col of res.columns) {
    const colName = col.name?.trim();
    if (!colName)
      continue; // should not happen
    col.setTag(DG.Tags.DbSchema, ciMapGet(columnSchemas, colName) ?? schemaName);
    const existingTableTag = col.getTag(DG.Tags.DbTable);
    if (!!existingTableTag && existingTableTag !== 'null') {
      // do nothing, the table is already set
      continue;
    }
    const tableForCol = ciMapGet(columnTables, colName);
    if (tableForCol) {
      col.setTag(DG.Tags.DbTable, tableForCol);
    } else if (joins
      .some(
        (tn) => /\(\d\)$/.test(colName) ||
          (colName.toLowerCase().startsWith(tn.tableName.toLowerCase()) && colName.length > tn.tableName.length + 1 && tn.select.some((s) => s.toLowerCase() === colName.substring(tn.tableName.length + 1).toLowerCase()))
      )
    ) {
      // not sure in this case, but try to avoid setting wrong table
      continue;
    } else {
      col.setTag(DG.Tags.DbTable, tableName); // main table, everything else should be there
    }
  }
  return res ?? DG.DataFrame.create(0);
}

function getNameBrackets(dataSource: string | undefined): [string, string] {
  // in future - substitute with direct api call
  switch (dataSource?.toLowerCase()) {
    case 'bigquery':
    case 'databricks':
    case 'mysql':
      return ['`', '`'];
    default:
      return ['"', '"']; // most common
  }
}

export function getConnectionNameBrackets(connection: DG.DataConnection): [string, string] {
  return getNameBrackets(connection.dataSource);
}


export async function queryDBMultiple(
  connection: DG.DataConnection | null,
  tableName: string,
  match: string,
  matchValues: (string | number)[]
): Promise<DG.DataFrame | string> {
  if (connection == null)
    return DG.DataFrame.create(0);
  if ((matchValues?.length ?? 0) < 1) return DG.DataFrame.create(0);

  // matchValues can be more than 1000, which is the limit for oracle and others
  // we need to split it in chunks of 500, to be on the safe side
  const chunkSize = 500;
  let outDataFrame: DG.DataFrame | null = null;
  for (let i = 0; i < matchValues.length; i += chunkSize) {
    const chunk = matchValues.slice(i, i + chunkSize);
    const q = connection.query(
      'getDBValueInfoMult',
      `
            --name: getDBValueInfoMult
            --output: dataframe result
            select * from ${tableName} where ${match} in (${chunk.map((v) => `'${v}'`).join(',')})
        `
    );
    const res = await q.apply({});
    if (res) {
      if (!outDataFrame) outDataFrame = res;
      else outDataFrame.append(res, true);
    }
  }
  return outDataFrame ?? DG.DataFrame.create(0);
}
