/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {QueryJoinOptions} from './types';

/**
 *
 * @param connection
 * @param tableName
 * @param match - Pass the @match and @matchValue as empty strings to get 1 random row (for example queries)
 * @param matchValue
 * @param schemaName
 * @param joinOptions
 * @param max1
 * @returns
 */
export async function queryDB(
  connection: DG.DataConnection | null,
  tableName: string,
  match: string,
  matchValue: string | number,
  schemaName: string,
  joinOptions: QueryJoinOptions[] = [],
  max1: boolean = false
): Promise<DG.DataFrame> {
  if (connection == null)
    return DG.DataFrame.create(0);
  const matchValueStr = typeof matchValue === 'string' ? `'${matchValue}'` : matchValue;
  const startCharCode = 98; // ascii b
  const [nbs, nbe] = getConnectionNameBrackets(connection);

  const applicableJoins = joinOptions.filter((opt) => opt.fromTable === tableName);

  const tableAliases = applicableJoins.map((_, i) => String.fromCharCode(startCharCode + i));
  const otherCols =
    applicableJoins
      .map((opt, i) => {
        const alias = tableAliases[i];
        return opt.select.map((col) => {
          // if the col is confugured with alias, make sure it has quotes
          // split cace insensitive way
          const parts = col.split(/ as /i);
          if (parts.length > 1)
            return `${alias}.${nbs}${parts[0].trim()}${nbe} as ${nbs}${parts[1].trim()}${nbe}`;

          return `${alias}.${nbs}${col}${nbe}`;
        }).join(', ');
      })
      .join(', ');
  const joinStr = applicableJoins
    .map((opt, i) => {
      const alias = tableAliases[i];
      return `left join ${nbs}${schemaName}${nbe}.${nbs}${opt.tableName}${nbe} ${alias} on a.${nbs}${opt.columnName}${nbe} = ${alias}.${nbs}${opt.onColumn}${nbe}`;
    })
    .join(' ');
  // let qb = connection.buildQuery(`"${schemaName}".${tableName}`).selectAll();
  // applicableJoins.forEach((opt, i) => {
  //   const alias = tableAliases[i];
  //   qb = qb.leftJoin(`${schemaName}.${opt.tableName}`, [opt.columnName], [opt.onColumn], alias);
  // });
  // qb = qb.where(match, ` = ${matchValueStr}`);
  // console.log(qb.build().query);
  const isExampleQuery = match === '' && matchValue === ''; // used for setup editor, will not happen otherwise if not intentionally

  const q = connection.query(
    'getDBValueInfo',
    `
        --name: getDBValueInfo
        --output: dataframe result
        select a.* ${(applicableJoins.length > 0 ? ',' : '')} ${otherCols} 
        from ${nbs}${schemaName}${nbe}.${nbs}${tableName}${nbe} a
        ${joinStr}
        ${isExampleQuery ? '' : `where a.${match} = ${matchValueStr}`}
        ${max1 || isExampleQuery ? 'limit 1' : ''}
    `
  );
  return (await q.apply({})) ?? DG.DataFrame.create(0);
}

export function getConnectionNameBrackets(connection: DG.DataConnection): [string, string] {
  // in future - substitute with direct api call
  const providerType = connection.dataSource?.toLowerCase();
  switch (providerType) {
    case 'bigquery':
    case 'databricks':
    case 'mysql':
      return ['`', '`'];
    default:
      return ['"', '"']; // most common
  }
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
