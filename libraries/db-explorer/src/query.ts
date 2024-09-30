import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {QueryJoinOptions} from './types';

export async function queryDB(
  connection: DG.DataConnection,
  tableName: string,
  match: string,
  matchValue: string | number,
  schemaName: string,
  joinOptions: QueryJoinOptions[] = []
): Promise<DG.DataFrame> {
  const matchValueStr = typeof matchValue === 'string' ? `'${matchValue}'` : matchValue;
  const startCharCode = 98; // ascii b

  const applicableJoins = joinOptions.filter((opt) => opt.fromTable === tableName);

  const tableAliases = applicableJoins.map((_, i) => String.fromCharCode(startCharCode + i));
  const otherCols =
    applicableJoins
      .map((opt, i) => {
        const alias = tableAliases[i];
        return opt.select.map((col) => `${alias}.${col}`).join(', ') + (i < applicableJoins.length - 1 ? ', ' : '');
      })
      .join(', ') + (applicableJoins.length > 0 ? ', ' : '');
  const joinStr = applicableJoins
    .map((opt, i) => {
      const alias = tableAliases[i];
      return `left join "${schemaName}".${opt.tableName} ${alias} on a.${opt.columnName} = ${alias}.${opt.onColumn}`;
    })
    .join(' ');
  const q = connection.query(
    'getDBValueInfo',
    `
        --name: getDBValueInfo
        --output: dataframe result
        select ${otherCols} a.* 
        from "${schemaName}".${tableName} a
        ${joinStr}
        where a.${match} = ${matchValueStr}
    `
  );
  return (await q.apply({})) ?? DG.DataFrame.create(0);
}

export async function queryDBMultiple(
  connection: DG.DataConnection,
  tableName: string,
  match: string,
  matchValues: (string | number)[]
): Promise<DG.DataFrame | string> {
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
