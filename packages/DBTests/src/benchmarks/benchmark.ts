import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataQuery, FuncCall} from 'datagrok-api/dg';
import dayjs from 'dayjs';

category('Benchmarks', () => {
  test('Sequential 100', async () => {
    const dataConnection = await grok.dapi.connections.filter(`name = "PostgreSQLDBTests"`).first();
    const times = [];
    for (let i = 0; i < 100; i++) {
      const dataQuery = await dataConnection.query('test', 'SELECT 1');
      const funcCall = dataQuery.prepare();
      funcCall.adHoc = true;
      times.push(await getCallTime(funcCall));
    }
    return getTestResult(times, times.length);
  }, {timeout: 120000});

  test('Parallel 200', async () => {
    const callCount: number = 200;
    const dataConnection = await grok.dapi.connections.filter(`name = "PostgreSQLDBTests"`).first();
    const calls = [];
    for (let i = 0; i < callCount; i++) {
      const dataQuery = await dataConnection.query('test', 'SELECT 1');
      const funcCall = dataQuery.prepare();
      funcCall.adHoc = true;
      calls.push(getCallTime(funcCall));
    }
    const times: number[] = [];
    const results = await Promise.allSettled(calls);
    results.forEach((result) => {
      if (result.status == 'fulfilled')
        times.push(result.value);
    });
    return getTestResult(times, callCount);
  });

  test('TestNormal', async () => {
    return `Execution time: ${await getDataQueryTime('PostgresqlPerfTestTableNormal')}`;
  });

  test('TestWide', async () => {
    return `Execution time: ${await getDataQueryTime('PostgresqlPerfTestTableWide')}`;
  });

  test('TestLong', async () => {
    return `Execution time: ${await getDataQueryTime('PostgresqlPerfTestTableLong')}`;
  });

  test('Compression int', async () => {
    const compressionOnTime = await getDataQueryTime('PostgresqlCompressionIntOn');
    const compressionOffTime = await getDataQueryTime('PostgresqlCompressionIntOff');
    expect(compressionOnTime < compressionOffTime * 2, true);
  }, {skipReason: 'Feature of compression in development'});

  test('Scalar float cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestFloat'));

  test('Scalar int cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestInt'));

  test('Scalar string cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestString'));

  test('TestWide table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableWide'));

  test('TestNormal table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableNormal'));

  test('Connection cache test', async () => await basicCacheTest('PostgresqlCachedConnTest'));

  test('Connection cache invalidation test', async () => {
    const connection = await grok.dapi.connections.filter(`name="PostgreSQLDBTestsCached"`).first();
    await invalidationCacheTest(connection.query('test1', 'SELECT * FROM MOCK_DATA;'), 1);
  });

  test('Query cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries.filter(`friendlyName="PostgresqlCacheInvalidateQueryTest"`).first();
    await invalidationCacheTest(dataQuery, 1);
  });
});

function getTestResult(times: number[], expectedCount: number): String {
  const totalTime = times.reduce((acc, currentValue) => acc + currentValue, 0);
  return `Total execution time of all calls, ms: ${totalTime}\n` +
      `Percentage of success: ${Math.round((times.length / expectedCount) * 100)}%\n` +
      `Average execution time, ms: ${totalTime / times.length}\n` +
      `Max execution time, ms: ${Math.max(...times)}\n` +
      `Min execution time, ms: ${Math.min(...times)}\n`;
}

async function getCallTime(call: FuncCall): Promise<number> {
  const start = Date.now();
  await call.call();
  return Date.now() - start;
}

async function getDataQueryTime(dataQueryName: string): Promise<number> {
  const startTime = Date.now();
  await grok.functions.eval('Dbtests:' + dataQueryName);
  return Date.now() - startTime;
}

async function invalidationCacheTest(dataQuery: DataQuery, minutes: number): Promise<void> {
  const funcCall1 = dataQuery.prepare();
  const firstExecutionTime = await getCallTime(funcCall1);
  funcCall1.started = dayjs().subtract(minutes, 'minute');
  await grok.dapi.functions.calls.save(funcCall1);
  const secondExecutionTime = await getCallTime(dataQuery.prepare());
  const isEqual: boolean = (secondExecutionTime <= firstExecutionTime + firstExecutionTime * 0.1) &&
      (secondExecutionTime >= firstExecutionTime - firstExecutionTime * 0.1);
  expect(isEqual, true,
    `The second execution time ${secondExecutionTime} ms
        is not approximately equals to the first execution time ${firstExecutionTime} ms`);
}

async function basicCacheTest(query: String): Promise<void> {
  const dataQuery = await grok.dapi.queries.filter(`friendlyName="${query}"`).first();
  const firstExecutionTime = await getCallTime(dataQuery.prepare());
  const secondExecutionTime = await getCallTime(dataQuery.prepare());
  expect(firstExecutionTime > secondExecutionTime * 2, true,
    `The first execution time ${firstExecutionTime} ms
        is no more than twice the second execution time ${secondExecutionTime} ms for ${query}`);
}
