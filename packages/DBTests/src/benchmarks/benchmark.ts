import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {FuncCall} from 'datagrok-api/dg';

category('Benchmarks', () => {
  test('Sequential 100', async () => {
    const dataConnection = await grok.dapi.connections.filter(`name = "PostgreSQLDBTests"`).first();
    const times = [];
    for (let i = 0; i < 100; i++) {
      const dataQuery = await dataConnection.query('test', 'SELECT 1');
      dataQuery.adHoc = true;
      times.push(await getCallTime(dataQuery.prepare()));
    }
    return getTestResult(times, times.length);
  }, {timeout: 120000});

  test('Parallel 200', async () => {
    const callCount: number = 200;
    const dataConnection = await grok.dapi.connections.filter(`name = "PostgreSQLDBTests"`).first();
    const calls = [];
    for (let i = 0; i < callCount; i++) {
      const dataQuery = await dataConnection.query('test', 'SELECT 1');
      dataQuery.adHoc = true;
      calls.push(getCallTime(dataQuery.prepare()));
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

  test('Scalar Cache test', async () => {
    const dataQueryName = 'PostgresqlScalarCacheTestTableLong';
    const firstExecutionTime = await getDataQueryTime(dataQueryName);
    const secondExecutionTime = await getDataQueryTime(dataQueryName);
    expect(firstExecutionTime > secondExecutionTime * 2, true);
  }, {skipReason: 'Feature of caching scalars in development'});

  test('Compression int', async () => {
    const compressionOnTime = await getDataQueryTime('PostgresqlCompressionIntOn');
    const compressionOffTime = await getDataQueryTime('PostgresqlCompressionIntOff');
    expect(compressionOnTime < compressionOffTime * 2, true);
  }, {skipReason: 'Feature of compression in development'});

  test('Cache test for Table_Wide', async () => {
    const dataQueryName = 'PostgresqlTestCacheTableWide';
    const firstExecutionTime = await getDataQueryTime(dataQueryName);
    const secondExecutionTime = await getDataQueryTime(dataQueryName);
    expect(firstExecutionTime > secondExecutionTime * 2, true);
  });

  test('Cache test for Table_Normal', async () => {
    const dataQueryName = 'PostgresqlTestCacheTableNormal';
    const firstExecutionTime = await getDataQueryTime(dataQueryName);
    const secondExecutionTime = await getDataQueryTime(dataQueryName);
    expect(firstExecutionTime > secondExecutionTime * 2, true);
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
