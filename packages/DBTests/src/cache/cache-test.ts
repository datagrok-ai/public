import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {getCallTime} from '../benchmarks/benchmark';

category('Client cache', () => {
  before(async () => {
    await grok.functions.clientCache.start();
  });

  test('Scalars cache test', async () => {
    await grok.functions.clientCache.clear();
    await basicCacheTest('PostgresqlScalarCacheTestClient');
  });

  test('TestNormal table cache test', async () => {
    await grok.functions.clientCache.clear();
    await basicCacheTest('PostgresqlTestCacheTableNormalClient');
  });

  test('Performance: 5 heavy DataQuery with no cache', async () => {
    return await runHeavy(false);
  }, {timeout: 300000});

  test('Performance: 5 heavy DataQuery with cache', async () => {
    return await runHeavy(true);
  }, {timeout: 300000});

  after(async () => {
    await grok.functions.clientCache.clear();
  });
});

category('Server cache', () => {
  const testConnections: String[] = ['PostgreSQLDBTests', 'PostgreSQLDBTestsCached'];

  before(async () => {
    await cleanCache(testConnections);
  });

  test('Scalars cache test', async () => await basicCacheTest('PostgresqlScalarCacheTest'));

  test('TestWide table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableWide'), {timeout: 120000});

  test('TestNormal table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableNormal'), {timeout: 120000});

  test('Connection cache test', async () => await basicCacheTest('PostgresqlCachedConnTest'), {timeout: 120000});

  test('Connection cache invalidation test', async () => {
    const dataQuery: DG.DataQuery = await grok.functions.eval(`Dbtests:TestConnCache`);
    await invalidationCacheTest(dataQuery, 2);
  }, {timeout: 120000});

  test('Query cache invalidation test', async () => {
    const dataQuery = await grok.functions.eval(`Dbtests:PostgresqlCacheInvalidateQueryTest`);
    await invalidationCacheTest(dataQuery, 2);
  });

  test('Scalars cache invalidation test', async () => {
    const dataQuery = await grok.functions.eval(`Dbtests:PostgresqlScalarCacheInvalidationTest`);
    await invalidationCacheTest(dataQuery, 2);
  });

  test('Cached conn DataFrame id diff', async () => {
    const connection = await grok.dapi.connections.filter(`name="${testConnections[1]}"`).first();
    const funcCall = await connection.query('test', 'SELECT * FROM MOCK_DATA').prepare().call();
    const firstId = funcCall.outputs['result'].id;
    const funcCall1 = await connection.query('test', 'SELECT * FROM MOCK_DATA').prepare().call();
    const secondId = funcCall1.outputs['result'].id;
    expect(firstId !== secondId, true, 'Ids are the same');
  });

  after(async () => {
    await cleanCache(testConnections);
  });
});

async function invalidationCacheTest(dataQuery: DG.DataQuery, days: number): Promise<void> {
  const start = Date.now();
  const funcCall1 = await dataQuery.prepare().call();
  const firstExecutionTime = Date.now() - start;
  await delay(500);
  //@ts-ignore
  funcCall1.started = dayjs().subtract(days, 'day');
  await grok.dapi.functions.calls.save(funcCall1);
  await grok.functions.clientCache.clear();
  await delay(500);
  const secondExecutionTime = await getCallTime(dataQuery.prepare());
  const isEqual: boolean = (secondExecutionTime <= firstExecutionTime + firstExecutionTime * 0.5) &&
        (secondExecutionTime >= firstExecutionTime - firstExecutionTime * 0.5);
  // eslint-disable-next-line max-len
  expect(isEqual, true, `The second execution time ${secondExecutionTime} ms is not approximately equals to the first execution time ${firstExecutionTime} ms`);
}

async function basicCacheTest(query: string): Promise<void> {
  const dataQuery = await grok.functions.eval(`Dbtests:${query}`);
  const firstExecutionTime = await getCallTime(dataQuery.prepare());
  await delay(1000);
  const secondExecutionTime = await getCallTime(dataQuery.prepare());
  // eslint-disable-next-line max-len
  expect(firstExecutionTime > secondExecutionTime * 2, true, `The first execution time ${firstExecutionTime} ms is no more than twice the second execution time ${secondExecutionTime} ms for ${query}`);
}

async function cleanCache(connections: String[]): Promise<void> {
  for (const conn of connections) {
    await grok.functions.call('DropConnectionCache',
      {connection: await grok.dapi.connections.filter(`name="${conn}"`).first()});
  }
}

async function runHeavy(cache: boolean) {
  try {
    let queryName;
    if (!cache) {
      grok.functions.clientCache.stop();
      queryName = 'NotCachedHeavy';
    } else {
      await grok.functions.clientCache.clear();
      queryName = 'CachedHeavy';
    }

    const results = [];
    const dataQuery = await grok.functions.eval(`Dbtests:${queryName}`);
    if (cache)
      await dataQuery.prepare().call();
    for (let i = 0; i < 5; i++) {
      const start = Date.now();
      await dataQuery.prepare().call();
      results.push(Date.now() - start);
    }

    const sum = results.reduce((p, c) => p + c, 0);
    return DG.toDart({'Average time': sum / results.length,
      'Min time': Math.min(...results), 'Max time': Math.max(...results)});
  } catch (e) {
    grok.log.error(e);
  }
}
