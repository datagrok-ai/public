import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataQuery} from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {getCallTime} from '../benchmarks/benchmark';

category('Cache', () => {
  const clientSideCache: boolean = grok.shell.settings.clientSideCache;
  const testConnections: String[] = ['PostgreSQLDBTests', 'PostgreSQLDBTestsCached'];

  before(async () => {
    await cleanCache(testConnections);
    grok.shell.settings.clientSideCache = true;
  });

  test('Client function cache scalar int', async () => {
    grok.functions.register({
      signature: 'int getEpochDate()',
      run: () => {
        return Date.now();
      },
      options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    });

    await compareTwoFunctionResults('getEpochDate');
  });

  test('Client function cache scalar double', async () => {
    grok.functions.register({
      signature: 'double getRandomDouble()',
      run: () => {
        return Math.random();
      },
      options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    });

    await compareTwoFunctionResults('getRandomDouble');
  });

  test('Client function cache scalar string', async () => {
    grok.functions.register({
      signature: 'string getRandomStringDate()',
      run: () => {
        return new Date().toISOString();
      },
      options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    });

    await compareTwoFunctionResults('getRandomStringDate');
  });

  test('Client function cache invalidation', async () => {
    const functionName = 'getEpochDateInv()';
    const functionCallName = 'getEpochDateInv';
    grok.functions.register({
      signature: `int ${functionName}`,
      run: () => {
        return Date.now();
      },
      options: {'cache': 'true', 'cache.invalidateOn': '0 * * * *'},
    });
    const first = await compareTwoFunctionResults(functionCallName);
    let transaction: any;
    let db: any;
    const cacheName = 'function_results_cache';
    try {
      const request = window.indexedDB.open(cacheName, 4);
      request.onsuccess = function() {
        db = request.result;
        transaction = db.transaction(cacheName, 'readwrite');
        const objectStore = transaction.objectStore(cacheName);
        let map;
        const getRequest = objectStore.get(functionName);
        getRequest.onsuccess = function() {
          map = getRequest.result;
          map['.expires'] = new Date(1970).toISOString();
          objectStore.delete(functionName);
          objectStore.add(map, functionName);
        };
      };
      if (transaction != null)
        transaction.commit();
    } catch (e) {
      if (transaction != null)
        transaction.abort();
      throw e;
    } finally {
      if (db != null)
        db.close;
    }
    const second = await grok.functions.call(functionCallName);
    expect(first == second, false);
  });

  test('Scalars cache test', async () => await basicCacheTest('PostgresqlScalarCacheTest'));

  test('TestWide table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableWide'));

  test('TestNormal table cache test', async () => await basicCacheTest('PostgresqlTestCacheTableNormal'));

  test('Connection cache test', async () => await basicCacheTest('PostgresqlCachedConnTest'), {timeout: 120000});

  test('Connection cache invalidation test', async () => {
    const connection = await grok.dapi.connections.filter(`name="${testConnections[1]}"`).first();
    await invalidationCacheTest(connection.query('test1', 'SELECT *, pg_sleep(0.1) FROM MOCK_DATA;'), 2);
  }, {timeout: 120000});

  test('Query cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries.filter(`friendlyName="PostgresqlCacheInvalidateQueryTest"`).first();
    await invalidationCacheTest(dataQuery, 2);
  });

  test('Scalars cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries
      .filter(`friendlyName="PostgresqlScalarCacheInvalidationTest"`).first();
    await invalidationCacheTest(dataQuery, 2);
  });

  after(async () => {
    await cleanCache(testConnections);
    grok.shell.settings.clientSideCache = clientSideCache;
  });
});

async function invalidationCacheTest(dataQuery: DataQuery, days: number): Promise<void> {
  const funcCall1 = dataQuery.prepare();
  const firstExecutionTime = await getCallTime(funcCall1);
  await delay(100);
  funcCall1.started = dayjs().subtract(days, 'day');
  await grok.dapi.functions.calls.save(funcCall1);
  const secondExecutionTime = await getCallTime(dataQuery.prepare());
  const isEqual: boolean = (secondExecutionTime <= firstExecutionTime + firstExecutionTime * 0.5) &&
        (secondExecutionTime >= firstExecutionTime - firstExecutionTime * 0.5);
  expect(isEqual, true,
    `The second execution time ${secondExecutionTime} ms
        is not approximately equals to the first execution time ${firstExecutionTime} ms`);
}

async function basicCacheTest(query: String): Promise<void> {
  const dataQuery = await grok.dapi.queries.filter(`friendlyName="${query}"`).first();
  const funcCall1 = dataQuery.prepare();
  const firstExecutionTime = await getCallTime(funcCall1);
  await delay(100);
  const funcCall2 = dataQuery.prepare();
  const secondExecutionTime = await getCallTime(funcCall2);
  expect(firstExecutionTime > secondExecutionTime * 2, true,
    `The first execution time ${firstExecutionTime} ms
        is no more than twice the second execution time ${secondExecutionTime} ms for ${query}`);
}

async function compareTwoFunctionResults(funcName: string): Promise<any> {
  const first = await grok.functions.call(funcName);
  await new Promise((res) => setTimeout(res, 25));
  const second = await grok.functions.call(funcName);
  expect(first, second);
  return second;
}

async function cleanCache(connections: String[]): Promise<void> {
  for (const conn of connections) {
    await grok.functions.call('DropConnectionCache',
      {connection: await grok.dapi.connections.filter(`name="${conn}"`).first()});
  }
}
