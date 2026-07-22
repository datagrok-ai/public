import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

// Browser-bound DB tests relocated from DBTests so that package runs fully headless
// under the grok test Node pass. They call DBTests queries cross-package, so DBTests
// must be published on the stand (always true on CI and dev).
//
// Why each needs the browser:
// - Data sync: streams the query into a table view (`processed: false` + grok.shell.tv).
// - Client cache: grok.functions.clientCache is the IndexedDB-backed browser cache.
// - Benchmarks: "client cached" uses the client cache; "server cached" trips a Node
//   runtime bug in batched cached-result delivery (DataFrame.appendMerge on a null
//   target on cache hit) — move it back once that is fixed.

category('DB: Data sync', () => {
  test('grok connect streaming and data sync', async () => {
    const func: DG.Func = await grok.functions.eval('DbTests:PostgresqlTableWide');
    await func.prepare().call(false, undefined, {processed: false});
    await awaitCheck(() => grok.shell.tv?.table?.name === 'PostgresqlTableWide', 'Query first batch timeout', 30000);
    grok.shell.closeTable(grok.shell.tv!.table!);
    grok.shell.tv?.close();
  }, {stressTest: true});
});

category('DB: Client cache', () => {
  before(async () => {
    await grok.functions.clientCache.start();
  });

  test('Scalars cache test', async () => {
    await grok.functions.clientCache.clear();
    await basicCacheTest('PostgresqlScalarCacheTestClient', 2);
  });

  test('TestNormal table cache test', async () => {
    await grok.functions.clientCache.clear();
    await basicCacheTest('PostgresqlTestCacheTableNormalClient', 3);
  });

  after(async () => {
    await grok.functions.clientCache.clear();
  });
});

category('DB: Benchmarks', () => {
  test('Performance: TestWide client cached', async () => {
    const count = DG.Test.isInBenchmark ? 5 : 2;
    return await benchmarkQuery('PostgresqlTableWideCachedClient', count);
  }, {benchmark: true, stressTest: true});

  test('Performance: TestWide server cached', async () => {
    const count = DG.Test.isInBenchmark ? 5 : 2;
    return await benchmarkQuery('PostgresqlTableWideCachedServer', count);
  }, {timeout: 120000, benchmark: true, stressTest: true});
});

async function getDataQueryTime(dataQueryName: string): Promise<number> {
  const func: DG.Func = await grok.functions.eval('Dbtests:' + dataQueryName);
  const start = Date.now();
  await func.prepare().call();
  return Date.now() - start;
}

async function basicCacheTest(query: string, times: number): Promise<void> {
  const firstExecutionTime = await getDataQueryTime(query);
  await delay(1000);
  const secondExecutionTime = await getDataQueryTime(query);
  // eslint-disable-next-line max-len
  expect(firstExecutionTime > secondExecutionTime * times, true, `The first execution time ${firstExecutionTime} ms is no more than twice the second execution time ${secondExecutionTime} ms for ${query}`);
}

async function benchmarkQuery(query: string, count: number): Promise<object> {
  const results: number[] = [];
  for (let i = 0; i < count; i++)
    results.push(await getDataQueryTime(query));
  const sum = results.reduce((p, c) => p + c, 0);
  return DG.toDart({'Average time': sum / results.length,
    'Min time': Math.min(...results), 'Max time': Math.max(...results)});
}
