import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {DataQuery} from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {getCallTime} from '../benchmarks/benchmark';

category('Cache', () => {
  test('Scalar float cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestFloat'));

  test('Scalar int cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestInt'));

  test('Scalar string cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestString'));

  test('Scalar datetime cache test', async () => await basicCacheTest('PostgresqlScalarCacheTestDate'));

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

  test('Scalar float cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries
      .filter(`friendlyName="PostgresqlScalarCacheInvalidationTestFloat"`).first();
    await invalidationCacheTest(dataQuery, 1);
  });

  test('Scalar int cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries
      .filter(`friendlyName="PostgresqlScalarCacheInvalidationTestInt"`).first();
    await invalidationCacheTest(dataQuery, 1);
  });

  test('Scalar string cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries
      .filter(`friendlyName="PostgresqlScalarCacheInvalidationTestString"`).first();
    await invalidationCacheTest(dataQuery, 1);
  });

  test('Scalar date cache invalidation test', async () => {
    const dataQuery = await grok.dapi.queries
      .filter(`friendlyName="PostgresqlScalarCacheInvalidationTestDate"`).first();
    await invalidationCacheTest(dataQuery, 1);
  });
});

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
  const funcCall1 = dataQuery.prepare();
  funcCall1.adHoc = true;
  const firstExecutionTime = await getCallTime(funcCall1);
  const funcCall2 = dataQuery.prepare();
  funcCall2.adHoc = true;
  const secondExecutionTime = await getCallTime(funcCall2);
  expect(firstExecutionTime > secondExecutionTime * 2, true,
    `The first execution time ${firstExecutionTime} ms
        is no more than twice the second execution time ${secondExecutionTime} ms for ${query}`);
}
