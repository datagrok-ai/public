import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {benchmarkQuery} from './benchmark';

/** Rows returned by each streaming query in queries/clickhouse-streaming-test.sql. */
const EXPECTED_ROWS: {[query: string]: number} = {
  'ClickHouseStreamingInts': 1000000,
  'ClickHouseStreamingLowCardinality': 1000000,
  'ClickHouseStreamingMixed': 1000000,
  'ClickHouseStreamingWideStrings': 250000,
  'ClickHouseStreamingUnsupportedTypes': 200000,
};

// Benchmark-only: every query here moves 100 MB+ and the connection targets a shared demo DB,
// so these stay out of the normal CI run (same convention as the Benchmarks category).
category('ClickHouse Streaming', () => {
  for (const query of Object.keys(EXPECTED_ROWS))
    test(query.replace('ClickHouseStreaming', ''), async () => {
      const df: DG.DataFrame = await grok.functions.call('Dbtests:' + query);
      expect(df.rowCount, EXPECTED_ROWS[query]);
      return await benchmarkQuery(query, DG.Test.isInBenchmark ? 3 : 1);
    }, {timeout: 600000, benchmark: true, stressTest: true});

  // The chunk cut is every 50_000 rows, so these three land either side of the first boundary.
  test('Chunk boundary', async () => {
    for (const rows of [49999, 50000, 50001]) {
      const df: DG.DataFrame = await grok.functions.call('Dbtests:ClickHouseStreamingRows', {rows});
      expect(df.rowCount, rows);
    }
  }, {timeout: 600000, benchmark: true});

  // Both connections point at the same server/db (db.datagrok.ai:18125, `test`); they differ
  // only in the engine that carries the rows back, so the same SQL is a fair head-to-head.
  test('ADBC vs JDBC on identical SQL', async () => {
    const sql = 'SELECT number AS id, ' +
      "toLowCardinality(['Fail', 'Pending', 'Success'][(number % 3) + 1]) AS status " +
      'FROM numbers(1000000)';
    const times: {[connectionName: string]: number} = {};
    for (const name of ['ClickHouseStreaming', 'ClickHouseDBTests']) {
      const connection = await grok.dapi.connections.filter(`name = "${name}"`).first();
      const startTime = Date.now();
      const call = await connection.query('adhoc', sql).prepare().call();
      expect((call.getOutputParamValue() as DG.DataFrame).rowCount, 1000000);
      times[`${name} (${connection.parameters['connectionType'] ?? 'JDBC'}), ms`] = Date.now() - startTime;
    }
    return times;
  }, {timeout: 600000, benchmark: true});
}, {node: true});
