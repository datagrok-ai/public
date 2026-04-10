import {category, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

const WARMUP = 3;
const MEASURED = 5;

category('Databricks Benchmarks', () => {
  test('BENCH_SMALL (1K rows, mixed types)', async () => {
    return await benchmarkQuery('DatabricksBenchSmall', WARMUP, MEASURED);
  }, {timeout: 300000, benchmark: true});

  test('BENCH_MEDIUM (100K rows, mixed types)', async () => {
    return await benchmarkQuery('DatabricksBenchMedium', WARMUP, MEASURED);
  }, {timeout: 600000, benchmark: true});

  test('BENCH_LARGE (1M rows, mixed types)', async () => {
    return await benchmarkQuery('DatabricksBenchLarge', WARMUP, MEASURED);
  }, {timeout: 600000, benchmark: true});

  test('BENCH_CATEGORICAL (1M rows, 10 string cols)', async () => {
    return await benchmarkQuery('DatabricksBenchCategorical', WARMUP, MEASURED);
  }, {timeout: 600000, benchmark: true});

  test('SELECT 1 latency', async () => {
    return await benchmarkQuery('DatabricksSelect1', WARMUP, 10);
  }, {timeout: 120000, benchmark: true});
});

async function benchmarkQuery(queryName: string, warmup: number, measured: number): Promise<object> {
  // Warmup runs
  for (let i = 0; i < warmup; i++) {
    console.log(`  warmup ${i + 1}/${warmup}: ${queryName}`);
    await runQuery(queryName);
  }

  // Measured runs
  const times: number[] = [];
  const ttfrs: number[] = [];
  let rows = 0;
  let cols = 0;

  for (let i = 0; i < measured; i++) {
    console.log(`  measured ${i + 1}/${measured}: ${queryName}`);
    const result = await runQueryWithTiming(queryName);
    times.push(result.ttc);
    ttfrs.push(result.ttfr);
    rows = result.rows;
    cols = result.cols;
  }

  return {
    'Query': queryName,
    'Rows': rows,
    'Columns': cols,
    'Warmup runs': warmup,
    'Measured runs': measured,
    'Avg TTC (ms)': Math.round(times.reduce((a, b) => a + b, 0) / times.length),
    'Min TTC (ms)': Math.min(...times),
    'Max TTC (ms)': Math.max(...times),
    'Avg TTFR (ms)': Math.round(ttfrs.reduce((a, b) => a + b, 0) / ttfrs.length),
    'Min TTFR (ms)': Math.min(...ttfrs),
    'Max TTFR (ms)': Math.max(...ttfrs),
  };
}

async function runQuery(queryName: string): Promise<void> {
  await grok.functions.call('Dbtests:' + queryName);
}

async function runQueryWithTiming(queryName: string): Promise<{ttc: number; ttfr: number; rows: number; cols: number}> {
  const startTime = Date.now();
  let ttfr = 0;
  let ttfrSet = false;

  const sub = grok.functions.onParamsUpdated.subscribe(() => {
    if (!ttfrSet) {
      ttfr = Date.now() - startTime;
      ttfrSet = true;
    }
  });

  try {
    const call = await grok.functions.call('Dbtests:' + queryName);
    const ttc = Date.now() - startTime;
    if (!ttfrSet)
      ttfr = ttc;

    let rows = 0;
    let cols = 0;
    if (call instanceof DG.DataFrame) {
      rows = call.rowCount;
      cols = call.columns.length;
    }

    return {ttc, ttfr, rows, cols};
  }
  finally {
    sub.unsubscribe();
  }
}
