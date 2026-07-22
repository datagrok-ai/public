import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';

import {category, test, expect} from '@datagrok-libraries/test/src/test';

// Benchmark of the celery Node worker path (meta.queue / meta.server functions).
// Runs only with `grok test --benchmark`. Timings are printed as `BENCH:` console
// lines (and each test's own duration in the report approximates the same number).
category('Celery: node worker benchmark', () => {
  const containerNames = ['cvm-tests-queue', 'cvm-tests-queue-celery'];

  async function queueContainer(): Promise<DG.DockerContainer> {
    for (const name of containerNames) {
      const c = await grok.dapi.docker.dockerContainers.filter(`name = "${name}"`).first();
      if (c) return c;
    }
    throw new Error('CVMTests queue worker container not found');
  }

  function log(name: string, ms: number, extra: string = ''): void {
    console.log(`BENCH: ${name}: ${Math.round(ms)} ms${extra ? ' ' + extra : ''}`);
  }

  function makeTable(rows: number): DG.DataFrame {
    const ints = new Int32Array(rows);
    const floats = new Float32Array(rows);
    const strs = new Array<string>(rows);
    const bools = new Array<boolean>(rows);
    const dates = new Array<dayjs.Dayjs>(rows);
    const t0 = dayjs('2020-01-01T00:00:00Z');
    for (let i = 0; i < rows; i++) {
      ints[i] = i;
      floats[i] = i * 0.5;
      strs[i] = `row ${i % 1000}`;
      bools[i] = i % 2 === 0;
      dates[i] = t0.add(i % 86400, 'second');
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('i', ints),
      DG.Column.fromFloat32Array('d', floats),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 's', strs),
      DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'b', bools),
      DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'dt', dates),
    ]);
  }

  async function timedCall(fn: string, params: object): Promise<number> {
    const t0 = performance.now();
    await grok.functions.call(fn, params);
    return performance.now() - t0;
  }

  test('Cold start: spin-up + first call', async () => {
    const container = await queueContainer();
    await grok.dapi.docker.dockerContainers.stop(container.id, true);
    const ms = await timedCall('CVMTests:jsCvmInt', {x: 42});
    log('cold start (container start + worker init + loadPackage + call)', ms);
  }, {benchmark: true, timeout: 600000});

  test('Warm scalar call', async () => {
    const ms = await timedCall('CVMTests:jsCvmInt', {x: 42});
    log('warm scalar call', ms);
  }, {benchmark: true, timeout: 90000});

  test('Scalars x20', async () => {
    const times: number[] = [];
    for (let i = 0; i < 20; i++)
      times.push(await timedCall('CVMTests:jsCvmInt', {x: i}));
    times.sort((a, b) => a - b);
    const avg = times.reduce((a, b) => a + b, 0) / times.length;
    log('scalar x20', avg, `avg (min ${Math.round(times[0])}, ` +
      `median ${Math.round(times[10])}, max ${Math.round(times[19])})`);
  }, {benchmark: true, timeout: 300000});

  test('Small table round trip (1k x 5)', async () => {
    const df = makeTable(1000);
    const payload = df.toByteArray().length;
    const times: number[] = [];
    for (let i = 0; i < 5; i++)
      times.push(await timedCall('CVMTests:jsCvmDataframe', {df}));
    times.sort((a, b) => a - b);
    const avg = times.reduce((a, b) => a + b, 0) / times.length;
    log('small table (1k x 5) round trip', avg,
      `avg of 5 (min ${Math.round(times[0])}, max ${Math.round(times[4])}), d42 payload ${(payload / 1024).toFixed(1)} KB`);
    expect(times.length, 5);
  }, {benchmark: true, timeout: 300000});

  test('Large table round trip (500k x 5)', async () => {
    const df = makeTable(500000);
    const payload = df.toByteArray().length;
    const ms = await timedCall('CVMTests:jsCvmDataframe', {df});
    log('large table (500k x 5) round trip', ms,
      `d42 payload ${(payload / 1024 / 1024).toFixed(1)} MB`);
    const ms2 = await timedCall('CVMTests:jsCvmDataframe', {df});
    log('large table (500k x 5) round trip, 2nd', ms2);
  }, {benchmark: true, timeout: 600000});
});
