import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {CruxModule, CruxWorker, cruxMl} from '../crux/crux-loader';
import {cruxFit, cruxPredict, cruxTerminate} from '../crux/crux-model';
import {copyCol, expectRejected, irisMeasurements, predictColumn, regressionDataset} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const TIMEOUT = 120000;
const FITS = 8;

function frame(cols: DG.Column[]): DG.DataFrame {
  return DG.DataFrame.fromColumns(cols.map(copyCol));
}

async function iris(): Promise<{features: DG.Column[], target: DG.Column}> {
  const df = await grok.dapi.files.readCsv(IRIS);
  return {features: irisMeasurements(df).columns.toList(), target: df.col('Species')!};
}

function regressionData(rows: number, cols: number): {features: DG.Column[], y: DG.Column} {
  const df = regressionDataset(rows, cols, 1);
  return {features: Array.from({length: cols}, (_, j) => copyCol(df.columns.byIndex(j))), y: copyCol(df.col('y0')!)};
}

async function warmUp(): Promise<void> {
  const {features, y} = regressionData(20, 2);
  await Promise.all(Array.from({length: cruxMl.size}, () => cruxFit('regression', 'linear', frame(features), y)));
}

category('Crux pool', () => {
  test('Parallel fits give the fits of a sequential run', async () => {
    const {features, target} = await iris();
    const fit = (k: number) => cruxFit('classification', 'random_forest', frame(features), target,
      {estimator: {n_estimators: 20, max_depth: 5}, seed: k});
    const sequential: Uint8Array[] = [];
    for (let k = 0; k < FITS; ++k)
      sequential.push(await fit(k));

    cruxTerminate();
    const pending = Promise.all(Array.from({length: FITS}, (_, k) => fit(k)));
    expect(cruxMl.workerCount, Math.min(FITS, cruxMl.size));
    const parallel = await pending;

    for (let k = 0; k < FITS; ++k) {
      expect(parallel[k].length, sequential[k].length, `fit ${k}: blob size`);
      const a = await predictColumn(frame(features), parallel[k]);
      const b = await predictColumn(frame(features), sequential[k]);
      for (let i = 0; i < a.length; ++i)
        expect(a.get(i), b.get(i), `fit ${k}, row ${i}`);
    }
  }, {timeout: TIMEOUT});

  test('A model stays on the worker that holds it', async () => {
    const {features, target} = await iris();
    const bytes = await cruxFit('classification', 'logistic', frame(features), target, {preprocess: 'standard'});
    const proto = CruxWorker.prototype as unknown as
      {post: (op: string, args: object, transfer?: Transferable[]) => Promise<any>};
    const post = proto.post;
    let registrations = 0;
    proto.post = function(op: string, args: object, transfer?: Transferable[]) {
      if (op === 'fnResource')
        ++registrations;
      return post.call(this, op, args, transfer);
    };
    try {
      cruxTerminate();
      await cruxPredict(frame(features), bytes);
      const unrelated = (k: number) => cruxFit('classification', 'decision_tree', frame(features), target, {seed: k});
      await Promise.all([unrelated(1), unrelated(2), cruxPredict(frame(features), bytes)]);
      expect(registrations, 1, 'the model must not be registered on another worker');
    } finally {
      proto.post = post;
    }
  }, {timeout: TIMEOUT});

  test('Terminate rejects every pending call and the next call works', async () => {
    const {features, y} = regressionData(200, 3);
    await warmUp();
    const calls = Array.from({length: 4}, (_, k) => cruxFit('regression', 'random_forest', frame(features), y,
      {estimator: {n_estimators: 200}, seed: k}));
    cruxTerminate();
    await expectRejected(calls.map((c, k): [string, () => Promise<unknown>] => [`call ${k}`, () => c]),
      (e) => e.message === 'Crux worker terminated');
    expect(cruxMl.workerCount, 0);

    await cruxFit('regression', 'linear', frame(features), y);
    expect(cruxMl.workerCount, 1);
  }, {timeout: TIMEOUT});

  test('Pool size', async () => {
    expect(cruxMl.size >= 1 && cruxMl.size <= 4, true, `size ${cruxMl.size}`);
    expect(new CruxModule('eda', 2).size, 2);
  });

  test('BENCH pool', async () => {
    const {features, y} = regressionData(5000, 10);
    const fit = () => cruxFit('regression', 'xgboost', frame(features), y, {estimator: {n_estimators: 100}});
    await warmUp();

    let t0 = performance.now();
    for (let k = 0; k < FITS; ++k)
      await fit();
    const sequential = performance.now() - t0;

    t0 = performance.now();
    await Promise.all(Array.from({length: FITS}, fit));
    const parallel = performance.now() - t0;

    console.log(`BENCH pool: sequential ${sequential.toFixed(0)} ms, parallel ${parallel.toFixed(0)} ms on ` +
      `${cruxMl.size} workers`);
    if (DG.Test.isInBenchmark && cruxMl.size > 1)
      expect(parallel < sequential, true, `parallel ${parallel} ms, sequential ${sequential} ms`);
  }, {timeout: 600000, benchmark: true});
});
