import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

// tests for dimensionality reduction

import {category, expect, test, uiSettled} from '@datagrok-libraries/test/src/test';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {KnownMetrics, NumberMetricsNames, StringMetricsNames, VectorMetricsNames}
  from '@datagrok-libraries/ml/src/typed-metrics';
import {multiColReduceDimensionality}
  from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {multiColWebGPUKNN, WEBGPUDISTANCE} from '@datagrok-libraries/math';
import {WEBGSLAGGREGATION} from '@datagrok-libraries/math/src/webGPU/multi-col-distances/webGPU-aggregation';
import {WebGPUUMAP} from '@datagrok-libraries/math/src/webGPU/umap';
import {getGPUDevice} from '@datagrok-libraries/math/src/webGPU/getGPUDevice';

const DEMOG_COLNAMES = {
  SUBJ: 'subj',
  STUDY: 'study',
  SITE: 'site',
  AGE: 'age',
  SEX: 'sex',
  RACE: 'race',
  DISEASE: 'disease',
  WEIGHT: 'weight',
  HEIGHT: 'height',
} as const;
category('Dimensionality reduction: UMAP', () => {
  test('Numeric column', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.AGE], DimReductionMethods.UMAP, [NumberMetricsNames.Difference]);
  }, {timeout: 30000});

  test('String column', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.SEX], DimReductionMethods.UMAP, [StringMetricsNames.Onehot]);
  }, {timeout: 30000});

  test('Numeric and string columns', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.SEX, DEMOG_COLNAMES.AGE], DimReductionMethods.UMAP,
      [StringMetricsNames.Onehot, NumberMetricsNames.Difference]);
  });

  test('All demog columns', async () => {
    const allDemogCols = grok.data.demo.demog(10).columns.toList()
      .filter((col) => Object.values(DEMOG_COLNAMES).includes(col.name as any)); ;
    const distFuncs = allDemogCols.map((col) => col.type === DG.COLUMN_TYPE.STRING ?
      StringMetricsNames.Onehot : NumberMetricsNames.Difference);
    const colNames = allDemogCols.map((col) => col.name);
    await testDimensionalityReductionUI( colNames, DimReductionMethods.UMAP, distFuncs);
  });
});

category('Dimensionality reduction: T-SNE', () => {
  test('Numeric column', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.AGE], DimReductionMethods.T_SNE, [NumberMetricsNames.Difference]);
  }, {timeout: 30000});

  test('String column', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.SEX], DimReductionMethods.T_SNE, [StringMetricsNames.Onehot]);
  }, {timeout: 30000});

  test('Numeric and string columns', async () => {
    await testDimensionalityReductionUI(
      [DEMOG_COLNAMES.SEX, DEMOG_COLNAMES.AGE], DimReductionMethods.T_SNE,
      [StringMetricsNames.Onehot, NumberMetricsNames.Difference]);
  });

  test('All demog columns', async () => {
    const allDemogCols = grok.data.demo.demog(10).columns.toList()
      .filter((col) => Object.values(DEMOG_COLNAMES).includes(col.name as any));
    const distFuncs = allDemogCols.map((col) => col.type === DG.COLUMN_TYPE.STRING ?
      StringMetricsNames.Onehot : NumberMetricsNames.Difference);
    const colNames = allDemogCols.map((col) => col.name);
    await testDimensionalityReductionUI(colNames, DimReductionMethods.T_SNE, distFuncs);
  });
});

category('Dimensionality reduction: UMAP random seed', () => {
  const NUMERIC_COLS = [DEMOG_COLNAMES.AGE, DEMOG_COLNAMES.HEIGHT, DEMOG_COLNAMES.WEIGHT];

  test('GPU fit: same seed, same coordinates', async () => {
    if (!await getGPUDevice())
      return 'No WebGPU device in this browser: nothing to check';
    const df = grok.data.demo.demog(1000);
    const cols = NUMERIC_COLS.map((name) => df.col(name)!.toList());
    const vectors = Array.from({length: df.rowCount}, (_, i) => Float32Array.from(cols, (col) => col[i]));
    const knn = (await multiColWebGPUKNN([vectors], 15, [WEBGPUDISTANCE.EUCLIDEAN], WEBGSLAGGREGATION.MANHATTAN,
      [1], [{}]))!;
    const fit = async (seed: number) => {
      const umap = new WebGPUUMAP(vectors.length, {random: seededRandom(seed)});
      umap.setPrecomputedKNN(knn.knnIndexes, knn.knnDistances);
      return (await umap.fit())!;
    };
    const first = await fit(1);
    expect(sameCoordinates(await fit(1), first), true, 'seed 1 gave other coordinates on the second run');
    expect(sameCoordinates(await fit(2), first), false, 'seeds 1 and 2 gave the same coordinates');
  }, {timeout: 60000});

  // what the Reduce Dimensionality dialog runs on numeric columns, in the dimensionality reduction worker
  test('Same seed, same coordinates, on the GPU and on the CPU', async () => {
    const df = grok.data.demo.demog(1000);
    const cols = NUMERIC_COLS.map((name) => df.col(name)!);
    const embed = async (randomSeed: string, useWebGPU: boolean) => {
      const before = df.columns.names();
      await multiColReduceDimensionality(df, cols, DimReductionMethods.UMAP,
        cols.map(() => NumberMetricsNames.Difference), cols.map(() => 1), cols.map(() => undefined), 'EUCLIDEAN',
        false, false, {preprocessingFuncArgs: cols.map(() => ({})), randomSeed, useWebGPU},
        {[BYPASS_LARGE_DATA_WARNING]: true}, null, {}, VectorMetricsNames.Euclidean);
      return df.columns.names().filter((name) => !before.includes(name)).map((name) => df.col(name)!.toList());
    };
    const gpu = await getGPUDevice() != null;
    const firstRuns: number[][][] = [];
    for (const useWebGPU of gpu ? [true, false] : [false]) {
      const device = useWebGPU ? 'GPU' : 'CPU';
      const first = await embed('42', useWebGPU);
      expect(first.length, 2, `${device}: no embedding columns added`);
      expect(sameCoordinates(await embed('42', useWebGPU), first), true,
        `${device}: seed 42 gave other coordinates on the second run`);
      expect(sameCoordinates(await embed('43', useWebGPU), first), false,
        `${device}: seeds 42 and 43 gave the same coordinates`);
      firstRuns.push(first);
    }
    // a failed GPU fit falls back to the CPU fit, which keeps to the seed too: the checks above would pass without it
    if (gpu)
      expect(sameCoordinates(firstRuns[0], firstRuns[1]), false, 'the GPU run gave the CPU coordinates');
  }, {timeout: 120000});
});

category('DBSCAN', () => {
  test('two dense groups and a far point', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromFloat32Array('x', new Float32Array([0, 0.01, 0, 0.01, 5, 5.01, 5, 5.01, 20])),
      DG.Column.fromFloat32Array('y', new Float32Array([0, 0, 0.01, 0.01, 5, 5, 5.01, 5.01, 20])),
    ]);
    const clusters: DG.Column = await grok.functions.call('EDA:DBSCAN',
      {df, xCol: df.col('x'), yCol: df.col('y'), epsilon: 0.02, minPts: 4});
    expect(clusters.name, 'Cluster (DBSCAN)');
    expect(clusters.toList().join(', '), '1, 1, 1, 1, 2, 2, 2, 2, -1');
  });
});

/** mulberry32: numbers in [0, 1) from a seed */
function seededRandom(seed: number): () => number {
  return () => {
    seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

function sameCoordinates(a: ArrayLike<number>[], b: ArrayLike<number>[]): boolean {
  if (a.length !== b.length)
    return false;
  for (let d = 0; d < a.length; d++) {
    if (a[d].length !== b[d].length)
      return false;
    for (let i = 0; i < a[d].length; i++) {
      if (a[d][i] !== b[d][i])
        return false;
    }
  }
  return true;
}

async function testDimensionalityReductionUI(
  columns: string[], methodName: DimReductionMethods, metrics: KnownMetrics[],
) {
  const df = grok.data.demo.demog(100);
  const _tv = grok.shell.addTableView(df);
  const dimRedResult = await multiColReduceDimensionality(
    df, columns.map((c) => df.col(c)!), methodName, metrics,
    columns.map(() => 1), columns.map(() => undefined),
    'EUCLIDEAN', true, true, {preprocessingFuncArgs: columns.map(() => ({}))});
  expect(!!dimRedResult, true, 'No scatterplot returned');
  const addedEmbeddingsCols = df.columns.names().filter((c) => c.toLowerCase().startsWith('embed'));
  expect(addedEmbeddingsCols.length, 2, 'Wrong number of embeddings added');
  const clusterColName = df.columns.names().find((c) => c.toLowerCase().startsWith('cluster'));
  expect(!!clusterColName, true, 'No cluster column added');
  for (const embedColName of addedEmbeddingsCols) {
    const c = df.col(embedColName)!;
    expect(new Array(c.length).fill(null).every((_, i) => !c.isNone(i) && !isNaN(c.get(i))), true,
      'Embedding column has null-ish values');
  }
  await uiSettled(`${methodName} of ${columns}`);
}
