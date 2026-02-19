import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

// tests for dimensionality reduction

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {KnownMetrics, NumberMetricsNames, StringMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {multiColReduceDimensionality}
  from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';

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
  await new Promise((resolve) => setTimeout(resolve, 500));
}
