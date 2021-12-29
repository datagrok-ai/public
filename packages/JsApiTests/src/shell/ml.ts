import {after, before, category, expect, test} from "../test";
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('ML', () => {
  // test('Apply Model', async () => {
  //   const data = grok.data.demo.demog();
  //   const resultDf = await grok.ml.applyModel(
  //     // 'Demo:PredictSexByBasicDemographics',
  //     'Donufriienko:PredictSEXByAGEHEIGHTWEIGHTUsingDistributedRandomForest',
  //     data,
  //     {'SEX': 'SEX', 'AGE': 'AGE', 'HEIGHT': 'HEIGHT', 'WEIGHT': 'WEIGHT'},
  //     true,
  //   );

  //   expect((resultDf.columns as DG.ColumnList).names().includes('outcome'), true);
  // });

  test('Cluster', async () => {
    const data = await grok.data.loadTable('https://public.datagrok.ai/demo/xclara.csv');
    const resultDf = await grok.ml.cluster(data, ["V1", "V2"], 3);

    expect(resultDf.getCol('clusters').categories.length, 3);
  });

  test('Missing Values Imputation', async () => {
    const data = await grok.data.loadTable('https://public.datagrok.ai/demo/demog.csv');
    const resultDf = await grok.ml.missingValuesImputation(
      data, ['age', 'height', 'weight'], ['age', 'height', 'weight'], 5,
    );

    expect(resultDf.getCol('age').isNone(1192), false);
    expect(resultDf.getCol('height').isNone(2220), false);
    expect(resultDf.getCol('weight').isNone(2221), false);
  });

  test('PCA', async () => {
    const data = await grok.data.loadTable('https://public.datagrok.ai/demo/cars.csv');
    const resultDf = await grok.ml.pca(
      data, ["wheel.base", "length", "width", "height", "city.mpg", "price"], 2, true, true,
    );
    
    expect((resultDf.columns as DG.ColumnList).names().includes('PCA0'), true);
    expect((resultDf.columns as DG.ColumnList).names().includes('PCA1'), true);
  });

  test('Random Data', async () => {
    const seed = 42;
    const data = grok.data.demo.demog();

    await grok.ml.randomData(data, 'normal', {sd: 3.0, mean: 1.0}, seed);
    await grok.ml.randomData(data, 'uniform', {min: 0.0, max: 1.0}, seed);
    await grok.ml.randomData(data, 'binomial', {size: 100, prob: 0.7}, seed);

    //too naive?
    const normalCol = data.getCol('normal');
    expect(0.9 <= normalCol.stats.avg && normalCol.stats.avg <= 1.1, true);
    expect(2.9 <= normalCol.stats.stdev && normalCol.stats.stdev <= 3.1, true);
    const uniformCol = data.getCol('uniform');
    expect(uniformCol.stats.min >= 0.0 && uniformCol.stats.max <= 1.0, true);
    const binomialCol = data.getCol('binomial');
    expect(binomialCol.stats.min >= 0 && binomialCol.stats.max <= 100, true);
    expect(68 <= binomialCol.stats.avg && binomialCol.stats.avg <= 72, true);
  });
});

