import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

category('ML', () => {
  test('Cluster', async () => {
    const data = await _package.files.readCsv('datasets/xclara.csv');
    const resultDf = await grok.ml.cluster(data, ['V1', 'V2'], 3);
    expect(resultDf.getCol('clusters').categories.length, 3);
  });

  test('PCA', async () => {
    const data = await _package.files.readCsv('datasets/cars.csv');
    const resultDf = await grok.ml.pca(
      data, ['wheel.base', 'length', 'width', 'height', 'city.mpg', 'price'], 2, true, true,
    );
    expect((resultDf.columns as DG.ColumnList).names().includes('PCA0'), true);
    expect((resultDf.columns as DG.ColumnList).names().includes('PCA1'), true);
  });
}, {owner: 'kamelichev@datagrok.ai'});
