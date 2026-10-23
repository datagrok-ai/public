import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {LEGACY_MODEL_ERROR, cruxPredict} from '../crux/crux-model';
import {numericFeatures} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const CARS = 'System:DemoFiles/cars.csv';
const FOLDER = 'System:AppData/Eda/legacy-models';
const TIMEOUT = 30000;

type Fixture = {file: string, dataset: string, target: string, exclude: string[]};

const CARS_REGRESSION = {dataset: CARS, target: 'price', exclude: []};
const IRIS_CLASSIFICATION = {dataset: IRIS, target: 'Species', exclude: ['col 1']};

const FIXTURES: Fixture[] = [
  {file: 'linear-regression', ...CARS_REGRESSION},
  {file: 'pls-regression', ...CARS_REGRESSION},
  {file: 'softmax', ...IRIS_CLASSIFICATION},
  {file: 'xgboost', ...IRIS_CLASSIFICATION},
  {file: 'svm', ...IRIS_CLASSIFICATION},
];

category('Legacy models', () => {
  for (const fixture of FIXTURES) {
    test(fixture.file, async () => {
      const bytes = await grok.dapi.files.readAsBytes(`${FOLDER}/${fixture.file}.bin`);
      const df = await grok.dapi.files.readCsv(fixture.dataset);
      const features = numericFeatures(df, fixture.target, fixture.exclude);
      await expectExceptionAsync(async () => {
        await cruxPredict(DG.DataFrame.fromColumns(features.toList()), bytes);
      }, (e) => e.message === LEGACY_MODEL_ERROR);
    }, {timeout: TIMEOUT});
  }
});
