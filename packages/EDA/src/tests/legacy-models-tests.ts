import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, expectExceptionAsync, test} from '@datagrok-libraries/test/src/test';

import {SoftmaxClassifier} from '../softmax-classifier';
import {getPredictionByLinearRegression, unpackLinearRegressionModel} from '../regression';
import {PlsModel} from '../pls/pls-ml';
import {SVM} from '../svm';
import {XGBooster} from '../xgbooster';
import {initXgboost} from '../../wasm/xgbooster';
import {initSvm} from '../../wasm/svm';
import {LEGACY_MODEL_ERROR, cruxPredict} from '../crux/crux-model';
import {accuracy, numericFeatures, r2} from './utils';

const IRIS = 'System:DemoFiles/iris.csv';
const CARS = 'System:DemoFiles/cars.csv';
const FOLDER = 'System:AppData/Eda/legacy-models';
const TIMEOUT = 30000;

type Fixture = {
  file: string,
  dataset: string,
  target: string,
  exclude: string[],
  floor: number,
  score: (target: DG.Column, pred: DG.Column) => number,
  predict: (bytes: Uint8Array, features: DG.ColumnList) => DG.Column,
};

const CARS_REGRESSION = {dataset: CARS, target: 'price', exclude: [], floor: 0.4, score: r2};
const IRIS_CLASSIFICATION = {dataset: IRIS, target: 'Species', exclude: ['col 1'], floor: 0.8, score: accuracy};

const FIXTURES: Fixture[] = [
  {file: 'linear-regression', ...CARS_REGRESSION, predict: (bytes, features) => {
    const {params, names} = unpackLinearRegressionModel(bytes);
    return getPredictionByLinearRegression(features, params, names);
  }},
  {file: 'pls-regression', ...CARS_REGRESSION, predict: (bytes, features) => new PlsModel(bytes).predict(features)},
  {file: 'softmax', ...IRIS_CLASSIFICATION,
    predict: (bytes, features) => new SoftmaxClassifier(undefined, bytes).predict(features)},
  {file: 'xgboost', ...IRIS_CLASSIFICATION, predict: (bytes, features) => new XGBooster(bytes).predict(features)},
  {file: 'svm', ...IRIS_CLASSIFICATION, predict: (bytes, features) => new SVM(bytes).predict(features)},
];

category('Legacy models', () => {
  before(async () => {
    await initXgboost();
    await initSvm();
  });

  for (const fixture of FIXTURES) {
    test(fixture.file, async () => {
      const bytes = await grok.dapi.files.readAsBytes(`${FOLDER}/${fixture.file}.bin`);
      const df = await grok.dapi.files.readCsv(fixture.dataset);
      const features = numericFeatures(df, fixture.target, fixture.exclude);
      const score = fixture.score(df.col(fixture.target)!, fixture.predict(bytes, features));
      expect(score > fixture.floor, true, `${fixture.file}: score ${score} (expected > ${fixture.floor})`);
      await expectExceptionAsync(async () => {
        await cruxPredict(DG.DataFrame.fromColumns(features.toList()), bytes);
      }, (e) => e.message === LEGACY_MODEL_ERROR);
    }, {timeout: TIMEOUT});
  }
});
