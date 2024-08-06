// Tests of numerical methods

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';

import {ODEs} from '../solver-tools/solver-defs';
import {mrt} from '../solver-tools/mrt-method';
import {ros3prw} from '../solver-tools/ros3prw-method';
import {ros34prw} from '../solver-tools/ros34prw-method';

import {evaluateMethod} from './testing-utils';

const TIMEOUT = 4000;
const TINY = 0.1;

const methods = new Map([
  ['MRT', mrt],
  ['ROSP3PRw', ros3prw],
  ['ROSP34PRw', ros34prw],
]);

category('Correctness', () => {
  methods.forEach((method, name) => test(`The ${name} method`, async () => {
    const error = evaluateMethod(method);
    console.log(`Error (${name}): ${error}`);

    expect(
      error < TINY,
      true,
      `The ${name} method failed, too big error: ${error}; expected: < ${TINY}`,
    );
  }, {timeout: TIMEOUT}));

  /*test('Correctness', async () => {
    // Prepare data
    const df = classificationDataset(ROWS_K, MIN_COLS, true);
    const features = df.columns;
    const target = features.byIndex(MIN_COLS);
    features.remove(target.name);

    // Fit & pack trained model
    const model = new SoftmaxClassifier({
      classesCount: target.categories.length,
      featuresCount: features.length,
    });

    await model.fit(features, target);
    const modelBytes = model.toBytes();

    // Unpack & apply model
    const unpackedModel = new SoftmaxClassifier(undefined, modelBytes);
    const prediction = unpackedModel.predict(features);

    // Evaluate accuracy
    const acc = accuracy(target, prediction);
    expect(
      acc > MIN_ACCURACY,
      true,
      `Softmax failed, too small accuracy: ${acc}; expected: <= ${MIN_ACCURACY}`,
    );
  }, {timeout: TIMEOUT});*/
}); // Softmax
