import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  category,
  test,
  expect,
  expectArray,
  assure,
} from '@datagrok-libraries/test/src/test';
import {DistanceMatrix} from '@datagrok-libraries/ml/src/distance-matrix';
import {mapToFixed} from './utils/array-utils';

const validDistanceMatrix5x5 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10];
const validDistanceMatrixSize = 6;
const invalidDistanceMatrix = [1, 2, 3, 4];
const invalidDistanceMatrixSize = 5;
const arrayToCalcDistance = [1, 2, 3, 4, 5];

const validDistanceMatrix5x5Squared = [0, 4, 9, 16, 25, 36, 49, 64, 81, 100];

const validDistanceMatrix5x5Normalized = [
  0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
];

const validDistanceMatrix5x5SquareRoot = [
  0, 1.4142135623730951, 1.7320508075688772, 2, 2.23606797749979,
  2.449489742783178, 2.6457513110645907, 2.8284271247461903, 3,
  3.1622776601683795,
];

const validDistanceMatrix5x5Double = [0, 4, 6, 8, 10, 12, 14, 16, 18, 20];

const calculatedDistanceMatrix = [1, 2, 3, 4, 1, 2, 3, 1, 2, 1];

category('DistanceMatrix', () => {
  test('initializeWithValidData', async () => {
    const dm = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    expectArray(dm.data, validDistanceMatrix5x5);
  });

  test('initializeWithValidSize', async () => {
    let error: any = null;
    try {
      const dm = new DistanceMatrix(undefined, validDistanceMatrixSize);
      expect(dm.size, validDistanceMatrixSize);
    } catch (e) {
      error = e;
    } finally {
      expect(error, null);
    }
  });

  test('initializeWithInvalidData', async () => {
    let error: any = null;
    try {
      const _dm = new DistanceMatrix(new Float32Array(invalidDistanceMatrix));
    } catch (e) {
      error = e;
    } finally {
      assure.notNull(error);
    }
  });

  test('initializeWithInvalidSize', async () => {
    let error: any = null;
    try {
      const _dm = new DistanceMatrix(undefined, invalidDistanceMatrixSize);
    } catch (e) {
      error = e;
    } finally {
      expect(error, null);
    }
  });

  test('squareDistanceMatrix', async () => {
    const dm = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    dm.square();
    expectArray(mapToFixed(dm.data), mapToFixed(validDistanceMatrix5x5Squared));
  });

  test('normalizeDistanceMatrix', async () => {
    const dm = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    dm.normalize();
    expectArray(
      mapToFixed(dm.data),
      mapToFixed(validDistanceMatrix5x5Normalized),
    );
  });

  test('squareRootDistanceMatrix', async () => {
    const dm = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    dm.sqrt();
    expectArray(
      mapToFixed(dm.data),
      mapToFixed(validDistanceMatrix5x5SquareRoot),
    );
  });

  test('addDistanceMatrix', async () => {
    const dm = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    const dm2 = new DistanceMatrix(new Float32Array(validDistanceMatrix5x5));
    dm.add(dm2);
    expectArray(mapToFixed(dm.data), mapToFixed(validDistanceMatrix5x5Double));
  });

  test('calcDistanceMatrixNumeric', async () => {
    const dm = DistanceMatrix.calc(arrayToCalcDistance, (a, b) =>
      Math.abs(a - b),
    );
    expectArray(mapToFixed(dm.data), mapToFixed(calculatedDistanceMatrix));
  });
});
