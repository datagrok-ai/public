// Tests for Pareto Front Computations
// Performance tests for the Pareto optimality algorithm

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {getParetoMask} from '../pareto-optimization/pareto-computations';
import {OPT_TYPE, NumericArray} from '../pareto-optimization/defs';

const TIMEOUT = 5000;

// Test dataset sizes 
const ROWS_COUNT = 1000000;
const M = 1000000;
const COLS_COUNT = 2;
const suffix = M < 1e6 ? 'K' : 'M';
const DATASET_SIZE_LABEL = `${ROWS_COUNT / M}${suffix} points, ${COLS_COUNT}D`;

/** Generates synthetic numeric data for Pareto front testing */
function generateSyntheticData(nPoints: number, nDims: number, seed: number = 42): NumericArray[] {
  const data: NumericArray[] = [];

  // Simple deterministic pseudo-random generator for reproducibility
  let rng = seed;
  const random = () => {
    rng = (rng * 1664525 + 1013904223) % 4294967296;
    return rng / 4294967296;
  };

  for (let d = 0; d < nDims; d++) {
    const column = new Float32Array(nPoints);
    for (let i = 0; i < nPoints; i++) {
      // Generate values with some correlation to create realistic Pareto fronts
      column[i] = random() * 100 + (d * 10);
    }
    data.push(column);
  }

  return data;
}

/** Generates optimization sense array */
function generateSense(nDims: number, pattern: 'all-min' | 'all-max' | 'mixed'): OPT_TYPE[] {
  const sense: OPT_TYPE[] = [];

  for (let d = 0; d < nDims; d++) {
    if (pattern === 'all-min') {
      sense.push(OPT_TYPE.MIN);
    } else if (pattern === 'all-max') {
      sense.push(OPT_TYPE.MAX);
    } else {
      // Mixed: alternate between MIN and MAX
      sense.push(d % 2 === 0 ? OPT_TYPE.MIN : OPT_TYPE.MAX);
    }
  }

  return sense;
}

/** Generates null indices set for testing missing value handling */
function generateNullIndices(nPoints: number, nullRatio: number): Set<number> {
  const nullCount = Math.floor(nPoints * nullRatio);
  const nullIndices = new Set<number>();

  // Distribute null indices evenly
  const step = Math.floor(nPoints / nullCount);
  for (let i = 0; i < nullCount; i++) {
    nullIndices.add(i * step);
  }

  return nullIndices;
}

/** Validates Pareto mask result */
function validateParetoMask(mask: boolean[], nPoints: number): void {
  if (mask.length !== nPoints) {
    throw new Error(`Invalid mask length: expected ${nPoints}, got ${mask.length}`);
  }

  const optimalCount = mask.filter(x => x).length;
  if (optimalCount === 0) {
    throw new Error('No optimal points found');
  }

  if (optimalCount === nPoints) {
    grok.shell.warning('All points are optimal - data may be degenerate');
  }
}

category('Pareto optimization', () => {

  test(`Performance: ${DATASET_SIZE_LABEL}`, async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data = generateSyntheticData(ROWS_COUNT, COLS_COUNT);
      const sense = generateSense(COLS_COUNT, 'mixed');
      mask = getParetoMask(data, sense, ROWS_COUNT);
      validateParetoMask(mask, ROWS_COUNT);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  // Tests for different optimization patterns
  test(`Performance: ${DATASET_SIZE_LABEL}, all minimize`, async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data = generateSyntheticData(ROWS_COUNT, COLS_COUNT);
      const sense = generateSense(COLS_COUNT, 'all-min');
      mask = getParetoMask(data, sense, ROWS_COUNT);
      validateParetoMask(mask, ROWS_COUNT);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  test(`Performance: ${DATASET_SIZE_LABEL}, all maximize`, async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data = generateSyntheticData(ROWS_COUNT, COLS_COUNT);
      const sense = generateSense(COLS_COUNT, 'all-max');
      mask = getParetoMask(data, sense, ROWS_COUNT);
      validateParetoMask(mask, ROWS_COUNT);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  // Tests with missing values
  test(`Performance: ${DATASET_SIZE_LABEL} with 10% null indices`, async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data = generateSyntheticData(ROWS_COUNT, COLS_COUNT);
      const sense = generateSense(COLS_COUNT, 'mixed');
      const nullIndices = generateNullIndices(ROWS_COUNT, 0.1);
      mask = getParetoMask(data, sense, ROWS_COUNT, nullIndices);
      validateParetoMask(mask, ROWS_COUNT);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  test(`Performance: ${DATASET_SIZE_LABEL} with 25% null indices`, async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data = generateSyntheticData(ROWS_COUNT, COLS_COUNT);
      const sense = generateSense(COLS_COUNT, 'mixed');
      const nullIndices = generateNullIndices(ROWS_COUNT, 0.25);
      mask = getParetoMask(data, sense, ROWS_COUNT, nullIndices);
      validateParetoMask(mask, ROWS_COUNT);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  // Edge cases
  test('Edge case: Empty dataset', async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data: NumericArray[] = [new Float32Array(0), new Float32Array(0)];
      const sense = generateSense(COLS_COUNT, 'mixed');
      mask = getParetoMask(data, sense, 0);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(mask!.length, 0, 'Empty dataset should return empty mask');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  test('Edge case: Single point', async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const data: NumericArray[] = [new Float32Array([1.0]), new Float32Array([2.0])];
      const sense = generateSense(COLS_COUNT, 'mixed');

      mask = getParetoMask(data, sense, 1);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(mask!.length, 1, 'Single point dataset should return mask with one element');
    expect(mask![0], true, 'Single point should be optimal');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});

  test('Edge case: All identical points', async () => {
    let mask: boolean[] | null = null;
    let error: Error | null = null;

    try {
      const nPoints = 100;
      const data: NumericArray[] = [
        new Float32Array(nPoints).fill(5.0),
        new Float32Array(nPoints).fill(10.0),
      ];
      const sense = generateSense(COLS_COUNT, 'mixed');

      mask = getParetoMask(data, sense, nPoints);
    } catch (e) {
      error = e as Error;
      grok.shell.error(error.message);
    }

    expect(mask !== null, true, 'Failed to compute Pareto mask');
    expect(mask!.length, 100, 'Should return mask with correct length');
    const optimalCount = mask!.filter(x => x).length;
    expect(optimalCount > 0, true, 'At least some identical points should be optimal');
    expect(error === null, true, error?.message ?? '');
  }, {timeout: TIMEOUT});
});
