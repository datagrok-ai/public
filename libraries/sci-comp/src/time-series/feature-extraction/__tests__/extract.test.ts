import {extractFeatures} from '../extract';
import {
  constantData, rampData, alternatingData, stepData, spikeData,
  negativeData, nearConstantOutlierData,
  constantExpected, rampExpected, alternatingExpected, stepExpected,
  spikeExpected, negativeExpected, nearConstantOutlierExpected,
  toleranceFor, assertClose, assertExact,
  makeDataFrame, makeMultiSampleDataFrame,
} from './helpers';

function assertFeatures(
  result: ReturnType<typeof extractFeatures>,
  sampleIdx: number,
  colName: string,
  expected: Record<string, number>,
) {
  for (const [featureKey, expectedVal] of Object.entries(expected)) {
    const fullName = `${colName}__${featureKey}`;
    const col = result.columns.find((c) => c.name === fullName);
    if (!col)
      throw new Error(`Missing feature column: ${fullName}`);

    const actual = col.data[sampleIdx];
    const tol = toleranceFor(fullName);
    if (tol === 'exact')
      assertExact(actual, expectedVal, fullName);
    else
      assertClose(actual, expectedVal, tol, fullName);
  }
}

describe('extractFeatures — single sample, all features', () => {
  const datasets: [string, Float64Array, Record<string, number>][] = [
    ['constant', constantData, constantExpected],
    ['ramp', rampData, rampExpected],
    ['alternating', alternatingData, alternatingExpected],
    ['step', stepData, stepExpected],
    ['spike', spikeData, spikeExpected],
    ['negative', negativeData, negativeExpected],
    ['near_constant_outlier', nearConstantOutlierData, nearConstantOutlierExpected],
  ];

  test.each(datasets)('%s', (_name, data, expected) => {
    const df = makeDataFrame(data, 'value');
    const result = extractFeatures(df);

    expect(result.nSamples).toBe(1);
    expect(result.sampleIds[0]).toBe(1);

    assertFeatures(result, 0, 'value', expected);
  });
});

describe('extractFeatures — multiple samples', () => {
  it('3 samples with different data', () => {
    const df = makeMultiSampleDataFrame([
      {id: 1, data: rampData},
      {id: 2, data: constantData},
      {id: 3, data: alternatingData},
    ]);
    const result = extractFeatures(df);

    expect(result.nSamples).toBe(3);
    expect(Array.from(result.sampleIds)).toEqual([1, 2, 3]);

    assertFeatures(result, 0, 'value', rampExpected);
    assertFeatures(result, 1, 'value', constantExpected);
    assertFeatures(result, 2, 'value', alternatingExpected);
  });
});

describe('extractFeatures — multiple columns', () => {
  it('2 columns produce correct prefixes and double feature count', () => {
    const n = 10;
    const df = {
      ids: new Uint32Array(n).fill(1),
      time: new Float64Array(Array.from({length: n}, (_, i) => i)),
      rowCount: n,
      columns: [
        {name: 'pH', data: rampData},
        {name: 'concentration', data: constantData},
      ],
    };
    const result = extractFeatures(df);

    expect(result.nSamples).toBe(1);

    const phCols = result.columns.filter((c) => c.name.startsWith('pH__'));
    const concCols = result.columns.filter((c) => c.name.startsWith('concentration__'));
    expect(phCols.length).toBeGreaterThan(0);
    expect(concCols.length).toBeGreaterThan(0);
    expect(phCols.length).toBe(concCols.length);

    assertFeatures(result, 0, 'pH', rampExpected);
    assertFeatures(result, 0, 'concentration', constantExpected);
  });
});

describe('extractFeatures — validation', () => {
  it('validates when options.validate=true', () => {
    const df = makeDataFrame(new Float64Array([1, 2, NaN, 4, 5]), 'value');
    expect(() => extractFeatures(df, {validate: true})).toThrow(/NaN/);
  });
});
