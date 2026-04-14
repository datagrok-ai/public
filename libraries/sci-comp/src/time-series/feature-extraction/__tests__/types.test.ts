import {extractFeatures} from '../extract';
import {TimeSeriesDataFrame} from '../types';

function makeTypedDf(
  data: Int32Array | Uint32Array | Float32Array | Float64Array,
): TimeSeriesDataFrame {
  const n = data.length;
  return {
    ids: new Uint32Array(n).fill(1),
    time: new Float64Array(Array.from({length: n}, (_, i) => i)),
    rowCount: n,
    columns: [{name: 'value', data}],
  };
}

describe('input type compatibility', () => {
  const intData = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
  const f64 = new Float64Array(intData);
  const refResult = extractFeatures(makeTypedDf(f64));

  it('Int32Array produces same features as Float64Array', () => {
    const result = extractFeatures(makeTypedDf(new Int32Array(intData)));
    for (let i = 0; i < refResult.columns.length; i++) {
      const diff = Math.abs(result.columns[i].data[0] - refResult.columns[i].data[0]);
      expect(diff).toBeLessThan(1e-10);
    }
  });

  it('Uint32Array produces same features as Float64Array', () => {
    const result = extractFeatures(makeTypedDf(new Uint32Array(intData)));
    for (let i = 0; i < refResult.columns.length; i++) {
      const diff = Math.abs(result.columns[i].data[0] - refResult.columns[i].data[0]);
      expect(diff).toBeLessThan(1e-10);
    }
  });

  it('Float32Array produces features within Float32 precision', () => {
    const f32Data = new Float32Array([1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.1]);
    const f64Data = new Float64Array(f32Data); // promotes to f64
    const r32 = extractFeatures(makeTypedDf(f32Data));
    const r64 = extractFeatures(makeTypedDf(f64Data));

    for (let i = 0; i < r64.columns.length; i++) {
      const diff = Math.abs(r32.columns[i].data[0] - r64.columns[i].data[0]);
      expect(diff).toBeLessThan(1e-5);
    }
  });
});
