import {validate} from '../dataframe';
import {extractFeatures} from '../extract';
import {TimeSeriesDataFrame} from '../types';

function cleanDf(): TimeSeriesDataFrame {
  return {
    ids: new Uint32Array([1, 1, 1, 1, 1]),
    time: new Float64Array([0, 1, 2, 3, 4]),
    rowCount: 5,
    columns: [{name: 'value', data: new Float64Array([1, 2, 3, 4, 5])}],
  };
}

describe('validate', () => {
  it('passes on clean data', () => {
    expect(() => validate(cleanDf())).not.toThrow();
  });

  it('throws on NaN in Float64Array column', () => {
    const df = cleanDf();
    (df.columns[0].data as Float64Array)[3] = NaN;
    expect(() => validate(df)).toThrow(/NaN.*index 3/);
  });

  it('throws on Inf', () => {
    const df = cleanDf();
    (df.columns[0].data as Float64Array)[0] = Infinity;
    expect(() => validate(df)).toThrow(/non-finite/);
  });

  it('skips NaN check for Int32Array', () => {
    const df: TimeSeriesDataFrame = {
      ids: new Uint32Array([1, 1, 1]),
      time: new Float64Array([0, 1, 2]),
      rowCount: 3,
      columns: [{name: 'counts', data: new Int32Array([1, 2, 3])}],
    };
    expect(() => validate(df)).not.toThrow();
  });

  it('throws on mismatched lengths', () => {
    const df = cleanDf();
    df.columns[0] = {name: 'value', data: new Float64Array([1, 2, 3, 4])};
    expect(() => validate(df)).toThrow(/data\.length.*!== rowCount/);
  });

  it('throws on non-contiguous ids', () => {
    const df: TimeSeriesDataFrame = {
      ids: new Uint32Array([1, 2, 1, 2]),
      time: new Float64Array([0, 0, 1, 1]),
      rowCount: 4,
      columns: [{name: 'value', data: new Float64Array([1, 2, 3, 4])}],
    };
    expect(() => validate(df)).toThrow(/contiguous/);
  });

  it('not called when options.validate=false', () => {
    const df = cleanDf();
    (df.columns[0].data as Float64Array)[2] = NaN;
    // Should not throw — NaN propagates to results
    expect(() => extractFeatures(df, {validate: false})).not.toThrow();
  });
});
