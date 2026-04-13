/**
 * Example: extracting features from a single time series (1 sample, 1 column).
 *
 * Run: npx tsx src/time-series/feature-extraction/examples/single-sample.ts
 */
import {extractFeatures} from '..';
import type {TimeSeriesDataFrame} from '..';

/* ================================================================== */
/*  Build the input DataFrame                                          */
/* ================================================================== */

const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
  time: new Float64Array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
  rowCount: 10,
  columns: [{
    name: 'temperature',
    data: new Float64Array([20.1, 20.5, 21.0, 20.8, 21.3, 21.7, 22.0, 21.5, 21.8, 22.2]),
  }],
};

/* ================================================================== */
/*  Print the input                                                    */
/* ================================================================== */

console.log('Input DataFrame:');
console.log(`  rowCount: ${df.rowCount}`);
console.log(`  columns:  [${df.columns.map((c) => c.name).join(', ')}]`);
console.log('');
console.log('  id  time  temperature');
console.log('  --  ----  -----------');
for (let i = 0; i < df.rowCount; i++)
  console.log(`  ${df.ids[i]}   ${df.time[i].toFixed(0).padStart(4)}  ${df.columns[0].data[i].toFixed(1)}`);

/* ================================================================== */
/*  Extract features                                                   */
/* ================================================================== */

const result = extractFeatures(df);

console.log(`\nFeatureMatrix: ${result.nSamples} sample × ${result.columns.length} features\n`);
for (const col of result.columns)
  console.log(`  ${col.name.padEnd(52)} ${col.data[0]}`);
