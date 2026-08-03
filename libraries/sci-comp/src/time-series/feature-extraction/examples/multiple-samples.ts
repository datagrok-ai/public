/**
 * Example: comparing features across multiple samples (experiments).
 *
 * Three pH time series are packed into one DataFrame with different ids.
 * Features are extracted for all samples at once and compared side-by-side.
 *
 * Run: npx tsx src/time-series/feature-extraction/examples/multiple-samples.ts
 */
import {extractFeatures} from '..';
import type {TimeSeriesDataFrame} from '..';

/* ================================================================== */
/*  Build the input DataFrame                                          */
/* ================================================================== */

const df: TimeSeriesDataFrame = {
  ids: new Uint32Array([1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]),
  time: new Float64Array([0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4]),
  rowCount: 15,
  columns: [{
    name: 'pH',
    data: new Float64Array([
      // Sample 1: steadily rising
      7.0, 7.1, 7.2, 7.3, 7.4,
      // Sample 2: constant
      7.0, 7.0, 7.0, 7.0, 7.0,
      // Sample 3: noisy fluctuations
      6.8, 7.3, 6.9, 7.5, 7.1,
    ]),
  }],
};

/* ================================================================== */
/*  Print the input                                                    */
/* ================================================================== */

console.log('Input DataFrame:');
console.log(`  rowCount: ${df.rowCount}`);
console.log(`  columns:  [${df.columns.map((c) => c.name).join(', ')}]`);
console.log('');
console.log('  id  time    pH');
console.log('  --  ----  ----');
for (let i = 0; i < df.rowCount; i++)
  console.log(`  ${df.ids[i]}   ${df.time[i].toFixed(0).padStart(4)}  ${df.columns[0].data[i].toFixed(1)}`);

/* ================================================================== */
/*  Extract features                                                   */
/* ================================================================== */

const result = extractFeatures(df);

console.log(`\nFeatureMatrix: ${result.nSamples} samples × ${result.columns.length} features`);
console.log(`Sample IDs: [${Array.from(result.sampleIds).join(', ')}]`);

/* ================================================================== */
/*  Compare key features across samples                                */
/* ================================================================== */

const keyFeatures = [
  'pH__mean',
  'pH__standard_deviation',
  'pH__linear_trend__attr_slope',
  'pH__linear_trend__attr_rvalue',
  'pH__has_duplicate',
  'pH__cid_ce__normalize_false',
  'pH__longest_strike_above_mean',
  'pH__number_crossing_m__m_0',
];

console.log('\nSide-by-side comparison (sample 1 = rising, 2 = constant, 3 = noisy):');
console.log('');
const header = '  ' + 'feature'.padEnd(36) +
  'sample 1'.padStart(10) + 'sample 2'.padStart(10) + 'sample 3'.padStart(10);
console.log(header);
console.log('  ' + '-'.repeat(header.length - 2));

for (const name of keyFeatures) {
  const col = result.columns.find((c) => c.name === name);
  if (!col) continue;
  const label = name.replace('pH__', '');
  const vals = Array.from(col.data).map((v) => v.toPrecision(4).padStart(10)).join('');
  console.log(`  ${label.padEnd(36)}${vals}`);
}
