/**
 * Example: extracting features from multiple value columns simultaneously.
 *
 * A single sample with two measurement channels (pH and concentration).
 * Each column produces its own set of features, prefixed with the column name.
 *
 * Run: npx tsx src/time-series/feature-extraction/examples/multiple-columns.ts
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
  columns: [
    {
      name: 'pH',
      data: new Float64Array([7.0, 7.1, 7.2, 7.3, 7.4, 7.3, 7.2, 7.1, 7.0, 6.9]),
    },
    {
      name: 'concentration',
      data: new Float64Array([0.5, 0.8, 1.2, 1.5, 2.0, 2.3, 2.8, 3.1, 3.5, 4.0]),
    },
  ],
};

/* ================================================================== */
/*  Print the input                                                    */
/* ================================================================== */

console.log('Input DataFrame:');
console.log(`  rowCount: ${df.rowCount}`);
console.log(`  columns:  [${df.columns.map((c) => c.name).join(', ')}]`);
console.log('');
console.log('  id  time    pH  concentration');
console.log('  --  ----  ----  -------------');
for (let i = 0; i < df.rowCount; i++) {
  const ph = df.columns[0].data[i].toFixed(1);
  const conc = df.columns[1].data[i].toFixed(1);
  console.log(`  ${df.ids[i]}   ${df.time[i].toFixed(0).padStart(4)}  ${ph.padStart(4)}  ${conc.padStart(13)}`);
}

/* ================================================================== */
/*  Extract features                                                   */
/* ================================================================== */

const result = extractFeatures(df);

const phCols = result.columns.filter((c) => c.name.startsWith('pH__'));
const concCols = result.columns.filter((c) => c.name.startsWith('concentration__'));

console.log(`\nFeatureMatrix: ${result.nSamples} sample × ${result.columns.length} total features`);
console.log(`  pH features:            ${phCols.length}`);
console.log(`  concentration features: ${concCols.length}`);

/* ================================================================== */
/*  Compare the same features across columns                           */
/* ================================================================== */

const compareFeatures = [
  'mean', 'standard_deviation', 'variance',
  'linear_trend__attr_slope', 'linear_trend__attr_rvalue',
  'cid_ce__normalize_false', 'mean_change',
  'count_above_mean', 'count_below_mean',
];

console.log('\nSide-by-side comparison:');
console.log('');
const header = '  ' + 'feature'.padEnd(30) + 'pH'.padStart(12) + 'concentration'.padStart(16);
console.log(header);
console.log('  ' + '-'.repeat(header.length - 2));

for (const feat of compareFeatures) {
  const phCol = result.columns.find((c) => c.name === `pH__${feat}`);
  const concCol = result.columns.find((c) => c.name === `concentration__${feat}`);
  if (!phCol || !concCol) continue;

  const phVal = phCol.data[0].toPrecision(5).padStart(12);
  const concVal = concCol.data[0].toPrecision(5).padStart(16);
  console.log(`  ${feat.padEnd(30)}${phVal}${concVal}`);
}
