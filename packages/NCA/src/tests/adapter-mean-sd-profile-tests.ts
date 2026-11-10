/**
 * Unit tests for `buildMeanSdProfile` (VR-03 mean ± SD profile data primitive).
 *
 * Behaviour-level fixture: a hand-built observations frame → assert mean /
 * sample-SD (Welford) / N / Upper / Lower / log variants, the N=1→SD=0 case,
 * non-finite-row skipping, and log10-of-non-positive → FLOAT_NULL.
 */

import {category, test, expect, expectFloat} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {buildMeanSdProfile} from '../adapter/mean-sd-profile';

/** Build a small (Group, Time, Conc) observations frame. */
function obsFrame(
  groups: ReadonlyArray<string>, times: ReadonlyArray<number>, concs: ReadonlyArray<number>,
): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Grp', Array.from(groups)),
    DG.Column.fromFloat64Array('T', Float64Array.from(times)),
    DG.Column.fromFloat64Array('C', Float64Array.from(concs)),
  ]);
}

function rowFor(df: DG.DataFrame, grp: string, t: number): number {
  const g = df.col('Group')!; const tm = df.col('Time')!;
  for (let i = 0; i < df.rowCount; i++)
    if (String(g.get(i)) === grp && (tm.get(i) as number) === t) return i;
  return -1;
}

category('Adapter: mean±SD profile (VR-03 primitive)', async () => {
  test('B1: arithmetic mean + sample SD (n−1) per group×time, with Upper/Lower', async () => {
    // Group A @ t=1: [2,4,6] → mean 4, sample SD sqrt((4+0+4)/2)=2.
    const df = obsFrame(['A', 'A', 'A'], [1, 1, 1], [2, 4, 6]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    expect(out.rowCount, 1, 'one (group×time) row');
    const r = rowFor(out, 'A', 1);
    expectFloat(out.col('Mean')!.get(r) as number, 4, 1e-9, 'mean = 4');
    expectFloat(out.col('SD')!.get(r) as number, 2, 1e-9, 'sample SD = 2');
    expect(out.col('N')!.get(r) as number, 3, 'N = 3');
    expectFloat(out.col('Upper')!.get(r) as number, 6, 1e-9, 'Upper = mean+SD = 6');
    expectFloat(out.col('Lower')!.get(r) as number, 2, 1e-9, 'Lower = mean−SD = 2');
  });

  test('B1: N=1 cell → SD=0 (band collapses to the mean)', async () => {
    const df = obsFrame(['B'], [1], [5]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'B', 1);
    expect(out.col('N')!.get(r) as number, 1, 'N = 1');
    expectFloat(out.col('SD')!.get(r) as number, 0, 1e-12, 'N=1 → SD = 0');
    expectFloat(out.col('Upper')!.get(r) as number, 5, 1e-9, 'Upper = mean (band collapsed)');
    expectFloat(out.col('Lower')!.get(r) as number, 5, 1e-9, 'Lower = mean (band collapsed)');
  });

  test('B1: identical multi-obs values → SD exactly 0 (Welford, no cancellation)', async () => {
    const df = obsFrame(['I', 'I', 'I', 'I'], [2, 2, 2, 2], [7, 7, 7, 7]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'I', 2);
    expect(out.col('N')!.get(r) as number, 4, 'N = 4');
    expect(out.col('SD')!.get(r) as number, 0, 'identical values → SD exactly 0');
  });

  test('B1: non-finite time or conc rows are skipped', async () => {
    // Two valid A@1 obs + one NaN-conc row (skipped) + one NaN-time row (skipped).
    const df = obsFrame(['A', 'A', 'A', 'A'], [1, 1, 1, NaN], [10, 20, NaN, 30]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    expect(out.rowCount, 1, 'only the valid A@1 cell survives');
    const r = rowFor(out, 'A', 1);
    expect(out.col('N')!.get(r) as number, 2, 'N = 2 (NaN rows skipped)');
    expectFloat(out.col('Mean')!.get(r) as number, 15, 1e-9, 'mean of [10,20] = 15');
  });

  test('B1: arithmetic Lower clamps to 0 when mean−SD < 0 (linear-axis band floor)', async () => {
    // P @ t=2: [1,9] → mean 5, SD sqrt(((1-5)^2+(9-5)^2)/1)=sqrt(32)≈5.657 → Lower=max(0,−0.657)=0.
    const df = obsFrame(['P', 'P'], [2, 2], [1, 9]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const p = rowFor(out, 'P', 2);
    expectFloat(out.col('Mean')!.get(p) as number, 5, 1e-9, 'mean = 5');
    expectFloat(out.col('Lower')!.get(p) as number, 0, 1e-9, 'Lower clamped to 0 (mean−SD < 0)');
  });

  test('B1-geo: geometric mean ×/÷ gSD over c>0 (2,4,8 → geoMean 4, geoLower 2, geoUpper 8)', async () => {
    // geo mean of [2,4,8] = (2·4·8)^(1/3) = 4; sd(ln) = ln 2 → gSD = 2 → band [2, 8].
    const df = obsFrame(['A', 'A', 'A'], [1, 1, 1], [2, 4, 8]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'A', 1);
    expect(out.col('nPos')!.get(r) as number, 3, 'nPos = 3');
    expectFloat(out.col('geoMean')!.get(r) as number, 4, 1e-9, 'geoMean = 4');
    expectFloat(out.col('geoLower')!.get(r) as number, 2, 1e-9, 'geoLower = geoMean / gSD = 2');
    expectFloat(out.col('geoUpper')!.get(r) as number, 8, 1e-9, 'geoUpper = geoMean · gSD = 8');
  });

  test('B1-geo: BLQ/zero excluded from geometric (nPos counts only c>0)', async () => {
    // [0,4,4] → arithmetic N=3 (zero included); geometric over the two 4s → geoMean 4, gSD 1.
    const df = obsFrame(['A', 'A', 'A'], [1, 1, 1], [0, 4, 4]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'A', 1);
    expect(out.col('N')!.get(r) as number, 3, 'arithmetic N = 3 (zero included)');
    expect(out.col('nPos')!.get(r) as number, 2, 'nPos = 2 (zero excluded)');
    expectFloat(out.col('geoMean')!.get(r) as number, 4, 1e-9, 'geoMean = 4');
    expectFloat(out.col('geoLower')!.get(r) as number, 4, 1e-9, 'identical c>0 → gSD 1 → band collapses');
    expectFloat(out.col('geoUpper')!.get(r) as number, 4, 1e-9, 'geoUpper = geoMean');
  });

  test('B1-geo: all-BLQ cell (nPos=0) → geometric columns FLOAT_NULL', async () => {
    const df = obsFrame(['Z', 'Z'], [1, 1], [0, 0]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'Z', 1);
    expect(out.col('nPos')!.get(r) as number, 0, 'nPos = 0');
    expect(out.col('geoMean')!.isNone(r), true, 'geoMean → FLOAT_NULL (isNone)');
    expect(out.col('geoLower')!.isNone(r), true, 'geoLower → FLOAT_NULL');
    expect(out.col('geoUpper')!.isNone(r), true, 'geoUpper → FLOAT_NULL');
  });

  test('B1-geo: nPos=1 → geoMean kept, BAND = FLOAT_NULL (SD undefined, no false zero-spread)', async () => {
    // One BLQ + one measurable → geometric mean over the single positive; SD undefined
    // at n=1 → geoLower/geoUpper null (a collapsed band would read as zero uncertainty).
    const df = obsFrame(['P', 'P'], [1, 1], [0, 5]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const r = rowFor(out, 'P', 1);
    expect(out.col('nPos')!.get(r) as number, 1, 'nPos = 1');
    expectFloat(out.col('geoMean')!.get(r) as number, 5, 1e-9, 'geoMean = 5');
    expect(out.col('geoLower')!.isNone(r), true, 'nPos=1 → geoLower FLOAT_NULL (no estimable spread)');
    expect(out.col('geoUpper')!.isNone(r), true, 'nPos=1 → geoUpper FLOAT_NULL');
  });

  test('B1-order: numeric dose labels get a NATURAL category order (10,30,100,300 not 10,100,30,300)', async () => {
    // The colour legend of the consuming scatterPlot follows the Group column's
    // `.category-order` tag; without it, numeric labels sort string-lexical.
    const df = obsFrame(['100', '10', '300', '30'], [1, 1, 1, 1], [4, 1, 8, 2]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const order = out.col('Group')!.getCategoryOrder();
    expect(order.join(','), '10,30,100,300', 'Group categories ordered by magnitude, not lexical');
  });

  test('B1-order: non-numeric groups degrade to plain lexical order (no regression)', async () => {
    // Sex F/M — numeric-aware compare falls back to lexical for non-numeric labels.
    const df = obsFrame(['M', 'F', 'M', 'F'], [1, 1, 1, 1], [3, 1, 4, 2]);
    const out = buildMeanSdProfile(df, {timeCol: 'T', concCol: 'C', groupCol: 'Grp'});
    const order = out.col('Group')!.getCategoryOrder();
    expect(order.join(','), 'F,M', 'non-numeric groups keep lexical order');
  });
});
