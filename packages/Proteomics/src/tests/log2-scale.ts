import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {SEMTYPE} from '../utils/proteomics-types';
import {log2TransformColumns, copyAsLog2Columns} from '../parsers/shared-utils';
import {getIntensityOriginals, detectCurrentLog2Applied, applyLog2Scale}
  from '../analysis/log2-scale';

/** DataFrame with two raw intensity samples (semType INTENSITY, no log2 columns
 * yet). Values are large so they unambiguously round-trip through log2. */
function makeRawDf(): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('Protein ID', ['P1', 'P2', 'P3']),
    DG.Column.fromFloat32Array('Intensity_S1', new Float32Array([1000, 4000, 16000])),
    DG.Column.fromFloat32Array('Intensity_S2', new Float32Array([2000, 8000, 32000])),
  ]);
  df.col('Protein ID')!.semType = SEMTYPE.PROTEIN_ID;
  df.col('Intensity_S1')!.semType = SEMTYPE.INTENSITY;
  df.col('Intensity_S2')!.semType = SEMTYPE.INTENSITY;
  return df;
}

category('Log2 Scale', () => {
  test('getIntensityOriginals: only the pristine (non-log2) intensity columns', async () => {
    const df = makeRawDf();
    log2TransformColumns(df, ['Intensity_S1', 'Intensity_S2']);
    const originals = getIntensityOriginals(df);
    expect(originals.length, 2, 'two originals');
    expect(originals.includes('Intensity_S1') && originals.includes('Intensity_S2'), true,
      'includes both raw columns');
    expect(originals.some((n) => n.startsWith('log2(')), false, 'excludes the log2 copies');
  });

  test('detectCurrentLog2Applied: transform vs copy is read from the values', async () => {
    const t = makeRawDf();
    log2TransformColumns(t, ['Intensity_S1', 'Intensity_S2']);
    expect(detectCurrentLog2Applied(t, getIntensityOriginals(t)), 'transform',
      'log2(x) columns read as transform');

    const c = makeRawDf();
    copyAsLog2Columns(c, ['Intensity_S1', 'Intensity_S2']);
    expect(detectCurrentLog2Applied(c, getIntensityOriginals(c)), 'copy',
      'copied columns read as copy');
  });

  test('applyLog2Scale: transform→copy rebuilds values, clears stale tags, is idempotent', async () => {
    const df = makeRawDf();
    log2TransformColumns(df, ['Intensity_S1', 'Intensity_S2']);
    // Simulate a pipeline that already ran on the (wrong-scale) data.
    df.setTag('proteomics.normalized', 'true');
    df.setTag('proteomics.de_complete', 'true');
    df.setTag('proteomics.de_method', 't-test');

    // log2 transform is currently applied; before flip the copy equals log2(1000).
    expect(Math.abs(df.col('log2(Intensity_S1)')!.get(0) - Math.log2(1000)) < 1e-3, true,
      'starts transformed');

    // Flip to "already log2" → copies the raw value straight through.
    const changed = applyLog2Scale(df, true);
    expect(changed, true, 'a real scale change returns true');
    expect(df.col('log2(Intensity_S1)')!.get(0), 1000, 'value is now the raw copy');
    expect(detectCurrentLog2Applied(df, getIntensityOriginals(df)), 'copy', 'now reads as copy');

    // Stale downstream state is cleared so viewers re-gate.
    expect(df.getTag('proteomics.normalized') || '', '', 'normalized tag cleared');
    expect(df.getTag('proteomics.de_complete') || '', '', 'de_complete tag cleared');
    expect(df.getTag('proteomics.de_method') || '', '', 'de_method tag cleared');

    // No original columns were duplicated by the rebuild.
    expect(getIntensityOriginals(df).length, 2, 'still exactly two originals');

    // Re-applying the same scale is a no-op.
    expect(applyLog2Scale(df, true), false, 'no-op when scale already matches returns false');
  });

  test('applyLog2Scale: copy→transform restores the log2 transform', async () => {
    const df = makeRawDf();
    copyAsLog2Columns(df, ['Intensity_S1', 'Intensity_S2']);
    const changed = applyLog2Scale(df, false);
    expect(changed, true, 'scale change returns true');
    expect(Math.abs(df.col('log2(Intensity_S2)')!.get(2) - Math.log2(32000)) < 1e-3, true,
      'transform re-applied to raw value');
    expect(detectCurrentLog2Applied(df, getIntensityOriginals(df)), 'transform', 'reads as transform');
  });

  test('getIntensityOriginals: empty when no intensity columns', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Protein ID', ['P1', 'P2']),
    ]);
    df.col('Protein ID')!.semType = SEMTYPE.PROTEIN_ID;
    expect(getIntensityOriginals(df).length, 0, 'no originals');
    expect(applyLog2Scale(df, true), false, 'applyLog2Scale is a no-op with nothing to rescale');
  });
});
