/**
 * F3 — λz `_lz_*` service-column exclusion (AC-3.1).
 *
 * @see src/adapter/lambdaz-manual-points.ts (isLzServiceColumn)
 * @see src/adapter/nca-repository.ts (computeSnapshotHash)
 *
 * The critical assertion: a λz override must NOT change the source frame's snapshot hash — a
 * changed hash marks the view dirty and stops the summary write-back (nca-ems-storage-migration
 * PR-2 / PR-13), so an override would spuriously read as an edit of the data.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {computeSnapshotHash, isLzServiceColumn} from '../adapter';
import {NcaController} from '../controller/nca-controller';
import {DEFAULT_RULES} from '../adapter/defaults';
import {makeSyntheticObs, singleSubjectContract, singleSubjectScalars} from './test-utils';

category('F3: λz service-column exclusion', () => {
  test('isLzServiceColumn matches only the _lz-prefixed columns', async () => {
    for (const n of ['_lz_action', '_lz_user', '_lz_time', '_lz_srchash', '_lz_future'])
      expect(isLzServiceColumn(n), true, `${n} is a service column`);
    for (const n of ['Subject', 'Time', 'conc', 'lz_action', 'LAMBDA_Z'])
      expect(isLzServiceColumn(n), false, `${n} is NOT a service column`);
  });

  test('a λz override never touches the source frame, so it cannot flag drift (AD-2 / AC-S2.1.4)', async () => {
    // Formerly asserted that computeContentHash IGNORED `_lz_*` columns on the source. Overrides
    // now live in the controller's keyed store (nca-workflow-stages AD-2): the source frame never
    // gains a service column, so the hash is stable by construction — asserted end to end.
    const df = makeSyntheticObs({
      subjects: ['S1', 'S1', 'S1', 'S1', 'S1'],
      times: Float64Array.from([0, 1, 2, 4, 8]),
      concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039]),
    });
    const h0 = computeSnapshotHash(df);
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: 0, action: 'force-out'}]);
    expect(df.columns.names().some(isLzServiceColumn), false, 'no service column on the source');
    expect(computeSnapshotHash(df), h0, 'an override leaves the source snapshot hash unchanged');

    // A genuine source edit still changes the hash.
    df.col('conc')!.set(0, 999);
    expect(computeSnapshotHash(df) === h0, false, 'a real conc edit still changes the hash');
  });

  test('SD-19: an edit in the middle of a long table changes the snapshot hash', async () => {
    // The retired sampled hash read only the first and last 10 rows of a table over 20 rows, so an
    // edit of row 15 went unseen — and a Save would then skip the new snapshot.
    const n = 30;
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('Subject', Array.from({length: n}, (_, i) => `S${i % 3}`)),
      DG.Column.fromFloat64Array('Time', Float64Array.from({length: n}, (_, i) => i)),
      DG.Column.fromFloat64Array('conc', Float64Array.from({length: n}, (_, i) => 100 / (i + 1))),
    ]);
    const h0 = computeSnapshotHash(df);
    df.col('conc')!.set(14, 7.5);
    expect(computeSnapshotHash(df) === h0, false, 'row 15 edited → a different hash');
    const h1 = computeSnapshotHash(df);
    df.col('conc')!.set(14, 100 / 15);
    expect(computeSnapshotHash(df), h0, 'the value restored → the original hash');
    expect(h1 === h0, false);
  });
});
