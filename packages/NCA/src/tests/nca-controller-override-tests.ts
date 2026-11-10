/**
 * Controller — NcaController interactive λz override category (Chunk B).
 *
 * @see src/controller/nca-controller.ts (NcaController, computeSingleProfile)
 *
 * End-to-end over the verified core: overrideLambdaZPoints writes the keyed store (AD-2),
 * the resolver builds manualPoints (augRowId + auto-window), computeNca runs on
 * the manual set, the strict ≥3 floor suppresses a 2-point fit to Not estimated,
 * and patchResultRow updates ONE row in the SAME frame instance.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {NcaController, runNcaPipeline} from '../controller/nca-controller';
import type {OverrideAction} from '../controller/nca-controller';
import {DEFAULT_RULES} from '../adapter/defaults';
import {AuditLog} from '../adapter/audit-log';
import {makeSyntheticObs, singleSubjectContract, singleSubjectScalars} from './test-utils';
import {LAMBDAZ_NONCONTIGUOUS, isLzServiceColumn} from '../adapter/lambdaz-manual-points';

const FORCE_S1_TO_2PTS: OverrideAction[] = [
  {rowIdx: 0, action: 'force-out'}, {rowIdx: 1, action: 'force-out'},
  {rowIdx: 2, action: 'force-out'}, {rowIdx: 3, action: 'force-in'},
  {rowIdx: 4, action: 'force-in'},
];

// Two IV-bolus subjects, clean log-linear decay; S1 rows 0-4, S2 rows 5-9.
function twoSubjectDf(): DG.DataFrame {
  return makeSyntheticObs({
    subjects: ['S1', 'S1', 'S1', 'S1', 'S1', 'S2', 'S2', 'S2', 'S2', 'S2'],
    times: Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8]),
    concs: Float64Array.from([10, 5, 2.5, 0.625, 0.039, 20, 10, 5, 1.25, 0.078]),
  });
}

function rowOf(df: DG.DataFrame, subject: string): number {
  const c = df.col('SUBJECT')!;
  for (let i = 0; i < c.length; i++) if (c.getString(i) === subject) return i;
  return -1;
}

category('Controller: NcaController λz override', () => {
  test('no override → computeAll matches the plain runNcaPipeline auto path (AC-3.2)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    const outcome = ctrl.computeAll();
    expect(outcome.status, 'completed', 'controller completed');

    const auto = runNcaPipeline(
      twoSubjectDf(), singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    expect(auto.status, 'completed', 'baseline completed');
    if (outcome.status !== 'completed' || auto.status !== 'completed') return;
    // Same λz for both subjects on the clean, override-free frame.
    for (const s of ['S1', 'S2']) {
      const a = auto.results.col('LAMBDA_Z')!.get(rowOf(auto.results, s));
      const b = outcome.results.col('LAMBDA_Z')!.get(rowOf(outcome.results, s));
      expect(b, a, `${s}: controller auto λz must equal runNcaPipeline`);
    }
  });

  test('force down to 2 points → Not estimated (AC-4.3); other subject + frame identity intact', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    const results = ctrl.results!;
    const s2LzBefore = results.col('LAMBDA_Z')!.get(rowOf(results, 'S2'));
    expect(results.col('STATUS')!.getString(rowOf(results, 'S1')), 'ok', 'S1 starts ok');

    // Force S1 down to exactly 2 terminal points (rows 3,4 in; 0,1,2 out).
    const actions: OverrideAction[] = [
      {rowIdx: 0, action: 'force-out'}, {rowIdx: 1, action: 'force-out'},
      {rowIdx: 2, action: 'force-out'}, {rowIdx: 3, action: 'force-in'},
      {rowIdx: 4, action: 'force-in'},
    ];
    ctrl.overrideLambdaZPoints('S1', actions);

    // Same frame instance (in-place patch, §9.4).
    expect(ctrl.results === results, true, 'results frame identity preserved (no rebuild)');
    // S1 → Not estimated: partial status, λz suppressed.
    const s1 = rowOf(results, 'S1');
    expect(results.col('STATUS')!.getString(s1), 'partial', 'S1 downgraded (Not estimated)');
    expect(results.col('LAMBDA_Z')!.isNone(s1), true, 'S1 λz suppressed');
    expect(results.col('HALF_LIFE')!.isNone(s1), true, 'S1 t½ suppressed');
    // S2 untouched.
    expect(results.col('LAMBDA_Z')!.get(rowOf(results, 'S2')), s2LzBefore, 'S2 λz unchanged');
    // Recorded in the keyed store; the source frame is NOT mutated (AC-S2.1.4).
    expect(ctrl.currentLzAction(3), 'force-in', 'S1 row3 marked force-in');
    expect(ctrl.currentLzAction(0), 'force-out', 'S1 row0 marked force-out');
    expect(df.columns.names().some(isLzServiceColumn), false, 'no _lz_* column on the source frame');
    expect(ctrl.getLzOverrides().length, 5, 'five keyed overrides');
  });

  test('manual 3-point selection stays estimable (ok, λz present)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    // force-in the 3 terminal points (rows 2,3,4); force-out 0,1.
    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'force-out'}, {rowIdx: 1, action: 'force-out'},
      {rowIdx: 2, action: 'force-in'}, {rowIdx: 3, action: 'force-in'},
      {rowIdx: 4, action: 'force-in'},
    ]);
    const s1 = rowOf(ctrl.results!, 'S1');
    expect(ctrl.results!.col('STATUS')!.getString(s1), 'ok', 'S1 still ok with 3 points');
    expect(ctrl.results!.col('LAMBDA_Z')!.isNone(s1), false, 'S1 λz estimated on 3 points');
    expect(ctrl.results!.col('LAMBDAZ_NPOINTS')!.get(s1), 3, 'exactly 3 points used');
  });

  test('reset to auto restores the auto fit (AC-4.4)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    const s1 = rowOf(ctrl.results!, 'S1');
    const autoLz = ctrl.results!.col('LAMBDA_Z')!.get(s1);

    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'force-out'}, {rowIdx: 1, action: 'force-out'},
      {rowIdx: 2, action: 'force-out'}, {rowIdx: 3, action: 'force-in'},
      {rowIdx: 4, action: 'force-in'},
    ]);
    expect(ctrl.results!.col('STATUS')!.getString(s1), 'partial', 'overridden to Not estimated');

    // Reset every S1 point to auto.
    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'auto'}, {rowIdx: 1, action: 'auto'}, {rowIdx: 2, action: 'auto'},
      {rowIdx: 3, action: 'auto'}, {rowIdx: 4, action: 'auto'},
    ]);
    expect(ctrl.results!.col('STATUS')!.getString(s1), 'ok', 'reset back to ok');
    expect(ctrl.results!.col('LAMBDA_Z')!.get(s1), autoLz, 'λz restored to the auto value');
  });

  test('sticky across setRules: override survives a global AUC-method change (AC-1.4)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'force-out'}, {rowIdx: 1, action: 'force-out'},
      {rowIdx: 2, action: 'force-out'}, {rowIdx: 3, action: 'force-in'},
      {rowIdx: 4, action: 'force-in'},
    ]);
    // Change a global rule and recompute the whole cohort.
    ctrl.setRules({...DEFAULT_RULES, aucMethod: 'linear'});
    const s1 = rowOf(ctrl.results!, 'S1');
    // The manual override is re-applied on top of the new rules → still Not estimated.
    expect(ctrl.results!.col('STATUS')!.getString(s1), 'partial',
      'override preserved across setRules (not reset to auto)');
  });

  test('reset logs reset-lambdaz, forward override logs override-lambdaz (F-4)', async () => {
    const df = twoSubjectDf();
    const log = new AuditLog();
    const ctrl = new NcaController(
      df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES, {log});
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints('S1', FORCE_S1_TO_2PTS);
    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'auto'}, {rowIdx: 1, action: 'auto'}, {rowIdx: 2, action: 'auto'},
      {rowIdx: 3, action: 'auto'}, {rowIdx: 4, action: 'auto'},
    ]);
    const actions = log.entries().map((e) => e.action);
    expect(actions.includes('override-lambdaz'), true, 'forward override audited');
    expect(actions[actions.length - 1], 'reset-lambdaz', 'a reset is audited distinctly (F-4)');
  });

  test('overriddenSubjects / overrideUserFor track the ✏️ override state (F7)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    // A clean Compute has no overrides.
    expect(ctrl.overriddenSubjects().size, 0, 'no overrides after a clean Compute');
    expect(ctrl.overrideUserFor('S1'), null, 'S1 not overridden');

    ctrl.overrideLambdaZPoints('S1', FORCE_S1_TO_2PTS);
    // Map MEMBERSHIP is the overridden test (the user label may be '' in a headless
    // context where no shell user resolves).
    expect(ctrl.overriddenSubjects().has('S1'), true, 'S1 reported overridden');
    expect(ctrl.overriddenSubjects().has('S2'), false, 'S2 not overridden');
    expect(ctrl.overrideUserFor('S1') !== null, true, 'S1 overrideUserFor non-null');
    expect(ctrl.overrideUserFor('S2'), null, 'S2 overrideUserFor null');

    // Reset clears the override → both accessors report it gone.
    ctrl.overrideLambdaZPoints('S1', [
      {rowIdx: 0, action: 'auto'}, {rowIdx: 1, action: 'auto'}, {rowIdx: 2, action: 'auto'},
      {rowIdx: 3, action: 'auto'}, {rowIdx: 4, action: 'auto'},
    ]);
    expect(ctrl.overrideUserFor('S1'), null, 'reset clears the override');
    expect(ctrl.overriddenSubjects().size, 0, 'no overrides after reset');
  });

  test('fitRowsBySubject is REAL membership — a gappy manual set excludes the interior point (F6)', async () => {
    // The honest-marker contract, and the exact case the time-window proxy got wrong:
    // an interior point that sits inside [tFirst,tLast] but is NOT in a gappy fit set.
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    // Data-driven off the ACTUAL auto fit (the terminal search drops Cmax, so it is a
    // subset of the rows, not all of them).
    const auto = ctrl.fitRowsBySubject().get('S1')!;
    expect(auto.size >= 3, true, 'the auto terminal fit uses at least 3 points');
    const members = [...auto].sort((a, b) => a - b);
    const lo = members[0];
    const hi = members[members.length - 1];
    const interior = members[1]; // strictly between lo and hi in time

    // Force the interior member OUT. Its time still lies inside [tFirst,tLast] — the
    // time-window rule would keep calling it "included"; real membership must not.
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: interior, action: 'force-out'}]);
    const gappy = ctrl.fitRowsBySubject().get('S1')!;
    expect(gappy.has(interior), false, 'the forced-out interior point is NOT in the fit (honest membership)');
    expect(gappy.has(lo) && gappy.has(hi), true, 'the points around the gap stay in the fit');

    // Reset it → back to auto ⇒ it rejoins the fit (the reset the user did on t=2).
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: interior, action: 'auto'}]);
    expect(ctrl.fitRowsBySubject().get('S1')!.has(interior), true,
      'a reset point returns to the auto fit — membership tracks the reported points, not a window');

    // S2 membership untouched throughout (single-subject cache patch).
    expect(ctrl.fitRowsBySubject().get('S2')!.size > 0, true, 'S2 keeps its own fit membership');
  });

  test('non-contiguous λz selection is flagged + warned (decision 2)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    // Auto-best-fit is contiguous by construction.
    expect(ctrl.nonContiguousSubjects().has('S1'), false, 'auto fit is contiguous');

    // Force an INTERIOR member of the auto fit out → gappy selection.
    const members = [...ctrl.fitRowsBySubject().get('S1')!].sort((a, b) => a - b);
    const interior = members[1];
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: interior, action: 'force-out'}]);
    expect(ctrl.nonContiguousSubjects().has('S1'), true, 'the gappy selection is flagged non-contiguous');

    // The advisory reached the warnings stream (WARNINGS_JSON on the patched row).
    const s1 = rowOf(ctrl.results!, 'S1');
    const warns = JSON.parse(ctrl.results!.col('WARNINGS_JSON')!.getString(s1)) as Array<{code: string}>;
    expect(warns.some((w) => w.code === LAMBDAZ_NONCONTIGUOUS), true,
      'the LAMBDAZ_NONCONTIGUOUS warning is in the stream');

    // Reset the interior point → contiguity restored, flag + warning gone.
    ctrl.overrideLambdaZPoints('S1', [{rowIdx: interior, action: 'auto'}]);
    expect(ctrl.nonContiguousSubjects().has('S1'), false, 'reset restores contiguity');
    const warns2 = JSON.parse(ctrl.results!.col('WARNINGS_JSON')!.getString(s1)) as Array<{code: string}>;
    expect(warns2.some((w) => w.code === LAMBDAZ_NONCONTIGUOUS), false,
      'the warning is cleared once the selection is contiguous again');
  });

  test('reconcile drops an override whose source datum changed (AC-3.4 / AC-S2.1.3)', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints('S1', FORCE_S1_TO_2PTS);

    // No source edit → nothing invalidated (no false positive from the canonical form).
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 0, 'no edit → override preserved');
    expect(ctrl.currentLzAction(3), 'force-in', 'row 3 still force-in');

    // Correct row 3's concentration → its override must drop on the next compute.
    df.col('conc')!.set(3, 0.7);
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 1, 'exactly the edited-row override is dropped');
    expect(ctrl.currentLzAction(3), 'auto', 'row 3 reset to auto (ALCOA+)');
    // An unedited overridden row is preserved.
    expect(ctrl.currentLzAction(4), 'force-in', 'unedited row 4 preserved');
    expect(ctrl.getLzOverrides().length, 4, 'store pruned by one');
  });

  test('D22 / AD-2: overrides are frame-independent — they survive a rebuilt, reordered frame', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    ctrl.overrideLambdaZPoints('S1', FORCE_S1_TO_2PTS);
    const saved = JSON.parse(JSON.stringify(ctrl.getLzOverrides()));
    const s1Status = ctrl.results!.col('STATUS')!.getString(rowOf(ctrl.results!, 'S1'));

    // A NEW frame (as a derived-mode Compute rebuilds it): S2 rows first, then S1 — every
    // S1 row index differs from the frame the overrides were made on.
    const rebuilt = makeSyntheticObs({
      subjects: ['S2', 'S2', 'S2', 'S2', 'S2', 'S1', 'S1', 'S1', 'S1', 'S1'],
      times: Float64Array.from([0, 1, 2, 4, 8, 0, 1, 2, 4, 8]),
      concs: Float64Array.from([20, 10, 5, 1.25, 0.078, 10, 5, 2.5, 0.625, 0.039]),
    });
    // (i) the same controller, frame swapped (a full Compute in a derived mode);
    ctrl.setInputs(rebuilt, singleSubjectContract(rebuilt), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();
    expect(ctrl.consumeReconciled(), 0, 'nothing dropped');
    expect(ctrl.currentLzAction(8), 'force-in', 'S1 t=4 (now row 8) still force-in');
    expect(ctrl.results!.col('STATUS')!.getString(rowOf(ctrl.results!, 'S1')), s1Status, 'same fit');
    // (ii) a fresh controller seeded from the persisted list (Save → reopen).
    const reopened = new NcaController(rebuilt, singleSubjectContract(rebuilt), singleSubjectScalars(), DEFAULT_RULES);
    reopened.setLzOverrides(saved);
    reopened.computeAll();
    expect(reopened.currentLzAction(5), 'force-out', 'S1 t=0 (now row 5) still force-out');
    expect(reopened.overriddenSubjects().has('S1'), true, '✎ source survives reopen');
    expect(reopened.results!.col('STATUS')!.getString(rowOf(reopened.results!, 'S1')), s1Status,
      'reopened fit equals the overridden fit');
  });

  test('who/when accessors: overrideTimeFor / overriddenSubjectTimes / lzProvenanceForRow', async () => {
    const df = twoSubjectDf();
    const ctrl = new NcaController(df, singleSubjectContract(df), singleSubjectScalars(), DEFAULT_RULES);
    ctrl.computeAll();

    // Before any override there is no who/when — the per-row accessor is the not-overridden test.
    expect(ctrl.overrideTimeFor('S1'), null, 'S1 not overridden → null time');
    expect(ctrl.overriddenSubjectTimes().size, 0, 'no overridden subjects yet');
    expect(ctrl.lzProvenanceForRow(3), null, 'an auto row has no provenance');

    ctrl.overrideLambdaZPoints('S1', FORCE_S1_TO_2PTS);

    const t = ctrl.overrideTimeFor('S1');
    expect(typeof t === 'string' && t !== '', true, 'S1 override carries a non-empty timestamp');
    expect(Number.isNaN(Date.parse(t!)), false, 'the timestamp parses as a date');
    expect(ctrl.overriddenSubjectTimes().has('S1'), true, 'S1 is in the overridden-times map');
    expect(ctrl.overriddenSubjectTimes().has('S2'), false, 'S2 (auto) is not');

    // who and when come from the SAME representative point → the same stamp across accessors
    // (one override batch stamps every forced point with one ISO time).
    const prov = ctrl.lzProvenanceForRow(3); // row 3 = force-in
    expect(prov !== null, true, 'a forced row has provenance');
    expect(prov!.time, t, 'row provenance time == subject-level overrideTimeFor (same stamp)');
    // S2's rows were never touched → still auto → null.
    expect(ctrl.lzProvenanceForRow(rowOf(df, 'S2')), null, 'S2 first row stays auto → null');
  });
});
