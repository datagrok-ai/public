/**
 * Column profile, route-based suppression, and the fragility tooltip for the
 * per-subject results grid (`nca-results-grid-columns` F1 / F4 / F5 / F6).
 *
 * Everything here runs against `resultsGridVisibleOrder`, which is deliberately
 * pure and grid-free so the partition can be asserted without a live render.
 * That shape is not incidental: the study grid's header logic once needed a real
 * grid to exercise, and shipped the literal text "undefined undefined" past
 * every automated test.
 *
 * The load-bearing test is the FIRST one. Everything else checks a column list;
 * that one checks that no future edit can hide a fragility cue by accident.
 */

import {category, test, expect} from '@datagrok-libraries/test/src/test';
import * as DG from 'datagrok-api/dg';

import {
  resultsGridVisibleOrder, colourBandProjections, readFragility,
  readAcceptance, writeAcceptanceColumns, ACCEPTANCE_COL, FLAG_COLS, ACCEPT_REASONS_COL,
  refreshAcceptanceColumnsInPlace, DEFAULT_RESULTS_THRESHOLDS, REASON_PILLS_TAG,
  problemSortRank, buildFacetSubjectOrder, PROBLEM_RANK_COL,
  writeSubjectOverrideFlags, setSubjectOverrideFlag, LZ_OVERRIDDEN_COL, LZ_OVERRIDDEN_TIME_COL,
  applyReasonPillColours,
} from '../views/results-view';
import type {AcceptanceReading} from '../views/results-view';
import {writeResultsToDataFrame, FIXED_RESULT_COLUMNS} from '../adapter/results-writer';
import {DEFAULT_AGG_PARAMS} from '../adapter/results-aggregator';
import {criterionShortLabel} from '../adapter/acceptance';
import {powerGridDrawsPills, PILLS_MIN_POWERGRID} from '../views/powergrid-support';
import {nca} from '@datagrok-libraries/sci-comp';
import type {ProfileIndex} from '../adapter/profile-index';

// ── Frame builders ──────────────────────────────────────────────────────────
// Built through the REAL writer rather than by hand, so the frames under test
// carry the true schema, the true null representation (DG.FLOAT_NULL, not raw
// NaN) and the true route strings. A hand-rolled DataFrame would let the tests
// agree with an idea of the schema instead of the schema.

function values(overrides: Partial<nca.ParameterValues> = {}): nca.ParameterValues {
  return {
    cmax: 10, tmax: 1, aucLast: 50, aucInf: 55, pctExtrap: 9,
    lambdaZ: 0.1, halfLife: 6.93, cl: 1.8, vz: 18,
    aumcLast: 200, aumcInf: 240, mrt: 4.4, vss: 8, tlag: 0.25,
    pctExtrapAumc: 16,
    ...overrides,
  } as nca.ParameterValues;
}

function lambdaZ(over: Partial<nca.LambdaZResult> = {}): nca.LambdaZResult {
  return {
    lambdaZ: 0.1, intercept: 2.3, rSquared: 0.999, adjRSquared: 0.998,
    pointsUsed: Int32Array.from([2, 3, 4]), tStart: 2, tEnd: 8,
    spanRatio: (8 - 2) / (Math.LN2 / 0.1),
    ...over,
  };
}

function result(
  status: nca.ComputeResult['status'], vals: nca.ParameterValues,
  lz: nca.LambdaZResult | null,
): nca.ComputeResult {
  return {
    status, values: vals,
    provenance: {
      lambdaZ: lz,
      blqApplied: {conc: new Float64Array(0), excluded: new Int32Array(0)},
      aucMethod: 'linear-up-log-down', compensated: false, warnings: [],
    },
  } as nca.ComputeResult;
}

function idx(subject: string): ProfileIndex {
  return {subjectKey: subject, observationRows: Int32Array.from([0])};
}

/** A results frame of `n` subjects, all on one route. */
function frame(route: nca.RouteCode, n = 3, opts: {
  failed?: boolean; vss?: number; tlag?: number;
} = {}): DG.DataFrame {
  const keys: ProfileIndex[] = []; const routes: nca.RouteCode[] = [];
  const results: nca.ComputeResult[] = [];
  for (let i = 0; i < n; i++) {
    keys.push(idx(`S${i + 1}`));
    routes.push(route);
    const v = values({
      vss: opts.vss ?? NaN, tlag: opts.tlag ?? NaN,
    });
    results.push(opts.failed === true ?
      result('failed', values({
        cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
        lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN, aumcLast: NaN,
        aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN, pctExtrapAumc: NaN,
      }), null) :
      result('ok', v, lambdaZ()));
  }
  // Production-shaped: configureResultsDataFrame adds the ACCEPTANCE + flag columns
  // before the grid renders, so the profile tests must run against a frame that
  // carries them (ACCEPTANCE is banded → in the compact floor).
  return writeAcceptanceColumns(writeResultsToDataFrame(keys, routes, results));
}

/** A mixed IV-bolus + PO frame — the case where ROUTE and CL_BASIS both vary. */
function mixedFrame(): DG.DataFrame {
  const keys = [idx('S1'), idx('S2')];
  const routes = [nca.ROUTE_IV_BOLUS, nca.ROUTE_PO];
  // S2 is 'partial', NOT 'ok': STATUS is banded and now auto-hides in compact when
  // EVERY row is ok (statusSuppressed). A non-ok row here keeps STATUS's band able
  // to fire, so this stays the canonical "every banded column can fire" frame the
  // AD-2 hard-floor test needs (values still present → the profile shape is unchanged).
  const results = [
    result('ok', values({vss: 8, tlag: NaN}), lambdaZ()),
    result('partial', values({vss: NaN, tlag: 0.25}), lambdaZ()),
  ];
  return writeAcceptanceColumns(
    writeResultsToDataFrame(keys, routes, results, Float64Array.from([5, NaN])));
}

category('Results grid: column profile (F1)', () => {
  test('AD-2 HARD FLOOR: every colour-banded column survives the compact profile', async () => {
    // THE load-bearing invariant. A banded column the default view hides is a
    // fragility cue the scientist never sees — honest-uncertainty (thesis #4)
    // made invisible by a column profile.
    //
    // Both sides are DERIVED: the floor comes from the same declaration that
    // paints the cells, and the visible set comes from the real order function.
    // Nothing here restates a list of banded columns, and that is the point —
    // three independent readers restated that set and produced three different
    // wrong answers (6, 4, and the then-correct 7). A fourth cannot ship.
    const df = mixedFrame();
    const visible = new Set(resultsGridVisibleOrder(df, false));
    const banded = [...new Set(
      colourBandProjections().filter((b) => b.grid === 'results').map((b) => b.column),
    )];
    expect(banded.length > 0, true, 'the results grid declares colour bands');
    for (const col of banded) {
      expect(visible.has(col), true,
        `${col} is colour-banded but the compact profile hides it. Its tint would ` +
        'be invisible by default. Add it to the compact set — or, if it is ' +
        'route-inapplicable on this frame, prove its band cannot fire here.');
    }
  });

  test('AC-5.5: a banded column may be hidden ONLY when its band cannot fire', async () => {
    // The precise form of the invariant. There are now TWO legal exceptions, both
    // under the SAME "band cannot fire" test: VSS (route-structural — banded AND
    // route-conditional, so on an all-PO frame every VSS cell is null and the
    // cross-column tint's `own === null` guard means it can never paint) and STATUS
    // (value-based — its partial/failed bands cannot fire when every row is `ok`;
    // asserted in the dedicated STATUS test below). The rule is "never hide a banded
    // column whose band CAN FIRE" — not "never hide a banded column".
    const po = frame(nca.ROUTE_PO, 3, {tlag: 0.25});
    const visible = new Set(resultsGridVisibleOrder(po, false));
    expect(visible.has('VSS'), false, 'VSS is hidden on an all-PO frame (never populated)');
    // ... and the premise that makes that legal is itself asserted, not assumed:
    const vss = po.col('VSS')!;
    let anyValue = false;
    for (let i = 0; i < vss.length; i++) if (!vss.isNone(i)) anyValue = true;
    expect(anyValue, false,
      'the VSS-hiding exception rests on every cell being null on this frame; ' +
      'if a value appeared here the band COULD fire and hiding would be a defect');
    // Every OTHER banded column whose band CAN fire here stays. STATUS is NOT in
    // this list: this all-PO frame is also all-'ok', so its band cannot fire either
    // and it is legitimately hidden too (see the dedicated STATUS test below).
    for (const col of ['PCT_EXTRAP', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS',
      'PCT_AUMC_EXTRAP', 'SPAN_RATIO', 'MRT', 'ACCEPTANCE'])
      expect(visible.has(col), true, `${col} stays visible on an all-PO frame`);
  });

  test('AC-1.5: the compact profile is the VR-20 scan set plus every fragility cue', async () => {
    const df = mixedFrame();
    const visible = resultsGridVisibleOrder(df, false);
    for (const col of ['SUBJECT', 'CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF',
      'HALF_LIFE', 'CL', 'VZ'])
      expect(visible.includes(col), true, `${col} (VR-20 core) is in the compact set`);

    expect(visible[0], 'SUBJECT', 'identity leads the grid');
    // And the secondary parameters are NOT — this is what makes it compact.
    for (const col of ['AUMC_LAST', 'AUMC_INF', 'LAMBDA_Z', 'LAMBDAZ_R2',
      'LAMBDAZ_TFIRST', 'LAMBDAZ_TLAST', 'LAMBDAZ_INTERCEPT', 'AUC_METHOD',
      'COMPENSATED', 'DOSE_NORM'])
      expect(visible.includes(col), false, `${col} is hidden in the compact profile`);
  });

  test('AC-1.3/1.7: JSON audit columns hidden in BOTH profiles; Full shows everything else', async () => {
    const df = mixedFrame();
    for (const extended of [false, true]) {
      const visible = new Set(resultsGridVisibleOrder(df, extended));
      expect(visible.has('WARNINGS_JSON'), false, `WARNINGS_JSON hidden (extended=${extended})`);
      expect(visible.has('PROVENANCE_JSON'), false, `PROVENANCE_JSON hidden (extended=${extended})`);
    }
    const full = resultsGridVisibleOrder(df, true);
    // Full means full: every schema column except the two JSON ones, PLUS the six
    // presentation-added λz acceptance columns (ACCEPTANCE + Violated + 4 flags) — those
    // are in PRESENTATION_COLUMN_ORDER but not in the writer's FIXED_RESULT_COLUMNS.
    expect(full.length, FIXED_RESULT_COLUMNS.length - 2 + 6,
      `Full profile shows every column but the two JSON ones, plus the 6 acceptance ` +
      `columns; got ${full.length}`);
    expect(full.length > resultsGridVisibleOrder(df, false).length, true,
      'Full is strictly wider than compact');
  });

  test('AC-1.6: compact column count per frame shape (measured, and it corrects the spec)', async () => {
    // The plan tabulated these counts as PO 17 / IV-bolus 17 / IV-infusion 16 /
    // mixed 20. Run against real writer output, PO is **16**, not 17: the table
    // added TLAG for an extravascular frame but did not subtract VSS, which
    // AC-5.2 hides on a frame with no IV rows. The other three reproduce exactly.
    //
    // This is the same error class the plan's own review caught once already (an
    // IV-bolus / IV-infusion conflation in a count) — which is why the counts are
    // pinned here against generated frames instead of transcribed from prose.
    const counts = {
      po: resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), false).length,
      bolus: resultsGridVisibleOrder(frame(nca.ROUTE_IV_BOLUS, 3, {vss: 8}), false).length,
      infusion: resultsGridVisibleOrder(frame(nca.ROUTE_IV_INFUSION, 3, {vss: 8}), false).length,
      mixed: resultsGridVisibleOrder(mixedFrame(), false).length,
    };
    // Net: +1 for ACCEPTANCE (banded, route-agnostic) but −1 for STATUS, which now
    // auto-hides on these all-'ok' single-route frames (statusSuppressed). The two
    // cancel. +1 for Violated (acceptance-reason-pills): the default fit (span 0.87,
    // 3 points) is Flagged on every one of these frames, so the pills column shows.
    expect(counts.po, 17,
      `single-route PO: 17 columns (16-column base − VSS + TLAG + ACCEPTANCE − STATUS + ` +
      `Violated, all rows ok and Flagged). Got ${counts.po}.`);
    expect(counts.bolus, 18,
      `IV-bolus: 16 + C0_EXTRAP + ACCEPTANCE − STATUS (all ok) + Violated = 18. Got ${counts.bolus}.`);
    expect(counts.infusion, 17,
      `IV-infusion: 16 + ACCEPTANCE − STATUS (all ok) + Violated — VSS applies (infusion IS ` +
      `intravascular) but C0_EXTRAP does not (IV-BOLUS only) and TLAG does not (IV). ` +
      `Got ${counts.infusion}.`);
    // 20 → 18 → 19. Parameters-only (2026-08-13) dropped both ROUTE and CL_BASIS
    // from the mixed frame; review restored CL_BASIS for its varying case
    // (see AC-5.4), so the mixed frame keeps one of the two. The
    // three single-route counts never moved: RG-J was already hiding both wherever
    // they did not vary.
    expect(counts.mixed, 21,
      `mixed IV-bolus + PO: 16 + C0_EXTRAP + TLAG + CL_BASIS + ACCEPTANCE + Violated = 21 (ROUTE ` +
      `stays Full-only; CL_BASIS returns because it varies; STATUS stays because the ` +
      `frame has a non-ok row). Got ${counts.mixed}.`);
    // The point of the profile, stated as a number rather than a claim.
    const fullWidth = resultsGridVisibleOrder(mixedFrame(), true).length;
    return `compact PO=${counts.po} bolus=${counts.bolus} infusion=${counts.infusion} ` +
      `mixed=${counts.mixed} vs Full=${fullWidth} ` +
      `(${Math.round((1 - counts.mixed / fullWidth) * 100)}% narrower on a mixed frame)`;
  });

  test('AC-1.8: carried metadata is FULL-only, kept (never dropped), right after Subject', async () => {
    // `appendGroupingColumns` adds compound / dose / sex / study id / matrix /
    // SMILES AFTER the grid is first configured, so these names are unknown to the
    // fixed schema. Two separate guarantees, and only the first one changed on
    // 2026-08-13:
    //   1. Compact shows SUBJECT + PARAMETERS, so they are Full-only. On a real rat
    //      cassette this block was nine columns before the first parameter, seven of
    //      them constant for the whole study.
    //   2. They are KEPT, never dropped — the silent-drop guard. Full lists every
    //      one, and the grid's column chooser reaches any single one individually.
    const df = mixedFrame();
    df.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Compound', ['DGK-1', 'DGK-2']));
    df.columns.add(DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'F']));

    const compact = resultsGridVisibleOrder(df, false);
    expect(compact.includes('Compound'), false, 'Compound is not in the compact profile');
    expect(compact.includes('Sex'), false, 'Sex is not in the compact profile');
    expect(compact[0], 'SUBJECT', 'identity still leads the grid');
    // The claim is "no CARRIED-METADATA block between identity and the results",
    // not "STATUS is second" — which it was only until CL_BASIS returned to the
    // compact profile on varying frames (AC-5.4). Assert the intent: everything
    // after SUBJECT belongs to the fixed schema.
    const carried = new Set(['Compound', 'Sex']);
    expect(compact.some((n) => carried.has(n)), false,
      'no carried metadata column appears anywhere in the compact profile');
    expect(compact[1], 'CL_BASIS',
      'a fixed-schema column follows the subject immediately (CL_BASIS on this ' +
      'mixed-route frame — it varies, so AC-5.4 keeps it visible)');

    const full = resultsGridVisibleOrder(df, true);
    expect(full.includes('Compound'), true, 'Compound is kept, on the Full profile');
    expect(full.includes('Sex'), true, 'Sex is kept, on the Full profile');
    expect(full.indexOf('Compound'), 1, 'Compound sits immediately after SUBJECT');
    expect(full.indexOf('Sex'), 2, 'then Sex, in frame order');
    expect(full.indexOf('CMAX') > full.indexOf('Sex'), true,
      'parameters follow the metadata columns');
  });
});

category('Results grid: inapplicable-column suppression (F5)', () => {
  test('AC-5.2: the three route-gated columns have three DIFFERENT contracts', async () => {
    // An earlier plan grouped VSS and C0_EXTRAP under one "no IV rows" rule.
    // They do not share a contract: C0_EXTRAP is null for non-IV-BOLUS, so an
    // IV-INFUSION-only study would have kept a permanently-null C0_EXTRAP column
    // on screen — reproducing the very dead-column defect F5 exists to remove.
    const po = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), false));
    expect(po.has('VSS'), false, 'PO: VSS hidden (not IV)');
    expect(po.has('C0_EXTRAP'), false, 'PO: C0_EXTRAP hidden (not IV-bolus)');
    expect(po.has('TLAG'), true, 'PO: TLAG shown (extravascular)');

    const infusion = new Set(resultsGridVisibleOrder(
      frame(nca.ROUTE_IV_INFUSION, 3, {vss: 8}), false));
    expect(infusion.has('VSS'), true, 'IV-infusion: VSS shown (infusion IS intravascular)');
    expect(infusion.has('C0_EXTRAP'), false,
      'IV-infusion: C0_EXTRAP hidden — the narrower IV-BOLUS-only contract. ' +
      'Sharing one discriminator with VSS would leave this column permanently null.');
    expect(infusion.has('TLAG'), false, 'IV-infusion: TLAG hidden (IV route)');

    const bolus = new Set(resultsGridVisibleOrder(
      frame(nca.ROUTE_IV_BOLUS, 3, {vss: 8}), false));
    expect(bolus.has('VSS'), true, 'IV-bolus: VSS shown');
    expect(bolus.has('C0_EXTRAP'), true, 'IV-bolus: C0_EXTRAP shown');
    expect(bolus.has('TLAG'), false, 'IV-bolus: TLAG hidden');
  });

  // nca-cl-bw-normalised-units AC-8.2: the absolute companions and the body weight are Full-only;
  // VSS_ABS follows VSS's route-structural rule; the compact floor stays the derived banded set.
  test('AC-8.2: CL_ABS / VZ_ABS / VSS_ABS / BW_KG are Full-only; VSS_ABS hides with VSS', async () => {
    const companions = ['CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG'];
    const bolus = frame(nca.ROUTE_IV_BOLUS, 3, {vss: 8});
    const compact = new Set(resultsGridVisibleOrder(bolus, false));
    const full = new Set(resultsGridVisibleOrder(bolus, true));
    for (const c of companions) {
      expect(bolus.col(c) !== null, true, `${c} is written to the frame`);
      expect(compact.has(c), false, `${c} hidden in Compact`);
      expect(full.has(c), true, `${c} shown in Full`);
    }
    expect(compact.has('CL') && compact.has('VZ') && compact.has('VSS'), true, 'the per-kg family stays in Compact');
    const poFull = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), true));
    expect(poFull.has('VSS_ABS'), true, 'Full suppresses nothing by route');
  });

  test('AC-5.3: an ALL-FAILED study keeps the entire λz block on screen', async () => {
    // The honesty core of F5, and the case a naive implementation gets wrong:
    // this asserts a REFUSAL to hide. Every λz column here is null, so an
    // emptiness-based rule would helpfully tidy them away — and would thereby
    // hide the finding. "The fits all failed" is the result, not an absence.
    const dead = frame(nca.ROUTE_PO, 3, {failed: true});
    const visible = new Set(resultsGridVisibleOrder(dead, false));
    for (const col of ['STATUS', 'HALF_LIFE', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS',
      'SPAN_RATIO', 'PCT_EXTRAP', 'PCT_AUMC_EXTRAP']) {
      expect(visible.has(col), true,
        `${col} must stay visible on an all-failed study — it holds explicit ` +
        'nulls, and that is the finding');
    }
    // Confirm the premise: these really are all-null on this frame, so the test
    // is exercising the refusal rather than passing because values exist.
    const adj = dead.col('LAMBDAZ_ADJ_R2')!;
    let any = false;
    for (let i = 0; i < adj.length; i++) if (!adj.isNone(i)) any = true;
    expect(any, false, 'the all-failed frame really does carry an empty λz block');
  });

  test('Computation (STATUS) auto-hides in compact when all rows ok; returns on non-ok; always in Full', async () => {
    // Banded-but-hideable, exactly like VSS: STATUS's partial/failed bands cannot
    // fire when every row is ok, so an all-ok Computation column carries no cue and
    // is decluttered from the default view. Any non-ok row brings it back; the
    // all-FAILED case (the one the scientist most needs) keeps it; Full always shows it.
    const allOk = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), false));
    expect(allOk.has('STATUS'), false, 'all rows ok → Computation hidden in compact');

    const failed = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {failed: true}), false));
    expect(failed.has('STATUS'), true, 'an all-failed frame shows Computation — the band fires');

    const mixed = new Set(resultsGridVisibleOrder(mixedFrame(), false));
    expect(mixed.has('STATUS'), true, 'a frame with a non-ok (partial) row shows Computation');

    const full = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), true));
    expect(full.has('STATUS'), true, 'Full shows Computation even when every row is ok');
  });

  test('AC-5.4: ROUTE is FULL-only; CL_BASIS returns to compact when it VARIES', async () => {
    // The parameters-only rule first subsumed RG-J wholesale, hiding CL_BASIS even
    // on a mixed-route frame. That was a science defect: `CL` / `Vz`
    // carry a NEUTRAL header on a mixed frame by construction, so with CL_BASIS
    // hidden nothing on screen separated true CL from apparent CL/F — a ~1/F
    // misread on the parameter the primary persona ranks compounds by. RG-J is
    // restored for CL_BASIS ONLY, and only in its varying case.
    const single = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), false));
    expect(single.has('ROUTE'), false, 'single-route frame: ROUTE hidden in compact');
    expect(single.has('CL_BASIS'), false,
      'single-route frame: CL_BASIS is uniform, so it says nothing and stays hidden — ' +
      'the varying gate is what keeps this from undoing the parameters-only rule');

    const mixed = new Set(resultsGridVisibleOrder(mixedFrame(), false));
    expect(mixed.has('CL_BASIS'), true,
      'mixed frame: CL_BASIS is VISIBLE in compact — it is the only per-row signal ' +
      'that separates true CL from apparent CL/F, and the header cannot carry it');
    expect(mixed.has('ROUTE'), false,
      'mixed frame: ROUTE stays hidden — it is a proxy for the basis, and the basis ' +
      'has its own column; showing both spends a column to repeat one fact');

    // Full keeps both, on either frame shape.
    const mixedFull = new Set(resultsGridVisibleOrder(mixedFrame(), true));
    expect(mixedFull.has('ROUTE'), true, 'Full shows ROUTE');
    expect(mixedFull.has('CL_BASIS'), true, 'Full shows CL_BASIS');
  });

  test('Full profile suppresses nothing by route — full means full', async () => {
    const po = new Set(resultsGridVisibleOrder(frame(nca.ROUTE_PO, 3, {tlag: 0.25}), true));
    expect(po.has('VSS'), true, 'Full shows VSS even on an all-PO frame');
    expect(po.has('C0_EXTRAP'), true, 'Full shows C0_EXTRAP even when never populated');
    expect(po.has('ROUTE'), true, 'Full shows ROUTE even when uniform');
  });

  test('no route information suppresses nothing (unknown is not evidence of absence)', async () => {
    const df = mixedFrame();
    df.columns.remove('ROUTE');
    const visible = new Set(resultsGridVisibleOrder(df, false));
    for (const col of ['VSS', 'C0_EXTRAP', 'TLAG']) {
      expect(visible.has(col), true,
        `${col} stays when the frame carries no route — an unknown route is not ` +
        'a reason to hide a parameter');
    }
  });
});

category('Results grid: scope boundaries (F3 / F6)', () => {
  test('AC-3.8: SPAN_RATIO is excluded from the aggregated "By group" parameters', async () => {
    // Decided on merit, not effort: the group MEAN of a diagnostic answers no
    // question. A cohort's mean span ratio is not a quantity anyone reports,
    // unlike a mean Cmax. Including it would also have added ~4 stat columns
    // (MEAN_/SD_/N_/CV_PCT_) to a grid already flagged for clutter.
    expect(DEFAULT_AGG_PARAMS.includes('SPAN_RATIO'), false,
      'SPAN_RATIO must not be aggregated — the mean of a diagnostic is not a statistic');
    // 12 → 15 when the absolute companions joined (nca-cl-bw-normalised-units AD-4); still no SPAN_RATIO.
    expect(DEFAULT_AGG_PARAMS.length, 15, 'the aggregated parameter set is 15');
  });

  test('AC-3.2: the new column is reserved against passthrough collision by derivation', async () => {
    // RESERVED_RESULT_COLS derives from FIXED_RESULT_COLUMNS, so a source column
    // literally named SPAN_RATIO stops being offerable as a passthrough with no
    // edit anywhere — the designed behaviour, asserted rather than assumed.
    expect(FIXED_RESULT_COLUMNS.includes('SPAN_RATIO'), true,
      'SPAN_RATIO is in the single source of truth for reserved result-column names');
  });

  test('F6: the frame keeps every column, so a full export survives a compact grid', async () => {
    // The decoupling premise, asserted at the level this module owns: the profile
    // is a GRID-side concern and never touches the frame. Measured on the
    // platform (F0 / AC-0.5): a grid-aware export drops hidden columns while a
    // frame-level one carries all of them — which is exactly why the results
    // panel's download button exports the frame, not the grid.
    const df = mixedFrame();
    const compact = resultsGridVisibleOrder(df, false);
    expect(compact.length < df.columns.length, true, 'the compact profile hides columns');
    const csvHeader = df.toCsv().split('\n')[0];
    for (const col of ['PCT_EXTRAP', 'PCT_AUMC_EXTRAP', 'LAMBDAZ_ADJ_R2', 'SPAN_RATIO']) {
      expect(csvHeader.includes(col), true,
        `${col} must survive a frame-level export even though the compact grid ` +
        'hides it — losing the fragility columns on the one path where data ' +
        'leaves the screen is the defect this guards');
    }
  });

  test('the thresholds interface carries no span field (AD-4 structural guard)', async () => {
    // The span band's two constants are module-local by design. Every field of
    // ResultsColorThresholds is a projection FROM NcaRules, so adding one here
    // would oblige the next maintainer — correctly following that rule — to
    // source it from ncaRules.lambdaZ.minSpanRatio, silently setting sci-comp's
    // gate. Keeping the interface free of a span field means there is no field
    // to project, so the invariant that would cause the defect now prevents it.
    const keys = Object.keys(DEFAULT_RESULTS_THRESHOLDS);
    for (const k of keys) {
      expect(/span/i.test(k), false,
        `ResultsColorThresholds gained "${k}". If a span threshold belongs to the ` +
        'Rules→presentation projection, it will be sourced from NcaRules and will ' +
        'start setting sci-comp minSpanRatio as a side effect. Keep it module-local.');
    }
    // 6 fields since lambdaz-adjr2-threshold-spec-v2: five are NcaRules projections;
    // lambdaAcceptanceRSquared (R_acc) is the ONE documented exception — an nca-studio
    // acceptance threshold with NO sci-comp counterpart, so it is sourced from its own
    // mapping-pane input, not projected. It is SAFE from the AD-4 wire-up defect precisely
    // because there is no ncaRules.lambdaZ.acceptanceRSquared to accidentally source. A span
    // field would NOT be safe (minSpanRatio DOES exist) — hence the /span/ guard stays.
    expect(keys.length, 6, 'the projection surface is 6 fields (5 projected + R_acc, the documented exception)');
    expect(keys.includes('lambdaAcceptanceRSquared'), true,
      'R_acc is the acceptance threshold field (applied verbatim, not derived from the gate)');
  });
});

category('Results grid: MRT/Vss fragility reading (F4)', () => {
  // AC-4.5 requires BOTH directions, and the spec is explicit that this test is
  // written net-new: a grep for a pre-existing cell-branch assertion returns
  // nothing. It is asserted against `readFragility` rather than through the grid
  // because `onCellTooltip` hands its callback to the platform and returns only a
  // subscription — the rendered text is unreachable from a test. The reading is
  // data, so its CONTENT is checkable; that is the same reason the header
  // tooltips became `RichTipSpec` instead of `() => HTMLElement`.

  /** The three sibling columns the tint reads, from a real writer-built frame. */
  function siblings(df: DG.DataFrame): [DG.Column | null, DG.Column | null, DG.Column | null] {
    return [df.col('PCT_AUMC_EXTRAP'), df.col('LAMBDAZ_ADJ_R2'), df.col('LAMBDAZ_NPOINTS')];
  }

  test('AC-4.5 (fires): a fragile row reports it AND names the condition with this row\'s value', async () => {
    // adj-R² 0.70 is below the 0.85 floor; %AUMCextrap 16 and n=3 are not the
    // trigger here, so exactly one condition must fire and it must carry the
    // LIVE value — the whole point of F4 over the generic threshold sentence.
    const keys = [idx('S1')];
    const results = [result('ok', values(), lambdaZ({adjRSquared: 0.70}))];
    const df = writeResultsToDataFrame(keys, [nca.ROUTE_IV_BOLUS], results);
    const [pct, adj, npt] = siblings(df);

    const r = readFragility(pct, adj, npt, 0, DEFAULT_RESULTS_THRESHOLDS);
    expect(r.fragile, true, 'a sub-floor adj. R² must read as fragile');
    const joined = r.reasons.join(' · ');
    expect(joined.includes('adj. R²'), true, `the fired condition is named: ${joined}`);
    expect(joined.includes('0.700'), true,
      `the row's OWN value appears, not just the threshold: ${joined}`);
  });

  test('AC-4.5 (does not fire): a clean row reads not-fragile with no reasons', async () => {
    // The negative direction. An untinted cell must keep DG's default tooltip
    // rather than gain one explaining a tint it has not got.
    //
    // The window is FOUR points on purpose. The shared frame builder's λz uses
    // three, and `fewPts` is `nPts <= minPoints` (3) — so the default fixture is
    // fragile BY DESIGN, not by accident. Writing this test against it asserted
    // the opposite of the truth and failed, which is the check working: a fit
    // sitting exactly on the floor is fragile, and the negative case needs a
    // profile that clears every one of the three conditions.
    const keys = [idx('S1')];
    const results = [result('ok', values({pctExtrapAumc: 16}), // < 20 warn
      lambdaZ({adjRSquared: 0.998, pointsUsed: Int32Array.from([2, 3, 4, 5])}))];
    const df = writeResultsToDataFrame(keys, [nca.ROUTE_IV_BOLUS], results);
    const [pct, adj, npt] = siblings(df);

    const r = readFragility(pct, adj, npt, 0, DEFAULT_RESULTS_THRESHOLDS);
    expect(r.fragile, false, 'a clean fit must not read as fragile');
    expect(r.reasons.length, 0, 'and must offer no explanation to show');
  });

  test('AC-4.1: a fit sitting exactly ON the n-points floor is fragile (<=, not <)', async () => {
    // Pinned because it is the surprise the negative test above tripped over, and
    // because the boundary is a deliberate choice the code comments defend: a fit
    // AT the algorithmic minimum is already fragile, and one somehow below it more
    // so. A future "tidy-up" to `<` would silently un-flag every 3-point fit —
    // the most common fragile case in the reference corpus.
    const keys = [idx('S1')];
    const results = [result('ok', values({pctExtrapAumc: 16}),
      lambdaZ({adjRSquared: 0.998, pointsUsed: Int32Array.from([2, 3, 4])}))];
    const df = writeResultsToDataFrame(keys, [nca.ROUTE_IV_BOLUS], results);
    const [pct, adj, npt] = siblings(df);

    const r = readFragility(pct, adj, npt, 0, DEFAULT_RESULTS_THRESHOLDS);
    expect(r.fragile, true, 'n == lambdaMinPoints must still read as fragile');
    expect(r.reasons.join(' ').includes('N points'), true,
      'and the reason names the point count that fired');
  });

  test('AC-4.3: a FAILED row is not reported fragile — the FLOAT_NULL sentinel is not a value', async () => {
    // The class-C trap this cycle exists downstream of. On a failed profile the
    // siblings hold DG.FLOAT_NULL (2.68e-34 — FINITE and POSITIVE) and INT_NULL,
    // so a bare Column.get() would make `adjR2 < 0.85` and `nPts <= 3` both TRUE
    // and report "fragile — adj. R² = 0.000" on a profile that never fit at all.
    // `numAt` is isNone-first, so the reading must come back empty instead.
    const df = frame(nca.ROUTE_IV_BOLUS, 1, {failed: true});
    const [pct, adj, npt] = siblings(df);
    expect(adj!.isNone(0), true, 'premise: the failed row really does hold a null adj. R²');

    const r = readFragility(pct, adj, npt, 0, DEFAULT_RESULTS_THRESHOLDS);
    expect(r.fragile, false,
      'a profile with NO fit must not be described as a FRAGILE fit — that is a ' +
      'statement about a regression that never happened');
    expect(r.reasons.length, 0, 'and it must invent no numbers to justify one');
  });
});

category('Results grid: λz acceptance verdict + flags (nca-grid-status-spec)', () => {
  // A single-subject frame with EXACT λz diagnostics, built through the real writer
  // so SPAN_RATIO / LAMBDAZ_ADJ_R2 / LAMBDAZ_NPOINTS / PCT_EXTRAP carry the true
  // null representation. spanRatio is set directly (the writer copies lz.spanRatio
  // verbatim); N points is the pointsUsed length.
  function oneRow(over: {
    span?: number; adjR2?: number; npts?: number; pctExtrap?: number; failed?: boolean;
  }): DG.DataFrame {
    if (over.failed === true) {
      const res = result('failed', values({
        cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
        lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN, aumcLast: NaN,
        aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN, pctExtrapAumc: NaN,
      }), null);
      return writeResultsToDataFrame([idx('S1')], [nca.ROUTE_IV_BOLUS], [res]);
    }
    const npts = over.npts ?? 6;
    const lz = lambdaZ({
      spanRatio: over.span ?? 5,
      adjRSquared: over.adjR2 ?? 0.99,
      pointsUsed: Int32Array.from(Array.from({length: npts}, (_, k) => k)),
    });
    const res = result('ok', values({pctExtrap: over.pctExtrap ?? 5}), lz);
    return writeResultsToDataFrame([idx('S1')], [nca.ROUTE_IV_BOLUS], [res]);
  }
  function reading(df: DG.DataFrame): AcceptanceReading {
    return readAcceptance(
      df.col('SPAN_RATIO'), df.col('LAMBDAZ_ADJ_R2'),
      df.col('LAMBDAZ_NPOINTS'), df.col('PCT_EXTRAP'), 0, DEFAULT_RESULTS_THRESHOLDS);
  }

  test('a clean profile is Accepted with every flag met and no reasons', async () => {
    const r = reading(oneRow({span: 5, adjR2: 0.99, npts: 6, pctExtrap: 5}));
    expect(r.verdict, 'Accepted', 'all four criteria met → Accepted');
    for (const s of [r.flags.span, r.flags.rsq, r.flags.next, r.flags.extrap])
      expect(s, 'met', 'every criterion reads met on a clean profile');
    expect(r.reasons.length, 0, 'and there is nothing to explain');
  });

  test('span < 2 flags the profile and names the threshold in SHORT form (no measured value)', async () => {
    const bad = reading(oneRow({span: 1.99}));
    expect(bad.flags.span, 'violated', 'span 1.99 < 2 violates');
    expect(bad.verdict, 'Flagged', 'any violation → Flagged');
    // Short form (user decision 2026-08-25): the ACCEPTANCE reason names the criterion +
    // threshold, NOT the row's own value — the colour band carries magnitude, and the
    // shorter phrase reads better in the facet caption + tooltip. The MRT/Vss readFragility
    // tooltip still carries the value (AC-4.5) — a deliberately different explanatory need.
    expect(bad.reasons.join(' ').includes('span ratio < 2'), true,
      `the reason names the threshold in short form: ${bad.reasons.join(' · ')}`);
    expect(bad.reasons.join(' ').includes('1.99'), false,
      'and does NOT restate the measured value (that is the readFragility contract, not this one)');
    // The boundary: exactly 2 is NOT below 2.
    expect(reading(oneRow({span: 2.0})).flags.span, 'met', 'span == 2 is met (< is strict)');
  });

  test('RSQ fires at the ACCEPTANCE threshold R_acc (default 0.90), boundary inclusive', async () => {
    // lambdaz-adjr2-threshold-spec-v2 §4: a finished fit fails acceptance iff
    // adjR² < R_acc; adjR² ≥ R_acc passes (boundary inclusive on the good side). R_acc is
    // the DEFAULT_RESULTS_THRESHOLDS value 0.90 (PKNCA min.hl.r.squared), a SEPARATE field
    // from the sci-comp point-selection gate minRSquared (0.85 by default) — the flag is
    // no longer derived from the gate (the removed goodR).
    const border = reading(oneRow({adjR2: 0.88}));
    expect(border.flags.rsq, 'violated', 'a 0.88 fit (< 0.90) fails acceptance');
    expect(border.verdict, 'Flagged', 'and the profile reads Flagged overall');
    expect(reading(oneRow({adjR2: 0.90})).flags.rsq, 'met', '0.90 == R_acc is met (boundary inclusive, < is strict)');
    expect(reading(oneRow({adjR2: 0.99})).flags.rsq, 'met', 'a clean 0.99 fit is met');
  });

  test('R_acc is applied VERBATIM — the entered value is the criterion, no derived recomputation', async () => {
    // The defect this spec fixes (§2) — the OLD goodR = min(minR + 0.05, 1.0) applied a
    // DERIVED threshold, so entering 0.95 flagged a 0.96 fit (goodR was 0.975/1.0). Now the
    // entered R_acc IS the line — at R_acc = 0.95, exactly 0.94 fails and 0.96 passes — the
    // number applied equals the number set. minR is irrelevant to acceptance (kept at its
    // default here; the R_acc ≥ minR constraint is enforced at the input, not here).
    const rAccAt = (rAcc: number, adjR2: number) => {
      const df = oneRow({adjR2});
      const t = {...DEFAULT_RESULTS_THRESHOLDS, lambdaAcceptanceRSquared: rAcc};
      return readAcceptance(df.col('SPAN_RATIO'), df.col('LAMBDAZ_ADJ_R2'),
        df.col('LAMBDAZ_NPOINTS'), df.col('PCT_EXTRAP'), 0, t).flags.rsq;
    };
    expect(rAccAt(0.95, 0.94), 'violated', 'entered 0.95 → 0.94 fails (applied at 0.95, not a derived line)');
    expect(rAccAt(0.95, 0.96), 'met', 'entered 0.95 → 0.96 passes (would have FAILED under the old goodR)');
    expect(rAccAt(0.95, 0.95), 'met', 'exactly the entered threshold passes (§4 boundary inclusive)');
    expect(rAccAt(0.80, 0.82), 'met', 'a soft 0.80 threshold accepts 0.82 (entered value applied verbatim)');
  });

  test('NEXT fires at <= minPoints (a fit AT the floor is fragile), not < ', async () => {
    // Mirrors the readFragility AC-4.1 boundary and the peach N==minP band. The
    // spec literal "< 3" would never fire (sci-comp rejects sub-minimum fits).
    expect(reading(oneRow({npts: 3})).flags.next, 'violated', 'N == 3 (the floor) violates');
    expect(reading(oneRow({npts: 4})).flags.next, 'met', 'N == 4 is met');
  });

  test('EXTRAP fires above the warn threshold 20%', async () => {
    expect(reading(oneRow({pctExtrap: 25})).flags.extrap, 'violated', '25% > 20% violates');
    expect(reading(oneRow({pctExtrap: 20})).flags.extrap, 'met', '20% is met (> is strict)');
  });

  test('a profile with NO λz fit is "Not estimated", every flag n/a, no reasons', async () => {
    // The honest third state: unevaluated criteria are neither passed nor failed.
    // Reads the FLOAT_NULL/INT_NULL sentinels through numAt — a bare get() would
    // report violated criteria on a fit that never happened.
    const r = reading(oneRow({failed: true}));
    expect(r.verdict, 'Not estimated', 'no λz → the verdict is Not estimated, not Accepted');
    for (const s of [r.flags.span, r.flags.rsq, r.flags.next, r.flags.extrap])
      expect(s, 'n/a', 'every criterion is n/a when λz was not estimated');
    expect(r.reasons.length, 0, 'and it invents no reasons');
  });

  test('multiple violations are all listed as reasons', async () => {
    const r = reading(oneRow({span: 1.5, adjR2: 0.80, npts: 3, pctExtrap: 30}));
    expect(r.verdict, 'Flagged', 'all four out of bounds → Flagged');
    expect(r.reasons.length, 4, 'every fired criterion contributes a reason');
  });

  test('the not-evaluated verdict is a VISIBLE red value "Not estimated" (present, banded)', async () => {
    // The ACCEPTANCE column stores the literal "Not estimated" (banded CRITICAL /
    // red) for an un-evaluated verdict — NOT a blank cell. λz-less means no
    // terminal-phase parameters at all, the most severe outcome, so the scientist
    // must see it rather than read an empty cell.
    const dead = writeAcceptanceColumns(frame(nca.ROUTE_PO, 2, {failed: true}));
    const col = dead.col(ACCEPTANCE_COL)!;
    expect(col.isNone(0), false, 'a not-evaluated verdict is a present, visible value');
    expect(col.get(0), 'Not estimated', 'and its label is "Not estimated"');
    // It is a banded category (so the red fill can attach), like Accepted/Flagged.
    const bandedCats = colourBandProjections()
      .filter((b) => b.grid === 'results' && b.column === ACCEPTANCE_COL)
      .map((b) => b.ruleKey);
    expect(bandedCats.includes('Not estimated'), true,
      '"Not estimated" is a declared ACCEPTANCE band (red), not an unpainted value');
  });

  test('ACCEPTANCE is compact-visible (banded floor); the 4 flags are Full-only', async () => {
    const compact = new Set(resultsGridVisibleOrder(mixedFrame(), false));
    expect(compact.has(ACCEPTANCE_COL), true,
      'the roll-up verdict rides the derived banded floor into the compact profile');
    for (const flag of FLAG_COLS)
      expect(compact.has(flag), false, `${flag} is Full-only (non-banded, filter/count column)`);

    const full = new Set(resultsGridVisibleOrder(mixedFrame(), true));
    expect(full.has(ACCEPTANCE_COL), true, 'Full keeps ACCEPTANCE');
    for (const flag of FLAG_COLS)
      expect(full.has(flag), true, `Full shows ${flag}`);
  });

  test('ACCEPTANCE joins the DERIVED banded set (so AD-2 keeps it visible with no restated list)', async () => {
    const banded = new Set(
      colourBandProjections().filter((b) => b.grid === 'results').map((b) => b.column));
    expect(banded.has(ACCEPTANCE_COL), true,
      'ACCEPTANCE is in buildColourRules, so the compact floor picks it up by derivation');
    for (const flag of FLAG_COLS)
      expect(banded.has(flag), false, `${flag} is deliberately NOT banded (Full-only)`);
  });

  test('writeAcceptanceColumns is idempotent (upsert, not append)', async () => {
    const df = frame(nca.ROUTE_IV_BOLUS, 3, {vss: 8}); // builder already added them once
    const before = df.columns.length;
    writeAcceptanceColumns(df, DEFAULT_RESULTS_THRESHOLDS);
    writeAcceptanceColumns(df, DEFAULT_RESULTS_THRESHOLDS);
    expect(df.columns.length, before, 're-running adds no duplicate columns');
    for (const name of [ACCEPTANCE_COL, ACCEPT_REASONS_COL, ...FLAG_COLS])
      expect(df.col(name) !== null, true, `${name} present exactly once`);
  });
});

// ── λz reason pills (acceptance-reason-pills) ───────────────────────────────
category('Results grid: λz reason pills (acceptance-reason-pills)', () => {
  /** One IV-bolus subject per spec, through the real writer + writeAcceptanceColumns. */
  function rows(specs: Array<{span: number; adjR2: number; npts: number; pctExtrap: number}>): DG.DataFrame {
    const results = specs.map((s) => result('ok', values({pctExtrap: s.pctExtrap}), lambdaZ({
      spanRatio: s.span, adjRSquared: s.adjR2,
      pointsUsed: Int32Array.from(Array.from({length: s.npts}, (_, k) => k)),
    })));
    return writeAcceptanceColumns(writeResultsToDataFrame(
      specs.map((_, i) => idx(`S${i + 1}`)), specs.map(() => nca.ROUTE_IV_BOLUS), results));
  }
  const CLEAN = {span: 5, adjR2: 0.99, npts: 6, pctExtrap: 5};

  test('AC-4: short labels per criterion, following the live thresholds', async () => {
    const t = DEFAULT_RESULTS_THRESHOLDS;
    expect(criterionShortLabel('SPAN_LT', t), 'span < 2');
    expect(criterionShortLabel('ADJR2_LT', t), 'R² < 0.9');
    expect(criterionShortLabel('NPTS_LE', t), 'N ≤ 3');
    expect(criterionShortLabel('EXTRAP_GT', t), '%extrap > 20%');
    const soft = {...t, lambdaAcceptanceRSquared: 0.85, lambdaMinPoints: 4, extrapWarnPct: 25};
    expect(criterionShortLabel('ADJR2_LT', soft), 'R² < 0.85', 'R_acc moves the label');
    expect(criterionShortLabel('NPTS_LE', soft), 'N ≤ 4', 'min points moves the label');
    expect(criterionShortLabel('EXTRAP_GT', soft), '%extrap > 25%', 'the warn threshold moves the label');
  });

  test('pill-colors AC-1/AC-2: each criterion keeps its own fill, also under edited thresholds', async () => {
    // Keyed by criterion CODE: a swapped mapping, or a fill lost when a threshold renames the pill,
    // fails here (a set-of-fills check would not catch either).
    const want = {SPAN_LT: '#BFE4F7', ADJR2_LT: '#F5DB98', NPTS_LE: '#EDC6DD', EXTRAP_GT: '#F3C6A6'} as const;
    const soft = {...DEFAULT_RESULTS_THRESHOLDS, lambdaAcceptanceRSquared: 0.85, lambdaMinPoints: 4, extrapWarnPct: 25};
    for (const t of [DEFAULT_RESULTS_THRESHOLDS, soft]) {
      const grid = DG.Viewer.grid(rows([{span: 1.5, adjR2: 0.80, npts: 3, pctExtrap: 30}, CLEAN]));
      applyReasonPillColours(grid, t);
      const colours = grid.col(ACCEPT_REASONS_COL)!.categoryColors;
      expect(Object.keys(colours).length, 4, 'one entry per criterion');
      for (const [code, hex] of Object.entries(want)) {
        const label = criterionShortLabel(code as keyof typeof want, t);
        const c = colours[label];
        expect(c !== undefined, true, `'${label}' has a fill`);
        expect(DG.Color.toHtml(typeof c === 'string' ? DG.Color.fromHtml(c) : c).toUpperCase(), hex,
          `'${label}' (${code}) is ${hex}`);
      }
    }
  });

  test('AC-2: the cell value is the violated labels in criterion order; the column is a Tags column', async () => {
    const df = rows([{span: 1.5, adjR2: 0.80, npts: 3, pctExtrap: 30}, {...CLEAN, npts: 3}, CLEAN]);
    const col = df.col(ACCEPT_REASONS_COL)!;
    expect(col.getString(0), 'span < 2,R² < 0.9,N ≤ 3,%extrap > 20%', 'all four, criterion order');
    expect(col.getString(1), 'N ≤ 3', 'one violated criterion → one pill');
    expect(col.getString(2), '', 'an Accepted profile has no pill');
    // `tagReasonPillsColumn` writes the tag names as literals (headless module) — pin them to DG's.
    expect(DG.TAGS.CELL_RENDERER, 'cell.renderer', 'DG tag name unchanged');
    expect(DG.TAGS.MULTI_VALUE_SEPARATOR, '.multi-value-separator', 'DG tag name unchanged');
    expect(col.getTag(DG.TAGS.CELL_RENDERER), 'Tags', 'drawn by the platform Tags renderer');
    expect(col.getTag(DG.TAGS.MULTI_VALUE_SEPARATOR), ',', 'split on the reasons separator (multi-value filter)');
    expect(df.col(ACCEPTANCE_COL)!.getString(0), 'Flagged', 'the verdict cell keeps the bare verdict');
  });

  test('AC-1: compact shows Violated right after Acceptance when a profile is Flagged', async () => {
    const order = resultsGridVisibleOrder(rows([{...CLEAN, npts: 3}, CLEAN]), false);
    const at = order.indexOf(ACCEPTANCE_COL);
    expect(order[at + 1], ACCEPT_REASONS_COL, `Violated follows Acceptance: ${order.join(', ')}`);
  });

  test('AC-3: no Flagged profile → Violated hidden in compact, shown in Full', async () => {
    const df = rows([CLEAN, CLEAN]);
    expect(resultsGridVisibleOrder(df, false).includes(ACCEPT_REASONS_COL), false, 'compact hides an empty column');
    expect(resultsGridVisibleOrder(df, true).includes(ACCEPT_REASONS_COL), true, 'Full always shows it');
  });

  test('Settings switch OFF hides Violated in BOTH profiles, even with Flagged rows', async () => {
    const df = rows([{...CLEAN, npts: 3}, CLEAN]);
    df.setTag(REASON_PILLS_TAG, 'false');
    expect(resultsGridVisibleOrder(df, false).includes(ACCEPT_REASONS_COL), false, 'compact');
    expect(resultsGridVisibleOrder(df, true).includes(ACCEPT_REASONS_COL), false, 'Full');
    df.setTag(REASON_PILLS_TAG, 'true');
    expect(resultsGridVisibleOrder(df, false).includes(ACCEPT_REASONS_COL), true, 'ON restores the rule');
  });

  test('PowerGrid version gate: a release after 1.8.4, or a debug build of 1.8.4+, draws the pills', async () => {
    expect(PILLS_MIN_POWERGRID, '1.8.5');
    for (const [v, want] of [
      ['1.8.4', false], ['1.7.22', false], ['1.4.3', false], ['1.8.3.X-abc', false], ['', false], ['dev', false],
      ['1.8.4.X-44ea6cd6', true], ['1.8.4.X', true], ['1.9.0.X-1f2e', true],
      ['1.8.5', true], ['1.8.10', true], ['1.9.0', true], ['2.0.0', true],
    ] as Array<[string, boolean]>)
      expect(powerGridDrawsPills(v), want, `version '${v}'`);
  });

  test('AC-5: an in-place refresh re-derives the pills from the patched λz numerics', async () => {
    const df = rows([CLEAN]);
    expect(df.col(ACCEPT_REASONS_COL)!.getString(0), '', 'clean before');
    df.col('LAMBDAZ_ADJ_R2')!.set(0, 0.8);
    refreshAcceptanceColumnsInPlace(df, DEFAULT_RESULTS_THRESHOLDS);
    expect(df.col(ACCEPT_REASONS_COL)!.getString(0), 'R² < 0.9', 'the override-driven change reaches the pill');
  });
});

// ── Problem-first sort (results-problem-first-sort) ─────────────────────────
/** One row of EACH disposition bucket, DELIBERATELY out of problem-first order
 *  (an Accepted row leads), so the rank must actually reorder. Buckets are
 *  deterministic against the default thresholds (minPoints 3, R_acc 0.90, warn 20). */
function mixedStatusFrame(): DG.DataFrame {
  const keys = [idx('S_acc'), idx('S_failed'), idx('S_partial'), idx('S_flag')];
  const routes = [nca.ROUTE_PO, nca.ROUTE_PO, nca.ROUTE_PO, nca.ROUTE_PO];
  const results = [
    // Strong fit — every acceptance criterion met → Accepted (rank 3).
    result('ok', values({pctExtrap: 5, pctExtrapAumc: 8}),
      lambdaZ({pointsUsed: Int32Array.from([2, 3, 4, 5, 6]),
        tStart: 1, tEnd: 40, rSquared: 0.999, adjRSquared: 0.999,
        spanRatio: (40 - 1) / (Math.LN2 / 0.1)})),
    // No measurable point → failed (rank 0). Verdict resolves to Not estimated.
    result('failed', values({
      cmax: NaN, tmax: NaN, aucLast: NaN, aucInf: NaN, pctExtrap: NaN,
      lambdaZ: NaN, halfLife: NaN, cl: NaN, vz: NaN, aumcLast: NaN,
      aumcInf: NaN, mrt: NaN, vss: NaN, tlag: NaN, pctExtrapAumc: NaN,
    }), null),
    // Cmax/AUClast obtained but λz not fit → partial (rank 1). Every λz-derived
    // field is NaN, so the verdict is Not estimated; STATUS='partial' fixes the rank.
    result('partial', values({
      aucInf: NaN, pctExtrap: NaN, lambdaZ: NaN, halfLife: NaN, cl: NaN,
      vz: NaN, aumcInf: NaN, mrt: NaN, vss: NaN, pctExtrapAumc: NaN,
    }), null),
    // Default fit — span 0.87 (< 2) and 3 points (≤ min) → Flagged (rank 2).
    result('ok', values(), lambdaZ()),
  ];
  return writeAcceptanceColumns(writeResultsToDataFrame(keys, routes, results));
}

category('Results grid: problem-first sort (results-problem-first-sort)', () => {
  test('problemSortRank pins the ratified order failed < Not estimated < Flagged < Accepted', async () => {
    // STATUS is primary: failed and partial BOTH carry the Not-estimated verdict,
    // so only STATUS separates an engine failure (0) from a λz-less profile (1).
    expect(problemSortRank('failed', 'Not estimated'), 0, 'failed leads');
    expect(problemSortRank('failed', 'Accepted'), 0, 'failed wins over any verdict');
    expect(problemSortRank('partial', 'Not estimated'), 1, 'partial is the Not-estimated block');
    expect(problemSortRank('ok', 'Not estimated'), 1, 'defensive: ok+NA ranks with Not estimated');
    expect(problemSortRank('ok', 'Flagged'), 2, 'Flagged');
    expect(problemSortRank('ok', 'Accepted'), 3, 'Accepted sinks last');
    expect(problemSortRank(null, 'Accepted'), 3, 'null status falls through to the verdict');
    expect(problemSortRank('ok', 'nonsense'), 4, 'unknown verdict sinks below everything');
  });

  test('_PROBLEM_RANK is written consistently with the verdict beside it, every row', async () => {
    // Read STATUS + ACCEPTANCE back off the frame and assert the rank column equals
    // the pure function — the rank can never drift from the verdict it encodes.
    const df = mixedStatusFrame();
    const statusCol = df.col('STATUS')!;
    const acceptCol = df.col(ACCEPTANCE_COL)!;
    const rankCol = df.col(PROBLEM_RANK_COL)!;
    for (let i = 0; i < df.rowCount; i++) {
      const st = statusCol.isNone(i) ? null : statusCol.getString(i);
      const verdict = acceptCol.getString(i);
      expect(rankCol.getString(i), String(problemSortRank(st, verdict)),
        `row ${i} (${st}/${verdict}) rank matches the pure function`);
    }
  });

  test('the four buckets reorder to failed → partial → Flagged → Accepted', async () => {
    const df = mixedStatusFrame();
    const subj = df.col('SUBJECT')!;
    const rankCol = df.col(PROBLEM_RANK_COL)!;
    const rankOf = new Map<string, number>();
    for (let i = 0; i < df.rowCount; i++) rankOf.set(subj.getString(i), Number(rankCol.getString(i)));
    expect(rankOf.get('S_failed'), 0, 'the failed profile is at the very top');
    expect(rankOf.get('S_partial'), 1, 'partial (Not estimated) sits directly below failed');
    expect(rankOf.get('S_flag'), 2, 'a flagged ok profile ranks above an accepted one');
    expect(rankOf.get('S_acc'), 3, 'a clean Accepted profile sinks to the bottom');
    const order = [...rankOf.entries()].sort((a, b) => a[1] - b[1]).map((e) => e[0]);
    expect(order.join(','), 'S_failed,S_partial,S_flag,S_acc', 'problem-first row order');
  });

  test('_PROBLEM_RANK is a hidden helper — invisible in BOTH grid profiles', async () => {
    const df = mixedStatusFrame();
    for (const extended of [false, true]) {
      const visible = new Set(resultsGridVisibleOrder(df, extended));
      expect(visible.has(PROBLEM_RANK_COL), false,
        `_PROBLEM_RANK never rides the grid (extended=${extended})`);
    }
  });

  test('buildFacetSubjectOrder mirrors the grid order and drops unknown subjects', async () => {
    const df = mixedStatusFrame();
    const order = buildFacetSubjectOrder(df);
    const bySortKey = [...order.entries()].sort((a, b) => a[1] - b[1]).map((e) => e[0]);
    expect(bySortKey.join(','), 'S_failed,S_partial,S_flag,S_acc',
      'facet order reproduces the grid problem-first order');
    expect(order.has('S_ghost'), false,
      'a subject absent from the results frame is not in the order map (facet leaves it in place)');
  });
});

// ── λz manual-override ✏️ flag (F7 / AC-7.1) ─────────────────────────────────
function subjRow(df: DG.DataFrame, subject: string): number {
  const c = df.col('SUBJECT')!;
  for (let i = 0; i < c.length; i++) if (c.getString(i) === subject) return i;
  return -1;
}

category('Results grid: λz override ✏️ flag (F7)', () => {
  test('writeSubjectOverrideFlags marks the overridden subject and clears the rest', async () => {
    const df = frame(nca.ROUTE_PO, 3);
    writeSubjectOverrideFlags(df, new Map([['S2', 'Dr. Smith']]));
    const col = df.col(LZ_OVERRIDDEN_COL)!;
    expect(col.getString(subjRow(df, 'S2')), 'Dr. Smith', 'S2 carries the analyst label');
    expect(col.getString(subjRow(df, 'S1')), '', 'S1 is not overridden → blank');
    expect(col.getString(subjRow(df, 'S3')), '', 'S3 is not overridden → blank');

    // A subject dropping out of the map (override reset away) is cleared to blank.
    writeSubjectOverrideFlags(df, new Map([['S1', 'Ana']]));
    expect(col.getString(subjRow(df, 'S1')), 'Ana', 'S1 now overridden');
    expect(col.getString(subjRow(df, 'S2')), '', 'S2 cleared (no longer in the map)');
  });

  test('an override with a blank analyst still shows the ✏️ (non-empty sentinel)', async () => {
    const df = frame(nca.ROUTE_PO, 2);
    // Map MEMBERSHIP means overridden even when the analyst label is '' (headless).
    writeSubjectOverrideFlags(df, new Map([['S1', '']]));
    const col = df.col(LZ_OVERRIDDEN_COL)!;
    expect(col.getString(subjRow(df, 'S1')) !== '', true,
      'a blank-user override still stores a non-empty cell so the ✏️ renders');
  });

  test('setSubjectOverrideFlag updates one cell; null clears it', async () => {
    const df = frame(nca.ROUTE_PO, 3);
    writeSubjectOverrideFlags(df, new Map()); // ensure the column exists
    const col = df.col(LZ_OVERRIDDEN_COL)!;
    setSubjectOverrideFlag(df, 'S3', 'Bo');
    expect(col.getString(subjRow(df, 'S3')), 'Bo', 'S3 flagged after a single-subject toggle');
    expect(col.getString(subjRow(df, 'S1')), '', 'other rows untouched');
    setSubjectOverrideFlag(df, 'S3', null);
    expect(col.getString(subjRow(df, 'S3')), '', 'null clears the flag (reverted to auto)');
  });

  test('the _lz_overridden service column is hidden in BOTH grid profiles', async () => {
    const df = frame(nca.ROUTE_PO, 2);
    writeSubjectOverrideFlags(df, new Map([['S1', 'x']]));
    expect(resultsGridVisibleOrder(df, true).includes(LZ_OVERRIDDEN_COL), false,
      'Full profile must not surface the _lz_overridden service column');
    expect(resultsGridVisibleOrder(df, false).includes(LZ_OVERRIDDEN_COL), false,
      'Compact profile must not surface it either');
  });
});

// ── λz-override who/when TIMESTAMP projection (_lz_overridden_time) ───────────
category('Results grid: λz override who/when timestamp', () => {
  test('writeSubjectOverrideFlags writes _lz_overridden_time when given a time map', async () => {
    const df = frame(nca.ROUTE_PO, 3);
    writeSubjectOverrideFlags(
      df, new Map([['S2', 'Dr. Smith']]), new Map([['S2', '2026-09-08T14:30:00.000Z']]));
    const tcol = df.col(LZ_OVERRIDDEN_TIME_COL)!;
    expect(tcol.getString(subjRow(df, 'S2')), '2026-09-08T14:30:00.000Z', 'S2 carries the raw ISO time');
    expect(tcol.getString(subjRow(df, 'S1')), '', 'S1 not overridden → blank time');
  });

  test('the 2-arg writeSubjectOverrideFlags is backward-compatible (no time column created)', async () => {
    const df = frame(nca.ROUTE_PO, 2);
    writeSubjectOverrideFlags(df, new Map([['S1', 'Ana']])); // legacy 2-arg call
    expect(df.col(LZ_OVERRIDDEN_TIME_COL), null, 'no _lz_overridden_time created by the 2-arg call');
    expect(df.col(LZ_OVERRIDDEN_COL)!.getString(subjRow(df, 'S1')), 'Ana', 'the flag is still written');
  });

  test('setSubjectOverrideFlag writes then clears _lz_overridden_time', async () => {
    const df = frame(nca.ROUTE_PO, 3);
    setSubjectOverrideFlag(df, 'S3', 'Bo', '2026-01-02T03:04:00.000Z');
    const tcol = df.col(LZ_OVERRIDDEN_TIME_COL)!;
    expect(tcol.getString(subjRow(df, 'S3')), '2026-01-02T03:04:00.000Z', 'time set on the toggle');
    setSubjectOverrideFlag(df, 'S3', null, null);
    expect(tcol.getString(subjRow(df, 'S3')), '', 'null user/time clears the timestamp');
  });

  test('the _lz_overridden_time service column is hidden in BOTH grid profiles', async () => {
    const df = frame(nca.ROUTE_PO, 2);
    writeSubjectOverrideFlags(
      df, new Map([['S1', 'x']]), new Map([['S1', '2026-09-08T00:00:00.000Z']]));
    expect(resultsGridVisibleOrder(df, true).includes(LZ_OVERRIDDEN_TIME_COL), false,
      'Full profile must not surface _lz_overridden_time');
    expect(resultsGridVisibleOrder(df, false).includes(LZ_OVERRIDDEN_TIME_COL), false,
      'Compact profile must not surface it either');
  });
});
