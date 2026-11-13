/**
 * Adapter — route-normalize category (per-row-route Slice A/B).
 *
 * @see src/adapter/route-normalize.ts (implementation under test)
 *
 * Covers the tolerant route normalizer and its column orchestrator: synonym /
 * case / punctuation matching, NCI C66729 SEND preferred terms (must map, NOT
 * block), bare-IV → bolus default + infusion-duration override (incl. per-row
 * AC-A8), IV-adjacency BLOCK vs non-IV `other` warning, the < 50% column
 * pre-scan + its exact-50% boundary, within-subject inconsistency, and the
 * defensive scalar/column conflict.
 *
 * **Footgun gate:** every assertion uses `expect(cond === expected, true, msg)`
 * form — never `expect(x, undefined, msg)`.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {
  normalizeRouteString,
  normalizeRouteColumn,
  routeScalarColumnConflict,
  BARE_IV_DEFAULTED_TO_BOLUS,
} from '../adapter/route-normalize';

/** Build a 2-column frame (Route + Subject) plus an optional duration column. */
function frame(
  routes: ReadonlyArray<string | null>,
  subjects: ReadonlyArray<string>,
  durations?: ReadonlyArray<number | null>,
): {route: DG.Column; subject: DG.Column; dur: DG.Column | null} {
  const cols: DG.Column[] = [
    DG.Column.fromStrings('Route', routes.map((r) => r ?? '')),
    DG.Column.fromStrings('Subject', [...subjects]),
  ];
  let dur: DG.Column | null = null;
  if (durations !== undefined) {
    const arr = new Float64Array(durations.length);
    for (let i = 0; i < durations.length; i++)
      arr[i] = durations[i] === null ? DG.FLOAT_NULL : (durations[i] as number);
    dur = DG.Column.fromFloat64Array('Dur', arr);
    cols.push(dur);
  }
  const df = DG.DataFrame.fromColumns(cols);
  return {route: df.col('Route')!, subject: df.col('Subject')!, dur};
}

category('Adapter: route-normalize (string)', () => {
  test('synonyms / case / punctuation map to the right RouteCode', async () => {
    const cases: Array<[string, nca.RouteCode]> = [
      ['PO', nca.ROUTE_PO], ['p.o.', nca.ROUTE_PO], ['Oral', nca.ROUTE_PO], ['gavage', nca.ROUTE_PO],
      ['SC', nca.ROUTE_SC], ['s.c.', nca.ROUTE_SC], ['subcutaneous', nca.ROUTE_SC],
      ['IM', nca.ROUTE_IM], ['intramuscular', nca.ROUTE_IM],
      ['IV-bolus', nca.ROUTE_IV_BOLUS], ['bolus', nca.ROUTE_IV_BOLUS],
      ['IV-infusion', nca.ROUTE_IV_INFUSION], ['  Intravenous   Infusion ', nca.ROUTE_IV_INFUSION],
    ];
    for (const [s, code] of cases) {
      const r = normalizeRouteString(s, false);
      expect(r.status === 'mapped' && r.code === code, true, `'${s}' → ${code} (got ${JSON.stringify(r)})`);
    }
  });

  test('NCI C66729 SEND preferred terms map (do NOT false-positive BLOCK)', async () => {
    // INTRAVENOUS BOLUS INJECTION contains "intravenous"+"bolus" — would hit the
    // IV-adjacency guard if synonym lookup did not precede it (R1-SF-4).
    const send: Array<[string, nca.RouteCode]> = [
      ['INTRAVENOUS BOLUS INJECTION', nca.ROUTE_IV_BOLUS],
      ['SUBCUTANEOUS INJECTION', nca.ROUTE_SC],
      ['INTRAMUSCULAR INJECTION', nca.ROUTE_IM],
      ['ORAL ADMINISTRATION', nca.ROUTE_PO],
    ];
    for (const [s, code] of send) {
      const r = normalizeRouteString(s, false);
      expect(r.status === 'mapped' && r.code === code, true, `SEND '${s}' must map to ${code}`);
    }
  });

  test('bare IV defaults to IV-bolus; infusion-positive flips to IV-infusion', async () => {
    const bolus = normalizeRouteString('IV', false);
    expect(bolus.status === 'mapped' && bolus.code === nca.ROUTE_IV_BOLUS, true, 'bare IV → bolus');
    expect(bolus.status === 'mapped' && bolus.bareIvDefaulted === true, true, 'bare IV sets defaulted flag');

    const inf = normalizeRouteString('intravenous', true);
    expect(inf.status === 'mapped' && inf.code === nca.ROUTE_IV_INFUSION, true, 'bare IV + duration → infusion');
    expect(inf.status === 'mapped' && inf.bareIvDefaulted === false, true, 'no default warning when infusion');
  });

  test('explicit IV-infusion with no positive duration flags infusionMissingDuration', async () => {
    const r = normalizeRouteString('IV-infusion', false);
    expect(r.status === 'mapped' && r.infusionMissingDuration === true, true, 'missing-duration flag set');
    const r2 = normalizeRouteString('IV-infusion', true);
    expect(r2.status === 'mapped' && r2.infusionMissingDuration === false, true, 'cleared with duration');
  });

  test('IV-adjacent unmappable BLOCKS; non-IV unknown → other; intraperitoneal is non-IV', async () => {
    expect(normalizeRouteString('IV-push-fast', false).status, 'block', 'IV-adjacent unmappable blocks');
    expect(normalizeRouteString('iv drip rapid', false).status, 'block', 'infusion-adjacent blocks');
    expect(normalizeRouteString('intraperitoneal', false).status, 'other', 'IP is non-IV unknown → other');
    expect(normalizeRouteString('topical', false).status, 'other', 'topical → other');
    expect(normalizeRouteString('', false).status, 'empty', 'blank → empty');
  });
});

category('Adapter: route-normalize (column)', () => {
  test('a STRING IV/PO/SC column encodes to RouteCodes without throwing (AC-A3)', async () => {
    const {route, subject} = frame(
      ['IV', 'PO', 'SC', 'po'], ['A', 'B', 'C', 'D']);
    const res = normalizeRouteColumn(route, subject, {});
    expect(res.errors.length, 0, 'no blocking errors');
    expect(res.codes[0], nca.ROUTE_IV_BOLUS, 'IV → bolus');
    expect(res.codes[1], nca.ROUTE_PO, 'PO');
    expect(res.codes[2], nca.ROUTE_SC, 'SC');
    expect(res.codes[3], nca.ROUTE_PO, 'po (lowercase)');
    // bare-IV default warning attached to subject A.
    const a = res.warningsBySubject.get('A') ?? [];
    expect(a.some((w) => w.code === BARE_IV_DEFAULTED_TO_BOLUS), true, 'bare-IV warning for A');
  });

  test('IV-adjacent unmappable row → UNMAPPABLE_ROUTE error (AC-A5)', async () => {
    const {route, subject} = frame(['IV-push-fast', 'PO'], ['A', 'B']);
    const res = normalizeRouteColumn(route, subject, {});
    expect(res.errors.some((e) => e.code === 'UNMAPPABLE_ROUTE'), true, 'IV-adjacent blocks');
  });

  test('non-IV unknown row → other + UNMAPPABLE_ROUTE warning (AC-A5)', async () => {
    const {route, subject} = frame(['intraperitoneal', 'PO'], ['A', 'B']);
    const res = normalizeRouteColumn(route, subject, {});
    expect(res.errors.length, 0, 'no blocking error for non-IV unknown');
    expect(res.codes[0], nca.ROUTE_OTHER, 'IP → other');
    const a = res.warningsBySubject.get('A') ?? [];
    expect(a.some((w) => w.code === 'UNMAPPABLE_ROUTE' && w.severity === 'warning'), true, 'IP → warning');
  });

  test('SEX-as-route column blocks with UNMAPPABLE_ROUTE_COLUMN (<50% map, AC-A6)', async () => {
    const {route, subject} = frame(['M', 'F', 'M', 'F'], ['A', 'B', 'C', 'D']);
    const res = normalizeRouteColumn(route, subject, {});
    expect(res.errors.some((e) => e.code === 'UNMAPPABLE_ROUTE_COLUMN'), true, 'SEX column blocked');
    // No per-row processing on a blocked column.
    expect(res.codes[0] === DG.INT_NULL, true, 'codes untouched on column block');
  });

  test('exactly 50% mapping is NOT blocked (strict < boundary, R2-N4)', async () => {
    // 2 of 4 unique values map (PO, IV) → fraction = 0.5 → not blocked.
    const {route, subject} = frame(['PO', 'IV', 'xyz', 'topical'], ['A', 'B', 'C', 'D']);
    const res = normalizeRouteColumn(route, subject, {});
    expect(res.errors.some((e) => e.code === 'UNMAPPABLE_ROUTE_COLUMN'), false,
      '50% map must NOT trigger the column block');
  });

  test('within-subject route change → ROUTE_INCONSISTENT_WITHIN_SUBJECT warning (AC-B2)', async () => {
    const {route, subject} = frame(['PO', 'SC'], ['A', 'A']);
    const res = normalizeRouteColumn(route, subject, {});
    const a = res.warningsBySubject.get('A') ?? [];
    expect(a.some((w) => w.code === 'ROUTE_INCONSISTENT_WITHIN_SUBJECT'), true, 'inconsistent route warned');
  });

  test('per-row bare-IV disambiguation from a mixed duration column (AC-A8)', async () => {
    // One IV subject with duration>0 → infusion; one with 0 → bolus + warning.
    const {route, subject, dur} = frame(['IV', 'IV'], ['A', 'B'], [2.0, 0]);
    const res = normalizeRouteColumn(route, subject, {infusionDurationCol: dur});
    expect(res.codes[0], nca.ROUTE_IV_INFUSION, 'A: duration>0 → infusion');
    expect(res.codes[1], nca.ROUTE_IV_BOLUS, 'B: duration 0 → bolus');
    const a = res.warningsBySubject.get('A') ?? [];
    const b = res.warningsBySubject.get('B') ?? [];
    expect(a.some((w) => w.code === BARE_IV_DEFAULTED_TO_BOLUS), false, 'A (infusion) not defaulted');
    expect(b.some((w) => w.code === BARE_IV_DEFAULTED_TO_BOLUS), true, 'B (bolus) defaulted + warned');
  });

  test('explicit IV-infusion + no duration → ROUTE_INFUSION_DURATION_MISSING warning (AC-B1a)', async () => {
    // A: no duration → warned + still computes; B: duration 2.0 → no warning.
    const {route, subject, dur} = frame(['IV-infusion', 'IV-infusion'], ['A', 'B'], [0, 2.0]);
    const res = normalizeRouteColumn(route, subject, {infusionDurationCol: dur});
    expect(res.codes[0], nca.ROUTE_IV_INFUSION, 'still resolves to IV-infusion (computes, not NaN)');
    const a = res.warningsBySubject.get('A') ?? [];
    const b = res.warningsBySubject.get('B') ?? [];
    expect(a.some((w) => w.code === 'ROUTE_INFUSION_DURATION_MISSING'), true, 'A (no duration) warned');
    expect(b.some((w) => w.code === 'ROUTE_INFUSION_DURATION_MISSING'), false, 'B (duration 2.0) not warned');
  });

  test('defensive scalar/column conflict warning (column wins)', async () => {
    const conflict = routeScalarColumnConflict(nca.ROUTE_PO, nca.ROUTE_IV_BOLUS);
    expect(conflict !== null && conflict.code === 'ROUTE_SCALAR_COLUMN_CONFLICT', true, 'conflict warned');
    expect(routeScalarColumnConflict(nca.ROUTE_PO, nca.ROUTE_PO) === null, true, 'no conflict when equal');
    expect(routeScalarColumnConflict(nca.ROUTE_PO, null) === null, true, 'no conflict when no scalar');
  });
});
