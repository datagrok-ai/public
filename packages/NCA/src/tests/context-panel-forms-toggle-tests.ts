/**
 * Study-overview context panel — two fixed defects (topic
 * context-panel-forms-toggle-dose-dup):
 *
 *  BUG A — the "NCA parameters" Forms viewer must MIRROR the grid's "Full" toggle. Covered by
 *    {@link studyParamsFormFields} over a REAL {@link buildStudyGrid} frame: compact → Mean-only
 *    SAR set, no %CV/n, no secondary params, no C–t / Warnings; extended → the full set; the C–t
 *    thumbnail (a DATA_FRAME column) is never a form field either way.
 *
 *  BUG B — {@link rowIdentityLabel} must not print a group value twice when two `groupValues`
 *    entries share a typographic canonical key (the `friendlyName`-collision case). The study grid
 *    merges those into ONE key column ({@link buildStudyGrid}); the label must dedup on the SAME
 *    `normalizeGroupKey`, so the panel/legend label == the grid's key columns.
 *
 * Real-fixture (NOT mirror-pattern, CLAUDE.md rule 16): the production `buildStudyGrid` /
 * `rowIdentityLabel` / `studyParamsFormFields` run on a constructed tier, and the assertions are
 * the cross-function INVARIANT (label group-segment count == grid key-column count).
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  buildStudyGrid, KEY_SEP, ANALYSIS_COL, ROUTE_COL, CURVE_COL, FLAGS_COL,
  ANALYSIS_ID_COL, COMPOUND_ID_COL, PARAM_TAG, PASSTHROUGH_TAG, dedupGroupValuesByCanonical,
} from '../adapter';
import type {NcaStudy, CompoundTierRow, MeanCv, GroupKeyValue} from '../adapter';
import {
  rowIdentityLabel, resolveClickedRow, studyParamsFormFields, studyGridVisibleOrder,
  STUDY_EXTENDED_TAG,
} from '../views/nca-study-overview';

function mv(mean: number, cvPct: number, n: number): MeanCv {
  return {mean, cvPct, n};
}

/** A tier row with the given per-key group values. `compoundId` is the composite key (the values
 *  joined by KEY_SEP), matching what `buildSummaryCache` emits (and what `resolveClickedRow`
 *  matches on). */
function tierRow(groupValues: GroupKeyValue[], route: string, over: Partial<CompoundTierRow> = {}): CompoundTierRow {
  const compoundId = groupValues.length === 0 ? 'WHOLE' : groupValues.map((g) => g.value).join(KEY_SEP);
  return {
    compoundId, isWholeAnalysis: groupValues.length === 0, groupValues,
    keyInventory: [], passthrough: [],
    n: 3, route, doseNorm: 'per_kg',
    cmax: mv(10, 5, 3), tmax: mv(1, 20, 3), aucInf: mv(100, 8, 3), thalf: mv(2, 12, 3),
    cl: mv(1.5, 9, 3), vz: mv(3, 11, 3), mrt: mv(4, 10, 3), vss: mv(6, 9, 3),
    fPct: null, pctExtrap: mv(15, 4, 3), curve: null, warnings: [], ...over,
  } as unknown as CompoundTierRow;
}

function studyWith(rows: CompoundTierRow[], groupByCols: string[], label = 'Screen A'): NcaStudy {
  return {
    id: 's1', label: 'Study', kind: 'study',
    createdAt: '2026-07-16T00:00:00.000Z', updatedAt: '2026-07-16T00:00:00.000Z',
    access: {canEdit: true, canDelete: true, canShare: true},
    analyses: [{
      id: 'a1', label, source: {kind: 'tableId', ref: 't'},
      columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols, passthroughCols: [],
      sourceSignature: {kind: 'content-hash', rowCount: 6, mtime: null, hash: 'h',
        computedAt: '2026-07-16T00:00:00.000Z'},
      summaryCache: {
        analysisTier: {nCompounds: rows.length, nSubjects: 3, route: 'IV-bolus', dose: '(per-row)',
          status: 'ok', flagSummary: 'No flags'},
        compoundTier: rows,
      },
    }],
  } as unknown as NcaStudy;
}

// ── BUG B: dose-doubling in the identity label ──────────────────────────────────────────────────

category('context-panel: rowIdentityLabel canonical dedup (BUG B)', () => {
  // Two groupValues keys that collapse under normalizeGroupKey (casefold): 'Dose_mg_kg' and
  // 'dose_mg_kg', both value '1' — the friendlyName-collision that produced "… · 1 · F · 1 · …".
  const dupRow = tierRow(
    [{key: 'Dose_mg_kg', value: '1'}, {key: 'Sex', value: 'F'}, {key: 'dose_mg_kg', value: '1'}],
    'IV-bolus');
  const study = studyWith([dupRow], ['Dose_mg_kg', 'Sex'], 'New NCA');

  test('label dedups the canonical-duplicate key — dose printed ONCE', async () => {
    const resolved = resolveClickedRow(study, 'a1', dupRow.compoundId, 'IV-bolus');
    expect(resolved !== null, true, 'row resolves');
    const label = rowIdentityLabel(resolved!);
    expect(label, 'New NCA · 1 · F · IV-bolus', 'no doubled dose');
    // Exactly one "1" segment (the pre-fix bug printed two).
    const ones = label.split(' · ').filter((s) => s === '1').length;
    expect(ones, 1, 'dose value appears once');
  });

  test('label group-segment count == grid key-column count (cross-function invariant)', async () => {
    const df = buildStudyGrid(study);
    // Grid key columns = the columns that are neither fixed provenance nor param/curve — the merged
    // per-key set. buildStudyGrid merged Dose_mg_kg + dose_mg_kg into ONE column.
    const fixed = new Set<string>([ANALYSIS_COL, ROUTE_COL, CURVE_COL, FLAGS_COL,
      ANALYSIS_ID_COL, COMPOUND_ID_COL]);
    const keyCols = df.columns.names().filter((n) =>
      !fixed.has(n) && (df.col(n)?.getTag(PARAM_TAG) ?? '') === '' &&
      (df.col(n)?.getTag(PASSTHROUGH_TAG) ?? '') === '');
    expect(keyCols.length, 2, 'grid merges to two key columns (Dose_mg_kg, Sex)');

    const resolved = resolveClickedRow(study, 'a1', dupRow.compoundId, 'IV-bolus')!;
    // Label group segment = everything between analysis label and route.
    const segs = rowIdentityLabel(resolved).split(' · ');
    const groupSegs = segs.slice(1, segs.length - 1); // drop analysis + route
    expect(groupSegs.length, keyCols.length, 'label group segments match grid key columns');
  });

  test('a genuinely-distinct multi-key row is unaffected (no over-dedup)', async () => {
    const row = tierRow([{key: 'Dose_mg_kg', value: '5'}, {key: 'Sex', value: 'M'}], 'PO');
    const st = studyWith([row], ['Dose_mg_kg', 'Sex'], 'New NCA');
    const resolved = resolveClickedRow(st, 'a1', row.compoundId, 'PO')!;
    expect(rowIdentityLabel(resolved), 'New NCA · 5 · M · PO', 'distinct keys both kept');
  });
});

// ── BUG B: shared canonical dedup helper (applied at groupValues CONSTRUCTION + display) ─────────

category('context-panel: dedupGroupValuesByCanonical (BUG B source)', () => {
  test('collapses canonical-duplicate keys, first-seen wins', async () => {
    const out = dedupGroupValuesByCanonical([
      {key: 'Dose_mg_kg', value: '1'}, {key: 'Sex', value: 'F'}, {key: 'dose_mg_kg', value: '1'},
    ]);
    expect(out.length, 2, 'three entries collapse to two canonical keys');
    expect(out.map((g) => g.value).join(' · '), '1 · F', 'first-seen value kept, order preserved');
  });

  test('distinct keys pass through unchanged (no over-dedup)', async () => {
    const out = dedupGroupValuesByCanonical([{key: 'Dose_mg_kg', value: '5'}, {key: 'Sex', value: 'M'}]);
    expect(out.map((g) => g.value).join(' · '), '5 · M', 'both kept');
  });

  test('whitespace-only collision merges (casefold + strip)', async () => {
    const out = dedupGroupValuesByCanonical([{key: 'Route ', value: 'PO'}, {key: 'route', value: 'PO'}]);
    expect(out.length, 1, 'trailing-space Route and route are one canonical key');
  });
});

// ── BUG A: Forms viewer field set mirrors the "Full" toggle ─────────────────────────────────────

category('context-panel: studyParamsFormFields mirrors Full toggle (BUG A)', () => {
  const row = tierRow([{key: 'Dose_mg_kg', value: '1'}, {key: 'Sex', value: 'F'}], 'IV-bolus');
  const study = studyWith([row], ['Dose_mg_kg', 'Sex'], 'New NCA');

  test('compact profile: Mean-only SAR set, no %CV/n, no secondary, no C–t/Warnings', async () => {
    const df = buildStudyGrid(study);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (df.temp as Record<string, any>)[STUDY_EXTENDED_TAG] = false;
    const fields = studyParamsFormFields(df);
    expect(fields.includes(ANALYSIS_COL), true, 'Analysis shown');
    expect(fields.includes(ROUTE_COL), true, 'Route shown');
    expect(fields.includes('AUC∞ Mean'), true, 'lead param Mean shown');
    expect(fields.includes('AUC∞ %CV'), false, '%CV hidden in compact');
    expect(fields.includes('AUC∞ n'), false, 'n hidden in compact');
    expect(fields.includes('Cmax Mean'), false, 'secondary param hidden in compact');
    expect(fields.includes(CURVE_COL), false, 'C–t thumbnail is never a form field');
    expect(fields.includes(FLAGS_COL), false, 'Warnings hidden in compact');
  });

  test('extended profile: full Mean|%CV|n + secondary + Warnings, still no C–t', async () => {
    const df = buildStudyGrid(study);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (df.temp as Record<string, any>)[STUDY_EXTENDED_TAG] = true;
    const fields = studyParamsFormFields(df);
    expect(fields.includes('AUC∞ %CV'), true, '%CV shown in extended');
    expect(fields.includes('AUC∞ n'), true, 'n shown in extended');
    expect(fields.includes('Cmax Mean'), true, 'secondary param shown in extended');
    expect(fields.includes(FLAGS_COL), true, 'Warnings shown in extended');
    expect(fields.includes(CURVE_COL), false, 'C–t thumbnail is never a form field');
  });

  test('field set == grid visible order minus C–t (mirrors the grid exactly)', async () => {
    const df = buildStudyGrid(study);
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (df.temp as Record<string, any>)[STUDY_EXTENDED_TAG] = false;
    const fields = studyParamsFormFields(df);
    const expected = studyGridVisibleOrder(df, false).filter((n) => n !== CURVE_COL);
    expect(fields.join('|'), expected.join('|'), 'form fields track the grid visible order');
  });

  test('missing profile tag defaults to compact', async () => {
    const df = buildStudyGrid(study); // no STUDY_EXTENDED_TAG stamped
    const fields = studyParamsFormFields(df);
    expect(fields.includes('AUC∞ %CV'), false, 'defaults to compact when unstamped');
  });
});
