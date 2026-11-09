/**
 * NCA Studies — study-grid builder category (study-unified-grid AD-3).
 *
 * @see src/adapter/study-grid.ts (buildStudyGrid under test)
 *
 * The builder is a pure `NcaStudy → DG.DataFrame` projection over the cached
 * summaries — it computes nothing. These tests pin the projection contract:
 *  - AC-3a: Σ compoundTier rows + one placeholder per empty analysis; an
 *    Analysis provenance column carries the source label.
 *  - AC-3b: basis-routing — an intravascular case's CL lands in `CL`, an
 *    extravascular case's in `CL/F`; the unused column is null.
 *  - AC-3c: no internal id (MEAN_ / CV_PCT_ / N_ / CL_BASIS / PCT_EXTRAP /
 *    DOSE_NORM) leaks into an output column name; NE/warnings surface as Flags.
 */

import {expectNoRawNan} from './test-utils';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  buildStudyGrid, studyCellExclusion,
  RECOMPUTE_FLAG, STALE_ENGINE_FLAG, isStale as isStaleEngine, PARAM_TAG, STAT_TAG, NO_GROUP,
  CURVE_COL, FLAGS_COL, normalizeGroupKey, buildCurveFrame, buildOverlayCurveFrame, setNcaRepository,
  RESULTS_BASIS, studyCellBasisStale,
} from '../adapter';
import type {NcaStudy, CompoundTierRow, MeanCv, CompoundCurve} from '../adapter';
import {MemoryNcaRepository} from './memory-nca-repository';
import {
  renderStudyOverview, applyStudyGridHeaders, studyParamGloss, studyStatSentence,
  studyGridVisibleOrder, resolveClickedRow, updateCompoundPanel, renderProfileChart, buildDisclosures,
} from '../views/nca-study-overview';
import type {AnalysisEntry} from '../adapter';

function mv(mean: number, cvPct: number, n: number): MeanCv {
  return {mean, cvPct, n};
}

function compoundRow(
  compoundId: string, route: string, doseNorm: string,
  over: Partial<CompoundTierRow> = {},
): CompoundTierRow {
  // A whole-analysis case (passed `{isWholeAnalysis: true}`) has EMPTY groupValues (its
  // per-key cells render `—` unless a keyInventory entry fills them); a real grouped row
  // carries a single-key `Compound` value (study-grid-grouping-columns: the grid always
  // emits one named column per key; a composite test overrides `groupValues` via `over`).
  const whole = over.isWholeAnalysis === true;
  return {
    compoundId, isWholeAnalysis: false,
    groupValues: whole ? [] : [{key: 'Compound', value: compoundId}],
    // v9: per-tier-row candidate-column inventory (empty by default; a join test supplies one).
    keyInventory: [],
    // v10: per-tier-row passthrough (annotation) inventory (empty by default; a passthrough test supplies one).
    passthrough: [],
    n: 3, route, doseNorm,
    cmax: mv(10, 5, 3), tmax: mv(1, 20, 3), aucInf: mv(100, 8, 3), thalf: mv(2, 12, 3),
    cl: mv(1.5, 9, 3), vz: mv(3, 11, 3),
    mrt: mv(4, 10, 3), vss: mv(6, 9, 3),
    // %F defaults to null (IV rows / no pair) so the column stays hidden unless a
    // test explicitly supplies an EV fPct — the auto-hide contract (Decision 2).
    fPct: null, pctExtrap: mv(15, 4, 3),
    // v5: curve defaults to null (thumbnail absent); a test supplies one to exercise
    // the DATA_FRAME thumbnail projection.
    curve: null,
    warnings: [], ...over,
  };
}

/** A study with `analysisTier`-only stubs and the given compound tiers. */
function makeStudy(tiers: ReadonlyArray<{
  id: string; label: string; route: string; rows: ReadonlyArray<CompoundTierRow>;
}>, studyId = 's1'): NcaStudy {
  const analyses = tiers.map((t) => ({
    id: t.id, label: t.label,
    // A cache computed under the current reporting basis (nca-cl-bw-normalised-units AD-8) — an
    // unstamped entry would read as basis-stale and the grid would blank its CL / Vz / Vss.
    resultsBasis: RESULTS_BASIS,
    source: {kind: 'tableId', ref: 't'},
    columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols: [],
    sourceSignature: {
      kind: 'content-hash', rowCount: 6, mtime: null, hash: 'h',
      computedAt: '2026-07-03T00:00:00.000Z',
    },
    summaryCache: {
      analysisTier: {
        nCompounds: t.rows.length, nSubjects: 3, route: t.route, dose: '1 mg/kg',
        status: 'ok', flagSummary: 'No flags',
      },
      compoundTier: t.rows,
    },
  }));
  return {
    id: studyId, label: 'Study 1', kind: 'study',
    createdAt: '2026-07-03T00:00:00.000Z', updatedAt: '2026-07-03T00:00:00.000Z',
    access: {canEdit: true, canDelete: true, canShare: true},
    analyses,
  } as unknown as NcaStudy;
}

category('NCA Studies: study-grid builder (AD-3)', () => {
  test('AC-3a: Σ compoundTier rows + placeholder per empty analysis; Analysis provenance', async () => {
    const study = makeStudy([
      {id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg'),
        compoundRow('DGK-2', 'PO', 'per_kg'),
      ]},
      {id: 'a2', label: 'Screen B', route: 'PO', rows: []}, // empty → 1 placeholder
    ]);
    const df = buildStudyGrid(study);
    // AC-13 (F4 Arm 1): this producer's REAL output must carry no raw IEEE NaN.
    // Asserted here rather than in a duplicate fixture — the frame already exists.
    expectNoRawNan(df, 'study-grid.buildStudyGrid');
    expect(df.rowCount, 3, '2 compound rows + 1 placeholder');
    const analysisCol = df.col('Analysis')!;
    expect(analysisCol.get(0), 'Screen A', 'provenance carried');
    expect(analysisCol.get(2), 'Screen B', 'placeholder carries its analysis label');
    // Placeholder row is flagged to recompute (visible, never a wrong row).
    expect(df.col(FLAGS_COL)!.get(2), RECOMPUTE_FLAG, 'empty analysis → recompute affordance');
  });

  test('engine-stale rows (v14 F-N.6): the coordinator predicate marks the ROW in the Warnings cell', async () => {
    // Two analyses: a1 computed by the installed engine, a2 by an older one (or unknown —
    // a migrated pre-v14 file). The grid consults the INJECTED predicate per row; the default
    // (no predicate) marks nothing, so the frame-shape callers stay unchanged.
    const study = makeStudy([
      {id: 'a1', label: 'Fresh', route: 'IV-bolus', rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg')]},
      {id: 'a2', label: 'Old engine', route: 'PO', rows: [
        compoundRow('DGK-2', 'PO', 'per_kg'),
        compoundRow('DGK-3', 'PO', 'per_kg', {warnings: [{code: 'AUC_EXTRAP_HIGH', message: '%AUCextrap high'}]}),
      ]},
      {id: 'a3', label: 'Old + empty', route: 'PO', rows: []},
    ]);
    const stamped = {a1: '0.11.0', a2: '0.10.0', a3: null} as const;
    // The ONE predicate (SD-9), with the installed engine pinned to 0.11.0.
    const isStale = (a: AnalysisEntry) =>
      isStaleEngine({engineVersion: stamped[a.id as keyof typeof stamped], resultsBasis: a.resultsBasis,
        hasCache: true}, '0.11.0');
    const df = buildStudyGrid(study, isStale);
    const flags = (i: number) => String(df.col(FLAGS_COL)!.get(i));
    expect(flags(0), '', 'fresh row: clean');
    expect(flags(1), STALE_ENGINE_FLAG, 'stale row with no other flag: the cue alone');
    expect(flags(2).split('\n')[0], STALE_ENGINE_FLAG, 'stale row with a warning: the cue LEADS');
    expect(flags(2).includes('%AUCextrap high'), true, '…and the warning is kept');
    expect(flags(3), `${STALE_ENGINE_FLAG}\n${RECOMPUTE_FLAG}`, 'stale placeholder: cue + recompute affordance');
    // Same wording as the overview's per-analysis badge, so the two cues read as one.
    expect(STALE_ENGINE_FLAG, 'stale (refreshes on open)');
    // Default predicate: nothing marked.
    expect(String(buildStudyGrid(study).col(FLAGS_COL)!.get(1)), '', 'no predicate → no stale marks');
  });

  // nca-cl-bw-normalised-units AD-8 (architect B3) / AC-7.3 / AC-7.4.
  test('AC-7.4: CL / Vz / Vss headers carry the fixed per-kg units (no doseNorm-derived unit)', async () => {
    const study = makeStudy([{id: 'a1', label: 'A', route: 'IV-bolus', rows: [
      compoundRow('DGK-1', 'IV-bolus', 'absolute'), compoundRow('DGK-2', 'IV-bolus', 'per_kg')]}]);
    const df = buildStudyGrid(study);
    expect(df.col('CL Mean')!.getTag('units'), 'mL/min/kg', 'CL per kg, whatever the engine dose basis');
    expect(df.col('Vz Mean')!.getTag('units'), 'L/kg', 'Vz per kg');
    expect(df.col('Vss Mean')!.getTag('units'), 'L/kg', 'Vss per kg');
  });

  test('AC-7.3 / B3: a basis-stale analysis keeps its row, its stale line and AUC/Cmax — CL/Vz/Vss blank', async () => {
    const study = makeStudy([
      {id: 'fresh', label: 'Fresh', route: 'IV-bolus', rows: [compoundRow('DGK-1', 'IV-bolus', 'absolute')]},
      {id: 'old', label: 'Old basis', route: 'IV-bolus', rows: [compoundRow('DGK-2', 'IV-bolus', 'absolute')]},
    ]);
    // The old analysis was cached before the per-kg basis: no stamp.
    (study.analyses[1] as {resultsBasis?: string | null}).resultsBasis = null;
    const df = buildStudyGrid(study, (a) => isStaleEngine(
      {engineVersion: '0.11.0', resultsBasis: a.resultsBasis, hasCache: true}, '0.11.0'));
    expect(df.col('CL Mean')!.isNone(0), false, 'fresh row: CL shown');
    expect(df.col('CL Mean')!.isNone(1), true, 'basis-stale row: CL blank (no relabelled absolute mean)');
    expect(df.col('Vz Mean')!.isNone(1), true, 'Vz blank');
    expect(df.col('Vss Mean')!.isNone(1), true, 'Vss blank');
    expect(df.col('AUC∞ Mean')!.isNone(1), false, 'basis-independent AUC∞ stays');
    expect(df.col('Cmax Mean')!.isNone(1), false, 'basis-independent Cmax stays');
    expect(String(df.col(FLAGS_COL)!.get(1)).split('\n')[0], STALE_ENGINE_FLAG, 'the ONE stale wording leads');
    expect(studyCellBasisStale(df, 'CL Mean', 1), true, 'the blanked cell is identifiable for its tooltip');
    expect(studyCellBasisStale(df, 'AUC∞ Mean', 1), false, 'not on a basis-independent cell');
    expect(studyCellBasisStale(df, 'CL Mean', 0), false, 'not on a fresh row');
  });

  test('B3: a study whose EVERY analysis is basis-stale keeps CL / Vz / Vss as blank columns', async () => {
    const study = makeStudy([
      {id: 'old1', label: 'Old 1', route: 'IV-bolus', rows: [compoundRow('DGK-1', 'IV-bolus', 'absolute')]},
      {id: 'old2', label: 'Old 2', route: 'IV-bolus', rows: [compoundRow('DGK-2', 'IV-bolus', 'absolute')]},
    ]);
    for (const a of study.analyses) (a as {resultsBasis?: string | null}).resultsBasis = null;
    const df = buildStudyGrid(study, (a) => isStaleEngine(
      {engineVersion: '0.11.0', resultsBasis: a.resultsBasis, hasCache: true}, '0.11.0'));
    for (const p of ['CL', 'Vz', 'Vss']) {
      expect(df.col(`${p} Mean`) !== null, true, `${p} column kept (measured at QA: it vanished)`);
      expect(df.col(`${p} Mean`)!.isNone(0) && df.col(`${p} Mean`)!.isNone(1), true, `${p} cells blank`);
      expect(studyCellBasisStale(df, `${p} Mean`, 0), true, `${p} blank cell carries the stale tooltip`);
    }
    expect(df.col('CL/F Mean'), null, 'a group with nothing to blank still drops (all-IV study: no CL/F)');
  });

  test('per-key column: real group-by value shown; whole-analysis case → "—" (not the analysis label)', async () => {
    // study-grid-grouping-columns E1: the grid ALWAYS emits one named column per grouping
    // key (the generic `Group` column is retired). "Batch" groups by Compound → a `Compound`
    // column with the compound ids; "Multi-route" is whole-analysis (empty groupValues) →
    // its `Compound` cell is `—` (no group, no inventory value), never the analysis label.
    const study = makeStudy([
      {id: 'a1', label: 'Batch', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg'),
        compoundRow('DGK-2', 'IV-bolus', 'per_kg'),
      ]},
      {id: 'a2', label: 'Multi-route', route: 'PO', rows: [
        compoundRow('Multi-route', 'PO', 'absolute', {isWholeAnalysis: true}),
        compoundRow('Multi-route', 'IV-bolus', 'absolute', {isWholeAnalysis: true}),
      ]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.col('Compound') !== null, true, 'a named per-key Compound column (not a generic Group)');
    const g = df.col('Compound')!;
    expect(g.get(0), 'DGK-1', 'real compound value shown');
    expect(g.get(1), 'DGK-2');
    expect(g.get(2), NO_GROUP, 'whole-analysis case shows "—", not the analysis label');
    expect(g.get(3), NO_GROUP);
    // The generic `Group` column is retired — it must never be emitted.
    expect(df.col('Group'), null, 'no generic Group column (retired, design-decisions §4)');
    // Provenance is still carried by the Analysis column (no duplication lost).
    expect(df.col('Analysis')!.get(2), 'Multi-route', 'Analysis column carries the source label');
  });

  test('per-key column: renamed whole-analysis stays "—" (regression, group-col-stale-on-rename)', async () => {
    // The bug: an ungrouped analysis computed as "New NCA" then RENAMED to "theoph"
    // showed the stale "New NCA" instead of "—". The tier's compoundId is frozen at
    // compute time ("New NCA"); the entry.label is now "theoph". With the v7
    // isWholeAnalysis flag + empty groupValues the per-key cell is label-independent, so
    // the mismatch is irrelevant and the (single) key column stays "—".
    const study = makeStudy([
      {id: 'a1', label: 'theoph', route: 'PO', rows: [
        // A grouped sibling analysis creates a `Compound` key column; the ungrouped
        // theoph row must show `—` there (empty groupValues, no inventory).
        compoundRow('DGK-9', 'IV-bolus', 'per_kg'),
        // compoundId frozen to the PRE-rename label; flag is the durable sentinel.
        compoundRow('New NCA', 'PO', 'absolute', {isWholeAnalysis: true}),
      ]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.col('Compound')!.get(1), NO_GROUP,
      'renamed whole-analysis stays "—" (no stale label leak)');
    // And the drill-in still resolves the row post-rename via the FROZEN compoundId
    // back-pointer (v8 — the grid carries the compute-time compoundId in a hidden
    // column; for a whole-analysis case that is the pre-rename label, which is exactly
    // what the tier row stores, so the match is label-independent).
    const hit = resolveClickedRow(study, 'a1', 'New NCA', 'PO');
    expect(hit !== null, true, 'drill-in resolves the renamed whole-analysis row');
    expect(hit!.row.compoundId, 'New NCA', 'resolves the frozen-compoundId tier row');
  });

  test('E2 typographic merge: Strain/STRAIN/"Strain " → ONE column (AC-E2a)', async () => {
    // Three analyses grouping by the SAME key with different case/whitespace must merge
    // into ONE column on the typographic canonical id; the header is the first-seen label.
    // Uses a covariate key (Strain) that does NOT collide with a reserved grid column —
    // `route` normalizes to the reserved `Route` column (case-insensitive), a different case.
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('x', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Strain', value: 'Wistar'}]})]},
      {id: 'a2', label: 'B', route: 'IV-bolus', rows: [
        compoundRow('y', 'IV-bolus', 'per_kg', {groupValues: [{key: 'STRAIN', value: 'SD'}]})]},
      {id: 'a3', label: 'C', route: 'IV-bolus', rows: [
        compoundRow('z', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Strain ', value: 'Wistar'}]})]},
    ]);
    const df = buildStudyGrid(study);
    const family = df.columns.names().filter((nm) => normalizeGroupKey(nm) === 'strain');
    expect(family.length, 1, 'Strain/STRAIN/"Strain " merge into ONE column (whitespace/case only)');
    expect(family[0], 'Strain', 'header = first-seen friendly label');
    const col = df.col(family[0])!;
    expect(col.get(0), 'Wistar', 'row 0 (analysis A) value');
    expect(col.get(1), 'SD', 'row 1 (analysis B) value under the merged column');
    expect(col.get(2), 'Wistar', 'row 2 (analysis C) value under the merged column');
  });

  test('E2 punctuation boundary: Dose(mg/kg)==Dose (mg/kg); ≠ Dose; Day -1 ≠ Day 1; Sex ≠ Gender', async () => {
    // Whitespace-only strip merges the unit-suffixed dose variants but KEEPS punctuation, so
    // bare Dose stays separate and signed study-days / semantic synonyms never over-merge.
    expect(normalizeGroupKey('Dose(mg/kg)'), normalizeGroupKey('Dose (mg/kg)'),
      'Dose(mg/kg) and Dose (mg/kg) normalize equal (differ only by a space)');
    expect(normalizeGroupKey('Dose') !== normalizeGroupKey('Dose (mg/kg)'), true,
      'bare Dose is a SEPARATE column (punctuation preserved)');
    expect(normalizeGroupKey('Day -1') !== normalizeGroupKey('Day 1'), true,
      'Day -1 ≠ Day 1 — a CDISC baseline-vs-dosing-day distinction must NOT merge (AC-E2d)');
    expect(normalizeGroupKey('Sex') !== normalizeGroupKey('Gender'), true,
      'Sex ≠ Gender — semantic synonyms are out of typographic scope (AC-E2c)');
    // End-to-end: Dose (mg/kg) and Dose(mg/kg) → one column; bare Dose → a second.
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('x', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Dose (mg/kg)', value: '10'}]})]},
      {id: 'a2', label: 'B', route: 'IV-bolus', rows: [
        compoundRow('y', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Dose(mg/kg)', value: '30'}]})]},
      {id: 'a3', label: 'C', route: 'IV-bolus', rows: [
        compoundRow('z', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Dose', value: '5'}]})]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.columns.names().filter((nm) => normalizeGroupKey(nm) === 'dose(mg/kg)').length, 1,
      'Dose (mg/kg)/Dose(mg/kg) → one merged column');
    expect(df.columns.names().filter((nm) => normalizeGroupKey(nm) === 'dose').length, 1,
      'bare Dose → its OWN separate column (two columns total)');
  });

  test('E4 join render: grouped → value; present-but-ungrouped → counted; absent → — (AC-E4a/b)', async () => {
    // Analysis A groups by Sex (creates the Sex column). Analysis B does NOT group by Sex
    // but its subpopulation carries Sex∈{F×15,M×15} in its keyInventory → counted cell.
    // Analysis C has no Sex → —.
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('x', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Sex', value: 'F'}]})]},
      {id: 'a2', label: 'B', route: 'IV-bolus', rows: [
        compoundRow('B', 'IV-bolus', 'per_kg', {isWholeAnalysis: true,
          keyInventory: [{key: 'Sex', values: [{value: 'F', count: 15}, {value: 'M', count: 15}]}]})]},
      {id: 'a3', label: 'C', route: 'IV-bolus', rows: [
        compoundRow('C', 'IV-bolus', 'per_kg', {isWholeAnalysis: true, keyInventory: []})]},
    ]);
    const df = buildStudyGrid(study);
    const sex = df.col('Sex')!;
    expect(sex.get(0), 'F', 'A grouped by Sex → its group value (precedence step 1)');
    expect(sex.get(1), 'F (15); M (15)', 'B present-but-ungrouped → COUNTED composition (AC-E4a)');
    expect(sex.get(2), NO_GROUP, 'C has no Sex → — (genuinely absent, AC-E4b)');
  });

  test('E4 pooled row discloses its imbalance (AC-E4e); single value+count never blank (AC-E4d)', async () => {
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('x', 'IV-bolus', 'per_kg', {groupValues: [{key: 'Sex', value: 'F'}]})]},
      // theoph: whole-analysis, imbalanced F:27 M:3 — the pooled thumbnail hides this; the
      // join cell must surface it (thesis #4), NOT collapse to `—`.
      {id: 'a2', label: 'theoph', route: 'PO', rows: [
        compoundRow('theoph', 'PO', 'absolute', {isWholeAnalysis: true,
          keyInventory: [{key: 'Sex', values: [{value: 'F', count: 27}, {value: 'M', count: 3}]}]})]},
      {id: 'a3', label: 'allM', route: 'IV-bolus', rows: [
        compoundRow('allM', 'IV-bolus', 'per_kg', {isWholeAnalysis: true,
          keyInventory: [{key: 'Sex', values: [{value: 'M', count: 30}]}]})]},
    ]);
    const df = buildStudyGrid(study);
    const sex = df.col('Sex')!;
    expect(sex.get(1), 'F (27); M (3)', 'pooled theoph discloses its imbalance (AC-E4e), not —');
    expect(sex.get(2), 'M (30)', 'all-M single value + count (AC-E4d), never blank');
  });

  test('E1 zero key columns when NO analysis groups (whole-analysis-only study)', async () => {
    // theoph-only, whole-analysis → the union of grouping keys is empty → NO key columns
    // at all (Analysis / Route / params only). The generic Group column is retired.
    const study = makeStudy([
      {id: 'a1', label: 'theoph', route: 'PO', rows: [
        compoundRow('theoph', 'PO', 'absolute', {isWholeAnalysis: true})]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.col('Group'), null, 'no generic Group column (retired)');
    expect(df.col('Compound'), null, 'no key column at all (nothing groups → zero key columns, AC-E1a)');
    expect(df.col('Analysis') !== null, true, 'Analysis provenance still present');
    expect(df.col('Route') !== null, true, 'Route still present');
  });

  test('E1 regression: a CompoundID key (case-collides with the hidden back-pointer) → clean column', async () => {
    // DG column names are CASE-INSENSITIVE: a user key 'CompoundID' must not collide with the
    // internal hidden compoundId back-pointer at DataFrame.fromColumns (the regression the live
    // grok test caught — fixture-05 groups by CompoundID). buildStudyGrid must not throw and
    // must produce a clean 'CompoundID' column.
    const study = makeStudy([
      {id: 'a1', label: 'Screen', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg', {groupValues: [{key: 'CompoundID', value: 'DGK-1'}]}),
        compoundRow('DGK-2', 'IV-bolus', 'per_kg', {groupValues: [{key: 'CompoundID', value: 'DGK-2'}]}),
      ]},
    ]);
    const df = buildStudyGrid(study); // must NOT throw "Column named 'compoundId' already exists"
    expect(df.col('CompoundID') !== null, true, 'a key named CompoundID gets a clean named column');
    expect(df.col('CompoundID')!.get(0), 'DGK-1', 'the CompoundID column carries the key value');
  });

  test('AC-3b: basis-routing — IV case fills CL, EV case fills CL/F (unused column null)', async () => {
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('IV-cpd', 'IV-bolus', 'per_kg', {cl: mv(2.0, 5, 3)}),
        compoundRow('PO-cpd', 'PO', 'per_kg', {cl: mv(8.0, 5, 3)}),
      ]},
    ]);
    const df = buildStudyGrid(study);
    const clMean = df.col('CL Mean')!;
    const clfMean = df.col('CL/F Mean')!;
    // Row 0 = IV: CL filled, CL/F null. Row 1 = PO: CL null, CL/F filled.
    expect(clMean.isNone(0), false, 'IV case fills CL');
    expect(clMean.get(0), 2.0, 'IV CL value');
    expect(clMean.isNone(1), true, 'PO case leaves true CL null');
    expect(clfMean.isNone(1), false, 'PO case fills CL/F');
    expect(clfMean.get(1), 8.0, 'PO CL/F value');
    expect(clfMean.isNone(0), true, 'IV case leaves CL/F null');
  });

  test('AC-4a scaffold: each multi-stat param spans exactly Mean|%CV|n; %AUCextrap standalone', async () => {
    // The native column groups (configureStudyGrid → df.temp[".columnGroups"]) are
    // built by grouping columns on PARAM_TAG where STAT_TAG is set. This pins that
    // tag contract on the builder output: a group per triple-param with exactly
    // {Mean,%CV,n}, and %AUCextrap tagged as a param but WITHOUT a stat (standalone).
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('IV-cpd', 'IV-bolus', 'per_kg'),
        compoundRow('PO-cpd', 'PO', 'per_kg'),
      ]},
    ]);
    const df = buildStudyGrid(study);
    const groups = new Map<string, Set<string>>();
    let extrapStandalone = false;
    for (const col of df.columns) {
      const param = col.getTag(PARAM_TAG);
      const stat = col.getTag(STAT_TAG);
      if (param === null || param === '') continue;
      if (stat !== null && stat !== '') {
        if (!groups.has(param)) groups.set(param, new Set());
        groups.get(param)!.add(stat);
      } else if (param === '%AUCextrap') extrapStandalone = true;
    }
    // Both CL (IV case) and CL/F (PO case) are present; each is a full triple.
    for (const p of ['Cmax', 'AUC∞', 't½', 'CL', 'CL/F', 'Vz', 'Vz/F', 'MRT', 'Vss']) {
      const stats = groups.get(p);
      expect(stats !== undefined, true, `param group '${p}' present`);
      expect(stats!.has('Mean') && stats!.has('%CV') && stats!.has('n'), true,
        `param '${p}' spans exactly Mean|%CV|n`);
      expect(stats!.size, 3, `param '${p}' has no extra stat columns`);
    }
    expect(extrapStandalone, true, '%AUCextrap is a standalone param column (no stat → not grouped)');
    expect(groups.has('%AUCextrap'), false, '%AUCextrap is NOT a Mean|%CV|n group (B3)');
  });

  test('Vss group drops on an all-extravascular screen (IV-only param, hide-empty); MRT stays', async () => {
    // Vss is IV-only (sci-comp NaN for EV — an apparent Vss/F would be absorption-
    // confounded). The grid routes it by r.iv, so an all-PO screen self-drops the Vss
    // group (addGroup's all-empty guard), exactly like Vz→Vz/F. MRT is route-uniform
    // and stays. Exercises the EV-drop path the IV+PO fixtures never hit.
    const study = makeStudy([
      {id: 'a1', label: 'PO-study', route: 'PO', rows: [
        compoundRow('PO-cpd-1', 'PO', 'per_kg'),
        compoundRow('PO-cpd-2', 'PO', 'per_kg'),
      ]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.col('Vss Mean'), null, 'Vss column absent on an all-EV screen (IV-only, hide-empty)');
    expect(df.col('Vz Mean'), null, 'Vz (IV basis) absent — EV rows route to Vz/F');
    expect(df.col('Vz/F Mean') !== null, true, 'Vz/F present (extravascular basis)');
    expect(df.col('MRT Mean') !== null, true, 'MRT present (route-uniform, not IV-gated)');
  });

  test('%F column auto-hides when no case carries a finite %F (Decision 2: unpaired → no column)', async () => {
    // Every row's fPct is null (the compoundRow default) — mirrors an all-IV screen
    // (fixture 05: 20 IV compounds, no EV arm anywhere). The addGroup all-empty guard
    // must omit the WHOLE %F column, exactly like an all-EV screen drops Vss.
    const study = makeStudy([
      {id: 'a1', label: 'All-IV', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg'),
        compoundRow('DGK-2', 'IV-bolus', 'per_kg'),
      ]},
    ]);
    const df = buildStudyGrid(study);
    expect(df.col('%F Mean'), null, '%F column absent when no case carries a finite %F');
    expect(df.col('%F %CV'), null, '%F %CV column absent too');
  });

  test('%F column appears when a paired EV row carries a finite %F; IV rows blank (null)', async () => {
    // A paired screen: one EV compound with a real %F (mean 42%) + one IV compound
    // (fPct null → F ≡ 1, never shown). The %F column must appear (Decision 2: a
    // paired arm exists), the EV row carries its value, the IV row is blank.
    const study = makeStudy([
      {id: 'a1', label: 'Paired', route: 'IV-bolus', rows: [
        compoundRow('IV-ref', 'IV-bolus', 'per_kg'), // fPct null (default)
        compoundRow('EV-cpd', 'PO', 'per_kg', {fPct: mv(42, 18, 3)}),
      ]},
    ]);
    const df = buildStudyGrid(study);
    const fMean = df.col('%F Mean');
    expect(fMean !== null, true, '%F column present once a paired EV compound carries a value');
    expect(fMean!.isNone(0), true, 'IV row leaves %F blank (null — F ≡ 1 by definition)');
    expect(fMean!.isNone(1), false, 'EV row fills %F');
    expect(fMean!.get(1), 42, 'EV %F value carried');
    // Gloss + header resolve (map-drift guard for the new PARAM_TAG).
    expect(df.col('%F Mean')!.getTag(PARAM_TAG), '%F', '%F tagged as its own param');
    const gloss = studyParamGloss('%F');
    expect(gloss !== null && gloss.length > 0, true, '%F resolves to a non-empty gloss (no map drift)');
  });

  test('headers: param columns get "stat param"; untagged columns keep their name (no "undefined")', async () => {
    // Regression: getTag returns `undefined` (not null/'') for an unset tag, so a
    // naive null-only guard writes "undefined undefined" onto the leading
    // Analysis/Compound/Route/Flags headers. This pins that they stay untouched
    // and the param columns get the qualified header.
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg')]},
    ]);
    const df = buildStudyGrid(study);
    applyStudyGridHeaders(df);
    // Untagged leading + per-key columns: no friendlyName written → NEVER "undefined ...".
    for (const name of ['Analysis', 'Compound', 'Route', FLAGS_COL]) {
      const fn = df.col(name)!.getTag('friendlyName');
      expect(fn === null || fn === undefined || fn === '', true,
        `'${name}' must keep its own header (got friendlyName='${String(fn)}')`);
    }
    // Param columns: qualified header.
    expect(df.col('Cmax Mean')!.getTag('friendlyName'), 'Mean Cmax', 'Cmax mean header');
    expect(df.col('CL n')!.getTag('friendlyName'), 'n CL', 'CL n header');
    // %AUCextrap standalone keeps its own name (no friendlyName).
    const extrapFn = df.col('%AUCextrap')!.getTag('friendlyName');
    expect(extrapFn === null || extrapFn === undefined || extrapFn === '', true,
      '%AUCextrap standalone keeps its column name');
  });

  test('header tooltips: every PARAM_TAG + STAT_TAG resolves (map-drift guard)', async () => {
    // Silent-drift guard (study-grid-header-tooltips): the header tooltip resolvers
    // (studyParamGloss / studyStatSentence) are backed by hand-maintained maps
    // (STUDY_PARAM_TO_AGG / BASIS_GLOSS) bridging DISPLAY labels to the shared
    // AGG_PARAM gloss. If buildStudyGrid ever emits a new PARAM_TAG not covered by
    // either map, that header would silently fall back to title-only. The IV+PO
    // fixture forces every param into existence (both CL & CL/F, Vz & Vz/F).
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        compoundRow('IV-cpd', 'IV-bolus', 'per_kg'),
        compoundRow('PO-cpd', 'PO', 'per_kg'),
      ]},
    ]);
    const df = buildStudyGrid(study);
    for (const col of df.columns) {
      const param = col.getTag(PARAM_TAG);
      const stat = col.getTag(STAT_TAG);
      if (param === null || param === undefined || param === '') continue;
      const gloss = studyParamGloss(param);
      expect(gloss !== null && gloss.length > 0, true,
        `param '${param}' must resolve to a non-empty gloss (map drift?)`);
      if (stat !== null && stat !== undefined && stat !== '') {
        const sentence = studyStatSentence(stat, param);
        expect(sentence !== null && sentence.length > 0, true,
          `stat '${stat}' (param '${param}') must resolve to a non-empty sentence`);
      }
    }
    // The `n` sentence must NOT cite the aggregated grid's cohort-N column, which
    // does not exist in this grid (decision-audit fix: no dangling reference).
    expect(studyStatSentence('n', 'Cmax')!.includes('group\'s N'), false,
      'study n sentence must not reference a nonexistent cohort-N column');
  });

  test('S3 chemist profile: compact = lead-trio + CL/t½ means only; Extended reveals the full set', async () => {
    // Paired IV+EV screen so every param group exists (CL & CL/F, Vz & Vz/F, %F).
    const study = makeStudy([
      {id: 'a1', label: 'SAR', route: 'IV-bolus', rows: [
        compoundRow('IV-ref', 'IV-bolus', 'per_kg'),
        compoundRow('EV-cpd', 'PO', 'per_kg', {fPct: mv(42, 18, 3)}),
      ]},
    ]);
    const df = buildStudyGrid(study);

    const compact = studyGridVisibleOrder(df, false);
    // Leading context always present; the per-key grouping column is headed by its key's
    // friendly label (`Compound` here) — the generic `Group` column is retired
    // (design-decisions §4, study-grid-grouping-columns).
    expect(compact.includes('Compound'), true, 'per-key grouping column shown (design-decisions §4)');
    expect(compact.includes('Group'), false, 'no generic Group column (retired)');
    expect(compact.includes('Analysis'), true, 'Analysis provenance shown');
    // Lead trio + CL/t½ MEAN columns present, in order (AUC∞ · MRT · %F │ CL · t½).
    const leadMeans = ['AUC∞ Mean', 'MRT Mean', '%F Mean', 'CL Mean', 't½ Mean'];
    for (const m of leadMeans)
      expect(compact.includes(m), true, `compact shows ${m}`);
    // Header casing: the %F column is %-prefixed (design-decisions §6).
    expect(compact.includes('%F Mean'), true, '%F header is %-prefixed (not F%)');
    // Compact HIDES: %CV / n of any param, Cmax, Vss, Vz, %AUCextrap.
    for (const hidden of ['AUC∞ %CV', 'AUC∞ n', '%F n', 'Cmax Mean', 'Vss Mean',
      'Vz Mean', '%AUCextrap'])
      expect(compact.includes(hidden), false, `compact hides ${hidden}`);
    // Ordering: lead trio precedes CL and t½.
    expect(compact.indexOf('AUC∞ Mean') < compact.indexOf('CL Mean'), true,
      'AUC∞ before CL (lead trio first)');
    expect(compact.indexOf('%F Mean') < compact.indexOf('t½ Mean'), true, '%F before t½');
    // Warnings is Full-only (user-directed): hidden in compact.
    expect(compact.includes(FLAGS_COL), false, 'compact hides the Warnings column');

    const extended = studyGridVisibleOrder(df, true);
    expect(extended.includes(FLAGS_COL), true, 'Full profile shows the Warnings column');
    // Extended reveals the full Mean|%CV|n triples + the secondary params.
    for (const shown of ['AUC∞ Mean', 'AUC∞ %CV', 'AUC∞ n', '%F %CV', '%F n',
      'Cmax Mean', 'Vss Mean', 'Vz Mean', '%AUCextrap'])
      expect(extended.includes(shown), true, `extended shows ${shown}`);
    // Extended is a strict superset of compact (a hidden-then-shown column, S3 DoD).
    for (const c of compact)
      expect(extended.includes(c), true, `extended keeps compact column ${c}`);
    expect(extended.length > compact.length, true, 'extended reveals more columns than compact');
    // Secondary-zone order per spec S3: `Cmax · Vss · Vz · %AUCextrap` (Vss before Vz).
    expect(extended.indexOf('Cmax Mean') < extended.indexOf('Vss Mean'), true, 'Cmax before Vss');
    expect(extended.indexOf('Vss Mean') < extended.indexOf('Vz Mean'), true,
      'Vss before Vz (spec secondary order)');
    expect(extended.indexOf('Vz Mean') < extended.indexOf('%AUCextrap'), true, 'Vz before %AUCextrap');
  });

  test('AC-3c: no internal id leaks into a column name; NE + warnings surface as Flags', async () => {
    const study = makeStudy([
      {id: 'a1', label: 'A', route: 'IV-bolus', rows: [
        // An unfit case: AUC∞ NaN → NE; plus a warning (v6: {code, message}).
        compoundRow('cpd', 'IV-bolus', 'per_kg', {
          aucInf: mv(NaN, NaN, 0),
          warnings: [{code: 'LOW_R_SQUARED', message: 'Terminal-phase R² below 0.90 (fragile λz fit).'}],
        }),
      ]},
    ]);
    const df = buildStudyGrid(study);
    const banned = ['MEAN_', 'CV_PCT_', 'N_', 'CL_BASIS', 'PCT_EXTRAP', 'DOSE_NORM'];
    for (const name of df.columns.names()) {
      for (const b of banned)
        expect(name.includes(b), false, `column '${name}' must not leak internal id '${b}'`);
    }
    // %-prefix display label present as a single column (design B3).
    expect(df.col('%AUCextrap') !== null, true, '%AUCextrap single display column present');
    expect(df.col('%AUCextrap Mean'), null, '%AUCextrap is NOT a Mean|%CV|n triple (B3)');
    // v6: the Flags cell carries the human-readable text IN the cell (exportable) —
    // the NE token AND the warning message, `; `-joined. No separate hidden Warnings
    // column (removed so table export stays clean, user-directed 2026-07-07).
    const flags = df.col(FLAGS_COL)!.get(0) as string;
    expect(flags.includes('NE'), true, 'unfit AUC∞ → NE token in the Warnings cell');
    expect(flags.includes('fragile λz fit'), true, 'the human-readable message is IN the cell (exportable)');
    expect(FLAGS_COL, 'Warnings', 'the column is named "Warnings" (renamed from "Flags")');
    expect(df.col('Flags'), null, 'no legacy "Flags" column');
  });

  test('S4b: inline C–t thumbnail is a DATA_FRAME column whose cell carries the compound curve points', async () => {
    // A compound row WITH a cached curve (2 dose series) → the projection emits a
    // DATA_FRAME column whose cell holds a small (Time/Conc/Series) frame with one
    // row per (series × point). Verified against the real projection output (rule 16).
    const curve = {
      series: [
        {group: '10', time: [0, 1, 2], geoMean: [100, 50, 25]},
        {group: '30', time: [0, 1], geoMean: [300, 150]},
      ],
    };
    const study = makeStudy([{
      id: 'an-1', label: 'Screen', route: 'IV-bolus',
      rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg', {curve})],
    }]);
    const df = buildStudyGrid(study);
    const curveCol = df.col(CURVE_COL);
    expect(curveCol !== null, true, 'C–t DATA_FRAME column present when a case carries a curve');
    expect(curveCol!.type, DG.COLUMN_TYPE.DATA_FRAME, 'C–t column is a DATA_FRAME column');
    const cell = curveCol!.get(0) as DG.DataFrame;
    expect(cell !== null && cell !== undefined, true, 'curve cell is a DataFrame (not null)');
    // 3 points (series 10) + 2 points (series 30) = 5 rows; one line per dose (BUG-08).
    expect(cell.rowCount, 5, 'cell frame carries every curve point across both series');
    for (const c of ['Time', 'Conc', 'Series'])
      expect(cell.col(c) !== null, true, `curve cell frame has a '${c}' column`);
    const seriesVals = new Set<string>();
    const sCol = cell.col('Series')!;
    for (let i = 0; i < cell.rowCount; i++) seriesVals.add(sCol.getString(i));
    expect(seriesVals.has('10') && seriesVals.has('30'), true, 'both dose series present in the cell frame');
  });

  test('S4b: C–t column is absent when no case carries a curve (hide-all-empty parity)', async () => {
    // Every row's curve is null (default) → the DATA_FRAME column is omitted, so a
    // pre-v5 / no-observations study shows no empty thumbnail column.
    const study = makeStudy([{
      id: 'an-1', label: 'Screen', route: 'IV-bolus',
      rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg')], // curve defaults to null
    }]);
    const df = buildStudyGrid(study);
    expect(df.col(CURVE_COL), null, 'no curve anywhere → C–t column dropped (hide-all-empty)');
  });

  test('S4b: C–t sits after Group/Route in the visible column order (compact and extended)', async () => {
    const curve = {series: [{group: '10', time: [0, 1], geoMean: [100, 50]}]};
    const study = makeStudy([{
      id: 'an-1', label: 'Screen', route: 'IV-bolus',
      rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg', {curve})],
    }]);
    const df = buildStudyGrid(study);
    for (const extended of [false, true]) {
      const order = studyGridVisibleOrder(df, extended);
      const ci = order.indexOf(CURVE_COL);
      expect(ci >= 0, true, `C–t is visible in the ${extended ? 'extended' : 'compact'} profile`);
      expect(ci > order.indexOf('Compound'), true, 'C–t comes after the per-key grouping column');
      // C–t leads the parameter zone: it precedes the first param column (AUC∞ Mean).
      expect(ci < order.indexOf('AUC∞ Mean'), true, 'C–t precedes the parameter columns');
    }
  });
});

category('NCA Studies: study overview grid render (AD-4)', () => {
  /** Render the overview of `study` from an in-memory store (decision 7: the grid is logic, not storage). */
  const renderFromMemory = async (study: NcaStudy): Promise<DG.View> => {
    const mem = new MemoryNcaRepository();
    mem.seedStudy(study);
    const prev = setNcaRepository(mem);
    try {
      return await renderStudyOverview(study.id);
    } finally {
      setNcaRepository(prev);
    }
  };

  test('AC-4: renderStudyOverview mounts the unified grid without error', async () => {
    // A study with a computed compound tier so the grid is non-empty.
    const study = makeStudy([
      {id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg'),
      ]},
    ], `test-grid-render-${Date.now()}`);
    const view = await renderFromMemory(study);
    // The non-empty grid renders its sized container (empty-state has none).
    const box = view.root.querySelector('.nca-study-grid');
    expect(box !== null, true, 'unified grid container rendered');
  });

  test('S4b: overview renders the inline C–t thumbnail column without error (native plot.line)', async () => {
    // A curve-bearing tier exercises the render path end-to-end: configureStudyGrid
    // sets the DATA_FRAME column to an html cell + widened row/column, and
    // onCellPrepare renders the native plot.line. Asserts the render mounts and the
    // C–t column survives the compact column profile (visible), never throwing.
    const curve = {series: [{group: '1', time: [0, 1, 2, 4], geoMean: [100, 60, 30, 8]}]};
    const study = makeStudy([
      {id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg', {curve}),
      ]},
    ], `test-grid-thumb-${Date.now()}`);
    const view = await renderFromMemory(study);
    const box = view.root.querySelector('.nca-study-grid');
    expect(box !== null, true, 'unified grid container rendered with a thumbnail column');
  });
});

// ── S4c: row-click drill-in (compound detail context panel) ──────────────────
category('NCA Studies: compound drill-in (chemist-sar-grid S4c)', () => {
  test('resolveClickedRow maps (analysisId, Group, Route) back to the cached tier row', async () => {
    const study = makeStudy([{
      id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
        compoundRow('DGK-1', 'IV-bolus', 'per_kg'),
        compoundRow('DGK-2', 'IV-bolus', 'per_kg'),
      ],
    }]);
    // A real group value → matched by compound + route.
    const hit = resolveClickedRow(study, 'a1', 'DGK-2', 'IV-bolus');
    expect(hit !== null, true, 'grouped row resolves');
    expect(hit!.row.compoundId, 'DGK-2', 'resolves the exact compound clicked');
    expect(hit!.entry.id, 'a1', 'resolves the owning analysis');
    // Unknown analysis / unknown route → null (never a wrong row).
    expect(resolveClickedRow(study, 'nope', 'DGK-1', 'IV-bolus'), null, 'unknown analysisId → null');
    expect(resolveClickedRow(study, 'a1', 'DGK-1', 'PO'), null, 'no matching route → null');
  });

  test('resolveClickedRow: a whole-analysis case resolves via the frozen compoundId back-pointer', async () => {
    // An ungrouped case carries isWholeAnalysis; the grid shows "—" but the hidden
    // compoundId back-pointer holds the FROZEN compute-time compoundId. resolveClickedRow
    // matches on that (v8), not on a reconstructed label (v7 stale-label root cause).
    const study = makeStudy([{
      id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
        compoundRow('Screen A', 'IV-bolus', 'per_kg', {isWholeAnalysis: true}),
      ],
    }]);
    const hit = resolveClickedRow(study, 'a1', 'Screen A', 'IV-bolus');
    expect(hit !== null, true, 'whole-analysis row resolves via the frozen compoundId');
    expect(hit!.row.compoundId, 'Screen A', 'resolves the whole-analysis tier row');
  });

  test('F5 (Tmax): the study grid carries a hidden Tmax param group so the Forms viewer can show it', async () => {
    // The "NCA parameters" Forms viewer binds to the STUDY GRID's frame (M1) and shows its columns.
    // Tmax was card-only in the schema; a hidden Tmax param group is added to the grid frame so the
    // form can display it (user-directed B, context-panel-ct-profile). It carries the PARAM_TAG so
    // the gloss/format resolvers cover it, but is NOT in the visible compact/extended order (hidden
    // in the grid — the SAR scan does not lead with Tmax). Real builder output (rule 16).
    const study = makeStudy([{id: 'a1', label: 'Screen', route: 'IV-bolus', rows: [
      compoundRow('DGK-1', 'IV-bolus', 'per_kg', {tmax: mv(1.25, 20, 3)}),
    ]}]);
    const df = buildStudyGrid(study);
    const tmaxMean = df.col('Tmax Mean');
    expect(tmaxMean !== null, true, 'the grid frame carries a Tmax Mean column (for the Forms viewer)');
    expect(tmaxMean!.get(0), 1.25, 'Tmax Mean carries the cached value');
    expect(tmaxMean!.getTag(PARAM_TAG), 'Tmax', 'Tmax is PARAM_TAG-tagged (gloss/format resolvers cover it)');
    // Hidden in the grid: NOT in the visible compact/extended order (card-only param, not a scan lead).
    expect(studyGridVisibleOrder(df, true).includes('Tmax Mean'), false,
      'Tmax is hidden in the grid (present in the frame for the form, absent from the visible order)');
    // The gloss resolves (map-drift guard — every PARAM_TAG must resolve).
    expect(studyParamGloss('Tmax') !== null, true, 'Tmax resolves to a gloss (no map drift)');
  });

  test('S4c: updateCompoundPanel populates the right context panel from cache (the row-click target)', async () => {
    // The PRODUCTION drill-in path (context-panel-ct-profile): the grid's row-click / selection
    // subscriptions call updateCompoundPanel(study, grid) directly. It sets grok.shell.o to the
    // two-accordion panel: "NCA parameters" (native Forms viewer bound to the grid, M1) +
    // "Concentration–time profile" (the cached curve). Cache-only + synchronous. Driven here through
    // the SAME entry point production uses (a real grid over buildStudyGrid + a current row) — no
    // test-only wrapper, so the test exercises the mechanism production actually runs.
    grok.shell.o = null;
    try {
      const curve: CompoundCurve = {series: [{group: '1', time: [0, 1, 2], geoMean: [100, 50, 25],
        geoLower: [90, 45, 22], geoUpper: [110, 55, 28]}]};
      const study = makeStudy([{
        id: 'a1', label: 'Screen A', route: 'IV-bolus', rows: [
          compoundRow('DGK-7', 'IV-bolus', 'per_kg', {curve, fPct: null}),
        ],
      }]);
      const grid = DG.Viewer.grid(buildStudyGrid(study));
      grid.dataFrame.currentRowIdx = 0; // the clicked row (no selection → single-row path)
      updateCompoundPanel(study, grid);
      const o = grok.shell.o as HTMLElement | null;
      expect(o !== null && o !== undefined, true, 'context panel object set');
      const text = o!.textContent ?? '';
      expect(text.includes('DGK-7'), true, 'panel title shows the clicked compound identity');
      // The two accordion pane headers (BQ-3 title + F4 restructure), both rendered.
      expect(text.includes('NCA parameters'), true, 'panel has the "NCA parameters" accordion pane');
      expect(text.includes('Concentration–time profile'), true,
        'panel has the "Concentration–time profile" accordion pane');
    } finally {
      grok.shell.o = null;
    }
  });

  test('F1: buildCurveFrame(row.curve) == the thumbnail frame (correspondence by construction)', async () => {
    // The core correspondence contract (context-panel-ct-profile AC-F1.1): the panel and the grid
    // thumbnail share ONE builder, so their Time/Conc pairs are byte-identical for the same row.
    const curve: CompoundCurve = {series: [{group: 'A', time: [0, 1, 2, 4], geoMean: [100, 60, 30, 8],
      geoLower: [90, 54, 27, 7], geoUpper: [110, 66, 33, 9]}]};
    const study = makeStudy([{id: 'a1', label: 'S', route: 'IV-bolus',
      rows: [compoundRow('DGK-1', 'IV-bolus', 'per_kg', {curve})]}]);
    const gridFrame = buildStudyGrid(study).col(CURVE_COL)!.get(0) as DG.DataFrame; // thumbnail cell
    const panelFrame = buildCurveFrame(curve); // panel source
    expect(panelFrame.rowCount, gridFrame.rowCount, 'same point count');
    const gt = gridFrame.col('Time')!; const gc = gridFrame.col('Conc')!;
    const pt = panelFrame.col('Time')!; const pc = panelFrame.col('Conc')!;
    for (let i = 0; i < panelFrame.rowCount; i++) {
      expect(pt.get(i), gt.get(i), `Time[${i}] matches the thumbnail`);
      expect(pc.get(i), gc.get(i), `Conc[${i}] matches the thumbnail`);
    }
    // AC-F2.4: the panel frame carries the ×/÷ gSD band columns (whiskers); the thumbnail ignores them.
    expect(panelFrame.col('geoLower') !== null && panelFrame.col('geoUpper') !== null, true,
      'panel frame carries the band columns');
    expect(panelFrame.col('geoLower')!.get(0), 90, 'band value carried onto the frame');
  });

  test('F2: a legacy band-absent curve renders plain (no band cols populated), fresh carries the band', async () => {
    // AC-F2.2/F2.3 distinguisher at the frame level: a legacy curve (no geoLower/geoUpper on the
    // series) → the band columns are all FLOAT_NULL (plain line, the panel adds the "recompute"
    // disclosure); a fresh curve (band present) → real band values.
    const legacy: CompoundCurve = {series: [{group: 'A', time: [0, 1], geoMean: [100, 50]}]}; // no band
    const fresh: CompoundCurve = {series: [{group: 'A', time: [0, 1], geoMean: [100, 50],
      geoLower: [90, 45], geoUpper: [110, 55]}]};
    const legacyFrame = buildCurveFrame(legacy);
    expect(legacyFrame.col('geoLower')!.isNone(0), true, 'legacy curve → band column FLOAT_NULL (plain line)');
    const freshFrame = buildCurveFrame(fresh);
    expect(freshFrame.col('geoLower')!.get(0), 90, 'fresh curve → real band value');
  });

  test('F3: buildOverlayCurveFrame overlays N rows as N identity-labelled series (never truncated)', async () => {
    // AC-F3.1/F3.2: N selected rows → one frame, each series a distinct row-identity category; all
    // rows overlaid (nothing dropped). 14 rows exceeds the legibility threshold but every one rides.
    const entries = Array.from({length: 14}, (_v, i) => ({
      curve: {series: [{group: `g${i}`, time: [0, 1], geoMean: [100 - i, 50 - i],
        geoLower: [90 - i, 45 - i], geoUpper: [110 - i, 55 - i]}]} as CompoundCurve,
      label: `Analysis · C${i} · IV-bolus`,
    }));
    const frame = buildOverlayCurveFrame(entries);
    expect(frame.rowCount, 14 * 2, 'every selected row overlaid (14 curves × 2 points — never truncated)');
    const series = frame.col('Series')!;
    const cats = new Set<string>();
    for (let i = 0; i < frame.rowCount; i++) cats.add(series.getString(i));
    expect(cats.size, 14, 'exactly 14 distinct identity-labelled series');
    expect(cats.has('Analysis · C0 · IV-bolus'), true, 'series keyed by row identity');
  });

  test('F3 (AC-F3.2): a series.length>1 anomaly overlays as ONE labelled series, never multiple', async () => {
    // Post-E3 each row's curve is exactly one series; a >1-series curve is stale migrated data (the
    // v8→v9 migration nulls those, so this cannot occur under current paths). The overlay assembler
    // must NOT render it as multiple series — pushCurve forces the row identity onto every point, so
    // all series collapse to ONE category (and a data-integrity warning is logged).
    const anomaly = {
      curve: {series: [
        {group: 'a', time: [0, 1], geoMean: [100, 50]},
        {group: 'b', time: [0, 1], geoMean: [80, 40]},
      ]} as CompoundCurve,
      label: 'Analysis · C1 · IV-bolus',
    };
    const frame = buildOverlayCurveFrame([anomaly]);
    const s = frame.col('Series')!;
    const cats = new Set<string>();
    for (let i = 0; i < frame.rowCount; i++) cats.add(s.getString(i));
    expect(cats.size, 1, 'a >1-series anomaly renders as ONE row-labelled series, never multiple (AC-F3.2)');
    expect(cats.has('Analysis · C1 · IV-bolus'), true, 'the single category is the row identity');
  });

  test('F7: renderProfileChart declines a NONMEM row honestly (caption, not a wrong chart)', async () => {
    // Event-driven (NONMEM): the decline is KEPT this cycle (BQ-5 conservative default until the
    // AC-F7.1 fixture check confirms cached-curve TAD correctness). renderProfileChart shows the
    // caption pointing to the analysis view — never a silently-wrong chart (rule 14).
    const nonmemEntry = {label: 'NM',
      modes: {nonmem: true, sparse: {on: false, nominalTimeCol: null}}} as unknown as AnalysisEntry;
    const curve: CompoundCurve = {series: [{group: '', time: [0, 1], geoMean: [100, 50]}]};
    const row = compoundRow('NM', 'IV-bolus', 'per_kg', {curve, isWholeAnalysis: true});
    const host = ui.div([]);
    renderProfileChart(host, [{entry: nonmemEntry, row}]);
    const text = host.textContent ?? '';
    expect(text.includes('event-driven (NONMEM)'), true, 'NONMEM row declined with the caption');
    expect(host.querySelector('canvas, svg') === null, true, 'no chart rendered for the NONMEM decline');
  });

  test('F6/F2 buildDisclosures: cross-dose → warning bar w/ n; legacy band → recompute bar; clean → none', async () => {
    // AC-F6.1 + AC-F2.2/F2.3 disclosure-appears (the disclosure bars the profile pane shows).
    const entry = {label: 'S',
      modes: {nonmem: false, sparse: {on: false, nominalTimeCol: null}}} as unknown as AnalysisEntry;
    const freshCurve: CompoundCurve = {series: [{group: 'A', time: [0, 1], geoMean: [100, 50],
      geoLower: [90, 45], geoUpper: [110, 55]}]};
    // Clean row (single-dose, fresh band) → no disclosures (AC-F6.1 negative).
    const clean = compoundRow('C1', 'IV-bolus', 'per_kg', {curve: freshCurve});
    expect(buildDisclosures([{entry, row: clean}]).length, 0, 'clean row → no disclosures (AC-F6.1 negative)');
    // Cross-dose row (CROSS_DOSE_POOL warning) → one non-muted warning bar naming n (AC-F6.1).
    const crossRow = compoundRow('C2', 'IV-bolus', 'per_kg', {curve: freshCurve, n: 16,
      warnings: [{code: 'CROSS_DOSE_POOL', message: 'pools across dose levels'}]});
    const crossBars = buildDisclosures([{entry, row: crossRow}]);
    expect(crossBars.length, 1, 'cross-dose row → one disclosure (AC-F6.1)');
    expect((crossBars[0].textContent ?? '').includes('16'), true, 'the disclosure names n = 16');
    expect((crossBars[0].textContent ?? '').includes('dose'), true, 'the disclosure explains the dose pool');
    // Legacy band-absent curve (no geoLower/geoUpper) → the "recompute" band disclosure (AC-F2.2/F2.3).
    const legacyCurve: CompoundCurve = {series: [{group: 'A', time: [0, 1], geoMean: [100, 50]}]};
    const legacyRow = compoundRow('C3', 'IV-bolus', 'per_kg', {curve: legacyCurve});
    const legacyBars = buildDisclosures([{entry, row: legacyRow}]);
    expect(legacyBars.length, 1, 'legacy band-absent row → one disclosure (AC-F2.2)');
    expect((legacyBars[0].textContent ?? '').toLowerCase().includes('recompute'), true,
      'legacy-band disclosure says "recompute" — distinct from a genuine n=1 present-FLOAT_NULL band (AC-F2.3)');

    // F-N.6 (GROK-20960 AC-N.6.4): the context panel is the SECOND reader of the cached
    // numbers the study grid marks stale. A curve is read as data, not as a quotation of
    // the row that was clicked, so it must disclose an engine-stale cache on its own.
    const staleBars = buildDisclosures([{entry, row: clean}], () => true);
    expect(staleBars.length, 1, 'engine-stale cache → one disclosure on the curve panel');
    const staleText = (staleBars[0].textContent ?? '').toLowerCase();
    expect(staleText.includes('engine version'), true, `disclosure names the cause; got: ${staleText}`);
    expect(staleText.includes('reopen'), true, 'disclosure says how to refresh it');
    // Default predicate (no coordinator) → unchanged behaviour for every other caller.
    expect(buildDisclosures([{entry, row: clean}]).length, 0, 'no predicate → no stale bar');
  });

  test('F1/F6: a standard cached-curve row renders the chart cache-only (no source resolve)', async () => {
    // A standard (non-NONMEM) row with a cached curve → renderProfileChart mounts the scatter
    // synchronously off row.curve (cache-only). No GridStubResolver call is possible (none passed).
    const entry = {label: 'S',
      modes: {nonmem: false, sparse: {on: false, nominalTimeCol: null}}} as unknown as AnalysisEntry;
    const curve: CompoundCurve = {series: [{group: 'A', time: [0, 1, 2], geoMean: [100, 50, 25],
      geoLower: [90, 45, 22], geoUpper: [110, 55, 28]}]};
    const row = compoundRow('DGK-1', 'IV-bolus', 'per_kg', {curve});
    const host = ui.div([]);
    renderProfileChart(host, [{entry, row}]);
    // A chart viewer node mounted (cache-only render succeeded); no decline caption.
    expect((host.textContent ?? '').includes('event-driven'), false, 'no NONMEM decline for a standard row');
    expect(host.childElementCount > 0, true, 'a chart node mounted from the cached curve');
  });
});

category('NCA Studies: study grid — Exclude λz disclosure (F-S5.2)', () => {
  test('n = a of b on a partly-excluded row; "not reportable" when every contributor was excluded', async () => {
    const study = makeStudy([{id: 'a1', label: 'Screen A', route: 'PO', rows: [
      // Row 0: every t½ contributor excluded (n 0, 3 excluded); AUC∞ partly (2 of 3).
      compoundRow('DGK-1', 'PO', 'per_kg', {
        thalf: {mean: NaN, cvPct: NaN, n: 0, nExcluded: 3},
        aucInf: {mean: 100, cvPct: 8, n: 2, nExcluded: 1},
      }),
      compoundRow('DGK-2', 'PO', 'per_kg'), // nothing excluded (keeps the t½ group non-empty)
    ]}]);
    const df = buildStudyGrid(study);
    const all = studyCellExclusion(df, 't½ Mean', 0);
    expect(all?.notReportable, true, 'every contributor excluded → not reportable');
    expect(all?.total, 3);
    const part = studyCellExclusion(df, 'AUC∞ n', 0);
    expect(part?.n, 2, 'n');
    expect(part?.total, 3, 'n = 2 of 3');
    expect(part?.notReportable, false, 'a partly-excluded cell keeps its value');
    expect(studyCellExclusion(df, 'AUC∞ Mean', 1), null, 'a row with nothing excluded discloses nothing');
    expect(studyCellExclusion(df, 'Cmax Mean', 0), null, 'a parameter with no exclusion discloses nothing');
  });
});
