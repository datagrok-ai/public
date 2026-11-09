/**
 * NCA Studies — passthrough (annotation) columns (study-grid-passthrough-columns, v10).
 *
 * @see src/adapter/passthrough-gate.ts (passthroughCandidate gate — F1)
 * @see src/views/nca-analysis-open.ts (buildPassthroughByTier — F3)
 * @see src/adapter/study-grid.ts (passthroughCellValue + emission — F4)
 *
 * Three surfaces:
 *  1. The CAP-FREE candidate gate (F1) — a pure predicate; the RG-D cap-free proof lives here
 *     (a 25-distinct-value column IS offerable; per-observation-varying / role / reserved excluded).
 *  2. `buildPassthroughByTier` (F3) against REAL generated pipeline output (rule 16) — the 05
 *     multicompound fixture with an INJECTED per-compound SMILES + a >20-distinct per-subject
 *     column: cap-free carriage (AC-F3a), pooled multiplicity (AC-F3b), isStructure from the
 *     recipe (AC-F3c), the SEPARATE-from-E4 gate (a >20 column E4 drops but passthrough keeps).
 *  3. `buildStudyGrid` render (F4) — bare-when-unique, counted-when-many-text, `— (N)`-when-many-
 *     structure, E2 typographic merge, E1 sparse, PASSTHROUGH_TAG + the overview zone order.
 *
 * Runtime: `grok test --host local` (fixtures load via `grok.dapi.files`).
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {passthroughCandidate, isConstantWithinSubject, recipeRoleNames} from '../adapter';
import {buildSummaryCache, recipeOf} from '../adapter';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {
  buildStudyGrid, PASSTHROUGH_TAG, PASSTHROUGH_STRUCT, PASSTHROUGH_TEXT, NO_GROUP,
} from '../adapter';
import type {
  AnalysisRecipe, AnalysisEntry, NcaStudy, CompoundTierRow, MeanCv, PassthroughEntry, KeyValueCount,
} from '../adapter';
import type {PerProfileResult} from '../controller';
import {DEFAULT_RULES} from '../adapter/defaults';
import {
  loadFixture, multicompoundContract, multicompoundScalars, runPipelineOrThrow,
} from './test-utils';

const ISO = '2026-07-14T00:00:00.000Z';

/** A recipe carrying the group-by set + a passthrough selection. Cast through unknown (the
 *  fields buildSummaryCache/gate read); mirrors summary-cache-tests' recipeWith convenience. */
function recipeWith(
  groupByCols: string[],
  passthroughCols: {name: string; isStructure: boolean}[] = [],
  columns: Record<string, unknown> = {doseSource: 'column', doseScalar: 0},
): AnalysisRecipe {
  return {columns, units: {dose: 'mg/kg'}, modes: {}, ncaRules: {},
    groupByCols, passthroughCols} as unknown as AnalysisRecipe;
}

// ─── 1. Candidate gate (F1) — the cap-free predicate ─────────────────────────

category('NCA passthrough: candidate gate (cap-free)', () => {
  /** A source frame: 25 subjects × 3 timepoints. `SMILES` = one distinct per subject (25 distinct,
   *  constant within subject → the cap-free proof); `PerObs` varies within a subject; `Sex` is a
   *  role; `CompoundID` is a group-by; `CMAX` collides with a reserved results-frame name. */
  function gateFrame(): DG.DataFrame {
    const subj: string[] = []; const time: number[] = []; const conc: number[] = [];
    const smiles: string[] = []; const perObs: string[] = []; const sex: string[] = [];
    const compound: string[] = []; const cmax: number[] = [];
    for (let s = 0; s < 25; s++) {
      const sid = `S${String(s).padStart(2, '0')}`;
      for (let t = 0; t < 3; t++) {
        subj.push(sid); time.push(t); conc.push(100 - t * 10 - s);
        smiles.push(`C${'C'.repeat(s)}O`); // 25 distinct SMILES, constant within a subject
        perObs.push(`${sid}-${t}`); // varies within a subject → NOT constant
        sex.push(s % 2 === 0 ? 'F' : 'M');
        compound.push(sid); cmax.push(100 - s);
      }
    }
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings('SubjectID', subj),
      DG.Column.fromFloat64Array('Time', Float64Array.from(time)),
      DG.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
      DG.Column.fromStrings('SMILES', smiles),
      DG.Column.fromStrings('PerObs', perObs),
      DG.Column.fromStrings('Sex', sex),
      DG.Column.fromStrings('CompoundID', compound),
      DG.Column.fromFloat64Array('CMAX', Float64Array.from(cmax)),
    ]);
  }

  const recipe = recipeWith(['CompoundID'], [], {
    subject: 'SubjectID', time: 'Time', concentration: 'Conc', sexCol: 'Sex',
    doseSource: 'column', doseScalar: 0,
  });

  test('AC-F1b: a 25-distinct constant-within-subject column IS offerable (NO cardinality cap)', async () => {
    const df = gateFrame();
    const subjectCol = df.col('SubjectID');
    expect(passthroughCandidate(df.col('SMILES')!, recipe, subjectCol), true,
      '25-distinct SMILES offerable — the cap-free proof (RG-D)');
  });

  test('AC-F1b: a per-observation column that varies within a subject is NOT offered', async () => {
    const df = gateFrame();
    expect(passthroughCandidate(df.col('PerObs')!, recipe, df.col('SubjectID')), false,
      'per-timepoint-varying column rejected (not a valid annotation)');
  });

  test('AC-F1e: a mapped role (sexCol) and a group-by (CompoundID) are excluded (LBC-1)', async () => {
    const df = gateFrame();
    const subjectCol = df.col('SubjectID');
    expect(passthroughCandidate(df.col('Sex')!, recipe, subjectCol), false,
      'mapped sexCol excluded — no double-render vs its E4 keyInventory disclosure');
    expect(passthroughCandidate(df.col('CompoundID')!, recipe, subjectCol), false,
      'the grouping key is excluded (already a key column)');
    // The role exclusion set is built from RecipeColumns directly (not inputColumnNames).
    expect(recipeRoleNames(recipe).has('Sex'), true, 'Sex is in the RecipeColumns role set (LBC-1)');
  });

  test('AC-F1f: a source column colliding with a reserved results-frame name is EXCLUDED (LBC-3/R2-4)', async () => {
    const df = gateFrame();
    // `CMAX` collides with the results-writer / DEFAULT_AGG_PARAMS schema → appendGroupingColumns
    // would silently skip it, so buildPassthroughByTier would read the RESERVED column. Excluded.
    expect(passthroughCandidate(df.col('CMAX')!, recipe, df.col('SubjectID')), false,
      'reserved-name collision excluded (never rename-on-carry — would break the E2 merge)');
  });

  test('isConstantWithinSubject: null subject → global-constancy fallback', async () => {
    const constant = DG.Column.fromStrings('X', ['a', 'a', 'a']);
    const varying = DG.Column.fromStrings('Y', ['a', 'b', 'a']);
    expect(isConstantWithinSubject(constant, null), true, 'globally-constant column ok when no subject');
    expect(isConstantWithinSubject(varying, null), false, 'globally-varying column rejected when no subject');
  });
});

// ─── 2. buildPassthroughByTier (F3) — real generated output ──────────────────

category('NCA passthrough: buildPassthroughByTier (real fixture)', () => {
  /** Inject a per-compound SMILES (20 distinct — the SAR structure column) and a per-SUBJECT
   *  unique barcode (60 distinct, > the E4 cap) into the 05 fixture, run the real pipeline, and
   *  carry both onto the results frame. Returns the post-append results + perProfile. */
  async function injectedRun(): Promise<{results: DG.DataFrame; perProfile: ReadonlyArray<PerProfileResult>}> {
    const df = await loadFixture('05_multicompound_batch_screen.csv');
    const compoundCol = df.col('CompoundID')!;
    const subjectCol = df.col('SubjectID')!;
    // Deterministic per-compound SMILES (constant within a subject, since a subject is one compound).
    const smilesOf = (c: string): string => `C${'C'.repeat((c.charCodeAt(c.length - 1)) % 12)}O`;
    const smiles: string[] = []; const barcode: string[] = [];
    for (let i = 0; i < df.rowCount; i++) {
      smiles.push(smilesOf(compoundCol.getString(i)));
      barcode.push(`BC-${subjectCol.getString(i)}`); // unique per subject → 60 distinct (> cap 20)
    }
    df.columns.add(DG.Column.fromStrings('Canonical_SMILES', smiles));
    df.columns.add(DG.Column.fromStrings('Batch_Barcode', barcode));

    const contract = multicompoundContract(df);
    const {results, outcome} = runPipelineOrThrow(df, contract, multicompoundScalars(), DEFAULT_RULES);
    const added = appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);
    // RG-H verified against REAL output: the passthrough source columns ride onto the results frame.
    expect(added.includes('Canonical_SMILES'), true, 'AC-F3d: SMILES carried onto resultsDf');
    expect(added.includes('Batch_Barcode'), true, 'AC-F3d: barcode carried onto resultsDf');
    return {results, perProfile: outcome.perProfile};
  }

  const passSel = [
    {name: 'Canonical_SMILES', isStructure: true},
    {name: 'Batch_Barcode', isStructure: false},
  ];

  test('AC-F3a/c: grouped by compound keeps each row SMILES (1 value, count N); isStructure from recipe', async () => {
    const {results, perProfile} = await injectedRun();
    const cache = buildSummaryCache(results, perProfile, recipeWith(['CompoundID'], passSel), 'S');
    expect(cache.compoundTier.length, 20, '20 compound rows');
    for (const r of cache.compoundTier) {
      const smi = r.passthrough.find((p) => p.key === 'Canonical_SMILES');
      expect(smi !== undefined, true, `compound ${r.compoundId} carries a SMILES passthrough entry`);
      expect(smi!.values.length, 1, `compound ${r.compoundId}: one distinct SMILES (compound-homogeneous)`);
      expect(smi!.values[0].count, 3, `compound ${r.compoundId}: SMILES subject count = 3`);
      // AC-F3c: isStructure comes from the recipe spec, NOT re-derived from resultsDf (which
      // stripped semType). A recompute of the SAME analysis yields the SAME designation.
      expect(smi!.isStructure, true, `compound ${r.compoundId}: SMILES isStructure from recipe`);
      const bc = r.passthrough.find((p) => p.key === 'Batch_Barcode');
      expect(bc !== undefined && bc.isStructure === false, true, `${r.compoundId}: barcode isStructure false`);
    }
  });

  test('BUG-11 stress: a duplicated passthrough spec does NOT double the subject counts', async () => {
    // The dangerous half of the BUG-11 pattern. `accumulateByTier` keys its counts by column
    // NAME, so two same-named candidates would increment the SAME counter — turning a real
    // "3 subjects" into "6". Unlike a doubled label, a doubled COUNT reads as real data
    // (thesis #4), so the count must be correct by construction even if a caller ever hands
    // over a duplicate. Against the real pipeline, not a mirror (rule 16): the same fixture the
    // clean test above pins at count 3.
    const {results, perProfile} = await injectedRun();
    const dupSel = [...passSel, {name: 'Canonical_SMILES', isStructure: true}];
    const cache = buildSummaryCache(results, perProfile, recipeWith(['CompoundID'], dupSel), 'S');
    for (const r of cache.compoundTier) {
      const entries = r.passthrough.filter((p) => p.key === 'Canonical_SMILES');
      expect(entries.length, 1, `compound ${r.compoundId}: ONE entry, not one per repeated spec`);
      expect(entries[0].values[0].count, 3,
        `compound ${r.compoundId}: subject count stays 3 — never doubled to 6`);
    }
  });

  test('AC-F3a/b: whole-analysis pooled row discloses ALL distinct values + counts (cap-free, >20)', async () => {
    const {results, perProfile} = await injectedRun();
    // No compound grouping → one pooled row per (route × dose-norm). The 05 screen is single-route
    // (IV) per-kg → ONE pooled row spanning all 20 compounds / 60 subjects.
    const cache = buildSummaryCache(results, perProfile, recipeWith([], passSel), 'S');
    const pooled = cache.compoundTier.find((r) => r.isWholeAnalysis) ?? cache.compoundTier[0];
    const smi = pooled.passthrough.find((p) => p.key === 'Canonical_SMILES');
    const bc = pooled.passthrough.find((p) => p.key === 'Batch_Barcode');
    // AC-F3b: every distinct value disclosed with a count — never one silently picked (thesis #4).
    expect(smi !== undefined && smi.values.length > 1, true, 'pooled row discloses MANY distinct SMILES');
    // AC-F3a: the barcode has 60 distinct values (> the E4 cap of 20) yet is KEPT — the cap-free proof.
    expect(bc !== undefined && bc.values.length > 20, true,
      'a >20-distinct passthrough column is KEPT (cap-free — E4 would drop it)');
    const totalSubjects = (bc!.values as readonly KeyValueCount[]).reduce((a, v) => a + v.count, 0);
    expect(totalSubjects, 60, 'barcode counts sum to the 60 pooled subjects');
  });

  test('AC-F3e (RG-D): a >20-distinct column E4 drops is KEPT by passthrough (separate gates)', async () => {
    const {results, perProfile} = await injectedRun();
    const cache = buildSummaryCache(results, perProfile, recipeWith([], passSel), 'S');
    const pooled = cache.compoundTier.find((r) => r.isWholeAnalysis) ?? cache.compoundTier[0];
    // E4 keyInventory applies the cardinality cap: the 60-distinct barcode is NOT a keyInventory
    // candidate. The passthrough channel (cap-free) KEEPS it — proving the two gates are separate.
    expect(pooled.keyInventory.some((e) => e.key === 'Batch_Barcode'), false,
      'barcode DROPPED from keyInventory (E4 cardinality cap)');
    expect(pooled.passthrough.some((p) => p.key === 'Batch_Barcode'), true,
      'barcode KEPT in passthrough (cap-free) — RG-D separate gates');
  });

  test('empty passthroughCols → no passthrough entries (pre-v10 parity)', async () => {
    const {results, perProfile} = await injectedRun();
    const cache = buildSummaryCache(results, perProfile, recipeWith(['CompoundID'], []), 'S');
    for (const r of cache.compoundTier)
      expect(r.passthrough.length, 0, `compound ${r.compoundId}: no passthrough when none selected`);
  });
});

// ─── 3. buildStudyGrid render (F4) — display honesty ─────────────────────────

category('NCA passthrough: study-grid render', () => {
  function mv(mean: number): MeanCv {return {mean, cvPct: 5, n: 3};}
  /** A minimal grouped tier row with the given passthrough entries. */
  function row(compoundId: string, passthrough: PassthroughEntry[]): CompoundTierRow {
    return {
      compoundId, isWholeAnalysis: false,
      groupValues: [{key: 'CompoundID', value: compoundId}], keyInventory: [], passthrough,
      n: 3, route: 'IV-bolus', doseNorm: 'per_kg',
      cmax: mv(10), tmax: mv(1), aucInf: mv(100), thalf: mv(2), cl: mv(1), vz: mv(3),
      mrt: mv(4), vss: mv(6), fPct: null, pctExtrap: mv(15), curve: null, warnings: [],
    };
  }
  function analysis(
    id: string, rows: CompoundTierRow[], passthroughCols: {name: string; isStructure: boolean}[],
  ): AnalysisEntry {
    return {
      id, label: id, source: {kind: 'tableId', ref: id},
      columns: {} as never, modes: {} as never, ncaRules: {} as never, units: {} as never,
      groupByCols: ['CompoundID'], passthroughCols,
      sourceSignature: {kind: 'content-hash', rowCount: 3, mtime: null, hash: 'h', computedAt: ISO},
      summaryCache: {
        analysisTier: {nCompounds: rows.length, nSubjects: 3 * rows.length, route: 'IV-bolus',
          dose: '(per-row)', status: 'ok', flagSummary: 'No flags'},
        compoundTier: rows,
      },
    };
  }
  function study(analyses: AnalysisEntry[]): NcaStudy {
    return {id: 's', label: 'S', kind: 'study', createdAt: ISO, updatedAt: ISO,
      access: {canEdit: true, canDelete: true, canShare: true}, analyses};
  }
  const pe = (key: string, isStructure: boolean, ...vals: [string, number][]): PassthroughEntry =>
    ({key, isStructure, values: vals.map(([value, count]) => ({value, count}))});

  test('AC-F4b/c: bare-when-unique (lone structure → molecule value); counted-join when many-text', async () => {
    // Structure column named `Canonical_SMILES` (NOT the literal `SMILES`) — so a PASS proves the
    // structure tag rides on the persisted isStructure designation, NOT the retired name-match guard.
    const a = analysis('A', [
      row('DGK-1', [pe('Canonical_SMILES', true, ['CCO', 3]), pe('Series', false, ['S1', 3])]),
      row('DGK-2', [pe('Canonical_SMILES', true, ['c1ccccc1', 3]), pe('Series', false, ['S1', 2], ['S2', 1])]),
    ], [{name: 'Canonical_SMILES', isStructure: true}, {name: 'Series', isStructure: false}]);
    const grid = buildStudyGrid(study([a]));
    const smi = grid.col('Canonical_SMILES')!; const series = grid.col('Series')!;
    expect(smi.getString(0), 'CCO', 'lone structure rendered BARE (so it renders as a molecule)');
    expect(smi.getString(1), 'c1ccccc1', 'lone structure bare on row 2');
    expect(series.getString(0), 'S1', 'single-value text passthrough rendered bare');
    expect(series.getString(1), 'S1 (2); S2 (1)', 'multi-value text passthrough → counted join');
    // PASSTHROUGH_TAG carries the structure designation for the overview's semType logic — the tag
    // is set from isStructure, name-independently (a `Canonical_SMILES` col still tags structure).
    expect(smi.getTag(PASSTHROUGH_TAG), PASSTHROUGH_STRUCT, 'structure tag by designation, not name');
    expect(series.getTag(PASSTHROUGH_TAG), PASSTHROUGH_TEXT, 'Series tagged text');
  });

  test('AC-F5b (recipeOf half): recipeOf re-injects entry.passthroughCols for replay-recompute', async () => {
    // The headless-runnable half of the replay science gate: recipeOf(entry) must carry the
    // persisted passthrough selection (incl. isStructure) so a reopen recompute rebuilds the SAR
    // columns from the authoritative dialog-time pick (LBC-2). The full mountAnalysisIntoView replay
    // is exercised live via grok test; this pins the recipeOf injection against a silent drop.
    const a = analysis('A', [row('DGK-1', [])],
      [{name: 'Canonical_SMILES', isStructure: true}, {name: 'Series', isStructure: false}]);
    const recipe = recipeOf(a);
    expect(recipe.passthroughCols.length, 2, 'recipeOf carries the persisted passthrough specs');
    expect(recipe.passthroughCols[0].name, 'Canonical_SMILES', 'spec name preserved');
    expect(recipe.passthroughCols[0].isStructure, true, 'structure designation preserved (no re-derive)');
    expect(recipe.passthroughCols[1].isStructure, false, 'text designation preserved');
  });

  test('recipeOf carries the persisted adj-R² acceptance threshold (R_acc §7 — no silent drop on reopen)', async () => {
    // Regression pin for the reopen bug (user QA 2026-08-25): recipeOf (open-side) is a SEPARATE
    // field enumeration from getCurrentRecipe, so a recipe field must be listed here too or it is
    // silently dropped — which reset R_acc to the default on reopen. (The save-side entry build in
    // openAddToStudyDialog is the sibling enumeration, fixed in the same change; its full path is
    // exercised live via grok test, same split as the passthrough gate above.)
    const base = analysis('A', [row('DGK-1', [])], []);
    expect(recipeOf({...base, acceptanceRSquared: 0.86}).acceptanceRSquared, 0.86,
      'recipeOf carries the SAVED R_acc back to the pane, not the default');
    // Boolean-equality form: the DG-test `expect(x, undefined)` overload collapses to a
    // truthy check, so assert the `=== undefined` predicate instead.
    expect(recipeOf(base).acceptanceRSquared === undefined, true,
      'a pre-v13 entry (no field) carries undefined → the pane falls back to its default');
  });

  test('AC-F4d: a many-structure pooled cell → bare "—" (count on hover), never a single depiction', async () => {
    // User-approved 2026-07-14: the cell shows a bare '—'; the distinct count is disclosed on the
    // cell TOOLTIP (configureStudyGrid, a live/visual concern). What buildStudyGrid must guarantee
    // is that it NEVER emits one of the many structures as the cell value (thesis #4).
    const a = analysis('A', [
      row('POOL', [pe('Canonical_SMILES', true, ['CCO', 10], ['c1ccccc1', 15])]),
    ], [{name: 'Canonical_SMILES', isStructure: true}]);
    const grid = buildStudyGrid(study([a]));
    expect(grid.col('Canonical_SMILES')!.getString(0), NO_GROUP,
      'many-structure cell = bare "—" (never one of the 2 structures; count on tooltip — thesis #4)');
  });

  test('AC-F4a: E2 typographic merge (SMILES/smiles → 1 column) + E1 sparse (— where absent)', async () => {
    // Analysis A carries `SMILES`; Analysis B carries `smiles` (case variant) + a `Series` column
    // A lacks. The canonical-id merge folds SMILES/smiles into ONE column; Series is sparse on A.
    const a = analysis('A', [row('DGK-1', [pe('SMILES', true, ['CCO', 3])])],
      [{name: 'SMILES', isStructure: true}]);
    const b = analysis('B', [row('DGK-2', [pe('smiles', true, ['c1ccccc1', 3]), pe('Series', false, ['S9', 3])])],
      [{name: 'smiles', isStructure: true}, {name: 'Series', isStructure: false}]);
    const grid = buildStudyGrid(study([a, b]));
    const names = grid.columns.names();
    const smilesCols = names.filter((n) => n.toLowerCase() === 'smiles');
    expect(smilesCols.length, 1, 'SMILES/smiles merged into ONE column (E2 typographic merge)');
    const smi = grid.col(smilesCols[0])!;
    expect(smi.getString(0), 'CCO', 'row from A carries its SMILES');
    expect(smi.getString(1), 'c1ccccc1', 'row from B merges into the same column');
    // E1 sparse: Series exists (from B) but A's row has no value → NO_GROUP.
    expect(grid.col('Series')!.getString(0), NO_GROUP, 'a passthrough column absent from a row → —');
    expect(grid.col('Series')!.getString(1), 'S9', 'present where the analysis carries it');
  });

  test('AC-F4e: a study with NO passthrough selection renders zero passthrough columns (pre-v10 parity)', async () => {
    const a = analysis('A', [row('DGK-1', [])], []);
    const grid = buildStudyGrid(study([a]));
    const anyPassthrough = grid.columns.toList().some((c) => (c.getTag(PASSTHROUGH_TAG) ?? '') !== '');
    expect(anyPassthrough, false, 'no passthrough-tagged columns when nothing was selected');
  });
});
