/**
 * Multi-column "Group by" (composite grouping) — the v8 widen `groupByCol` →
 * `groupByCols` and its load-bearing seams.
 *
 * Coverage (pure-function + real-fixture, NOT mirror-pattern):
 *  - F3b `compositeGroupKey`: U+0001-joined, byte-identical, null → '' (no collision).
 *  - F3b/SF-2 `resolveAggregateGroupKeys`: composite-tuple mixed-basis bucketing —
 *    fires per CELL, not per dimension (AC-F3.6 in-unit form).
 *  - F2 `countDoseClusters`: log-space clusterer suppresses theoph (real fixture,
 *    1 nominal level) but fires on a jittered ≥2-arm escalation (DATA-GAP-31).
 *  - F3c `buildStudyGrid`: one named column per key (sparse `—`) for BOTH single-key and
 *    composite grouping — the generic `Group` column is retired (study-grid-grouping-columns).
 *  - F5 composite ordering: natural sort on the joined key is numeric-within-join.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {
  compositeGroupKey, KEY_SEP, resolveAggregateGroupKeys, buildStudyGrid, buildMeanSdProfile,
  aggregateNcaResults, NO_GROUP, countDoseClusters,
} from '../adapter';
import type {NcaStudy, CompoundTierRow, MeanCv, GroupKeyValue} from '../adapter';
import {
  fragilityNoticeText, dedupeGroupBySelection, FRAGILE_N, FRAGILE_LIST_MAX,
} from '../views/nca-mapping-pane';
import type {FragileCell} from '../views/nca-mapping-pane';
import {loadFixture} from './test-utils';

function mv(mean: number, cvPct: number, n: number): MeanCv {
  return {mean, cvPct, n};
}

/** A tier row with the given per-key group values (v8). `compoundId` is the composite
 *  key (compositeGroupKey over the values), matching what `buildSummaryCache` emits. */
function tierRow(groupValues: GroupKeyValue[], route: string, over: Partial<CompoundTierRow> = {}): CompoundTierRow {
  const compoundId = groupValues.length === 0 ? 'WHOLE' :
    groupValues.map((g) => g.value).join(KEY_SEP);
  return {
    compoundId, isWholeAnalysis: groupValues.length === 0, groupValues,
    keyInventory: [], // v9: per-tier-row candidate inventory (empty for these projection tests)
    passthrough: [], // v10: per-tier-row passthrough inventory (empty for these projection tests)
    n: 3, route, doseNorm: 'per_kg',
    cmax: mv(10, 5, 3), tmax: mv(1, 20, 3), aucInf: mv(100, 8, 3), thalf: mv(2, 12, 3),
    cl: mv(1.5, 9, 3), vz: mv(3, 11, 3), mrt: mv(4, 10, 3), vss: mv(6, 9, 3),
    fPct: null, pctExtrap: mv(15, 4, 3), curve: null, warnings: [], ...over,
  };
}

function studyWith(rows: CompoundTierRow[], groupByCols: string[]): NcaStudy {
  return {
    id: 's1', label: 'Study', kind: 'study',
    createdAt: '2026-07-09T00:00:00.000Z', updatedAt: '2026-07-09T00:00:00.000Z',
    access: {canEdit: true, canDelete: true, canShare: true},
    analyses: [{
      id: 'a1', label: 'Screen A', source: {kind: 'tableId', ref: 't'},
      columns: {}, modes: {}, ncaRules: {}, units: {}, groupByCols,
      sourceSignature: {kind: 'content-hash', rowCount: 6, mtime: null, hash: 'h',
        computedAt: '2026-07-09T00:00:00.000Z'},
      summaryCache: {
        analysisTier: {nCompounds: rows.length, nSubjects: 3, route: 'IV-bolus', dose: '(per-row)',
          status: 'ok', flagSummary: 'No flags'},
        compoundTier: rows,
      },
    }],
  } as unknown as NcaStudy;
}

category('multi-column-group-by: compositeGroupKey (F3b)', () => {
  test('U+0001-joined, byte-identical, null → empty segment (no collision)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('A', ['IV', 'IVbolus', 'M']),
      DG.Column.fromStrings('B', ['bolus', '', '10']),
    ]);
    const a = df.col('A')!; const b = df.col('B')!;
    // ['IV','bolus'] and ['IVbolus',''] must NOT collide (the .join('') hazard).
    expect(compositeGroupKey([a, b], 0) === compositeGroupKey([a, b], 1), false,
      '[\'IV\',\'bolus\'] and [\'IVbolus\',\'\'] are distinct keys (U+0001 separator)');
    expect(compositeGroupKey([a, b], 0), `IV${KEY_SEP}bolus`, 'joins with U+0001');
    // A single-column key equals that column's value (no separator).
    expect(compositeGroupKey([a], 2), 'M', 'single-column key = the value');
    // A null cell contributes an empty segment (not the string "null").
    const withNull = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'S', ['M', null] as unknown as string[]),
    ]);
    expect(compositeGroupKey([withNull.col('S')!], 1), '', 'null cell → empty segment');
  });
});

category('multi-column-group-by: resolveAggregateGroupKeys composite bucketing (SF-2)', () => {
  test('mixed-basis fires per composite CELL, not per single dimension', async () => {
    // SEX=M spans both CL (its DOSE=10 rows, IV) and CL/F (its DOSE=30 rows, PO), so
    // grouping by SEX ALONE mixes basis; but every (SEX,DOSE) CELL is single-basis, so
    // the composite (SEX,DOSE) grouping does NOT (AC-F3.6: no per-key over-trigger).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('SEX', ['M', 'M', 'M', 'M', 'F', 'F']),
      DG.Column.fromStrings('DOSE', ['10', '10', '30', '30', '10', '10']),
      DG.Column.fromStrings('CL_BASIS', ['CL', 'CL', 'CL/F', 'CL/F', 'CL', 'CL']),
    ]);
    const rSex = resolveAggregateGroupKeys(df, ['SEX']);
    expect(rSex.mixedBasis, true, 'SEX alone mixes CL basis (spans CL and CL/F)');
    // group-by-cl-basis-explicit-control F1: SEX-alone is mixed, but the resolver NO LONGER
    // appends CL_BASIS — it honours the user's key and reports the mixed cell (the caller
    // withholds the family instead of auto-stratifying).
    expect(rSex.keys.join(','), 'SEX', 'F1: no CL_BASIS auto-append even when SEX mixes basis');
    expect(rSex.mixedBasisKeys.has('M'), true, 'the mixed M cell is reported for the blank pass');
    expect(resolveAggregateGroupKeys(df, ['SEX', 'DOSE']).mixedBasis, false,
      'each (SEX,DOSE) cell is single-basis → composite grouping does not over-trigger');
    // The resolved keys carry the whole tuple; CL_BASIS is NOT appended for the composite.
    const {keys} = resolveAggregateGroupKeys(df, ['SEX', 'DOSE']);
    expect(keys.join(','), 'SEX,DOSE', 'composite keys = the tuple, no spurious CL_BASIS split');
  });
});

category('multi-column-group-by: countDoseClusters (F2 / DATA-GAP-31)', () => {
  test('theoph (real fixture, 1 nominal level) → 1 cluster; suppresses the cross-dose notice', async () => {
    const df = await loadFixture('theoph.csv');
    // theoph Dose = 10 distinct per-kg values in [3.1 … 5.86] (one nominal level, a
    // covariate). Verified against the REAL fixture, not a mirror (rule 16).
    expect(countDoseClusters(df.col('Dose')!), 1,
      'theoph per-kg dose distribution is one cluster (notice must NOT fire)');
  });

  test('a jittered ≥2-arm escalation → ≥2 clusters (notice fires)', async () => {
    // Linear 3-arm around 10 / 30 / 100 mg/kg with ±10% per-kg jitter.
    const lin3 = DG.Column.fromFloat64Array('Dose',
      Float64Array.from([9, 9.7, 10.4, 11, 28, 29.5, 31, 32, 95, 98, 102, 105]));
    expect(countDoseClusters(lin3), 3, 'linear 3-arm → 3 clusters');
    // Log-spaced (geometric) arms — the modal tox convention; must not over-split.
    const log = DG.Column.fromFloat64Array('Dose',
      Float64Array.from([0.9, 1, 1.1, 2.7, 3, 3.3, 9, 10, 11, 27, 30, 33]));
    expect(countDoseClusters(log) >= 2, true, 'log-spaced arms → ≥2 clusters (no within-arm over-split)');
  });
});

category('multi-column-group-by: buildStudyGrid per-key columns (F3c)', () => {
  test('≥2 keys → one column per key (sparse —); no generic Group column', async () => {
    const study = studyWith([
      tierRow([{key: 'Sex', value: 'M'}, {key: 'Dose', value: '10'}], 'IV-bolus'),
      tierRow([{key: 'Sex', value: 'F'}, {key: 'Dose', value: '30'}], 'IV-bolus'),
    ], ['SEX', 'DOSE']);
    const df = buildStudyGrid(study);
    expect(df.col('Sex') !== null, true, 'per-key Sex column present');
    expect(df.col('Dose') !== null, true, 'per-key Dose column present');
    expect(df.col('Group'), null, 'no generic Group column (retired, study-grid-grouping-columns)');
    // A per-key column of all-numeric values auto-types to numeric (so Dose sorts
    // numerically — a Variant-2 benefit); compare via String() to be type-robust.
    expect(df.col('Sex')!.get(0), 'M', 'Sex cell');
    expect(String(df.col('Dose')!.get(1)), '30', 'Dose cell');
    // Per-dimension filter/sort is possible because each key is its own column.
  });

  test('union of keys across analyses is sparse — a row missing a key shows —', async () => {
    const study = studyWith([
      tierRow([{key: 'Sex', value: 'M'}, {key: 'Dose', value: '10'}], 'IV-bolus'),
      // A second row that only carries Dose (no Sex) — the Sex cell must be `—`.
      tierRow([{key: 'Dose', value: '30'}], 'IV-bolus'),
    ], ['SEX', 'DOSE']);
    const df = buildStudyGrid(study);
    expect(df.col('Sex')!.get(1), NO_GROUP, 'a row missing the Sex key renders — (sparse union)');
    expect(String(df.col('Dose')!.get(1)), '30', 'its Dose key still shows');
  });

  test('a single key gets its OWN named column (Group retired, study-grid-grouping-columns)', async () => {
    const study = studyWith([
      tierRow([{key: 'Compound', value: 'DGK-1'}], 'IV-bolus'),
      tierRow([{key: 'Compound', value: 'DGK-2'}], 'IV-bolus'),
    ], ['CompoundID']);
    const df = buildStudyGrid(study);
    expect(df.col('Group'), null, 'no generic Group column even for single-key grouping (retired)');
    expect(df.col('Compound') !== null, true, 'single-key grouping gets a named Compound column');
    expect(df.col('Compound')!.get(0), 'DGK-1', 'the named column shows the single key value');
  });
});

category('multi-column-group-by: fragility banner text (F4 / DATA-GAP-30)', () => {
  const cell = (group: string, time: number, n: number): FragileCell => ({group, time, n});

  test('AC-F4.1: names each fragile cell, worst (smallest N) first', async () => {
    // Composite groups use U+0001 internally; the banner humanizes to ` · `.
    const cells = [
      cell(`M${KEY_SEP}10`, 24, 3), cell(`M${KEY_SEP}100`, 24, 1), cell(`F${KEY_SEP}30`, 12, 2)];
    const txt = fragilityNoticeText(cells, /* nominalUsed */ true);
    expect(txt.includes('N=1 for M · 100 at 24 h'), true, 'names the N=1 cell (humanized composite key)');
    expect(txt.includes('N=2 for F · 30 at 12 h'), true, 'names the N=2 cell');
    expect(txt.includes('N=3 for M · 10 at 24 h'), true, 'names the N=3 cell (N<FRAGILE_N=4)');
    // Worst-first: the N=1 cell precedes the N=2 and N=3 cells in the text.
    expect(txt.indexOf('N=1') < txt.indexOf('N=2') && txt.indexOf('N=2') < txt.indexOf('N=3'), true,
      'cells ordered worst (smallest N) first');
    expect(txt.includes('cells have N <'), false, 'enumerates (does not collapse) at 3 ≤ FRAGILE_LIST_MAX cells');
  });

  test('AC-F4.3: >FRAGILE_LIST_MAX fragile cells → collapse to count + worst (no enumeration)', async () => {
    // FRAGILE_LIST_MAX + 1 fragile cells → collapse. One is the worst (N=1).
    const cells: FragileCell[] = [];
    for (let i = 0; i < FRAGILE_LIST_MAX + 1; i++) cells.push(cell(`G${KEY_SEP}${i}`, i, i === 0 ? 1 : 2));
    const txt = fragilityNoticeText(cells, true);
    expect(txt.includes(`${FRAGILE_LIST_MAX + 1} (group × time) cells have N < ${FRAGILE_N}`), true,
      'collapsed summary states the count');
    expect(txt.includes('worst N=1 for G · 0 at 0 h'), true, 'collapsed summary names the worst cell');
    // Collapsed form must NOT enumerate every cell (the gap-18 flood it prevents).
    expect((txt.match(/N=2 for/g) ?? []).length, 0, 'does not enumerate the non-worst cells');
  });

  test('AC-F4.3 (DATA-GAP-30c / R2 NEW-3): low cardinality (≤MAX) enumerates at high fragile fraction', async () => {
    // 4 fragile cells (e.g. "3 of 4 fragile" = 75%) but ≤ FRAGILE_LIST_MAX → MUST enumerate,
    // NOT collapse — the regression pin against a fraction trigger with no floor.
    const cells = [
      cell(`M${KEY_SEP}2`, 1, 3), cell(`M${KEY_SEP}10`, 1, 2),
      cell(`F${KEY_SEP}2`, 1, 1), cell(`F${KEY_SEP}10`, 1, 3)];
    const txt = fragilityNoticeText(cells, true);
    expect(txt.includes('cells have N <'), false,
      'does NOT collapse at 4 ≤ FRAGILE_LIST_MAX (no unbounded fraction trigger)');
    expect(txt.includes('N=1 for F · 2 at 1 h'), true, 'enumerates the worst cell');
    expect((txt.match(/N=\d+ for/g) ?? []).length, 4, 'enumerates all four fragile cells');
  });

  test('observed-time pooling appends the nominal-time advice; nominal-time binning omits it', async () => {
    const cells = [cell('All', 24, 1)];
    expect(fragilityNoticeText(cells, /* nominalUsed */ false).includes('Map a nominal-time column'),
      true, 'observed-time (nominalUsed=false) appends the "Map a nominal-time column" fix');
    expect(fragilityNoticeText(cells, /* nominalUsed */ true).includes('nominal-time column'), false,
      'nominal-time binning (nominalUsed=true) omits the advice');
  });
});

category('multi-column-group-by: composite legend ordering (F5)', () => {
  test('natural sort on the joined composite key is numeric-within-join', async () => {
    // Synthesize a composite Group column (Sex × Dose) exactly as the view does, then
    // pool it through buildMeanSdProfile — its setCategoryOrder must order the doses
    // numerically within each Sex (M·2, M·10, M·100), not lexically (M·10, M·100, M·2).
    const sex = ['M', 'M', 'M', 'F']; const dose = ['2', '10', '100', '2'];
    const group = sex.map((s, i) => `${s}${KEY_SEP}${dose[i]}`);
    const obs = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('G', [...group, ...group]),
      DG.Column.fromFloat64Array('T', Float64Array.from([1, 1, 1, 1, 2, 2, 2, 2])),
      DG.Column.fromFloat64Array('C', Float64Array.from([5, 5, 5, 5, 3, 3, 3, 3])),
    ]);
    const profile = buildMeanSdProfile(obs, {timeCol: 'T', concCol: 'C', groupCol: 'G'});
    const order = profile.col('Group')!.categories;
    const idx = (s: string, d: string): number => order.indexOf(`${s}${KEY_SEP}${d}`);
    expect(idx('M', '2') < idx('M', '10') && idx('M', '10') < idx('M', '100'), true,
      'doses order numerically within the joined key (2 < 10 < 100), not lexically');
  });
});

category('group-by-duplicate-restore (BUG-11): the selection is a SET', () => {
  /** A stand-in results frame carrying the candidate keys of the reported study. */
  const resultsFrame = (): DG.DataFrame => DG.DataFrame.fromColumns([
    DG.Column.fromStrings('ROUTE', ['IV-bolus', 'PO']),
    DG.Column.fromInt32Array('Dose_mg_kg', Int32Array.from([1, 5])),
    DG.Column.fromStrings('Sex', ['F', 'M']),
    DG.Column.fromStrings('CL_BASIS', ['CL', 'CL/F']),
  ]);

  test('the exact reported value collapses to its 3 distinct keys, in pick order', async () => {
    // The user's ribbon read "(5) ROUTE, Dose_mg_kg, ROUTE, Sex, Dose_mg_kg" on reopen.
    const reported = ['ROUTE', 'Dose_mg_kg', 'ROUTE', 'Sex', 'Dose_mg_kg'];
    const cols = dedupeGroupBySelection(reported, resultsFrame());
    expect(cols.map((c) => c.name).join(','), 'ROUTE,Dose_mg_kg,Sex',
      'collapses to the distinct names, FIRST-seen order (the user\'s pick order) preserved');
  });

  test('the picker\'s growth step is a fixed point, not a ratchet (5 → 8 → 3)', async () => {
    // The measured defect: a no-op open→OK re-added every already-selected NAME (5 → 8, +3 =
    // the 3 distinct keys) because the dialog ticks `value ++ tableColumns` by name. Simulate
    // that exact commit on top of the already-duplicated value and pin that it converges.
    const frame = resultsFrame();
    const grown = ['ROUTE', 'Dose_mg_kg', 'ROUTE', 'Sex', 'Dose_mg_kg', 'ROUTE', 'Sex', 'Dose_mg_kg'];
    const once = dedupeGroupBySelection(grown, frame);
    expect(once.length, 3, 'the 8-entry grown value collapses back to 3');
    // Idempotent: re-applying to its own output must not drift (the write-back re-enters).
    const twice = dedupeGroupBySelection(once.map((c) => c.name), frame);
    expect(twice.map((c) => c.name).join(','), once.map((c) => c.name).join(','),
      'idempotent — the write-back cannot loop or reorder');
  });

  test('entries resolve to LIVE columns of the given frame, and unknown names drop', async () => {
    const frame = resultsFrame();
    // A stale value column (same name, different frame) is what seeds the duplicate: the
    // dialog lists it as an extra row beside the table's own. Resolution must return the
    // frame's OWN column object, not the one handed in.
    const stale = DG.DataFrame.fromColumns([DG.Column.fromStrings('ROUTE', ['IV-bolus', 'PO'])]);
    const cols = dedupeGroupBySelection(['ROUTE'], frame);
    expect(cols[0] === frame.col('ROUTE'), true, 'returns the live frame column');
    expect(cols[0] === stale.col('ROUTE'), false, 'never the stale same-named column');
    // A name the frame no longer carries (mapping changed) drops rather than throwing.
    expect(dedupeGroupBySelection(['Sex', 'Vanished', 'ROUTE'], frame).map((c) => c.name).join(','),
      'Sex,ROUTE', 'a name absent from the frame is dropped, the rest survive in order');
    expect(dedupeGroupBySelection([], frame).length, 0, 'empty selection stays empty (ungrouped)');
  });

  test('the invariant is NOT optional: aggregateNcaResults rejects a duplicated key list', async () => {
    // Why the set invariant has to hold at the SOURCE (this is the finding that makes the
    // pane-side dedup load-bearing rather than cosmetic): the aggregator emits one output
    // column per key, so a repeated key is not a harmless no-op partition — it THROWS.
    // Every caller survives today only because it dedupes first (`resolveAggregateGroupKeys`
    // → `Array.from(new Set(…))`; `buildSummaryCache` likewise). Pinned on the real
    // aggregator, not inferred (rule 16).
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromStrings('Sex', ['F', 'M', 'F', 'M']),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([10, 20, 30, 40])),
    ]);
    let threw: unknown = null;
    try {aggregateNcaResults(df, ['ROUTE', 'Sex', 'ROUTE'], {params: ['AUC_INF']});} catch (e) {threw = e;}
    expect(threw !== null, true, 'a duplicated key list is rejected, not silently double-grouped');
  });

  test('deduped keys bucket exactly like the clean pick (no partition changed by the repair)', async () => {
    // The other half: once collapsed, the grouping IS the 2-key grouping the user picked —
    // so repairing a poisoned selection cannot move a single number. This is what licenses
    // MIGRATIONS[11] to rewrite `groupByCols` without re-projecting the cached tier.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('ROUTE', ['IV-bolus', 'IV-bolus', 'PO', 'PO']),
      DG.Column.fromStrings('Sex', ['F', 'M', 'F', 'M']),
      DG.Column.fromFloat64Array('AUC_INF', Float64Array.from([10, 20, 30, 40])),
    ]);
    const viaDup = resolveAggregateGroupKeys(df, ['ROUTE', 'Sex', 'ROUTE']);
    expect(viaDup.keys.join(','), 'ROUTE,Sex', 'the poisoned tuple resolves to the clean key set');
    const clean = aggregateNcaResults(df, ['ROUTE', 'Sex'], {params: ['AUC_INF']});
    const repaired = aggregateNcaResults(df, viaDup.keys, {params: ['AUC_INF']});
    expect(repaired.rowCount, clean.rowCount, 'same partition — same row count');
    const means = (d: DG.DataFrame): string => {
      const c = d.col('MEAN_AUC_INF')!;
      return Array.from({length: d.rowCount}, (_v, i) => c.get(i) as number).sort((x, y) => x - y).join(',');
    };
    expect(means(repaired), means(clean), 'identical group means — the repair moves no number');
  });
});
