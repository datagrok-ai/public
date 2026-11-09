/**
 * NCA Studies — stored shape + snapshot hash + engine staleness.
 *
 * @see src/adapter/study-schema.ts, src/adapter/nca-repository.ts (under test)
 *
 * Coverage:
 * - A study with a populated passthrough selection survives the JSON round-trip its json columns
 *   make (AC-F5b) — the designation is never re-detected (LBC-2).
 * - FNV-1a is deterministic and order-sensitive.
 * - `computeSnapshotHash` is EXACT: equal frames hash equal; any change of a value (anywhere — not a
 *   head/tail sample), of the row count, of the column order or of a cell's type changes it; a
 *   missing cell is not a NaN.
 * - `isStale` — the ONE engine-staleness predicate (the design-decisions "Engine-staleness cue"):
 *   no cache → never stale; an unknown or a different engine → stale.
 *
 * The file-format migration ladder (v1..v16) and its tests were removed with the AppData store.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';

import {fnv1a} from '../adapter/study-schema';
import type {NcaStudy} from '../adapter/study-schema';
import {computeSnapshotHash, isStale, isBasisStale, snapshotFrame} from '../adapter/nca-repository';
import {RESULTS_BASIS} from '../adapter/unit-conversion';

category('NCA Studies: stored shape', () => {
  test('a study with a populated passthrough set survives the JSON round-trip (AC-F5b)', async () => {
    // The persist→reload half of the science gate: a real passthrough selection (name +
    // isStructure) + per-row passthrough entries round-trip through the json columns with the
    // designation intact (LBC-2: no re-detect).
    const study: NcaStudy = {
      id: 'p-1', label: 'Demo', kind: 'study',
      createdAt: '2026-06-28T00:00:00.000Z', updatedAt: '2026-06-28T00:00:00.000Z',
      access: {canEdit: true, canDelete: true, canShare: true},
      analyses: [{
        id: 'an-1', label: 'A1', source: {kind: 'stored', ref: 'System:DomainFiles/nca/analysis/x/source.d42'},
        columns: {} as never, modes: {} as never, ncaRules: {} as never, units: {} as never,
        groupByCols: ['CompoundID'],
        passthroughCols: [{name: 'SMILES', isStructure: true}, {name: 'Series', isStructure: false}],
        sourceSignature: {kind: 'content-hash', rowCount: 6, mtime: null, hash: 'abc',
          computedAt: '2026-07-13T00:00:00.000Z'},
        summaryCache: {
          analysisTier: {nCompounds: 1, nSubjects: 3, route: 'IV-bolus', dose: '(per-row)',
            status: 'ok', flagSummary: 'No flags'},
          compoundTier: [{
            compoundId: 'DGK-1', isWholeAnalysis: false,
            groupValues: [{key: 'CompoundID', value: 'DGK-1'}], keyInventory: [],
            passthrough: [{key: 'SMILES', values: [{value: 'CCO', count: 3}], isStructure: true}],
            n: 3, route: 'IV-bolus', doseNorm: 'per_kg',
            cmax: {mean: 10, cvPct: 5, n: 3}, tmax: {mean: 1, cvPct: 20, n: 3},
            aucInf: {mean: 100, cvPct: 8, n: 3}, thalf: {mean: 2, cvPct: 12, n: 3},
            cl: {mean: 1, cvPct: 9, n: 3}, vz: {mean: 3, cvPct: 11, n: 3},
            mrt: {mean: 4, cvPct: 10, n: 3}, vss: {mean: 6, cvPct: 9, n: 3},
            fPct: null, pctExtrap: {mean: 15, cvPct: 4, n: 3}, curve: null, warnings: [],
          }],
        },
      }],
    };
    const round = JSON.parse(JSON.stringify(study)) as NcaStudy;
    expect(JSON.stringify(round), JSON.stringify(study), 'byte-identical after the round-trip');
    const a = round.analyses[0];
    expect(a.passthroughCols.length, 2, 'both passthrough specs survive');
    expect(a.passthroughCols[0].isStructure, true, 'structure designation survives (LBC-2 — never re-detected)');
    expect(a.passthroughCols[1].isStructure, false, 'text designation survives');
    const pe = a.summaryCache!.compoundTier[0].passthrough[0];
    expect(pe.key, 'SMILES');
    expect(pe.values[0].count, 3, 'passthrough subject count survives');
    expect(pe.isStructure, true, 'passthrough entry structure flag survives');
  });
});

category('NCA Studies: snapshot hash', () => {
  test('fnv1a is deterministic and order-sensitive', async () => {
    expect(fnv1a('abc') === fnv1a('abc'), true, 'same input → same hash');
    expect(fnv1a('abc') !== fnv1a('acb'), true, 'reordered input → different hash');
    expect(fnv1a('') === fnv1a(''), true, 'empty stable');
  });

  test('equal frames ⇒ equal snapshot hash (no spurious "unsaved data")', async () => {
    const csv = 'SUBJ,TIME,CONC\nA,0,0\nA,1,5.5\nA,2,3.2\nB,0,0\nB,1,6.1\nB,2,2.8\n';
    expect(computeSnapshotHash(DG.DataFrame.fromCsv(csv)), computeSnapshotHash(DG.DataFrame.fromCsv(csv)));
  });

  test('row count, column order, a mid-table value and a cell type each change the hash (exact, SD-19)', async () => {
    const rows: string[] = ['SUBJ,TIME,CONC'];
    for (let i = 0; i < 25; i++) rows.push(`A,${i},${i + 1}`);
    const full = DG.DataFrame.fromCsv(rows.join('\n'));
    const h = computeSnapshotHash(full);
    expect(computeSnapshotHash(DG.DataFrame.fromCsv(rows.slice(0, 25).join('\n'))) !== h, true, 'a dropped row');
    const mid = rows.slice();
    mid[1 + 12] = 'A,12,999'; // a middle row — the retired sampled hash did not see it
    expect(computeSnapshotHash(DG.DataFrame.fromCsv(mid.join('\n'))) !== h, true, 'a mid-table value');
    const reordered = DG.DataFrame.fromCsv(['CONC,SUBJ,TIME', ...rows.slice(1).map((r) => {
      const [s, t, c] = r.split(',');
      return `${c},${s},${t}`;
    })].join('\n'));
    expect(computeSnapshotHash(reordered) !== h, true, 'column order is part of the snapshot');
    const asFloat = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('X', Float64Array.from([1, 2]))]);
    const asInt = DG.DataFrame.fromColumns([DG.Column.fromInt32Array('X', Int32Array.from([1, 2]))]);
    expect(computeSnapshotHash(asFloat) !== computeSnapshotHash(asInt), true, 'the cell type is part of it');
  });

  // PR-13: a compute adds derived `~` columns (`~RouteCode`, …) to the source frame. They are not
  // data — a hash over them marked every reopened analysis with a per-row route dirty.
  test('a derived `~` column does not change the hash; snapshotFrame drops it', async () => {
    const csv = 'SUBJ,TIME,CONC\nA,0,0\nA,1,5.5\nA,2,3.2\n';
    const df = DG.DataFrame.fromCsv(csv);
    const h = computeSnapshotHash(df);
    df.columns.add(DG.Column.fromInt32Array('~RouteCode', Int32Array.from([1, 1, 1])));
    expect(computeSnapshotHash(df), h, 'a derived column is not part of the snapshot');
    const stored = snapshotFrame(df);
    expect(stored.columns.names().join(','), 'SUBJ,TIME,CONC', 'the snapshot holds the source columns only');
    expect(df.col('~RouteCode') !== null, true, 'the view\'s own frame keeps its derived column');
    expect(snapshotFrame(DG.DataFrame.fromCsv(csv)).columns.length, 3);
  });

  test('a missing cell is not a NaN', async () => {
    const a = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('X', Float64Array.from([1, NaN]))]);
    const b = DG.DataFrame.fromColumns([DG.Column.fromFloat64Array('X', Float64Array.from([1, 0]))]);
    b.col('X')!.set(1, null);
    expect(computeSnapshotHash(a) !== computeSnapshotHash(b), true);
  });
});

category('NCA Studies: engine staleness (F-N.6)', () => {
  test('isStale: no cache → never stale; unknown or another engine → stale; the installed one → fresh', async () => {
    const installed = '0.11.0';
    expect(isStale({engineVersion: null, hasCache: false}, installed), false, 'never computed → never stale');
    expect(isStale({engineVersion: undefined, hasCache: true}, installed), true, 'unknown engine → stale');
    expect(isStale({engineVersion: null, hasCache: true}, installed), true, 'null engine → stale');
    expect(isStale({engineVersion: '0.10.0', hasCache: true, resultsBasis: RESULTS_BASIS}, installed), true,
      'another engine → stale');
    expect(isStale({engineVersion: installed, hasCache: true, resultsBasis: RESULTS_BASIS}, installed), false,
      'the installed engine + the current basis → fresh');
  });

  // nca-cl-bw-normalised-units AD-8 / AC-7.2: the ONE predicate also reads the reporting basis.
  test('AC-7.2: a cache with no / another results basis is stale; both stamps current → fresh', async () => {
    const installed = '0.11.0';
    const at = (resultsBasis: string | null | undefined) =>
      isStale({engineVersion: installed, resultsBasis, hasCache: true}, installed);
    expect(at(undefined), true, 'no basis stamp (written before the per-kg basis) → stale');
    expect(at(null), true, 'null basis → stale');
    expect(at('absolute-0'), true, 'another basis → stale');
    expect(at(RESULTS_BASIS), false, 'current basis + installed engine → fresh');
    expect(isBasisStale({resultsBasis: null, hasCache: false}), false, 'no cache → never basis-stale');
    expect(isBasisStale({resultsBasis: RESULTS_BASIS, hasCache: true}), false, 'current basis');
  });
});
