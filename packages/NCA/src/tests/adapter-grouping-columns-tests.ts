/**
 * Adapter — grouping-column propagation category (NOTE-11 gap 1).
 *
 * @see src/adapter/grouping-columns.ts (implementation under test)
 *
 * Coverage:
 * - Constant-per-subject metadata (CompoundID / Sex / Dose) is carried into the
 *   results frame, aligned by profile order; values read at observationRows[0].
 * - Per-row-varying columns are rejected (not valid group keys).
 * - Analysis inputs (subjectId / time / concentration) are excluded by name —
 *   even in the sparse 1-row-per-subject case where they are trivially constant.
 * - Type preservation (STRING / FLOAT) and name-collision skip.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {nca} from '@datagrok-libraries/sci-comp';

import {buildProfileIndices} from '../adapter/profile-index';
import {appendGroupingColumns} from '../adapter/grouping-columns';
import {aggregateNcaResults} from '../adapter/results-aggregator';
import {DEFAULT_RULES} from '../adapter/defaults';
import {runPipelineOrThrow} from './test-utils';
import type {ProfileScalars} from '../adapter/profile-inputs';
import type {ColumnContract} from '../adapter/column-contract';

function baseContract(): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: 'string'},
      time: {name: 'Time', required: true, rawType: 'float64', units: 'h'},
      concentration: {name: 'conc', required: true, rawType: 'float64', units: 'mg/L'},
    },
  };
}

/** Minimal results frame: a SUBJECT column in profile order (lex-sorted keys). */
function resultsStub(subjects: ReadonlyArray<string>): DG.DataFrame {
  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SUBJECT', Array.from(subjects)),
  ]);
}

category('Adapter: grouping-column propagation', () => {
  test('carries constant-per-subject metadata, skips varying + inputs', async () => {
    // 3 subjects, 2 rows each. CompoundID/Sex/Dose constant per subject;
    // Note varies within a subject (n1≠n2) → must be rejected.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S1', 'S2', 'S2', 'S3', 'S3']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 8, 4, 12, 6])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'A', 'A', 'B', 'B']),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Sex', ['M', 'M', 'F', 'F', 'M', 'M']),
      DG.Column.fromFloat64Array('Dose', Float64Array.from([10, 10, 10, 10, 50, 50])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Note', ['n1', 'n2', 'n1', 'n2', 'n1', 'n2']),
    ]);
    const contract = baseContract();
    const profiles = buildProfileIndices(df, contract);
    const results = resultsStub(profiles.map((p) => p.subjectKey));

    const added = appendGroupingColumns(results, df, profiles, contract);

    expect(added.includes('CompoundID'), true, 'CompoundID carried');
    expect(added.includes('Sex'), true, 'Sex carried');
    expect(added.includes('Dose'), true, 'Dose carried (valid group key)');
    expect(added.includes('Note'), false, 'varying Note rejected');
    expect(added.includes('Time'), false, 'Time (input) excluded');
    expect(added.includes('conc'), false, 'conc (input) excluded');
    expect(added.includes('Subject'), false, 'Subject (subjectId) excluded');

    // Values aligned by profile order (S1, S2, S3 lex-sorted).
    expect(results.col('CompoundID')!.get(0), 'A', 'S1 → A');
    expect(results.col('CompoundID')!.get(2), 'B', 'S3 → B');
    expect(results.col('Sex')!.get(1), 'F', 'S2 → F');
    expect(results.col('Dose')!.get(2), 50, 'S3 → 50');
    // Type preservation.
    expect(results.col('Dose')!.type, DG.COLUMN_TYPE.FLOAT, 'Dose stays FLOAT');
    expect(results.col('CompoundID')!.type, DG.COLUMN_TYPE.STRING, 'CompoundID stays STRING');
  });

  test('sparse 1-row-per-subject: inputs still excluded by name', async () => {
    // One row per subject → Time/conc are trivially constant, but must NOT be
    // carried (they are analysis inputs). CompoundID still carried.
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2', 'S3']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([1, 2, 3])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 8, 6])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', ['A', 'A', 'B']),
    ]);
    const contract = baseContract();
    const profiles = buildProfileIndices(df, contract);
    const results = resultsStub(profiles.map((p) => p.subjectKey));

    const added = appendGroupingColumns(results, df, profiles, contract);
    expect(added.includes('CompoundID'), true, 'CompoundID carried');
    expect(added.includes('Time'), false, 'Time excluded even when 1 row/subject');
    expect(added.includes('conc'), false, 'conc excluded even when 1 row/subject');
  });

  test('name collision with existing results column is skipped', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 0])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 8])),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'GROUP', ['X', 'Y']),
    ]);
    const contract = baseContract();
    const profiles = buildProfileIndices(df, contract);
    // Results frame already has a GROUP column → source GROUP must be skipped.
    const results = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'SUBJECT', profiles.map((p) => p.subjectKey)),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'GROUP', ['pre', 'existing']),
    ]);

    const added = appendGroupingColumns(results, df, profiles, contract);
    expect(added.includes('GROUP'), false, 'colliding GROUP skipped');
    expect(results.col('GROUP')!.get(0), 'pre', 'existing GROUP value untouched');
  });

  test('round-trip: pipeline → append → aggregate by carried CompoundID (unblocks #11)', async () => {
    // 2 compounds × 2 subjects, clean 5-point IV-bolus decay. CompoundID is the
    // group key that does not exist in the results frame until propagated.
    const subjects = ['C1S1', 'C1S2', 'C2S1', 'C2S2'];
    const times = [0, 1, 2, 4, 8];
    const sCol: string[] = [];
    const tCol: number[] = [];
    const cCol: number[] = [];
    const compound: string[] = [];
    for (const s of subjects) {
      for (const t of times) {
        sCol.push(s);
        tCol.push(t);
        cCol.push(10 * Math.pow(0.5, t)); // log-linear decay, λz ≈ ln2
        compound.push(s.startsWith('C1') ? 'CPD-1' : 'CPD-2');
      }
    }
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', sCol),
      DG.Column.fromFloat64Array('Time', Float64Array.from(tCol)),
      DG.Column.fromFloat64Array('conc', Float64Array.from(cCol)),
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'CompoundID', compound),
    ]);
    const contract = baseContract();
    const scalars: ProfileScalars = {
      dose: 100, route: nca.ROUTE_IV_BOLUS, doseUnits: 'mg', concentrationUnits: 'mg/L',
      timeUnits: 'h', bodyWeight: null, bodyWeightUnit: 'kg', infusionDuration: null, lloq: 0,
    };

    const {results, outcome} = runPipelineOrThrow(df, contract, scalars, DEFAULT_RULES);
    expect(results.col('CompoundID'), null, 'CompoundID absent from raw results');

    const added = appendGroupingColumns(results, df, outcome.perProfile.map((p) => p.index), contract);
    expect(added.includes('CompoundID'), true, 'CompoundID propagated into results');

    const agg = aggregateNcaResults(results, 'CompoundID');
    expect(agg.rowCount, 2, 'two compound groups');
    expect(agg.col('MEAN_CMAX') !== null, true, 'MEAN_CMAX present');
    expect(agg.col('N') !== null, true, 'N present');
    // Each group has 2 subjects.
    for (let i = 0; i < agg.rowCount; i++)
      expect(agg.col('N')!.get(i), 2, `group ${agg.col('CompoundID')!.get(i)} has N=2`);
  });

  test('throws on row/profile mismatch', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['S1', 'S2']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 0])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 8])),
    ]);
    const contract = baseContract();
    const profiles = buildProfileIndices(df, contract);
    const results = resultsStub(['only-one-row']); // rowCount 1 ≠ 2 profiles
    let threw = false;
    try {
      appendGroupingColumns(results, df, profiles, contract);
    } catch {
      threw = true;
    }
    expect(threw, true, 'mismatch must throw');
  });
});
