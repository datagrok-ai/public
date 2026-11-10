/**
 * Adapter — profile-index category.
 *
 * @see src/adapter/profile-index.ts (implementation under test)
 *
 * Profiles are returned sorted lexicographically by `subjectKey`; within a
 * profile the row order matches the DataFrame's natural order. UC-01 has no
 * separate dosing DataFrame — the index is `{subjectKey, observationRows}`.
 */

import * as DG from 'datagrok-api/dg';
import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';

import {buildProfileIndices} from '../adapter/profile-index';
import type {ColumnContract} from '../adapter/column-contract';

function strContract(): ColumnContract {
  return {
    observations: {
      subjectId: {name: 'Subject', required: true, rawType: 'string'},
      time: {name: 'Time', required: true, rawType: 'float64'},
      concentration: {name: 'conc', required: true, rawType: 'float64'},
    },
  };
}

category('Adapter: profile index', () => {
  test('three subjects → three profiles, lexicographically sorted', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject',
        ['C', 'C', 'A', 'A', 'B', 'B']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 8, 4, 12, 6])),
    ]);
    const profiles = buildProfileIndices(df, strContract());
    expect(profiles.length, 3, 'expected 3 profiles');
    expect(profiles[0].subjectKey, 'A', 'profile[0] should be A (lex sort)');
    expect(profiles[1].subjectKey, 'B', 'profile[1] should be B');
    expect(profiles[2].subjectKey, 'C', 'profile[2] should be C');
  });

  test('observationRows preserves original DataFrame row order within profile', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['B', 'A', 'A', 'B', 'A']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 0, 1, 1, 2])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([5, 10, 5, 2, 1])),
    ]);
    const profiles = buildProfileIndices(df, strContract());
    const a = profiles.find((p) => p.subjectKey === 'A')!;
    const b = profiles.find((p) => p.subjectKey === 'B')!;
    expectArray(Array.from(a.observationRows), [1, 2, 4]);
    expectArray(Array.from(b.observationRows), [0, 3]);
  });

  test('integer-id subject column groups by stringified value', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromInt32Array('Subject', Int32Array.from([2, 2, 1, 1])),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1, 0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5, 8, 4])),
    ]);
    const contract: ColumnContract = {
      observations: {
        subjectId: {name: 'Subject', required: true, rawType: 'int32'},
        time: {name: 'Time', required: true, rawType: 'float64'},
        concentration: {name: 'conc', required: true, rawType: 'float64'},
      },
    };
    const profiles = buildProfileIndices(df, contract);
    expect(profiles.length, 2, 'expected 2 profiles');
    expect(profiles[0].subjectKey, '1', 'lex-sort puts "1" before "2"');
    expect(profiles[1].subjectKey, '2');
  });

  test('observationRows is an Int32Array (typed-array invariant)', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'Subject', ['X', 'X']),
      DG.Column.fromFloat64Array('Time', Float64Array.from([0, 1])),
      DG.Column.fromFloat64Array('conc', Float64Array.from([10, 5])),
    ]);
    const [p] = buildProfileIndices(df, strContract());
    expect(p.observationRows instanceof Int32Array, true,
      `observationRows must be Int32Array; got ${p.observationRows.constructor.name}`);
  });
});
