import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';
import {findMutations} from '../utils/algorithms';
import * as type from '../utils/types';

category('Algorithms', () => {
  let activityCol: type.RawData;
  let monomerColumns: type.RawColumn[];
  let settings: type.PeptidesSettings;

  before(async () => {
    activityCol = DG.Column.fromList('int', 'test', [1, 2, 5]).getRawData();
    monomerColumns = [
      DG.Column.fromList('string', '1', 'AAA'.split('')),
      DG.Column.fromList('string', '2', 'BCC'.split('')),
      DG.Column.fromList('string', '3', 'CCD'.split('')),
    ].map((col) => ({
      name: col.name,
      rawData: col.getRawData(),
      cat: col.categories,
    }));
    settings = {maxMutations: 1, minActivityDelta: 2};
  });

  test('MutationCliffs', async () => {
    const substInfo: type.SubstitutionsInfo = findMutations(activityCol, monomerColumns, settings);
    expect(substInfo.has('C'), true);
    expect(substInfo.has('D'), true);
    expect(substInfo.has('A'), false);

    const c = substInfo.get('C')!;
    const d = substInfo.get('D')!;
    expect(c.has('3'), true);
    expect(d.has('3'), true);

    const c3 = c.get('3')!;
    const d3 = d.get('3')!;
    expect(c3.has(1), true);
    expect(d3.has(2), true);

    const c31 = c3.get(1)!;
    const d32 = d3.get(2)!;
    expect(c31[0], 2);
    expect(d32[0], 1);
  });

  test('MutationCliffs - Benchmark 5k', async () => {
    const df = (await _package.files.readBinaryDataFrames('tests/aligned_5k.d42'))[0];
    const activityCol: type.RawData = df.getCol('Activity').getRawData();
    const monomerCols: type.RawColumn[] = [];
    for (let i = 1; i < 16; ++i) {
      const col = df.getCol(i.toString());
      monomerCols.push({name: col.name, rawData: col.getRawData(), cat: col.categories});
    }
    DG.time('MutationCliffs', () => findMutations(activityCol, monomerCols));
  }, {skipReason: 'Benchmark'});
});
