import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, delay, before} from '@datagrok-libraries/utils/src/test';

import {_package} from '../package-test';
import {startAnalysis} from '../widgets/peptides';
import {PeptidesModel} from '../model';
import * as C from '../utils/constants';
import {scaleActivity} from '../utils/misc';
import {ALPHABET, TAGS, NOTATION, ALIGNMENT} from '@datagrok-libraries/bio';
import {findMutations} from '../utils/algorithms';
import * as type from '../utils/types';

category('Algorithms', () => {
  let activityCol: DG.Column<number>;
  let monomerColumns: DG.Column<string>[];
  let settings: type.PeptidesSettings;

  before(async () => {
    activityCol = DG.Column.fromList('int', 'test', [1, 2, 5]);
    monomerColumns = [
      DG.Column.fromList('string', '1', 'ABC'.split('')),
      DG.Column.fromList('string', '2', 'ACC'.split('')),
      DG.Column.fromList('string', '3', 'ACD'.split('')),
    ];
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
    expect(c3.has(2), true);
    expect(d3.has(3), true);

    const c32 = c3.get(2)!;
    const d33 = d3.get(3)!;
    expect(c32[0], 3);
    expect(d33[0], 2);
  });

  test('MutationCliffs - Benchmark 5k', async () => {
    const df = (await _package.files.readBinaryDataFrames('tests/aligned_5k.d42'))[0];
    const activityCol: DG.Column<number> = df.getCol('Activity');
    const monomerCols: DG.Column<string>[] = [];
    for (let i = 1; i < 16; ++i)
      monomerCols.push(df.getCol(i.toString()));
    DG.time('MutationCliffs', () => findMutations(activityCol, monomerCols));
  }, {skipReason: 'Benchmark'});
});
