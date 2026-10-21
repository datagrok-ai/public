/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {combineSequences} from '../polytool/pt-combine-dialog';

import {_package} from '../package-test';

category('PolyTool: Combine Sequences', () => {
  test('lists in turn', async () => {
    expect(combineSequences([['A', 'B'], ['C'], ['D', 'E']], '.').join(', '), 'A.C.D, A.C.E, B.C.D, B.C.E');
  });

  test('cyclized.csv seqs with itself', async () => {
    const seqs: string[] = DG.DataFrame.fromCsv(await _package.files.readAsText('samples/cyclized.csv')).getCol('seqs').toList();
    const res = combineSequences([seqs, seqs], '-');
    expect(res.length, 196);
    expect(res[0], 'R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI-R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI');
    expect(res[1], `${seqs[0]}-${seqs[1]}`);
    expect(res[195], `${seqs[13]}-${seqs[13]}`);
  });
});
