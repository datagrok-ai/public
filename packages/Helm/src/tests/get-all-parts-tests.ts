import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';

import {getAllParts} from '../helm-monomer-placer';
import {initHelmMainPackage} from './utils';

category('getAllParts', () => {
  /** All parts of Helm string are required for fallback rendering. */

  const tests: { [testName: string]: { tgt: string[]; seq: string } } = {
    'peptideSimpleSingleChar': {
      seq: 'PEPTIDE1{L.V.A}$$$$',
      tgt: ['PEPTIDE1{', 'L', '.', 'V', '.', 'A', '}$$$$'],
    },
    'peptideSimpleMultiChar': {
      seq: 'PEPTIDE1{[Ac].V.A}$$$$',
      /* doubt joining square brackets to separator or notation is correct */
      tgt: ['PEPTIDE1{[', 'Ac', '].', 'V', '.', 'A', '}$$$$'],
    }
  };

  before(async () => {
    await initHelmMainPackage();
  });

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      await _testAllParts(testData.seq, testData.tgt);
    });
  }
});

async function _testAllParts(seq: string | null, tgt: string[]): Promise<void> {
  const res: string[] = getAllParts(seq);
  expectArray(res, tgt);
}
