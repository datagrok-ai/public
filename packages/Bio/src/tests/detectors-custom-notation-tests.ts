import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test} from '@datagrok-libraries/utils/src/test';
import {ALIGNMENT, ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {DetectorTestData, DfReaderFunc, PosCol} from './utils/detectors-utils';
import {_testNeg, _testPos} from './detectors-tests';

category('detectors: custom', () => {
  const tests: DetectorTestData = {
    'cyclized1': {
      csv: `n,seq
1,R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI
2,C(1)-T-G-H-F-Y-P-C(1)-meI
3,R-F-C(1)-T-G-H-F-Y-P-C(1)
4,C(1)-T-G-H-F-H-P-C(1)
5,R-F-D(2)-T-G-H-F-Y-P-NH2(2)
6,R-F-aG(3)-T-G-H-F-Y-P-azG(3)-meI`,
      pos: {'seq': new PosCol(NOTATION.CUSTOM, ALIGNMENT.SEQ, ALPHABET.UN, 13, true, '-')}
    },
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const reader: DfReaderFunc = async (): Promise<DG.DataFrame> => {
        return DG.DataFrame.fromCsv(testData.csv);
      };
      for (const negColName of testData.neg ?? [])
        await _testNeg(reader, negColName);
      for (const [posColName, posCol] of Object.entries(testData.pos ?? {})) {
        await _testPos(reader, posColName, posCol.units, posCol.aligned,
          posCol.alphabet, posCol.alphabetSize, posCol.alphabetIsMultichar, posCol.separator);
      }
    });
  }
});
