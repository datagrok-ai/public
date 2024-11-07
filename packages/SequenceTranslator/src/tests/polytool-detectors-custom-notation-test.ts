import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, before} from '@datagrok-libraries/utils/src/test';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {ALIGNMENT, ALPHABET, NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {_testNeg, _testPos, DetectorTestData, DfReaderFunc, PosCol} from './utils/detect-macromolecule-utils';

category('PolyTool: detectors', () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  const tests: DetectorTestData = {
    'cyclized1': {
      csv: `n,seq
1,R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI
2,C(1)-T-G-H-F-Y-P-C(1)-meI
3,R-F-C(1)-T-G-H-F-Y-P-C(1)
4,C(1)-T-G-H-F-H-P-C(1)
5,R-F-D(2)-T-G-H-F-Y-P-NH2(2)
6,R-F-aG(4)-T-G-H-F-Y-P-azG(4)-meI`,
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
        await _testPos(reader, posColName, seqHelper, posCol.units, posCol.aligned,
          posCol.alphabet, posCol.alphabetSize, posCol.alphabetIsMultichar, posCol.separator);
      }
    });
  }
});
