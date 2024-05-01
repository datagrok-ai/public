import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {category, test} from '@datagrok-libraries/utils/src/test';
import {MonomerPlacer} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {MonomerLibManager} from '../utils/monomer-lib/lib-manager';

import {_package} from '../package-test';

category('renderers: monomerPlacer', () => {
  const tests = {
    splitter: {
      /**
       0 = Array(10) [0, 26, 45, 71, 97, 123, 142, 161, 187, 213],
       1 = Array(9) [0, 26, 52, 78, 104, 130, 156, 182, 208],
       2 = Array(8) [0, 26, 45, 71, 97, 123, 149, 175],
       * */
      csv: 'id,seq\n' +
        'id1,m1-M-m3-mon4-mon5-N-T-MON8-N9\n' + //Array(10) [0, 26, 52, 78, 104, 130, 156, 175, 201, 227]
        'id2,m1-mon2-m3-mon4-mon5-Num-MON8-N9\n' + //
        'id3,mon1-M-mon3-mon4-mon5-MON8-N9\n', //
      testList: [
        {src: {row: 0, x: -1}, tgt: {pos: null}},
        {src: {row: 1, x: 0}, tgt: {pos: null}},
        {src: {row: 1, x: 5}, tgt: {pos: 0}},
        {src: {row: 1, x: 6}, tgt: {pos: 0}},
        {src: {row: 1, x: 26}, tgt: {pos: 1}},
        {src: {row: 1, x: 160}, tgt: {pos: 6}},
        {src: {row: 1, x: 185}, tgt: {pos: 7}},
        {src: {row: 2, x: 140}, tgt: {pos: 5}},
        {src: {row: 2, x: 145}, tgt: {pos: 5}},
      ]
    },
    splitterMsa: {
      /** For charWidth=7 and sepWidth=12, MSA
       * Array(10) [5, 38, 71, 104, 137, 170, 203, 222, 255, 281]
       */
      csv: 'id,seq\n' +
        'id1,m1-M-m3-mon4-mon5-N-T-MON8-N9\n' +
        'id2,m1-mon2-m3-mon4-mon5-Num--MON8-N9\n' +
        'id3,\n' + // empty
        'id4,mon1-M-mon3-mon4-mon5---MON8-N9\n',
      testList: [
        {src: {row: 0, x: -1}, tgt: {pos: null}},
        {src: {row: 1, x: 0}, tgt: {pos: null}},
        {src: {row: 1, x: 1}, tgt: {pos: null}},
        {src: {row: 1, x: 4}, tgt: {pos: null}},
        {src: {row: 1, x: 5}, tgt: {pos: 0}},
        {src: {row: 1, x: 37}, tgt: {pos: 0}},
        {src: {row: 1, x: 38}, tgt: {pos: 1}},
        {src: {row: 1, x: 170}, tgt: {pos: 4}},
        {src: {row: 1, x: 200}, tgt: {pos: 5}},
        {src: {row: 2, x: 20}, tgt: {pos: null}}, // empty value
        {src: {row: 3, x: 170}, tgt: {pos: 4}},
        {src: {row: 3, x: 200}, tgt: {pos: 5}},
        {src: {row: 3, x: 282}, tgt: {pos: null}},
      ]
    },
    fastaMsa: {
      /** For charWidth=7 and sepWidth=12, MSA
       * Array(10) [0, 19, 38, 57, 76, 95, 114, 133, 152, 171]
       */
      csv: `id,seq
id1,QQYNIYPLT
id2,QQWSSFPYT
id3,
id3,QHIRE--LT
`,
      testList: [
        {src: {row: 1, x: -1}, tgt: {pos: null}},
        {src: {row: 1, x: 0}, tgt: {pos: null}},
        {src: {row: 1, x: 1}, tgt: {pos: null}},
        {src: {row: 1, x: 19}, tgt: {pos: 0}},
        {src: {row: 1, x: 170}, tgt: {pos: 8}},
        {src: {row: 1, x: 171}, tgt: {pos: 8}},
        {src: {row: 2, x: 5}, tgt: {pos: null}}, // empty value
        {src: {row: 3, x: 170}, tgt: {pos: 8}},
        {src: {row: 3, x: 181}, tgt: {pos: null}},
      ]
    },
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`getPosition-${testName}`, async () => {
      const df: DG.DataFrame = DG.DataFrame.fromCsv(testData.csv);
      await grok.data.detectSemanticTypes(df);
      const seqCol: DG.Column = df.getCol('seq');

      const monLength: number = 3;
      const charWidth: number = 7;
      const sepWidth: number = 12;
      const colTemp: MonomerPlacer = new MonomerPlacer(null, seqCol, _package.logger, () => {
        const sh = SeqHandler.forColumn(seqCol);
        return {
          seqHandler: sh,
          monomerCharWidth: charWidth,
          separatorWidth: sepWidth,
          monomerToShort: monomerToShort,
          monomerLengthLimit: monLength,
          monomerLib: MonomerLibManager.instance.getBioLib(),
        };
      });

      const testList = testData.testList;
      // simulate rendering
      for (let rowI: number = 0; rowI < seqCol.length; ++rowI)
        colTemp.getCellMonomerLengths(rowI);

      const errorList: string[] = [];
      for (const [test, _testI] of wu.enumerate(testList)) {
        const res = {pos: colTemp.getPosition(test.src.row, test.src.x)};
        if (test.tgt.pos != res.pos) {
          errorList.push(`Test src ${JSON.stringify(test.src)} expected tgt ${JSON.stringify(test.tgt)},` +
            ` but get ${JSON.stringify({res})}`);
        }
      }
      if (errorList.length > 0)
        throw new Error('Test failed error(s):\n' + errorList.join(', n'));
    });
  }
});
