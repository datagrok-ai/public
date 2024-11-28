/* eslint-disable max-lines-per-function */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {MonomerPlacer, hitBounds} from '@datagrok-libraries/bio/src/utils/cell-renderer-monomer-placer';
import {monomerToShort} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {_package} from '../package-test';

category('renderers: monomerPlacer', () => {
  let libHelper: IMonomerLibHelper;
  let libSettings: UserLibSettings;

  before(async () => {
    libHelper = await getMonomerLibHelper();
    libSettings = await getUserLibSettings();

    await libHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(libSettings);
    await libHelper.loadMonomerLib(true);
  });

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
        {src: {row: 1, x: 190}, tgt: {pos: 7}},
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
        'id4,mon1-M-mon3-mon4-mon5---MON8-N9\n', // [ 5, 38, 71, 104, 137, 170, 203, 236, 269, 295 ]
      testList: [
        {src: {row: 0, x: -1}, tgt: {pos: null}},
        {src: {row: 1, x: 0}, tgt: {pos: null}},
        {src: {row: 1, x: 1}, tgt: {pos: null}},
        {src: {row: 1, x: 4}, tgt: {pos: null}},
        {src: {row: 1, x: 5}, tgt: {pos: 0}},
        {src: {row: 1, x: 37}, tgt: {pos: 0}},
        {src: {row: 1, x: 38}, tgt: {pos: 1}},
        {src: {row: 1, x: 170}, tgt: {pos: 5}},
        {src: {row: 1, x: 200}, tgt: {pos: 5}},
        {src: {row: 2, x: 20}, tgt: {pos: null}}, // empty value
        {src: {row: 3, x: 170}, tgt: {pos: 5}},
        {src: {row: 3, x: 200}, tgt: {pos: 5}},
        {src: {row: 3, x: 297}, tgt: {pos: null}},
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
      const colTemp: MonomerPlacer = new MonomerPlacer(null, seqCol, _package.logger, monLength,
        () => {
          return {
            font: '12px monospace',
            fontCharWidth: charWidth,
            separatorWidth: sepWidth,
            monomerToShort: monomerToShort,
          };
        });
      await colTemp.init();

      const width: number = 10000;
      const testList = testData.testList;
      // simulate rendering
      for (let rowI: number = 0; rowI < seqCol.length; ++rowI)
        colTemp.getCellMonomerLengths(rowI, 10000);

      const errorList: string[] = [];
      for (const [test, _testI] of wu.enumerate(testList)) {
        const res = {pos: colTemp.getPosition(test.src.row, test.src.x, width)};
        if (test.tgt.pos != res.pos) {
          errorList.push(`Test src ${JSON.stringify(test.src)} expected tgt ${JSON.stringify(test.tgt)},` +
            ` but get ${JSON.stringify({res})}`);
        }
      }
      if (errorList.length > 0)
        throw new Error('Test failed error(s):\n' + errorList.join(', \n'));
    });
  }

  const boundsTestData = {
    bounds: [10, 20, 30, 40, 50, 60],
    tests: {
      left: {x: 3, tgt: null},
      c0left: {x: 10, tgt: 0},
      c0mid: {x: 12, tgt: 0},
      c0right: {x: 19, tgt: 0},
      c1left: {x: 20, tgt: 1},
      c2right: {x: 39, tgt: 2},
      c4left: {x: 50, tgt: 4},
      c4right: {x: 59, tgt: 4},
      max: {x: 60, tgt: null},
      right: {x: 65, tgt: null},
    }
  };

  for (const [testName, testData] of Object.entries(boundsTestData.tests)) {
    test('hitBounds-' + testName, async () => {
      const res = hitBounds(boundsTestData.bounds, testData.x);
      expect(res, testData.tgt);
    });
  }

  const lengthsTests = {
    mono1: {
      src: {
        csv: 'seq' + '\n' +
          'm1/m2/m3/m4/m5/m6/m7/m8/m9' + '\n' +
          'n1/m2/n3/m4/n5/m6/n7/m8/n9' + '\n' +
          'm1/n2/m3/n4/m5/n6/m7/n8/m9' + '\n',
      },
      tgt: {
        lengths: [5, 31, 57, 83, 109, 135, 161, 187, 213, 239],
      }
    },
    monoWithGaps: {
      src: {
        csv: 'seq' + '\n' +
          'm1/m2/m3/m4/m5/m6//m8/m9' + '\n' +
          'n1/m2/n3/m4/n5/m6//m8/n9' + '\n' +
          'm1/n2/m3/n4/m5/n6/m7/n8/m9' + '\n',
      },
      tgt: {
        lengths: [5, 31, 57, 83, 109, 135, 161, 187, 213, 239],
      }
    },
    monoWithGapColumn: {
      src: {
        csv: 'seq' + '\n' +
          'm1/m2/m3/m4/m5/m6//m8/m9' + '\n' +
          'n1/m2/n3/m4/n5/m6//m8/n9' + '\n' +
          'm1/n2/m3/n4/m5///n8/m9' + '\n',
      },
      tgt: {lengths: [5, 31, 57, 83, 109, 135, 161, 180, 206, 232],}
    },
  };

  for (const [testName, testData] of Object.entries(lengthsTests)) {
    test(`getCellMonomerLengths-${testName}`, async () => {
      const df: DG.DataFrame = DG.DataFrame.fromCsv(testData.src.csv);
      await grok.data.detectSemanticTypes(df);
      const seqCol = df.getCol('seq');

      const monLengthLimit: number = 3;
      const charWidth: number = 7;
      const sepWidth: number = 12;
      const colTemp = new MonomerPlacer(null, seqCol, _package.logger, monLengthLimit, () => {
        return {
          fontCharWidth: charWidth,
          font: '12px monospace',
          separatorWidth: sepWidth,
          monomerToShort: monomerToShort,
        };
      });
      await colTemp.init();
      const resLengths = colTemp.getCellMonomerLengths(0, 1000)[1];
      expectArray(resLengths, testData.tgt.lengths);
    });
  }
});
