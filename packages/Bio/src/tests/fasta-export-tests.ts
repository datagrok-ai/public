import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {saveAsFastaDo, wrapSequence} from '../utils/save-as-fasta';
import {splitterAsFasta} from '@datagrok-libraries/bio';

type SaveAsFastaTestArgs = { srcCsv: string, idCols: string [], seqCol: string, lineWidth: number, tgtFasta: string };

category('fastaExport', () => {
  enum WrapDataTest {
    single = 'single',
    multi = 'multi'
  }

  const wrapData: { [key: string]: { src: string, tgt: string[] } } = {
    [WrapDataTest.single]: {
      src: 'MDYKETLLMPKTDFPMRGGLP',
      tgt: ['MDYKETLLMP', 'KTDFPMRGGL', 'P'],
    },
    [WrapDataTest.multi]: {
      src: 'M[MeI]YKETLL[MeF]PKTDFPMRGGL[MeA]',
      tgt: ['M[MeI]YKETLL[MeF]P', 'KTDFPMRGGL', '[MeA]'],
    },
  };

  enum SaveAsFastaTests {
    test1 = 'test1',
    test2 = 'test2'
  }

  const saveAsFastaData: {
    [key: string]: SaveAsFastaTestArgs
  } = {
    [SaveAsFastaTests.test1]: {
      srcCsv: `id,seq
1,MDYKETLLMP
2,KTDFPMRGGL
3,P`,
      idCols: ['id'],
      seqCol: 'seq',
      lineWidth: 10,
      tgtFasta: `>1
MDYKETLLMP
>2
KTDFPMRGGL
>3
P
`
    },
    [SaveAsFastaTests.test2]: {
      srcCsv: `id,id2,seq
1,seqA,M[MeI]YKETLL[MeF]P
2,seqB,KTDFPMRGGL
3,seqC,[MeA]
`,
      idCols: ['id2', 'id'],
      seqCol: 'seq',
      lineWidth: 5,
      tgtFasta: `>seqA|1
M[MeI]YKE
TLL[MeF]P
>seqB|2
KTDFP
MRGGL
>seqC|3
[MeA]
`
    }
  };

  test('wrapSequenceSingle', async () => {
    _testWrapSequence(WrapDataTest.single, 10);
  });

  test('wrapSequenceMulti', async () => {
    _testWrapSequence(WrapDataTest.multi, 10);
  });

  test('saveAsFastaTest1', async () => {
    _testSaveAsFasta(saveAsFastaData[SaveAsFastaTests.test1]);
  });

  test('saveAsFastaTest2', async () => {
    _testSaveAsFasta(saveAsFastaData[SaveAsFastaTests.test2]);
  });

  function _testWrapSequence(testKey: string, lineWidth: number = 10) {
    const splitter = splitterAsFasta;

    const srcSeq: string = wrapData[testKey].src;
    const wrapRes: string[] = wrapSequence(srcSeq, splitter, lineWidth);
    const wrapTgt: string[] = wrapData[testKey].tgt;

    expectArray(wrapRes, wrapTgt);
  }

  async function _testSaveAsFasta(args: SaveAsFastaTestArgs) {
    const df: DG.DataFrame = DG.DataFrame.fromCsv(args.srcCsv);

    const seqCol: DG.Column = df.getCol(args.seqCol);
    const idCols: DG.Column[] = args.idCols.map((colName) => df.getCol(colName));

    const fastaRes: string = saveAsFastaDo(idCols, seqCol, args.lineWidth);
    expect(fastaRes, args.tgtFasta);
  }
});

