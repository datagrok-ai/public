import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {ISeqHandler} from '@datagrok-libraries/bio/src/utils/macromolecule/seq-handler';
import {ISeqHelper, getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {MmDistanceFunctionsNames, mmDistanceFunctions}
  from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

category('Distance', async () => {
  let seqHelper: ISeqHelper;

  before(async () => {
    seqHelper = await getSeqHelper();
  });

  const scoringMatrix = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
  ];

  const alphabetIndexes = {'F': 0, 'W': 1, 'R': 2, 'Y': 3};

  const prot1 = 'FWRWY';
  const prot2 = 'FWRWW';

  const prot3 = 'FWY';
  const prot4 = 'FWRWY';

  const prot5 = 'FWY';
  const prot6 = 'FWRRRRY';

  const protTable = `seq
FWRWYVKHPFWRWYVKHP
YNRWYVKHPYNRWYVKHP
MWRSWYCKHPMWRSWYCKHP`;

  const DNATable = `seq
ATAACGATAACG
ATCGAATCGA
ATCGAATCGA`;

  const MSATable = `seq
ATAACATAAC
ATCGAATCGA
ATCGAATCGA`;

  test('protein-distance-function', async () => {
    const sh = await _initMacromoleculeColumn(protTable, seqHelper);
    const distFunc = sh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.LEVENSHTEIN);
  });

  test('DNA-distance-function', async () => {
    const sh = await _initMacromoleculeColumn(DNATable, seqHelper);
    const distFunc = sh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.LEVENSHTEIN);
  });

  test('MSA-distance-function', async () => {
    const sh = await _initMacromoleculeColumn(MSATable, seqHelper);
    const distFunc = sh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.HAMMING);
  });

  test('levenstein-sub', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN]();
    _testDistance(prot1, prot2, df, 0.2);
  });
  test('levenstein-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN]();
    _testDistance(prot3, prot4, df, 0.4);
  });

  test('hamming', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.HAMMING]();
    _testDistance(prot3, prot4, df, 0.6);
  });

  // Note that here the result is actually an inverted value of alignment score, which is coorelated with distance
  // tests using default BLOSUM62 matrix are in agreement with the results of the online tool
  test('needleman-blosum62', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]({gapOpen: 8, gapExtend: 2});
    _testDistance(prot1, prot2, df, -6);
  });

  test('needleman-blosum62-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]({gapOpen: 8, gapExtend: 2});
    _testDistance(prot3, prot4, df, -3.667);
  });

  test('needleman-custom-sub', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1},
    );
    _testDistance(prot1, prot2, df, 0.2);
  });

  test('needleman-custom-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1},
    );
    _testDistance(prot3, prot4, df, 0.667);
  });

  test('needleman-custom-zero-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 0},
    );
    _testDistance(prot5, prot6, df, 0.333);
  });

  test('needleman-custom-half-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 2, gapExtend: 1},
    );
    _testDistance(prot5, prot6, df, 1.667);
  });

  test('needleman-custom-same-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1},
    );
    if (DG.Test.isInBenchmark) {
      // max length 10000 for Needleman-Wunsch
      const seq1 = /* length 3000 */ Array(1000).fill(prot5).join('');
      const seq2 = /* length 7000 */Array(1000).fill(prot6).join('');
      _testDistance(seq1, seq2, df, 1.333);
    } else
      _testDistance(prot5, prot6, df, 1.333);
  }, {benchmark: true});
});

async function _initMacromoleculeColumn(csv: string, seqHelper: ISeqHelper): Promise<ISeqHandler> {
  const srcDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  const seqCol = srcDf.col('seq')!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: seqCol}) as unknown as string;
  if (semType)
    seqCol.semType = semType;
  await grok.data.detectSemanticTypes(srcDf);
  const sh = seqHelper.getSeqHandler(seqCol);
  return sh;
}

function _testDistance(seq1: string, seq2: string, df: (a: string, b: string) => number, expected: number) {
  const d = df(seq1, seq2);
  expect(Number(d.toFixed(3)), Number(expected.toFixed(3)));
}

export function mapToFixed(ar: Float32Array | number[]) {
  return Array.from(ar).map((d) => Number(d.toFixed(3)));
}
