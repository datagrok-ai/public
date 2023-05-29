import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/utils/src/test';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {MmDistanceFunctionsNames, mmDistanceFunctions}
  from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

category('Distance', async () => {
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
      FWRWYVKHP
      YNRWYVKHP
      MWRSWYCKHP`;

  const DNATable = `seq
      ATAACG
      ATCGA
      ATCGA`;

  const MSATable = `seq
      ATAAC
      ATCGA
      ATCGA`;
  test('protein-distance-function', async () => {
    const uh = await _initMacromoleculeColumn(protTable);
    const distFunc = uh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH);
  });

  test('DNA-distance-function', async () => {
    const uh = await _initMacromoleculeColumn(DNATable);
    const distFunc = uh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.LEVENSHTEIN);
  });

  test('MSA-distance-function', async () => {
    const uh = await _initMacromoleculeColumn(MSATable);
    const distFunc = uh.getDistanceFunctionName();
    expect(distFunc, MmDistanceFunctionsNames.HAMMING);
  });

  test('levenstein-sub', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN]();
    _testDistance(prot1, prot2, df, 1);
  });
  test('levenstein-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN]();
    _testDistance(prot3, prot4, df, 2);
  });

  test('hamming', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.HAMMING]();
    _testDistance(prot3, prot4, df, 3);
  });

  // Note that here the result is actually an inverted value of alignment score, which is coorelated with distance
  // tests using default BLOSUM62 matrix are in agreement with the results of the online tool
  test('needleman-blosum62', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]();
    _testDistance(prot1, prot2, df, -35);
  });

  test('needleman-blosum62-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]();
    _testDistance(prot3, prot4, df, -14);
  });

  test('needleman-custom-sub', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1}
    );
    _testDistance(prot1, prot2, df, -4);
  });

  test('needleman-custom-del', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1}
    );
    _testDistance(prot3, prot4, df, -1);
  });

  test('needleman-custom-zero-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 0}
    );
    _testDistance(prot5, prot6, df, -2);
  });

  test('needleman-custom-half-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 2, gapExtend: 1}
    );
    _testDistance(prot5, prot6, df, 2);
  });

  test('needleman-custom-same-extend', async () => {
    const df = mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH](
      {scoringMatrix, alphabetIndexes, gapOpen: 1, gapExtend: 1}
    );
    if (DG.Test.isInBenchmark) {
      const seq1 = Array(10000).fill('FWRY').join('');
      const seq2 = Array(10000).fill('FYWRRY').join('');
      _testDistance(seq1, seq2, df, -20000);
    } else _testDistance(prot5, prot6, df, 1);
  });
});

async function _initMacromoleculeColumn(csv: string): Promise<UnitsHandler> {
  const srcDf: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  const seqCol = srcDf.col('seq')!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: seqCol}) as unknown as string;
  if (semType)
    seqCol.semType = semType;
  await grok.data.detectSemanticTypes(srcDf);
  const uh = new UnitsHandler(seqCol);
  return uh;
}

function _testDistance(seq1: string, seq2: string, df: (a: string, b: string) => number, expected: number) {
  const d = df(seq1, seq2);
  expect(d, expected);
}
