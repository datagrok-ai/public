import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('detectors', () => {
  //   const csvDf1: string = `col1
  // 1
  // 2
  // 3`;
  //
  //   const csvDf2: string = `col1
  // 4
  // 5
  // 6
  // 7`;
  //
  //   const csvDf3: string = `col1
  // 8
  // 9
  // 10
  // 11
  // 12`;

  const csvDfN1: string = `seq
ACGTC
CAGTGT
TTCAAC
`;

  /** Pure amino acids sequence */
  const csvDfAA1: string = `seq
FWPHEY
YNRQWYV
MKPSEYV
`;

  const csvDfSepNt: string = `seq
A*C*G*T*C
C*A*G*T*G*T
T*T*C*A*A*C
`;

  const csvDfSepPt: string = `seq
F-W-P-H-E-Y
Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V
`;

  const csvDfSepUn1: string = `seq
abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg
rut12-rty-her2-abc-cfr3-wert-rut12
`;

  const csvDfSepUn2: string = `seq
abc/dfgg/abc1/cfr3/rty/wert
rut12/her2/rty/wert//abc/abc1/dfgg
rut12/rty/her2/abc/cfr3//wert/rut12
`;

  const csvDfSepMsaN1: string = `seq
A-C--G-T--C-T
C-A-C--T--G-T
A-C-C-G-T-A-C-T
`;

  const csvDfMsaN1: string = `seq
AC-GT-CT
CAC-T-GT
ACCGTACT
`;

  const csvDfMsaAA1: string = `seq
FWR-WYV-KHP
YNR-WYV-KHP
MWRSWY-CKHP
`;

  // test('testDetectors1', async () => { await _testDetectors(csvDf1); });
  // test('testDetectors2', async () => { await _testDetectors(csvDf2); });
  // test('testDetectors3', async () => { await _testDetectors(csvDf3); });

  test('testDetectorsN1', async () => { await _testDetectorsN1(csvDfN1); });
  test('testDetectorsAA1', async () => { await _testDetectorsAA1(csvDfAA1); });
  test('testDetectorsMsaN1', async () => { await _testDetectorsMsaN1(csvDfMsaN1); });
  test('testDetectorsMsaAA1', async () => { await _testDetectorsMsaAA1(csvDfMsaAA1); });

  test('testDetectorsSepUn1', async () => { await _testDetectorsSepUn1(csvDfSepUn1); });
  test('testDetectorsSepUn2', async () => { await _testDetectorsSepUn2(csvDfSepUn2); });

  test('testDetectorsSepMsaN1', async () => { await _testDetectorsSepMsaN1(csvDfSepMsaN1); });
});

// export async function _testDetectors(csvDf: string) {
//   const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf);
//   await grok.data.detectSemanticTypes(df);
//
//   const col1: DG.Column = df.col('col1')!;
//   expect(col1.semType, null);
//   expect(col1.getTag(DG.TAGS.UNITS), null);
// }

export async function _testDetectorsN1(csvDfN1: string) {
  const dfN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1);
  await grok.data.detectSemanticTypes(dfN1);

  const col: DG.Column = dfN1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:NT');
}

export async function _testDetectorsAA1(csvDfAA1: string) {
  const dfAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfAA1);
  await grok.data.detectSemanticTypes(dfAA1);

  const col: DG.Column = dfAA1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
}

export async function _testDetectorsMsaN1(csvDfMsaN1: string) {
  const dfMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaN1);
  await grok.data.detectSemanticTypes(dfMsaN1);

  const col: DG.Column = dfMsaN1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:NT');
}

export async function _testDetectorsMsaAA1(csvDfMsaAA1: string) {
  const dfMsaAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaAA1);
  await grok.data.detectSemanticTypes(dfMsaAA1);

  const col: DG.Column = dfMsaAA1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:PT');
}

export async function _testDetectorsSepUn1(csvDfSepUn1: string) {
  const dfSepUn1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfSepUn1);
  await grok.data.detectSemanticTypes(dfSepUn1);

  const col: DG.Column = dfSepUn1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:UN');
}

export async function _testDetectorsSepUn2(csvDfSepUn2: string) {
  const dfSepUn2: DG.DataFrame = DG.DataFrame.fromCsv(csvDfSepUn2);
  await grok.data.detectSemanticTypes(dfSepUn2);

  const col: DG.Column = dfSepUn2.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:UN');
}

export async function _testDetectorsSepMsaN1(csvDfSepMsaN1: string) {
  const dfSepMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfSepMsaN1);
  await grok.data.detectSemanticTypes(dfSepMsaN1);

  const col: DG.Column = dfSepMsaN1.col('seq')!;
  expect(col.semType, 'MACROMOLECULE');
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ.MSA:NT');
}
