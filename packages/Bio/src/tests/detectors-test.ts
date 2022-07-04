import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {mmSemType} from '../const';

category('detectors', () => {
  const csvDf1: string = `col1
1
2
3`;

  const csvDf2: string = `col1
4
5
6
7`;

  const csvDf3: string = `col1
8
9
10
11
12`;

  const csvDfSmiles: string = `col1
CCCCN1C(=O)CN=C(c2cc(F)ccc12)C3CCCCC3
C1CCCCC1
CCCCCC
`;

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

  test('testDetectorsNegative1', async () => { await _testDetectorsNegative(csvDf1); });
  test('testDetectorsNegative2', async () => { await _testDetectorsNegative(csvDf2); });
  test('testDetectorsNegative3', async () => { await _testDetectorsNegative(csvDf3); });
  test('testDetectorsNegativeSmiles', async () => { await _testDetectorsNegative(csvDfSmiles); });

  test('testDetectorsN1', async () => { await _testDetectorsN1(csvDfN1); });
  test('testDetectorsAA1', async () => { await _testDetectorsAA1(csvDfAA1); });
  test('testDetectorsMsaN1', async () => { await _testDetectorsMsaN1(csvDfMsaN1); });
  test('testDetectorsMsaAA1', async () => { await _testDetectorsMsaAA1(csvDfMsaAA1); });

  test('testDetectorsSepNt', async () => { await _testDetectorsSepNt(csvDfSepNt, '*'); });
  test('testDetectorsSepPt', async () => { await _testDetectorsSepPt(csvDfSepPt, '-'); });
  test('testDetectorsSepUn1', async () => { await _testDetectorsSepUn(csvDfSepUn1, '-'); });
  test('testDetectorsSepUn2', async () => { await _testDetectorsSepUn(csvDfSepUn2, '/'); });

  test('testDetectorsSepMsaN1', async () => { await _testDetectorsSepMsaN1(csvDfSepMsaN1); });
});

export async function _testDetectorsNegative(csvDf: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csvDf);
  await grok.data.detectSemanticTypes(df);

  const col1: DG.Column = df.col('col1')!;
  expect(col1.semType == mmSemType, false);
}

export async function _testDetectorsN1(csvDfN1: string) {
  const dfN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1);
  await grok.data.detectSemanticTypes(dfN1);

  const col: DG.Column = dfN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:NT');
}

export async function _testDetectorsAA1(csvDfAA1: string) {
  const dfAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfAA1);
  await grok.data.detectSemanticTypes(dfAA1);

  const col: DG.Column = dfAA1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
}

export async function _testDetectorsMsaN1(csvDfMsaN1: string) {
  const dfMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaN1);
  await grok.data.detectSemanticTypes(dfMsaN1);

  const col: DG.Column = dfMsaN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:NT');
}

export async function _testDetectorsMsaAA1(csvDfMsaAA1: string) {
  const dfMsaAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaAA1);
  await grok.data.detectSemanticTypes(dfMsaAA1);

  const col: DG.Column = dfMsaAA1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:PT');
}

export async function _testDetectorsSepNt(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:NT');
  expect(col.getTag('separator'), separator);
}

export async function _testDetectorsSepPt(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:PT');
  expect(col.getTag('separator'), separator);
}

export async function _testDetectorsSepUn(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:UN');
  expect(col.getTag('separator'), separator);
}


export async function _testDetectorsSepMsaN1(csvDfSepMsaN1: string) {
  const dfSepMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfSepMsaN1);
  await grok.data.detectSemanticTypes(dfSepMsaN1);

  const col: DG.Column = dfSepMsaN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ.MSA:NT');
}
