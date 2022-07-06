import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {mmSemType} from '../const';
import {importFasta} from '../package';

type DfReaderFunc = () => Promise<DG.DataFrame>;

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

  const enum Samples {
    peptidesComplex = 'PeptidesComplex',
    fastaCsv = 'FastaCsv',
    msaComplex = 'MsaComplex',
    idCsv = 'IdCsv',
  }

  const samples: { [key: string]: string } = {
    'PeptidesComplex': 'System:AppData/Bio/samples/peptides_complex_msa.csv',
    'FastaCsv': 'System:AppData/Bio/samples/sample_FASTA.csv',
    'MsaComplex': 'System:AppData/Bio/samples/sample_MSA.csv',
    'IdCsv': 'System:AppData/Bio/samples/id.csv',
  };

  const _samplesDfs: { [key: string]: Promise<DG.DataFrame> } = {};
  const readSamplesCsv: (key: string) => DfReaderFunc = (key: string) => {
    return async () => {
      if (!(key in _samplesDfs)) {
        _samplesDfs[key] = (async (): Promise<DG.DataFrame> => {
          const csv: string = await grok.dapi.files.readAsText(samples[key]);
          const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
          await grok.data.detectSemanticTypes(df);
          return df;
        })();
      }
      return _samplesDfs[key];
    };
  };

  const _csvDfs: { [key: string]: Promise<DG.DataFrame> } = {};
  const readCsv: (key: string, csv: string) => DfReaderFunc = (key: string, csv: string) => {
    return async () => {
      if (!(key in _csvDfs)) {
        _csvDfs[key] = (async (): Promise<DG.DataFrame> => {
          const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
          await grok.data.detectSemanticTypes(df);
          return df;
        })();
      }
      return _csvDfs[key];
    };
  };

  test('Negative1', async () => { await _testNeg(readCsv('csvDf1', csvDf1), 'col1'); });
  test('Negative2', async () => { await _testNeg(readCsv('csvDf2', csvDf2), 'col1'); });
  test('Negative3', async () => { await _testNeg(readCsv('csvDf3', csvDf3), 'col1'); });
  test('NegativeSmiles', async () => { await _testNeg(readCsv('csvDfSmiles', csvDfSmiles), 'col1'); });

  test('N1', async () => { await _testN1(csvDfN1); });
  test('AA1', async () => { await _testAA1(csvDfAA1); });
  test('MsaN1', async () => { await _testMsaN1(csvDfMsaN1); });
  test('MsaAA1', async () => { await _testMsaAA1(csvDfMsaAA1); });

  test('SepNt', async () => { await _testSepNt(csvDfSepNt, '*'); });
  test('SepPt', async () => { await _testSepPt(csvDfSepPt, '-'); });
  test('SepUn1', async () => { await _testSepUn(csvDfSepUn1, '-'); });
  test('SepUn2', async () => { await _testSepUn(csvDfSepUn2, '/'); });

  test('SepMsaN1', async () => { await _testSepMsaN1(csvDfSepMsaN1); });

  test('SamplesFastaCsvPt', async () => {
    await _testSamplesFastaCsvPt();
  });
  test('SamplesFastaCsvNegativeEntry', async () => {
    await _testNeg(readSamplesCsv(Samples.fastaCsv), 'Entry');
  });
  test('SamplesFastaCsvNegativeLength', async () => {
    await _testNeg(readSamplesCsv(Samples.fastaCsv), 'Length');
  });
  test('SamplesFastaCsvNegativeUniProtKB', async () => {
    await _testNeg(readSamplesCsv(Samples.fastaCsv), 'UniProtKB');
  });

  test('SamplesFastaFastaPt', async () => { await _testSamplesFastaFastaPt(); });

  // System:AppData/Bio/samples/peptides_complex_align.csv contains monomers with spaces
  // test('SamplesPeptidesComplexUn', async () => {
  //   await _testSamplesPeptidesComplexUn();
  // });

  test('samplesPeptidesComplexNegativeID', async () => {
    await _testNeg(readSamplesCsv(Samples.peptidesComplex), 'ID');
  });
  test('SamplesPeptidesComplexNegativeMeasured', async () => {
    await _testNeg(readSamplesCsv(Samples.peptidesComplex), 'Measured');
  });
  test('SamplesPeptidesComplexNegativeValue', async () => {
    await _testNeg(readSamplesCsv(Samples.peptidesComplex), 'Value');
  });

  test('samplesMsaComplexUn', async () => {
    await _testPos(readSamplesCsv(Samples.msaComplex), 'MSA', 'separator:SEQ.MSA:UN', '/');
  });
  test('samplesMsaComplexNegativeActivity', async () => {
    await _testNeg(readSamplesCsv(Samples.msaComplex), 'Activity');
  });

  test('samplesIdCsvNegativeID', async () => {
    await _testNeg(readSamplesCsv(Samples.idCsv), 'ID');
  });
});

export async function _testNeg(readDf: DfReaderFunc, colName: string) {
  const df: DG.DataFrame = await readDf();

  const col: DG.Column = df.col(colName)!;
  expect(col.semType === mmSemType, false);
}

export async function _testPos(readDf: DfReaderFunc, colName: string, units: string, separator: string) {
  const df: DG.DataFrame = await readDf();

  const col: DG.Column = df.col(colName)!;
  expect(col.semType === mmSemType, true);
  expect(col.getTag(DG.TAGS.UNITS), units);
  if (separator)
    expect(col.getTag('separator'), separator);
}

export async function _testN1(csvDfN1: string) {
  const dfN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfN1);
  await grok.data.detectSemanticTypes(dfN1);

  const col: DG.Column = dfN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:NT');
}

export async function _testAA1(csvDfAA1: string) {
  const dfAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfAA1);
  await grok.data.detectSemanticTypes(dfAA1);

  const col: DG.Column = dfAA1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
}

export async function _testMsaN1(csvDfMsaN1: string) {
  const dfMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaN1);
  await grok.data.detectSemanticTypes(dfMsaN1);

  const col: DG.Column = dfMsaN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:NT');
}

export async function _testMsaAA1(csvDfMsaAA1: string) {
  const dfMsaAA1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfMsaAA1);
  await grok.data.detectSemanticTypes(dfMsaAA1);

  const col: DG.Column = dfMsaAA1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ.MSA:PT');
}

export async function _testSepNt(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:NT');
  expect(col.getTag('separator'), separator);
}

export async function _testSepPt(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:PT');
  expect(col.getTag('separator'), separator);
}

export async function _testSepUn(csv: string, separator: string) {
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ:UN');
  expect(col.getTag('separator'), separator);
}

export async function _testSepMsaN1(csvDfSepMsaN1: string) {
  const dfSepMsaN1: DG.DataFrame = DG.DataFrame.fromCsv(csvDfSepMsaN1);
  await grok.data.detectSemanticTypes(dfSepMsaN1);

  const col: DG.Column = dfSepMsaN1.col('seq')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ.MSA:NT');
}

export async function _testSamplesFastaCsvPt() {
  const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.csv');
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('sequence')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
  expect(col.getTag('separator'), null);
}

export async function _testSamplesFastaFastaPt() {
  const fasta: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/sample_FASTA.fasta');
  const df: DG.DataFrame = importFasta(fasta)[0];

  const col: DG.Column = df.col('sequence')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'fasta:SEQ:PT');
  expect(col.getTag('separator'), null);
}

export async function _testSamplesPeptidesComplexUn() {
  const csv: string = await grok.dapi.files.readAsText('System:AppData/Bio/samples/peptides_complex_aligned.csv');
  const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
  await grok.data.detectSemanticTypes(df);

  const col: DG.Column = df.col('AlignedSequence')!;
  expect(col.semType, mmSemType);
  expect(col.getTag(DG.TAGS.UNITS), 'separator:SEQ.MSA:UN');
  expect(col.getTag('separator'), '-');
}

