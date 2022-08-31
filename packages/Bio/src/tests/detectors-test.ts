import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {importFasta} from '../package';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

type DfReaderFunc = () => Promise<DG.DataFrame>;

category('detectors', () => {
  const csvDf1: string = `col1
1
2
3`;

  const csvDfEmpty: string = `id,col1
1,
2,
3,
4,
5,`;

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

  const csvDfDna1: string = `seq
ACGTC
CAGTGT
TTCAAC
`;

  const csvDfRna1: string = `seq
ACGUC
CAGUGU
UUCAAC
`;

  /** Pure amino acids sequence */
  const csvDfPt1: string = `seq
FWPHEY
YNRQWYV
MKPSEYV
`;

  const csvDfSepDna: string = `seq
A*C*G*T*C
C*A*G*T*G*T
T*T*C*A*A*C
`;

  const csvDfSepRna: string = `seq
A*C*G*U*C
C*A*G*U*G*U
U*U*C*A*A*C
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

  const csvDfSepMsaDna1: string = `seq
A-C--G-T--C-T
C-A-C--T--G-T
A-C-C-G-T-A-C-T
`;

  const csvDfMsaDna1: string = `seq
AC-GT-CT
CAC-T-GT
ACCGTACT
`;

  const csvDfMsaPt1: string = `seq
FWR-WYV-KHP
YNR-WYV-KHP
MWRSWY-CKHP
`;

  const enum Samples {
    peptidesComplex = 'peptidesComplex',
    peptidesSimple = 'peptidesSimple',
    fastaCsv = 'fastaCsv',
    fastaFasta = 'fastaFasta',
    fastaPtCsv = 'fastaPtCsv',
    msaComplex = 'msaComplex',
    helmCsv = 'helmCsv',
    testDemogCsv = 'testDemogCsv',
    testHelmCsv = 'testHelmCsv',
    testIdCsv = 'testIdCsv',
    testSmilesCsv = 'testSmilesCsv',
    testSmiles2Csv = 'testSmiles2Csv',
    testCerealCsv = 'testCerealCsv',
    testActivityCliffsCsv = 'testActivityCliffsCsv',
    testSpgi100 = 'testSpgi100',
    testUnichemSources = 'testUnichemSources',
    testDmvOffices = 'testDmvOffices',
    testAlertCollection = 'testAlertCollection',
  }

  const samples: { [key: string]: string } = {
    'fastaCsv': 'System:AppData/Bio/samples/sample_FASTA.csv',
    'fastaFasta': 'System:AppData/Bio/samples/sample_FASTA.fasta',
    'fastaPtCsv': 'System:AppData/Bio/samples/sample_FASTA_PT.csv',
    'msaComplex': 'System:AppData/Bio/samples/sample_MSA.csv',
    'helmCsv': 'System:AppData/Bio/samples/sample_HELM.csv',
    'peptidesComplex': 'System:AppData/Bio/tests/peptides_complex_msa.csv',
    'peptidesSimple': 'System:AppData/Bio/tests/peptides_simple_msa.csv',
    'testDemogCsv': 'System:AppData/Bio/tests/testDemog.csv',
    'testHelmCsv': 'System:AppData/Bio/tests/testHelm.csv',
    'testIdCsv': 'System:AppData/Bio/tests/testId.csv',
    'testSmilesCsv': 'System:AppData/Bio/tests/testSmiles.csv',
    'testSmiles2Csv': 'System:AppData/Bio/tests/testSmiles2.csv',
    'testActivityCliffsCsv': 'System:AppData/Bio/tests/testActivityCliffs.csv', // smiles
    'testCerealCsv': 'System:AppData/Bio/tests/testCereal.csv',
    'testSpgi100': 'System:AppData/Bio/tests/testSpgi100.csv',
    'testUnichemSources': 'System:AppData/Bio/tests/testUnichemSources.csv',
    'testDmvOffices': 'System:AppData/Bio/tests/testDmvOffices.csv',
    'testAlertCollection': 'System:AppData/Bio/tests/testAlertCollection.csv',
  };

  const _samplesDfs: { [key: string]: Promise<DG.DataFrame> } = {};

  function readSamples(key: string, readFile: (file: string) => Promise<DG.DataFrame> = readFileCsv): DfReaderFunc {
    return async () => {
      if (!(key in _samplesDfs)) {
        _samplesDfs[key] = (async (): Promise<DG.DataFrame> => {
          const df: DG.DataFrame = await readFile(samples[key]);
          // await grok.data.detectSemanticTypes(df);
          return df;
        })();
      }
      return _samplesDfs[key];
    };
  };

  async function readFileCsv(file: string): Promise<DG.DataFrame> {
    const csv: string = await grok.dapi.files.readAsText(file);
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    return df;
  }

  async function readFileFasta(file: string): Promise<DG.DataFrame> {
    const txt: string = await grok.dapi.files.readAsText(file);
    const df: DG.DataFrame = importFasta(txt)[0];
    return df;
  }

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


  test('NegativeEmpty', async () => { await _testNeg(readCsv('csvDfEmpty', csvDfEmpty), 'col1'); });
  test('Negative1', async () => { await _testNeg(readCsv('csvDf1', csvDf1), 'col1'); });
  test('Negative2', async () => { await _testNeg(readCsv('csvDf2', csvDf2), 'col1'); });
  test('Negative3', async () => { await _testNeg(readCsv('csvDf3', csvDf3), 'col1'); });
  test('NegativeSmiles', async () => { await _testNeg(readCsv('csvDfSmiles', csvDfSmiles), 'col1'); });

  test('Dna1', async () => {
    await _testPos(readCsv('csvDfDna1', csvDfDna1), 'seq', 'fasta', 'SEQ', 'DNA');
  });
  test('Rna1', async () => {
    await _testPos(readCsv('csvDfRna1', csvDfRna1), 'seq', 'fasta', 'SEQ', 'RNA');
  });
  test('AA1', async () => {
    await _testPos(readCsv('csvDfPt1', csvDfPt1), 'seq', 'fasta', 'SEQ', 'PT');
  });
  test('MsaDna1', async () => {
    await _testPos(readCsv('csvDfMsaDna1', csvDfMsaDna1), 'seq', 'fasta', 'SEQ.MSA', 'DNA');
  });

  test('MsaAA1', async () => {
    await _testPos(readCsv('csvDfMsaPt1', csvDfMsaPt1), 'seq', 'fasta', 'SEQ.MSA', 'PT');
  });

  test('SepDna', async () => {
    await _testPos(readCsv('csvDfSepDna', csvDfSepDna), 'seq', 'separator', 'SEQ', 'DNA', '*');
  });
  test('SepRna', async () => {
    await _testPos(readCsv('csvDfSepRna', csvDfSepRna), 'seq', 'separator', 'SEQ', 'RNA', '*');
  });
  test('SepPt', async () => {
    await _testPos(readCsv('csvDfSepPt', csvDfSepPt), 'seq', 'separator', 'SEQ', 'PT', '-');
  });
  test('SepUn1', async () => {
    await _testPos(readCsv('csvDfSepUn1', csvDfSepUn1), 'seq', 'separator', 'SEQ', 'UN', '-');
  });
  test('SepUn2', async () => {
    await _testPos(readCsv('csvDfSepUn2', csvDfSepUn2), 'seq', 'separator', 'SEQ', 'UN', '/');
  });

  test('SepMsaN1', async () => {
    await _testPos(readCsv('csvDfSepMsaDna1', csvDfSepMsaDna1), 'seq', 'separator', 'SEQ.MSA', 'DNA', '-');
  });

  test('SamplesFastaCsvPt', async () => {
    await _testPos(readSamples(Samples.fastaCsv), 'sequence', 'fasta', 'SEQ', 'PT');
  });
  test('SamplesFastaCsvNegativeEntry', async () => {
    await _testNeg(readSamples(Samples.fastaCsv), 'Entry');
  });
  test('SamplesFastaCsvNegativeLength', async () => {
    await _testNeg(readSamples(Samples.fastaCsv), 'Length');
  });
  test('SamplesFastaCsvNegativeUniProtKB', async () => {
    await _testNeg(readSamples(Samples.fastaCsv), 'UniProtKB');
  });

  test('SamplesFastaFastaPt', async () => {
    await _testPos(readSamples(Samples.fastaFasta, readFileFasta), 'sequence', 'fasta', 'SEQ', 'PT');
  });

  // peptidesComplex contains monomers with spaces in AlignedSequence columns, which are forbidden
  // test('samplesPeptidesComplexPositiveAlignedSequence', async () => {
  //   await _testPos(readSamples(Samples.peptidesComplex), 'AlignedSequence', 'separator:SEQ:UN', '-');
  // });
  test('samplesPeptidesComplexNegativeID', async () => {
    await _testNeg(readSamples(Samples.peptidesComplex), 'ID');
  });
  test('SamplesPeptidesComplexNegativeMeasured', async () => {
    await _testNeg(readSamples(Samples.peptidesComplex), 'Measured');
  });
  test('SamplesPeptidesComplexNegativeValue', async () => {
    await _testNeg(readSamples(Samples.peptidesComplex), 'Value');
  });

  test('samplesMsaComplexUn', async () => {
    await _testPos(readSamples(Samples.msaComplex), 'MSA', 'separator', 'SEQ.MSA', 'UN', '/');
  });
  test('samplesMsaComplexNegativeActivity', async () => {
    await _testNeg(readSamples(Samples.msaComplex), 'Activity');
  });

  test('samplesIdCsvNegativeID', async () => {
    await _testNeg(readSamples(Samples.testIdCsv), 'ID');
  });

  test('samplesSarSmallCsvNegativeSmiles', async () => {
    await _testNeg(readSamples(Samples.testSmilesCsv), 'smiles');
  });

  test('samplesHelmCsvHELM', async () => {
    await _testPos(readSamples(Samples.helmCsv), 'HELM', 'helm', null, null, null);
  });

  test('samplesHelmCsvNegativeActivity', async () => {
    await _testNeg(readSamples(Samples.helmCsv), 'Activity');
  });

  // sample_testHelm.csb
  // columns: ID,Test type,HELM string,Valid?,Mol Weight,Mol Formula,SMILES
  test('samplesTestHelmNegativeID', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'ID');
  });
  test('samplesTestHelmNegativeTestType', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'Test type');
  });
  test('samplesTestHelmPositiveHelmString', async () => {
    await _testPos(readSamples(Samples.testHelmCsv), 'HELM string', 'helm', null, null, null);
  });
  test('samplesTestHelmNegativeValid', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'Valid?');
  });
  test('samplesTestHelmNegativeMolWeight', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'Mol Weight');
  });
  test('samplesTestHelmNegativeMolFormula', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'Mol Formula');
  });
  test('samplesTestHelmNegativeSmiles', async () => {
    await _testNeg(readSamples(Samples.testHelmCsv), 'Smiles');
  });

  test('samplesTestDemogNegativeAll', async () => {
    const dfFunc: DfReaderFunc = readSamples(Samples.testDemogCsv);
    const df: DG.DataFrame = await dfFunc();

    for (const col of df.columns.toList())
      await _testNeg(dfFunc, col.name);
  });

  test('samplesTestSmiles2NegativeSmiles', async () => {
    await _testNeg(readSamples(Samples.testSmiles2Csv), 'SMILES');
  });

  test('samplesTestActivityCliffsNegativeSmiles', async () => {
    await _testNeg(readSamples(Samples.testActivityCliffsCsv), 'smiles');
  });

  test('samplesFastaPtPosSequence', async () => {
    await _testPos(readSamples(Samples.fastaPtCsv), 'sequence', 'fasta', 'SEQ', 'PT');
  });

  test('samplesTestCerealNegativeCerealName', async () => {
    await _testNeg(readSamples(Samples.testCerealCsv), 'cereal_name');
  });

  test('samplesTestSpgi100NegativeStereoCategory', async () => {
    await _testNeg(readSamples(Samples.testSpgi100), 'Stereo Category');
  });
  test('samplesTestSpgi100NegativeScaffoldNames', async () => {
    await _testNeg(readSamples(Samples.testSpgi100), 'Scaffold Names');
  });
  test('samplesTestSpgi100NegativePrimaryScaffoldName', async () => {
    await _testNeg(readSamples(Samples.testSpgi100), 'Primary Scaffold Name');
  });
  test('samplesTestSpgi100NegativeSampleName', async () => {
    await _testNeg(readSamples(Samples.testSpgi100), 'Sample Name');
  });

  test('samplesTestUnichemSourcesNegativeSrcUrl', async () => {
    await _testNeg(readSamples(Samples.testUnichemSources), 'src_url');
  });
  test('samplesTestUnichemSourcesNegativeBaseIdUrl', async () => {
    await _testNeg(readSamples(Samples.testUnichemSources), 'base_id_url');
  });

  test('samplesTestDmvOfficesNegativeOfficeName', async () => {
    await _testNeg(readSamples(Samples.testDmvOffices), 'Office Name');
  });
  test('samplesTestDmvOfficesNegativeCity', async () => {
    await _testNeg(readSamples(Samples.testDmvOffices), 'City');
  });

  test('samplesTestAlertCollectionNegativeSmarts', async () => {
    await _testNeg(readSamples(Samples.testAlertCollection), 'smarts');
  });
});

export async function _testNeg(readDf: DfReaderFunc, colName: string) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.col(colName)!;
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  if (col.semType === DG.SEMTYPE.MACROMOLECULE) {
    const msg = `Negative test detected semType='${col.semType}', units='${col.getTag(DG.TAGS.UNITS)}'.`;
    throw new Error(msg);
    // col.semType = '';
    // col.setTag(DG.TAGS.UNITS, '');
    // col.setTag('separator', '');
  }
}

export async function _testPos(
  readDf: DfReaderFunc, colName: string, units: string,
  aligned: string | null, alphabet: string | null, separator: string | null = null
) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.col(colName)!;
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  expect(col.semType === DG.SEMTYPE.MACROMOLECULE, true);
  expect(col.getTag(DG.TAGS.UNITS), units);
  expect(col.getTag('aligned'), aligned);
  expect(col.getTag('alphabet'), alphabet);
  if (separator)
    expect(col.getTag('separator'), separator);

  const uh = new UnitsHandler(col);
  if (!uh.isHelm()) {
    expect(uh.aligned, aligned);
    expect(uh.alphabet, alphabet);
  }
}

