import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectObject} from '@datagrok-libraries/utils/src/test';

import {importFasta} from '../package';
import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';

/*
// snippet to list df columns of semType='Macromolecule' (false positive)
const df = grok.shell.tableByName('SPGI');
for (let i = 0; i < df.columns.length; i++) {
  const col = df.columns.byIndex(i);
  if (col.semType == 'Macromolecule') {
  console.log( i + ' - ' + col.name + ' - ' + col.semType);
  }
}
 */

type DfReaderFunc = () => Promise<DG.DataFrame>;

category('detectors', () => {
  const enum csvTests {
    negEmpty = 'negEmpty',
    neg1 = 'neg1',
    neg2 = 'neg2',
    neg3 = 'neg3',
    negSmiles = 'negSmiles',
    fastaDna1 = 'csvFastaDna1',
    fastaRna1 = 'fastaRna1',
    fastaPt1 = 'fastaPt1',
    fastaUn = 'fastaUn',
    sepDna = 'sepDna',
    sepRna = 'sepRna',
    sepPt = 'sepPt',
    sepUn1 = 'sepUn1',
    sepUn2 = 'sepUn2',
    sepMsaDna1 = 'sepMsaDna1',
    fastaMsaDna1 = 'fastaMsaDna1',
    fastaMsaPt1 = 'fastaMsaPt1',
  }

  const csvData = new class {
    [csvTests.negEmpty]: string = `id,col1
1,
2,
3,
4,
5,`;
    [csvTests.neg1]: string = `col1
1
2
3`;
    [csvTests.neg2]: string = `col1
4
5
6
7`;
    [csvTests.neg3]: string = `col1
8
9
10
11
12`;
    [csvTests.negSmiles]: string = `col1
CCCCN1C(=O)CN=C(c2cc(F)ccc12)C3CCCCC3
C1CCCCC1
CCCCCC
`;
    [csvTests.fastaDna1]: string = `seq
ACGTC
CAGTGT
TTCAAC
`;
    [csvTests.fastaRna1]: string = `seq
ACGUC
CAGUGU
UUCAAC
`;
    /** Pure amino acids sequence */
    [csvTests.fastaPt1]: string = `seq
FWPHEY
YNRQWYV
MKPSEYV
`;
    [csvTests.fastaUn]: string = `seq
[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]D
[meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2]
[Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2]
`;
    [csvTests.sepDna]: string = `seq
A*C*G*T*C
C*A*G*T*G*T
T*T*C*A*A*C
`;
    [csvTests.sepRna]: string = `seq
A*C*G*U*C
C*A*G*U*G*U
U*U*C*A*A*C
`;
    [csvTests.sepPt]: string = `seq
F-W-P-H-E-Y
Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V
`;
    [csvTests.sepUn1]: string = `seq
abc-dfgg-abc1-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg
rut12-rty-her2-abc-cfr3-wert-rut12
`;
    [csvTests.sepUn2]: string = `seq
abc/dfgg/abc1/cfr3/rty/wert
rut12/her2/rty/wert//abc/abc1/dfgg
rut12/rty/her2/abc/cfr3//wert/rut12
`;
    [csvTests.sepMsaDna1]: string = `seq
A-C--G-T--C-T
C-A-C--T--G-T
A-C-C-G-T-A-C-T
`;
    [csvTests.fastaMsaDna1]: string = `seq
AC-GT-CT
CAC-T-GT
ACCGTACT
`;
    [csvTests.fastaMsaPt1]: string = `seq
FWR-WYV-KHP
YNR-WYV-KHP
MWRSWY-CKHP
`;
  }();

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
    testSmilesShort = 'testSmilesShort',
    testCerealCsv = 'testCerealCsv',
    testActivityCliffsCsv = 'testActivityCliffsCsv',
    testUnichemSources = 'testUnichemSources',
    testDmvOffices = 'testDmvOffices',
    testAlertCollection = 'testAlertCollection',
    testSpgi = 'testSpgi',
    testSpgi100 = 'testSpgi100',
    testUrl = 'testUrl',
  }

  const samples: { [key: string]: string } = {
    [Samples.fastaFasta]: 'System:AppData/Bio/data/sample_FASTA.fasta',
    [Samples.fastaPtCsv]: 'System:AppData/Bio/data/sample_FASTA_PT.csv',
    [Samples.msaComplex]: 'System:AppData/Bio/samples/sample_MSA.csv',
    [Samples.fastaCsv]: 'System:AppData/Bio/samples/sample_FASTA.csv',
    [Samples.helmCsv]: 'System:AppData/Bio/samples/sample_HELM.csv',
    [Samples.peptidesComplex]: 'System:AppData/Bio/tests/peptides_complex_msa.csv',
    [Samples.peptidesSimple]: 'System:AppData/Bio/tests/peptides_simple_msa.csv',
    [Samples.testDemogCsv]: 'System:AppData/Bio/tests/testDemog.csv',
    [Samples.testHelmCsv]: 'System:AppData/Bio/tests/testHelm.csv',
    [Samples.testIdCsv]: 'System:AppData/Bio/tests/testId.csv',
    [Samples.testSmilesCsv]: 'System:AppData/Bio/tests/testSmiles.csv',
    [Samples.testSmiles2Csv]: 'System:AppData/Bio/tests/testSmiles2.csv',
    [Samples.testSmilesShort]: 'System:AppData/Bio/tests/testSmilesShort.csv',
    [Samples.testActivityCliffsCsv]: 'System:AppData/Bio/tests/testActivityCliffs.csv', // smiles
    [Samples.testCerealCsv]: 'System:AppData/Bio/tests/testCereal.csv',
    [Samples.testUnichemSources]: 'System:AppData/Bio/tests/testUnichemSources.csv',
    [Samples.testDmvOffices]: 'System:AppData/Bio/tests/testDmvOffices.csv',
    [Samples.testAlertCollection]: 'System:AppData/Bio/tests/testAlertCollection.csv',
    [Samples.testSpgi100]: 'System:AppData/Bio/tests/testSpgi100.csv',
    [Samples.testSpgi]: 'System:AppData/Bio/tests/SPGI-derived.csv',
    [Samples.testUrl]: 'System:AppData/Bio/tests/testUrl.csv', // 10 rows
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

  const readCsv: (key: csvTests) => DfReaderFunc = (key: keyof typeof csvData) => {
    return async () => {
      // Always recreate test data frame from CSV for reproducible detector behavior in tests.
      const csv: string = csvData[key];
      const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
      await grok.data.detectSemanticTypes(df);
      return df;
    };
  };


  test('NegativeEmpty', async () => { await _testNeg(readCsv(csvTests.negEmpty), 'col1'); });
  test('Negative1', async () => { await _testNeg(readCsv(csvTests.neg1), 'col1'); });
  test('Negative2', async () => { await _testNeg(readCsv(csvTests.neg2), 'col1'); });
  test('Negative3', async () => { await _testNeg(readCsv(csvTests.neg3), 'col1'); });
  test('NegativeSmiles', async () => { await _testNeg(readCsv(csvTests.negSmiles), 'col1'); });

  test('FastaDna1', async () => {
    await _testPos(readCsv(csvTests.fastaDna1), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.DNA, 4, false);
  });
  test('FastaRna1', async () => {
    await _testPos(readCsv(csvTests.fastaRna1), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.RNA, 4, false);
  });
  test('FastaPt1', async () => {
    await _testPos(readCsv(csvTests.fastaPt1), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false);
  });
  test('FastaUn', async () => {
    await _testPos(readCsv(csvTests.fastaUn), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 12, true);
  });
  test('FastaMsaDna1', async () => {
    await _testPos(readCsv(csvTests.fastaMsaDna1), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.DNA, 4, false);
  });

  test('FastaMsaPt1', async () => {
    await _testPos(readCsv(csvTests.fastaMsaPt1), 'seq',
      NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.PT, 20, false);
  });

  test('SepDna', async () => {
    await _testPos(readCsv(csvTests.sepDna), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.DNA, 4, false, '*');
  });
  test('SepRna', async () => {
    await _testPos(readCsv(csvTests.sepRna), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.RNA, 4, false, '*');
  });
  test('SepPt', async () => {
    await _testPos(readCsv(csvTests.sepPt), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.PT, 20, false, '-');
  });
  test('SepUn1', async () => {
    await _testPos(readCsv(csvTests.sepUn1), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.UN, 8, true, '-');
  });
  test('SepUn2', async () => {
    await _testPos(readCsv(csvTests.sepUn2), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.UN, 9, true, '/');
  });

  test('SepMsaN1', async () => {
    await _testPos(readCsv(csvTests.sepMsaDna1), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ_MSA, ALPHABET.DNA, 4, false, '-');
  });

  test('samplesFastaCsv', async () => {
    await _testDf(readSamples(Samples.fastaCsv), {
      'Sequence': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false),
    });
  });

  test('samplesFastaFasta', async () => {
    await _testDf(readSamples(Samples.fastaFasta), {
      'sequence': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false),
    });
  });

  // peptidesComplex contains monomers with spaces in AlignedSequence columns, which are forbidden
  // test('samplesPeptidesComplexPositiveAlignedSequence', async () => {
  //   await _testPos(readSamples(Samples.peptidesComplex), 'AlignedSequence', 'separator:SEQ:UN', '-');
  // });
  test('samplesPeptidesComplex', async () => {
    await _testDf(readSamples(Samples.peptidesComplex), {} /* no positive */);
  });

  test('samplesMsaComplex', async () => {
    await _testDf(readSamples(Samples.msaComplex), {
      'MSA': new PosCol(NOTATION.SEPARATOR, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 161, true, '/'),
    });
  });

  test('samplesIdCsv', async () => {
    await _testDf(readSamples(Samples.testIdCsv), {} /* no positive */);
  });

  test('samplesSarSmallCsv', async () => {
    await _testDf(readSamples(Samples.testSmilesCsv), {} /* nopositive */);
  });

  test('samplesHelmCsv', async () => {
    await _testDf(readSamples(Samples.helmCsv), {
      'HELM': new PosCol(NOTATION.HELM, null, null, 160, true),
    });
  });

  // sample_testHelm.csv
  // columns: ID,Test type,HELM string,Valid?,Mol Weight,Mol Formula,SMILES
  test('samplesTestHelmCsv', async () => {
    await _testDf(readSamples(Samples.testHelmCsv), {
      'HELM string': new PosCol(NOTATION.HELM, null, null, 9, true),
    });
  });

  test('samplesTestDemogCsv', async () => {
    await _testDf(readSamples(Samples.testDemogCsv), {} /* no positive */);
  });

  test('samplesTestSmiles2Csv', async () => {
    await _testDf(readSamples(Samples.testSmiles2Csv), {} /* no positive */);
  });

  test('samplesTestSmilesShort', async () => {
    await _testDf(readSamples(Samples.testSmilesShort), {} /* no positive */);
  });

  test('samplesTestActivityCliffsNegativeSmiles', async () => {
    await _testDf(readSamples(Samples.testActivityCliffsCsv), {} /* no positive */);
  });

  test('samplesFastaPtCsv', async () => {
    await _testDf(readSamples(Samples.fastaPtCsv), {
      'sequence': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false),
    });
  });

  test('samplesTestCerealCsv', async () => {
    await _testDf(readSamples(Samples.testCerealCsv), {} /* no positive */);
  });

  test('samplesTestUnichemSources', async () => {
    await _testDf(readSamples(Samples.testUnichemSources), {} /* no positive */);
  });

  test('samplesTestDmvOffices', async () => {
    await _testDf(readSamples(Samples.testDmvOffices), {} /* no positive */);
  });

  test('samplesTestAlertCollection', async () => {
    await _testDf(readSamples(Samples.testAlertCollection), {} /* no positive */);
  });

  test('samplesTestSpgi', async () => {
    await _testDf(readSamples(Samples.testSpgi), {} /* no positive */);
  });

  test('samplesTestSpgi100', async () => {
    await _testDf(readSamples(Samples.testSpgi100), {} /* no positive */);
  });

  test('samplesTestUrl', async () => {
    await _testDf(readSamples(Samples.testUrl), {} /* no positive */);
  });
});

export async function _testNeg(readDf: DfReaderFunc, colName: string) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.getCol(colName)!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: col}) as unknown as string;
  if (semType)
    col.semType = semType;

  if (col.semType === DG.SEMTYPE.MACROMOLECULE) {
    const msg = `Negative test detected semType='${col.semType}', units='${col.getTag(DG.TAGS.UNITS)}'.`;
    throw new Error(msg);
    // col.semType = '';
    // col.setTag(DG.TAGS.UNITS, '');
    // col.setTag(NOTATION.SEPARATOR, '');
  }
}

export async function _testPos(
  readDf: DfReaderFunc, colName: string, units: string,
  aligned: string | null, alphabet: string | null, alphabetSize: number, alphabetIsMultichar: boolean,
  separator: string | null = null
) {
  const df: DG.DataFrame = await readDf();
  const col: DG.Column = df.col(colName)!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: col}) as unknown as string;
  if (semType)
    col.semType = semType;

  expect(col.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(col.getTag(DG.TAGS.UNITS), units);
  expect(col.getTag(bioTAGS.aligned), aligned);
  expect(col.getTag(bioTAGS.alphabet), alphabet);
  if (separator)
    expect(col.getTag(bioTAGS.separator), separator);

  const uh = new UnitsHandler(col);
  expect(uh.getAlphabetSize(), alphabetSize);
  expect(uh.getAlphabetIsMultichar(), alphabetIsMultichar);
  if (!uh.isHelm()) {
    expect(uh.aligned, aligned);
    expect(uh.alphabet, alphabet);
  }
}

class PosCol {
  constructor(
    public readonly units: string,
    public readonly aligned: string | null,
    public readonly alphabet: string | null,
    public readonly alphabetSize: number,
    public readonly alphabetIsMultichar: boolean,
    public readonly separator?: string
  ) { };
};

export async function _testDf(readDf: DfReaderFunc, posCols: { [colName: string]: PosCol }): Promise<void> {
  const df: DG.DataFrame = await readDf();
  const errList: string[] = [];
  for (const colName of df.columns.names()) {
    if (colName in posCols) {
      const p = posCols[colName];
      try {
        await _testPos(readDf, colName, p.units, p.aligned, p.alphabet,
          p.alphabetSize, p.alphabetIsMultichar, p.separator);
      } catch (err: any) {
        const errMsg: string = err.toString();
        errList.push(`Positive col '${colName}' failed: ${errMsg}`);
      }
    } else {
      try {
        await _testNeg(readDf, colName);
      } catch (err: any) {
        const errMsg: string = err.toString();
        errList.push(`Negative col '${colName}' failed: ${errMsg}`);
      }
    }
  }

  if (errList.length > 0)
    throw new Error(errList.join('\n'));
}
