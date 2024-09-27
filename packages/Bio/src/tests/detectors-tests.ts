import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect} from '@datagrok-libraries/utils/src/test';

import {ALIGNMENT, ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SeqHandler} from '@datagrok-libraries/bio/src/utils/seq-handler';

import {_testNeg, _testPos, DetectorTestData, DfReaderFunc, PosCol} from './utils/detectors-utils';
import {importFasta} from '../package';

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

category('detectors', () => {
  const enum csvTests {
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
    sepMsaUnWEmpty = 'sepMsaUnWEmpty',
    sepComplex = 'sepComplex',
    fastaMsaDna1 = 'fastaMsaDna1',
    fastaMsaPt1 = 'fastaMsaPt1',
    fastaMsaSameLength = 'fastaMsaSameLength',
    fastaExtSameLength = 'fastaExtSameLength',
    fastaMsaExtSameLength = 'fastaMsaExtSameLength',
    sepSameLength = 'sepSameLength',
    sepMsaSameLength = 'sepMsaSameLength',
    helmSameLength = 'helmSameLength',
  }

  const csvData2: DetectorTestData = {
    'negEmpty': {
      csv: `id,col1
1,
2,
3,
4,
5,`,
      neg: ['col1']
    },
    'negNum1': {
      csv: `col1
1
2
3`,
      neg: ['col1'],
    },
    'negNum2': {
      csv: `col1
4
5
6
7`,
      neg: ['col1'],
    },
    'negNum3': {
      csv: `col1
8
9
10
11
12`,
      neg: ['col1'],
    },

    'negSmiles': {
      csv: `col1
CCCCN1C(=O)CN=C(c2cc(F)ccc12)C3CCCCC3
C1CCCCC1
CCCCCC`,
      neg: ['col1'],
    },
    'negSmilesWithSquareBrackets': {
      csv: `col1
Cl.c1ccc2nc3ccccc3cc2c1
Oc1cccc2cc3ccccc3cc12
[SeH]c1ccc2ccccc2c1`,
      neg: ['col1'],
    },
    'negFastaUnSingleChar': {
      csv: `col1
Alanine
Cysteine
Aspartic acid
Glutamic acid
Phenylalanine`,
      neg: ['col1']
    },

    // Same length
    'fastaMsaSameLength': {
      csv: `seq
FWPHEYFWPHEYYV
YNRQWYVYNRQWYV
MKPSEYVMKPSEYV`,
      pos: {'seq': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.PT, 20, false, undefined)}
    },
    'fastaExtSameLength': {
      csv: `seq
FW[Ac]PHEYFWPH
YN[Re]VYNRQWYV
[Me]EYVMPS[Et]`,
      pos: {'seq': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.UN, 16, true, undefined)},
    },
    'fastaMsaExtSameLength': {
      csv: `seq
FW[Ac]PHEY[Re]WPH
YN[Re]VYNR[Ac]WYV
[Me]EYVMPSFW[Me]H`,
      pos: {'seq': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 14, true, undefined)},
    },
    'fastaMsaExtManyMinus': {
      csv: `seq
[D-Tic]-------[D-Tyr_Et][Tyr_ab-dehydroMe][dV][Cys_SEt]N[D-Orn][D-aThr]-[Phe_4Me]
[Phe_2F]--------[Tyr_ab-dehydroMe][dV][Aca]N[D-Orn][D-aThr]-[Phe_4Me]
[D-Tic]-[Hcy]QTWQ[Phe_4NH2][D-Tyr_Et][Tyr_ab-dehydroMe][dV][Cys_SEt]----[Phe_4Me]`,
      pos: {'seq': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 17, true, undefined)}
    },
    'sepSameLength': {
      csv: `seq
Aca-A-A-A-A-A-A-A-A-A-A-A-A-A-C-G-NH2
Aca-A-A-A-A-A-A-A-A-A-A-A-A-A-C-G-NH2
Aca-A-A-A-A-A-A-A-A-A-A-A-A-A-C-G-NH2`,
      pos: {'seq': new PosCol(NOTATION.SEPARATOR, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 5, true, '-')}
    },
    'sepMsaSameLength': {
      csv: `seq
Aca-A-A-A-A-A-A-A-A-A-A-A-A-A-Aca-G-NH2
Aca-A-Aca-A-A-A-meI-A-A-A-A-A-Aca-G-NH2
Aca-A-A-A-A-A-A-A-A-A-A-A-A-A-Aca-G-NH2`,
      pos: {'seq': new PosCol(NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.UN, 5, true, '-')}
    },
    'helmSameLength': {
      csv: `seq
PEPTIDE1{Ac(1).A.A.A.A.A.A.A.A.A.A.A.A.A.C(1).G.NH2}$$$$
PEPTIDE1{Ab(1).Y.V.K.H.P.F.W.R.W.Y.A.A.A.C(1).G.NH2}$$$$
PEPTIDE1{Ad(1).S.W.Y.C.K.H.P.M.W.A.A.A.A.C(1)-G-NH2}$$$$`,
      pos: {'seq': new PosCol(NOTATION.HELM, null, null, 19, undefined, undefined)}
    },
    'fastaNonDigitAlphabet': {
      csv: `flagC
"NMe-pyridazineH"
"Pyrrolo[2,3-c]pyridazineH"`,
      neg: ['flagC']
    }
  };

  const readCsv2: (key: keyof typeof csvData2) => DfReaderFunc = (key: keyof typeof csvData2) => {
    return async () => {
      const csv: string = csvData2[key].csv;
      const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
      await grok.data.detectSemanticTypes(df);
      return df;
    };
  };

  for (const [testName, testData] of Object.entries(csvData2)) {
    test(`csvData2-${testName}`, async () => {
      const reader = readCsv2(testName as csvTests);
      for (const negColName of testData.neg ?? [])
        await _testNeg(reader, negColName);
      for (const [posColName, posCol] of Object.entries(testData.pos ?? {})) {
        await _testPos(reader, posColName, posCol.units, posCol.aligned,
          posCol.alphabet, posCol.alphabetSize, posCol.alphabetIsMultichar, posCol.separator);
      }
    });
  }

  const csvData = new class {
    [csvTests.fastaDna1]: string = `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`;
    [csvTests.fastaRna1]: string = `seq
ACGUCACGUC
CAGUGUCAGUGU
UUCAACUUCAAC`;
    /** Pure amino acids sequence */
    [csvTests.fastaPt1]: string = `seq
FWPHEY
YNRQWYV
MKPSEYV`;
    [csvTests.fastaUn]: string = `seq
[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]DN
[meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][Aca]
[Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][Aca]`;
    [csvTests.sepDna]: string = `seq
A*C*G*T*C*A*C*G*T*C
C*A*G*T*G*T*C*A*G*T*G*T
T*T*C*A*A*C*T*T*C*A*A*C`;
    [csvTests.sepRna]: string = `seq
A*C*G*U*C*A*C*G*U*C
C*A*G*U*G*U*C*A*G*U*G*U
U*U*C*A*A*C*U*U*C*A*A*C`;
    [csvTests.sepPt]: string = `seq
F-W-P-H-E-Y-F-W-P-H-E-Y
Y-N-R-Q-W-Y-V-Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V-M-K-P-S-E-Y-V`;
    [csvTests.sepUn1]: string = `seq
abc-dfgg-abc1-cfr3-rty-wert-cfr3-rty-wert
rut12-her2-rty-wert-abc-abc1-dfgg-abc-abc1-dfgg
rut12-rty-her2-abc-cfr3-wert-rut12-cfr3-wert-rut12`;
    [csvTests.sepUn2]: string = `seq
abc/dfgg/abc1/cfr3/rty/wert/abc/dfgg/abc1/cfr3/rty/wert
rut12/her2/rty/wert//abc/abc1/dfgg/rut12/her2/rty/wert//abc/abc1/dfgg
rut12/rty/her2/abc/cfr3//wert/rut12/rut12/rty/her2/abc/cfr3//wert/rut12`;
    [csvTests.sepMsaDna1]: string = `seq
A-C--G-T--C-T-A-C--G-T--C-T
C-A-C--T--G-T-C-A-C--T--G-T
A-C-C-G-T-A-C-T-A-C-C-G-T-A-C-T`;
    [csvTests.sepMsaUnWEmpty]: string = `seq
m1-M-m3-mon4-mon5-N-T-MON8-N9-m1-M-m3-mon4-mon5-N-T-MON8-N9
m1-mon2-m3-mon4-mon5-Num--MON8-N9-m1-mon2-m3-mon4-mon5-Num--MON8-N9

mon1-M-mon3-mon4-mon5---MON8-N9-mon1-M-mon3-mon4-mon5---MON8-N9`;
    [csvTests.sepComplex]: string = `seq
Aca-F-K(AEEA-AEEA-R-Ac)-L-mF-V-Y-mNle-D-W-N-mF-Aca-G-NH2
Aca-F-K(AEEA-ARRA-W-Ac)-L-mF-V-Y-mNle-D-W-N-mF-Aca-G-NH2
Aca-F-K(AEEA-AEEA-Ac)-L-mF-V-Y-mNle-D-W-N-mF-Aca-G-NH2`;
    [csvTests.fastaMsaDna1]: string = `seq
AC-GT-CTAC-GT-CT
CAC-T-GTCAC-T-GT
ACCGTACTACCGTACT`;
    [csvTests.fastaMsaPt1]: string = `seq
FWR-WYV-KHPFWR-WYV-KHP
YNR-WYV-KHPYNR-WYV-KHP
MWRSWY-CKHPMWRSWY-CKHP`;
  }();

  const enum Samples {
    peptidesComplex = 'peptidesComplex',
    peptidesSimple = 'peptidesSimple',
    fastaCsv = 'fastaCsv',
    // fastaFasta = 'fastaFasta',
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
    // [Samples.fastaFasta]: 'System:AppData/Bio/samples/FASTA.fasta',
    [Samples.fastaPtCsv]: 'System:AppData/Bio/samples/FASTA_PT.csv',
    [Samples.msaComplex]: 'System:AppData/Bio/samples/MSA.csv',
    [Samples.fastaCsv]: 'System:AppData/Bio/samples/FASTA.csv',
    [Samples.helmCsv]: 'System:AppData/Bio/samples/HELM.csv',
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
    [Samples.testUrl]: 'System:AppData/Bio/tests/testUrl.csv',
  };

  const _samplesDfs: { [key: string]: Promise<DG.DataFrame> } = {};

  function readSamples(key: string, readFile: (file: string) => Promise<DG.DataFrame> = readFileCsv): DfReaderFunc {
    return async () => {
      if (!(key in _samplesDfs)) {
        _samplesDfs[key] = (async (): Promise<DG.DataFrame> => {
          const df: DG.DataFrame = await readFile(samples[key]);
          // await grok.data.detectSemanticTypes(df);
          return df;
        })().catch((err: any) => {
          delete _samplesDfs[key];
          throw err;
        });
      }
      return _samplesDfs[key];
    };
  };

  async function readFileCsv(file: string): Promise<DG.DataFrame> {
    const csv: string = await grok.dapi.files.readAsText(file);
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
    return df;
  }

  async function _readFileFasta(file: string): Promise<DG.DataFrame> {
    const txt: string = await grok.dapi.files.readAsText(file);
    const df: DG.DataFrame = importFasta(txt)[0];
    return df;
  }

  const readCsv: (key: keyof typeof csvData) => DfReaderFunc = (key: keyof typeof csvData) => {
    return async () => {
      // Always recreate test data frame from CSV for reproducible detector behavior in tests.
      const csv: string = csvData[key];
      const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
      await grok.data.detectSemanticTypes(df);
      return df;
    };
  };

  test('NegativeStartEnd', async () => { await _testNegList(['START', 'END']); });
  test('NegativeStartEndIntermediate', async () => { await _testNegList(['START', 'END', 'INTERMEDIATE']); });

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
  test('FastaPtGaps', () => _testPosList(['FW-PH-EYY', 'FYNRQWYV-', 'FKP-Q-SEYV'],
    NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false));
  test('FastaPtGapsMsa', () => _testPosList(['FW-PH-EYY', 'FYNRQWYV-', 'FKP-Q-SEY'],
    NOTATION.FASTA, ALIGNMENT.SEQ_MSA, ALPHABET.PT, 20, false));

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

  test('SepMsaUnWEmpty', async () => {
    await _testPos(readCsv(csvTests.sepMsaUnWEmpty), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ_MSA, ALPHABET.UN, 14, true);
  });

  test('SepComplex', async () => {
    await _testPos(readCsv(csvTests.sepComplex), 'seq',
      NOTATION.SEPARATOR, ALIGNMENT.SEQ, ALPHABET.UN, 17, true);
  });

  test('samplesFastaCsv', async () => {
    await _testDf(readSamples(Samples.fastaCsv), {
      'Sequence': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false),
    });
  });

  // test('samplesFastaFasta', async () => {
  //   await _testDf(readSamples(Samples.fastaFasta), {
  //     'sequence': new PosCol(NOTATION.FASTA, ALIGNMENT.SEQ, ALPHABET.PT, 20, false),
  //   });
  // });

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

export async function _testNegList(list: string[]): Promise<void> {
  const col: DG.Column = DG.Column.fromList(DG.TYPE.STRING, 'col1', list);
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (col.semType === DG.SEMTYPE.MACROMOLECULE) {
    const msg = `Negative test detected semType='${col.semType}', units='${col.meta.units}'.`;
    throw new Error(msg);
  }
}

export async function _testPosList(list: string[], units: NOTATION,
  aligned: ALIGNMENT, alphabet: ALPHABET, alphabetSize: number, alphabetIsMultichar: boolean,
  separator: string | null = null
): Promise<void> {
  const col: DG.Column = DG.Column.fromList(DG.TYPE.STRING, 'seq', list);
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  expect(col.semType, DG.SEMTYPE.MACROMOLECULE);
  expect(col.meta.units, units);
  expect(col.getTag(bioTAGS.aligned), aligned);
  expect(col.getTag(bioTAGS.alphabet), alphabet);
  if (separator)
    expect(col.getTag(bioTAGS.separator), separator);

  const sh = SeqHandler.forColumn(col);
  expect(sh.getAlphabetSize(), alphabetSize);
  expect(sh.getAlphabetIsMultichar(), alphabetIsMultichar);
  if (!sh.isHelm()) {
    expect(sh.aligned, aligned);
    expect(sh.alphabet, alphabet);
  }
}

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
