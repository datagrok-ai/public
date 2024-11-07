/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {before, after, category, test, expectArray, expect} from '@datagrok-libraries/utils/src/test';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {_package} from '../package-test';

const appPath = 'System:AppData/Bio';
const fileSource = new DG.FileSource(appPath);

const enum Tests {
  PT = 'peptides-fasta',
  DNA = 'dna-fasta',
  MSA_SEPARATOR = 'msa-separator',
  MSA_FASTA = 'msa-fasta',
}

const TestsData: { [testName: string]: { inPath: string, outPath: string } } = {
  [Tests.PT]: {
    inPath: 'tests/to-atomic-level-peptides-fasta-input.csv',
    outPath: 'tests/to-atomic-level-peptides-fasta-output.csv'
  },
  [Tests.DNA]: {
    inPath: 'tests/to-atomic-level-dna-fasta-input.csv',
    outPath: 'tests/to-atomic-level-dna-fasta-output.csv'
  },
  [Tests.MSA_SEPARATOR]: {
    inPath: 'tests/to-atomic-level-msa-separator-input.csv',
    outPath: 'tests/to-atomic-level-msa-separator-output.csv'
  },
  [Tests.MSA_FASTA]: {
    inPath: 'tests/to-atomic-level-msa-fasta-input.csv',
    outPath: 'tests/to-atomic-level-msa-fasta-output.csv'
  },
};

const inputColName = 'sequence';
const outputColName = 'molfile(sequence)';

category('toAtomicLevel', async () => {
  const sourceDf: { [key: string]: DG.DataFrame } = {};
  const targetDf: { [key: string]: DG.DataFrame } = {};

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  let seqHelper: ISeqHelper;
  let monomerLib: IMonomerLib;
  let rdKitModule: RDModule;

  before(async () => {
    rdKitModule = await getRdKitModule();
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    // Clear settings to test default
    await monomerLibHelper.loadMonomerLibForTests();

    monomerLib = monomerLibHelper.getMonomerLib();

    for (const [testName, testData] of Object.entries(TestsData)) {
      const inputPath = testData.inPath;

      sourceDf[testName] = DG.DataFrame.fromCsv((await fileSource.readAsText(testData.inPath)).replace(/\n$/, ''));
      await grok.data.detectSemanticTypes(sourceDf[testName]);
      targetDf[testName] = DG.DataFrame.fromCsv((await fileSource.readAsText(testData.outPath)).replace(/\n$/, ''));
    }
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  async function getTestResult(source: DG.DataFrame, target: DG.DataFrame): Promise<void> {
    const inputCol = source.getCol(inputColName);
    // await toAtomicLevel(source, inputCol, false);
    await grok.functions.call('Bio:toAtomicLevel', {table: source, seqCol: inputCol, nonlinear: false});
    const obtainedCol = source.getCol(outputColName);
    const expectedCol = target.getCol(outputColName);
    const obtainedArray: string[] = wu(obtainedCol.values()).map((mol) => polishMolfile(mol)).toArray();
    const expectedArray: string[] = wu(expectedCol.values()).map((mol) => polishMolfile(mol)).toArray();
    expectArray(obtainedArray, expectedArray);
  }

  for (const [testName, testData] of Object.entries(TestsData)) {
    test(`${testName}`, async () => {
      await getTestResult(sourceDf[testName], targetDf[testName]);
    });
  }

  enum csvTests {
    fastaDna = 'fastaDna',
    fastaRna = 'fastaRna',
    fastaPt = 'fastaPt',
    fastaUn = 'fastaUn',

    separatorDna = 'separatorDna',
    separatorRna = 'separatorRna',
    separatorPt = 'separatorPt',
    separatorUn = 'separatorUn',

    helm = 'helm',
  }

  const csvData: { [key in csvTests]: string } = {
    [csvTests.fastaDna]: `seq
ACGTCACGTC
CAGTGTCAGTGT
TTCAACTTCAAC`,
    [csvTests.fastaRna]: `seq
ACGUCACGUC
CAGUGUCAGUGU
UUCAACUUCAAC`,
    [csvTests.fastaPt]: `seq
FWPHEYFWPHEY
YNRQWYVYNRQWYV
MKPSEYVMKPSEYV`,
    [csvTests.fastaUn]: `seq
[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]D[meI][hHis][Aca]NT[dE][Thr_PO3H2][Aca]D
[meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][meI][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2]
[Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2][Lys_Boc][hHis][Aca][Cys_SEt]T[dK][Thr_PO3H2][Aca][Tyr_PO3H2]`,
    [csvTests.separatorDna]: `seq
A/C/G/T/C/A/C/G/T/C
C/A/G/T/G/T/C/A/G/T/G/T
T/T/C/A/A/C/T/T/C/A/A/C`,
    [csvTests.separatorRna]: `seq
A*C*G*U*C*A*C*G*U*C
C*A*G*U*G*U*C*A*G*U*G*U
U*U*C*A*A*C*U*U*C*A*A*C`,
    [csvTests.separatorPt]: `seq
F-W-P-H-E-Y-F-W-P-H-E-Y
Y-N-R-Q-W-Y-V-Y-N-R-Q-W-Y-V
M-K-P-S-E-Y-V-M-K-P-S-E-Y-V`,
    [csvTests.separatorUn]: `seq
meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-D-meI-hHis-Aca-N-T-dE-Thr_PO3H2-Aca-D
meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-meI-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2
Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2-Lys_Boc-hHis-Aca-Cys_SEt-T-dK-Thr_PO3H2-Aca-Tyr_PO3H2`,

    [csvTests.helm]: `seq
PEPTIDE1{meI.D-gGlu.Aca.N.T.dE.Thr_PO3H2.Aca.D.Thr_PO3H2.Aca.D}$$$
PEPTIDE1{meI.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.Thr_PO3H2.Aca.Tyr_PO3H2}$$$
PEPTIDE1{Lys_Boc.hHis.Aca.Cys_SEt.T.dK.Thr_PO3H2.Aca.Tyr_PO3H2.Thr_PO3H2.Aca.Tyr_PO3H2}$$$`,
  };

  /** Also detects semantic types
   * @param {string} key
   * @return {Promise<DG.DataFrame>}
   */
  async function readCsv(key: csvTests): Promise<DG.DataFrame> {
    // Always recreate test data frame from CSV for reproducible detector behavior in tests.
    const csv: string = csvData[key];
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv.replace(/\n$/, ''));
    await grok.data.detectSemanticTypes(df);
    return df;
  }

  test('fastaDna', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.fastaDna), 'seq', monomerLibHelper);
  });

  test('fastaRna', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.fastaRna), 'seq', monomerLibHelper);
  });

  test('fastaPt', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.fastaPt), 'seq', monomerLibHelper);
  });

  test('fastaUn', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.fastaUn), 'seq', monomerLibHelper);
  });

  test('separatorDna', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.separatorDna), 'seq', monomerLibHelper);
  });

  test('separatorDna', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.separatorRna), 'seq', monomerLibHelper);
  });

  test('separatorPt', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.separatorPt), 'seq', monomerLibHelper);
  });

  test('separatorUn', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.separatorUn), 'seq', monomerLibHelper);
  });

  test('helm', async () => {
    await _testToAtomicLevel(await readCsv(csvTests.helm), 'seq', monomerLibHelper);
  });

  test('ptFasta2', async () => {
    const srcCsv: string = `seq\nAR`;
    const tgtMol: string = await _package.files.readAsText('tests/to-atomic-level-pt-fasta-2.mol');

    const srcDf = DG.DataFrame.fromCsv(srcCsv);
    const seqCol = srcDf.getCol('seq');
    seqCol.semType = DG.SEMTYPE.MACROMOLECULE;
    seqCol.meta.units = NOTATION.FASTA;
    seqCol.setTag(bioTAGS.alphabet, ALPHABET.PT);
    const sh = seqHelper.getSeqHandler(seqCol);
    const resCol = (await _testToAtomicLevel(srcDf, 'seq', monomerLibHelper))!;
    expect(polishMolfile(resCol.get(0)), polishMolfile(tgtMol));
  });

  async function _testToAtomicLevel(
    df: DG.DataFrame, seqColName: string = 'seq', monomerLibHelper: IMonomerLibHelper
  ): Promise<DG.Column | null> {
    const seqCol: DG.Column<string> = df.getCol(seqColName);
    const res = await _toAtomicLevel(df, seqCol, monomerLib, seqHelper, rdKitModule);
    if (res.warnings.length > 0)
      _package.logger.warning(`_toAtomicLevel() warnings ${res.warnings.join('\n')}`);
    return res.molCol;
  }
});


function polishMolfile(mol: string): string {
  return mol.replaceAll('\r\n', '\n')
    .replace(/\n$/, '')
    .split('\n').map((l) => l.trimEnd()).join('\n');
}
