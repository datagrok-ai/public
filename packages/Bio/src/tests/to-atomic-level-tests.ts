/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, test, expectArray} from '@datagrok-libraries/utils/src/test';

import {toAtomicLevel} from '../package';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types/index';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, LibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';


const appPath = 'System:AppData/Bio';
const fileSource = new DG.FileSource(appPath);

const testNames: { [k: string]: string } = {
  PT: 'peptides fasta',
  DNA: 'dna fasta',
  MSA: 'msa separator',
};

const inputPath: { [k: string]: string } = {
  PT: 'tests/to-atomic-level-peptides-fasta-input.csv',
  DNA: 'tests/to-atomic-level-dna-fasta-input.csv',
  MSA: 'tests/to-atomic-level-msa-separator-input.csv',
};

const outputPath: { [k: string]: string } = {
  PT: 'tests/to-atomic-level-peptides-output.csv',
  DNA: 'tests/to-atomic-level-dna-output.csv',
  MSA: 'tests/to-atomic-level-msa-output.csv',
};

const inputColName = 'sequence';
const outputColName = 'molfile(sequence)';

category('toAtomicLevel', async () => {
  const sourceDf: { [key: string]: DG.DataFrame } = {};
  const targetDf: { [key: string]: DG.DataFrame } = {};

  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: LibSettings;

  before(async () => {
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    // Clear settings to test default
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadLibraries(true);

    for (const key in testNames) {
      sourceDf[key] = await fileSource.readCsv(inputPath[key]);
      await grok.data.detectSemanticTypes(sourceDf[key]);
      targetDf[key] = await fileSource.readCsv(outputPath[key]);
    }
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true);
  });

  async function getTestResult(source: DG.DataFrame, target: DG.DataFrame): Promise<void> {
    const inputCol = source.getCol(inputColName);
    await toAtomicLevel(source, inputCol, false);
    const obtainedCol = source.getCol(outputColName);
    const expectedCol = target.getCol(outputColName);
    const obtainedArray = [...obtainedCol.values()];
    const expectedArray = [...expectedCol.values()];
    expectArray(obtainedArray, expectedArray);
  }

  for (const key in testNames) {
    test(`${testNames[key]}`, async () => {
      await getTestResult(sourceDf[key], targetDf[key]);
    });
  }

  enum csvTests {
    fastaDna = 'fastaDna',
    fastaRna = 'fastaRna',
    fastaPt = 'fastaPt',

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
    const df: DG.DataFrame = DG.DataFrame.fromCsv(csv);
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
});

async function _testToAtomicLevel(df: DG.DataFrame, seqColName: string = 'seq', monomerLibHelper: IMonomerLibHelper) {
  const seqCol: DG.Column<string> = df.getCol(seqColName);
  const monomerLib: IMonomerLib = monomerLibHelper.getBioLib();
  const _resCol = await _toAtomicLevel(df, seqCol, monomerLib);
}
