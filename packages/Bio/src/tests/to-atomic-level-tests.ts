/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';

import {before, after, category, test, expectArray, expect} from '@datagrok-libraries/test/src/test';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {_toAtomicLevel} from '@datagrok-libraries/bio/src/monomer-works/to-atomic-level';
import {IMonomerLib, Monomer} from '@datagrok-libraries/bio/src/types/monomer-library';
import {ALPHABET, NOTATION, TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';

import {_package} from '../package-test';

const appPath = 'System:AppData/Bio';
const fileSource = new DG.FileSource(appPath);

const complexMonomerAllylRgroup: Monomer = {
  'symbol': 'allyl_mon',
  'name': 'monomer with Allyl R group',
  'molfile': '\n     RDKit          2D\n\n  9  8  0  0  0  0  0  0  0  0999 V2000\n    1.4434   -2.1667    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4434   -0.6667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.1443    0.0833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1547   -0.6667    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.4537    0.0833    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.7528   -0.6667    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    0.1443    1.5833    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4434    2.3333    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    2.7424    0.0833    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\n  2  3  1  0\n  3  4  1  6\n  4  5  1  0\n  5  6  1  0\n  3  7  1  0\n  7  8  1  0\n  2  9  1  0\nM  RGP  3   6   3   8   1   9   2\nM  END\n',
  'smiles': 'O=C([C@H](CS[*:3])N[*:1])[*:2]',
  'polymerType': 'PEPTIDE',
  'monomerType': 'Backbone',
  'naturalAnalog': 'C',
  'id': 16,
  'rgroups': [
    {
      'alternateId': 'R1-H',
      'capGroupName': 'H',
      'capGroupSmiles': '[H][*:1]',
      'label': 'R1'
    },
    {
      'alternateId': 'R2-OH',
      'capGroupName': 'OH',
      'capGroupSmiles': 'O[*:2]',
      'label': 'R2'
    },
    {
      'alternateId': 'R3-Allyl',
      'capGroupName': 'Allyl',
      'capGroupSmiles': 'C=C[*:3]',
      'label': 'R3'
    }
  ],
  'author': 'Admin',
  'createDate': '2026-02-18T14:48:41.723Z',
  'meta': {}
};

const complexMonomerWithComplexRgroup: Monomer = {
  'symbol': 'SomeComplex',
  'name': 'Some complex monomer with complex R group',
  'molfile': '\n     RDKit          2D\n\n 10  9  0  0  0  0  0  0  0  0999 V2000\n   -1.4289   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.7280    0.3750    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.0270   -0.3750    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1299    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.1299    1.8750    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1691    2.6250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    1.1691   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1691   -1.8750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4682   -2.6250    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n    2.4682    0.3750    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  4  1  1  6\n  4  5  1  0\n  5  6  1  0\n  4  7  1  0\n  7  8  1  0\n  8  9  1  0\n  7 10  1  0\nM  RGP  4   3   3   6   1   9   4  10   2\nM  END\n',
  'smiles': '[*:4]OC([C@H](CS[*:3])N[*:1])[*:2]',
  'polymerType': 'PEPTIDE',
  'monomerType': 'Backbone',
  'naturalAnalog': 'C',
  'id': 16,
  'rgroups': [
    {
      'alternateId': 'R1-H',
      'capGroupName': 'H',
      'capGroupSmiles': '[H][*:1]',
      'label': 'R1'
    },
    {
      'alternateId': 'R2-OH',
      'capGroupName': 'OH',
      'capGroupSmiles': 'O[*:2]',
      'label': 'R2'
    },
    {
      'alternateId': 'R3-Something',
      'capGroupName': 'Something',
      'capGroupSmiles': 'C=CC([*:3])=C',
      'label': 'R3'
    },
    {
      'alternateId': 'R4-SomethingElse',
      'capGroupName': 'SomethingElse',
      'capGroupSmiles': 'ClCCCC=CC([*:4])=CCC',
      'label': 'R4'
    }
  ],
  'author': 'Admin',
  'createDate': '2026-02-18T14:48:41.723Z',
  'meta': {}
};

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
      sourceDf[testName].name = testData.inPath.split('/').pop()!;
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
    // DG.Utils.download(source.name.endsWith('.csv') ? source.name : source.name + '.csv', source.toCsv());
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
    // DG.Utils.download('molfile.mol', polishMolfile(resCol.get(0)));

    expect(polishMolfile(resCol.get(0)), polishMolfile(tgtMol));
  });

  async function _testToAtomicLevelWithCustomMonomer(srcHelm: string, expectedSmiles: string): Promise<void> {
    let error: any = null;
    // first, patch the monomer library with a custom monomers
    const monomerLib = monomerLibHelper.getMonomerLib();
    // @ts-ignore
    monomerLib._monomers['PEPTIDE'][complexMonomerAllylRgroup.symbol] = complexMonomerAllylRgroup;
    // @ts-ignore
    monomerLib._monomers['PEPTIDE'][complexMonomerWithComplexRgroup.symbol] = complexMonomerWithComplexRgroup;

    try {
      const converter = await seqHelper.getHelmToMolfileConverter(monomerLib);
      const resMolFile = seqHelper.helmToAtomicLevelSingle(srcHelm, converter, true, true);
      const resSmiles = grok.chem.convert(resMolFile.molfile, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
      expect(resSmiles, expectedSmiles);
    } catch (err) {
      error = err;
    }
    // restore the monomer library to avoid affecting other tests
    // @ts-ignore
    delete monomerLib._monomers['PEPTIDE'][complexMonomerAllylRgroup.symbol];
    // @ts-ignore
    delete monomerLib._monomers['PEPTIDE'][complexMonomerWithComplexRgroup.symbol];

    if (error)
      throw error;
  }

  test('SingleHelmMonomerWithAllylGroups', async () => {
    const srcHelm = `PEPTIDE1{[${complexMonomerAllylRgroup.symbol}]}$$$$V2.0`;
    const expectedSmiles = 'C=CSC[C@H](N)C(=O)O';
    await _testToAtomicLevelWithCustomMonomer(srcHelm, expectedSmiles);
  });

  test('SingleHelmMonomerWithComplexRGroups', async () => {
    const srcHelm = `PEPTIDE1{[${complexMonomerWithComplexRgroup.symbol}]}$$$$V2.0`;
    const expectedSmiles = 'C=CC(=C)SC[C@H](N)C(O)OC(C=CCCCCl)=CCC';
    await _testToAtomicLevelWithCustomMonomer(srcHelm, expectedSmiles);
  });

  test('HelmPolymerWithComplexRGroups', async () => {
    const srcHelm = `PEPTIDE1{[dI].[Trp_Ome].[Asp_OMe].[D-Cit].[meG].[Phe_4NH2].[Phe_34diCl].[meY].[Pro_4Me3OH].[Met_O].[NMe2Abz].[Tyr_Ph4OH].[3Pal].[xiIle].[Tyr_35diI].[Ala_tBu]}|PEPTIDE2{[${complexMonomerAllylRgroup.symbol}].[${complexMonomerWithComplexRgroup.symbol}]}$PEPTIDE1,PEPTIDE1,16:R2-1:R1|PEPTIDE1,PEPTIDE2,1:R3-1:R1$$$V2.0`;
    const expectedSmiles = 'C=CSC[C@H](NCC[C@@H](C)[C@H]1NC(=O)[C@H](C(C)(C)C)NC(=O)[C@H](Cc2cc(I)c(O)c(I)c2)NC(=O)[C@H](C(C)CC)NC(=O)[C@H](Cc2cccnc2)NC(=O)[C@H](Cc2ccc(Oc3ccc(O)cc3)cc2)NC(=O)c2ccccc2N(C)C(=O)[C@H](CCS(C)=O)NC(=O)[C@@H]2C(O)C(C)CN2C(=O)[C@H](Cc2ccc(O)cc2)N(C)C(=O)[C@H](Cc2ccc(Cl)c(Cl)c2)NC(=O)[C@H](Cc2ccc(N)cc2)NC(=O)CN(C)C(=O)[C@@H](CCCNC(N)=O)NC(=O)[C@H](CC(=O)OC)NC(=O)[C@H](Cc2cn(OC)c3ccccc23)NC1=O)C(=O)N[C@@H](CSC(=C)C=C)C(O)OC(C=CCCCCl)=CCC';
    await _testToAtomicLevelWithCustomMonomer(srcHelm, expectedSmiles);
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
