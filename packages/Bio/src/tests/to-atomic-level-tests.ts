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

  // Explicit/inline SMILES monomers (e.g. `[*N(C)[C@H](C(=O)*)CCC |$_R1;...$|]`)
  // must round-trip through HELM → atomic level. They are resolved on-spot by
  // Datagrok's monomer functions (assigned `#P{N}` symbols) so the pseudo-molfile
  // carries a clean, space-free, resolvable symbol. Regression: post-migration
  // the raw SMILES (with its CXSMILES space) leaked into the pseudo-molfile and
  // got truncated, so the converter could not resolve the monomer.
  test('explicitSmilesMonomers', async () => {
    const cases: {name: string, helm: string, expected: string}[] = [
      {
        name: 'single-chain macrocycle',
        helm: `PEPTIDE1{[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@@H](C(=O)*)C |$_R1;;;;;_R2;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@H](C(=O)*)C |$_R1;;;;;_R2;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|]}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$`,
        expected: `CCC[C@H]1C(=O)O[C@H](Cc2ccccc2)C(=O)N(C)[C@@H](CCC)C(=O)O[C@H](C)C(=O)N(C)[C@@H](CCC)C(=O)O[C@H](Cc2ccccc2)C(=O)N(C)[C@@H](CCC)C(=O)O[C@@H](C)C(=O)N1C`,
      },
      {
        name: 'two-chain (named + inline SMILES monomers)',
        helm: `PEPTIDE1{[meL].[*O[C@@H](C(=O)*)C |$_R1;;;;;_R2;$|].[meL].[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@H](C(=O)*)C |$_R1;;;;;_R2;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|]}|PEPTIDE2{E.[*O[C@@H](C(=O)*)C |$_R1;;;;;_R2;$|].E.[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@H](C(=O)*)C |$_R1;;;;;_R2;$|].[*N(C)[C@H](C(=O)*)CCC |$_R1;;;;;;_R2;;;$|].[*O[C@@H](C(=O)*)Cc1ccccc1 |$_R1;;;;;_R2;;;;;;;$|]}$PEPTIDE2,PEPTIDE2,8:R2-1:R1|PEPTIDE1,PEPTIDE2,1:R1-3:R3|PEPTIDE2,PEPTIDE1,1:R3-8:R2$$$V2.0`,
        expected: `CCC[C@H]1C(=O)O[C@H](Cc2ccccc2)C(=O)N[C@H]2CCC(=O)C(=O)[C@@H](Cc3ccccc3)OC(=O)[C@H](CCC)N(C)C(=O)[C@H](C)OC(=O)[C@H](CCC)N(C)C(=O)[C@@H](Cc3ccccc3)OC(=O)[C@H](CC(C)C)N(C)C(=O)[C@@H](C)OC(=O)[C@H](CC(C)C)N(C)C(=O)CC[C@H](NC(=O)[C@@H](C)OC2=O)C(=O)O[C@H](Cc2ccccc2)C(=O)N(C)[C@@H](CCC)C(=O)O[C@@H](C)C(=O)N1C`,
      },
    ];
    const canonical = (smiles: string): string => {
      const mol = rdKitModule.get_mol(smiles);
      try { return mol.get_smiles(); } finally { mol.delete(); }
    };
    const converter = await seqHelper.getHelmToMolfileConverter(monomerLib);
    for (const {name, helm, expected} of cases) {
      const resMolFile = seqHelper.helmToAtomicLevelSingle(helm, converter, true, true);
      const resSmiles = grok.chem.convert(resMolFile.molfile, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
      expect(canonical(resSmiles), canonical(expected), `${name}: SMILES mismatch`);
    }
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

/** Tests for the linear HELM-RNA path: must preserve modified sugars,
 * phosphates, and bases per nucleotide. The non-linear (HELM via POM)
 * path is the reference; the linear path is expected to match it on
 * canonical SMILES for these inputs. */
category('toAtomicLevelHelmRna', async () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let seqHelper: ISeqHelper;
  let monomerLib: IMonomerLib;
  let rdKitModule: RDModule;

  before(async () => {
    rdKitModule = await getRdKitModule();
    seqHelper = await getSeqHelper();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await monomerLibHelper.loadMonomerLibForTests();
    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  // ---------- helpers --------------------------------------------------------

  /** Run the linear converter on a single HELM, returning both the molfile
   * and canonical SMILES. The molfile is the source of truth for structural
   * checks (atom indices, coordinates); the SMILES is kept for legacy /
   * presence-style assertions. */
  async function helmRnaLinear(srcHelm: string): Promise<{molfile: string; smiles: string}> {
    const srcCsv = `seq\n${srcHelm}`;
    const df = DG.DataFrame.fromCsv(srcCsv);
    await grok.data.detectSemanticTypes(df);
    const seqCol = df.getCol('seq');
    expect(seqCol.semType, DG.SEMTYPE.MACROMOLECULE);

    const res = await _toAtomicLevel(df, seqCol, monomerLib, seqHelper, rdKitModule);
    if (!res.molCol)
      throw new Error(`_toAtomicLevel returned no molCol for HELM '${srcHelm}'. ` +
        `Warnings: ${(res.warnings ?? []).join(' / ')}`);

    const molfile: string | null = res.molCol.get(0);
    if (!molfile)
      throw new Error(`_toAtomicLevel produced an empty molfile for HELM '${srcHelm}'`);
    let smiles: string;
    try {
      smiles = grok.chem.convert(molfile, grok.chem.Notation.Unknown, grok.chem.Notation.Smiles);
    } catch (err: any) {
      throw new Error(`SMILES conversion threw for HELM '${srcHelm}': ${err?.message ?? err}\n` +
        `--- MOLFILE START ---\n${molfile}\n--- MOLFILE END ---`);
    }
    if (smiles === 'MALFORMED_INPUT_VALUE' || /^MALFORMED/.test(smiles)) {
      throw new Error(`RDKit could not parse molfile produced for HELM '${srcHelm}'.\n` +
        `--- MOLFILE START ---\n${molfile}\n--- MOLFILE END ---`);
    }
    return {molfile, smiles};
  }

  /** Build an RDKit `RDMol` from the molfile, run `fn`, and free the mol.
   * Always pass the produced molfile (not its SMILES round-trip) — atom
   * indices and coordinates here are the same ones we want to assert on. */
  function withMol<T>(molfile: string, fn: (mol: any) => T): T {
    const mol = rdKitModule.get_mol(molfile);
    if (!mol || !mol.is_valid())
      throw new Error(`RDKit refused the produced molfile:\n${molfile}`);
    try {
      return fn(mol);
    } finally {
      mol.delete();
    }
  }

  /** True iff the molecule contains at least one match of the SMARTS query. */
  function hasSmarts(mol: any, smarts: string): boolean {
    const qmol = rdKitModule.get_qmol(smarts);
    try {
      const raw = mol.get_substruct_match(qmol);
      // RDKit JS returns the literal '{}' when there is no match.
      return !!raw && raw !== '{}';
    } finally {
      qmol.delete();
    }
  }

  /** Number of distinct matches of the SMARTS query in the molecule.
   * `get_substruct_matches` returns either '{}' (no match), a JSON array
   * of `{atoms,bonds}` objects, or — depending on the build — a single
   * match object. Normalise all three. */
  function countSmarts(mol: any, smarts: string): number {
    const qmol = rdKitModule.get_qmol(smarts);
    try {
      const raw = mol.get_substruct_matches(qmol);
      if (!raw || raw === '{}') return 0;
      const parsed = JSON.parse(raw);
      if (Array.isArray(parsed)) return parsed.length;
      // Single-match object
      if (parsed && typeof parsed === 'object' && Array.isArray(parsed.atoms))
        return parsed.atoms.length > 0 ? 1 : 0;
      return 0;
    } finally {
      qmol.delete();
    }
  }

  /** Atoms-by-element via a single-atom SMARTS — strictly counts the heavy
   * element (no false positives from `[Pa]`, `Si`, etc. that plain regex
   * on SMILES would produce). */
  function countAtoms(mol: any, atomicNumber: number): number {
    return countSmarts(mol, `[#${atomicNumber}]`);
  }

  /** SMARTS shortcuts used by several tests below. Bracketed atom specs are
   * deliberately permissive — the produced SMILES may render an atom
   * aromatic or kekulised depending on context. */
  const SMARTS = {
    // Generic phosphodiester backbone: C-O-P(=O)(X)-O-C with both bridging
    // oxygens present. X covers OH / O- (canonical p), SH / S- (sp), etc.
    PHOSPHODIESTER:
      '[#6][OX2][PX4](=[OX1])([OX2,SX2,OX1H,SX1H,OX1-,SX1-])[OX2][#6]',
    // Same but the non-bridging substituent is sulfur — phosphorothioate.
    PHOSPHOROTHIOATE_DIESTER:
      '[#6][OX2][PX4](=[OX1])([SX2,SX1H,SX1-])[OX2][#6]',
    // Direct sp3 C-P bond — appears ONLY when a bridging O on the linker
    // R-side has been (incorrectly) removed.
    DIRECT_C_P: '[CX4][PX4]',
    // P-O-P bridge — appears between consecutive phosphate units (a 5'/3'
    // di-/tri-phosphate or a run of phosphate linkers in the chain). Each
    // additional phosphate in a row adds one such bridge.
    P_O_P: '[PX4][OX2][PX4]',
    // Five-membered ring with exactly one ring oxygen — furanose.
    FURANOSE: '[#6;R]1[#6;R][#6;R][#6;R][O;R]1',
    // Adenine bicyclic core (aromatic Kekule-tolerant).
    ADENINE_RING: 'n1cnc2c1ncnc2N',
    // Cytosine 4-amino-pyrimidone.
    CYTOSINE_RING: 'Nc1ccn[cH0](=O)n1',
    // m5C: cytosine with a methyl at position 5.
    METHYL_CYTOSINE: '[CH3]c1cn([!#1])c(=O)nc1N',
    // 2'-fluoro on a sugar ring carbon (fl2r marker). Just `F` on a ring
    // sp3 C — no other monomer in our tests has fluorine, so this is
    // unambiguous; ring-position-specific patterns are too brittle to ring
    // traversal direction.
    FLUORO_ON_FURANOSE: '[F][CX4;R]',
    // Acetamide N-C(=O)-CH3 — GalNAc / N-acetyl marker.
    N_ACETYL: '[NX3]C(=O)[CH3]',
    // LNA-only marker: an sp3 carbon shared between two rings (R2). Plain
    // riboses have no such atom; LNA's bicyclic core puts C2', C3', C4'
    // each in two rings.
    LNA_BRIDGEHEAD: '[#6;R2]',
    // Methyl ether on a ring carbon (2'-OMe, the `m` ribose marker).
    TWO_PRIME_OME: '[CH3][OX2][#6;R]',
    // Biotin's cyclic urea (ureido) — 5-mem ring with N-C(=O)-N-C-C
    // pattern. The two C ring atoms are also bridgeheads to biotin's
    // thiolane ring (containing S), but we check that with a separate
    // ring-S query so this SMARTS stays robust to atom-order variations.
    BIOTIN_UREIDO: '[#7;R]1[#6;R](=[OX1])[#7;R][#6;R][#6;R]1',
    // Cholesterol gonane: four fused rings (3 cyclohexane + 1 cyclopentane).
    // Tested via two ring-counting heuristics rather than one rigid pattern,
    // see `looksLikeSteroid` below.
  } as const;

  /** Cholesterol detection: gonane has 4 fused rings; the D ring is a
   * cyclopentane (5-mem all-carbon) and the rest are cyclohexanes. None
   * of the other monomers we test against — sugars (always have a ring O),
   * nucleobases (always have N), biotin (5-mem rings have N or S), LNA
   * (5-mem rings have O) — produce an all-carbon 5-mem ring, so this
   * SMARTS is unique to steroids. We additionally require ≥ 4 ring
   * carbons in two rings (R2) to confirm a fused polycyclic system, not
   * an isolated cyclopentane. */
  function looksLikeSteroid(mol: any): boolean {
    const cyclopentane = hasSmarts(mol, '[#6]1[#6][#6][#6][#6]1');
    const fusedRingAtoms = countSmarts(mol, '[#6;R2]');
    return cyclopentane && fusedRingAtoms >= 4;
  }

  /** Parse a V3K molblock atom block into 0-indexed coordinate records.
   * The element symbol and x/y are sufficient for layout assertions; we
   * deliberately ignore z, charges, isotopes, etc. */
  function parseV3KAtoms(molfile: string): { element: string; x: number; y: number }[] {
    const atoms: { element: string; x: number; y: number }[] = [];
    const begin = molfile.indexOf('M  V30 BEGIN ATOM');
    if (begin < 0) return atoms;
    const end = molfile.indexOf('M  V30 END ATOM', begin);
    const block = molfile.substring(begin, end >= 0 ? end : molfile.length);
    // Capture the x/y coordinate TOKENS as `\S+` (not a strict numeric pattern)
    // so a non-finite coordinate ('NaN' / 'Infinity' / '-Infinity') is still
    // matched and surfaces as a NaN/Infinity from parseFloat — that is what
    // lets `expectNoNaN`'s finite-check below actually fire. A strict numeric
    // capture would silently skip the bad atom line and hide the defect.
    const lineRe = /^M\s+V30\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)/gm;
    let m: RegExpExecArray | null;
    while ((m = lineRe.exec(block))) {
      const idx = parseInt(m[1]) - 1;
      // Atoms are emitted in order; sanity check.
      if (idx !== atoms.length) continue;
      atoms.push({element: m[2], x: parseFloat(m[3]), y: parseFloat(m[4])});
    }
    return atoms;
  }

  /** Guard against the exact failure mode this suite was extended to fix:
   * NaN / Infinity atom coordinates. RDKit rejects such a molfile on the OCL
   * step, but we assert it directly on the produced molblock so the failure
   * message points at the layout bug rather than a downstream parse error.
   * Scans the raw V30 atom block for any non-finite coordinate token. */
  function expectNoNaN(molfile: string): void {
    const begin = molfile.indexOf('M  V30 BEGIN ATOM');
    const end = molfile.indexOf('M  V30 END ATOM', begin >= 0 ? begin : 0);
    const block = begin >= 0 ? molfile.substring(begin, end >= 0 ? end : molfile.length) : molfile;
    const bad = /\b(nan|-?inf(?:inity)?)\b/i.test(block);
    expect(bad, false, `produced molfile contains a non-finite (NaN/Infinity) coordinate:\n${molfile}`);
    // Also confirm every parsed coordinate is a finite number.
    for (const a of parseV3KAtoms(molfile)) {
      expect(Number.isFinite(a.x) && Number.isFinite(a.y), true,
        `non-finite coordinate on a ${a.element} atom (x=${a.x}, y=${a.y})`);
    }
  }

  /** Run a SMARTS against the molecule and collect every atom index that
   * appears in any match. Used to bin atoms by role (sugar / base / etc.). */
  function collectMatchedAtoms(mol: any, smarts: string): Set<number> {
    const set = new Set<number>();
    const qmol = rdKitModule.get_qmol(smarts);
    try {
      const raw = mol.get_substruct_matches(qmol);
      if (!raw || raw === '{}') return set;
      const parsed = JSON.parse(raw);
      const list = Array.isArray(parsed) ? parsed : [parsed];
      for (const m of list)
        for (const a of (m?.atoms ?? [])) set.add(a as number);
    } finally {
      qmol.delete();
    }
    return set;
  }

  /** Layout assertion: every atom in any nucleobase ring sits at a higher
   * Y than every sugar (furanose) ring atom. With the abnormal-sugar
   * override, the base is placed above the topmost atom of the sugar
   * cluster — including LNA's 2',4'-bridge oxygen / CH2. Without the
   * override the LNA bridge sits ABOVE the base attachment point and
   * this assertion fails. */
  function expectBaseAboveSugar(molfile: string): void {
    const atoms = parseV3KAtoms(molfile);
    if (atoms.length === 0) throw new Error(`failed to parse molblock atoms`);
    withMol(molfile, (mol) => {
      const sugarIdx = collectMatchedAtoms(mol, SMARTS.FURANOSE);
      // Base atoms = aromatic ring atoms (purines and pyrimidines aromatize
      // in RDKit's perception). Sugars are sp3, won't match `[a]`.
      const baseIdx = collectMatchedAtoms(mol, '[a]');
      if (sugarIdx.size === 0)
        throw new Error('no furanose ring atoms found — cannot verify layout');
      if (baseIdx.size === 0)
        throw new Error('no aromatic base atoms found — cannot verify layout');
      let maxSugarY = -Infinity;
      for (const i of sugarIdx) maxSugarY = Math.max(maxSugarY, atoms[i].y);
      let minBaseY = Infinity;
      for (const i of baseIdx) minBaseY = Math.min(minBaseY, atoms[i].y);
      expect(minBaseY > maxSugarY, true,
        `expected base atoms above sugar (minBaseY=${minBaseY.toFixed(3)}, ` +
        `maxSugarY=${maxSugarY.toFixed(3)})`);
    });
  }

  // Unmodified RNA HELM — regression baseline. The linear path must produce
  // a real RNA backbone: a furanose ring per nucleotide, two inter-nucleotide
  // phosphodiester linkers (C-O-P(=O)(O)-O-C) for three nucleotides, and
  // recognisable purine / pyrimidine bases attached to the sugars.
  test('rna-canonical', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{r(A)p.r(C)p.r(G)p}$$$$`);
    withMol(molfile, (mol) => {
      // 3 ribose furanose rings (one per nucleotide).
      const furanoses = countSmarts(mol, SMARTS.FURANOSE);
      expect(furanoses >= 3, true, `expected ≥ 3 furanose rings, got ${furanoses}`);
      // Inter-nucleotide phosphodiesters: r-r and r-r joints, so ≥ 2.
      // (The 3'-trailing P is a monoester and won't match the diester SMARTS.)
      const diesters = countSmarts(mol, SMARTS.PHOSPHODIESTER);
      expect(diesters >= 2, true,
        `expected ≥ 2 inter-nucleotide phosphodiester linkers, got ${diesters}`);
      // No direct sp3 C–P bond (would mean a bridging O was lost).
      const directCP = countSmarts(mol, SMARTS.DIRECT_C_P);
      expect(directCP, 0,
        `expected 0 direct C-P bonds (chain must use C-O-P-O-C), got ${directCP}`);
      // Purine ring (A and G are purines).
      const purines = countSmarts(mol, SMARTS.ADENINE_RING);
      expect(purines >= 1, true, `expected ≥ 1 purine ring, got ${purines}`);
      // Total phosphorus count: 3 (one per nucleotide as written).
      expect(countAtoms(mol, 15), 3, 'expected 3 phosphorus atoms');
    });
  });

  // Modified base — 5-methylcytosine. The methyl must end up at C5 of a
  // cytosine ring (not just any methyl on any ring), and only one m5C
  // appears in this row.
  test('rna-modified-base', async () => {
    const {molfile: plain} = await helmRnaLinear(`RNA1{r(C)p.r(A)p}$$$$`);
    const {molfile: mod} = await helmRnaLinear(`RNA1{r([m5C])p.r(A)p}$$$$`);
    withMol(plain, (mol) => {
      // No 5-methyl-cytosine in the plain version.
      expect(countSmarts(mol, SMARTS.METHYL_CYTOSINE), 0,
        'plain r(C) must not contain 5-methylcytosine');
    });
    withMol(mod, (mol) => {
      // Exactly one m5C ring; cytosine ring still present.
      expect(countSmarts(mol, SMARTS.METHYL_CYTOSINE), 1,
        'r([m5C]) must contain exactly one 5-methylcytosine ring');
    });
  });

  // Modified phosphate — phosphorothioate (Rsp). The S MUST be on the
  // phosphorus of the linker between positions 0 and 1 (not just somewhere
  // in the molecule), the linker must remain a diester (both bridging O
  // preserved), and the unmodified `p` at position 1 must stay unchanged.
  test('rna-modified-phosphate', async () => {
    const {molfile: plain} = await helmRnaLinear(`RNA1{r(A)p.r(C)p}$$$$`);
    const {molfile: mod} = await helmRnaLinear(`RNA1{r(A)[Rsp].r(C)p}$$$$`);
    withMol(plain, (mol) => {
      expect(countAtoms(mol, 16), 0, 'plain RNA must contain no sulfur');
      expect(countAtoms(mol, 15), 2, 'expected 2 phosphates in plain');
      expect(countSmarts(mol, SMARTS.PHOSPHODIESTER) >= 1, true,
        'plain inter-nucleotide diester must be present');
    });
    withMol(mod, (mol) => {
      // Sulfur is on phosphorus, not somewhere else.
      expect(hasSmarts(mol, '[PX4]=S') || hasSmarts(mol, '[PX4][SX2,SX1H,SX1-]'),
        true, 'sulfur must be bonded to a phosphorus atom');
      // Phosphorothioate diester has both bridging oxygens around the P.
      expect(countSmarts(mol, SMARTS.PHOSPHOROTHIOATE_DIESTER), 1,
        'expected exactly one phosphorothioate diester linker');
      // 2 phosphates total (Rsp + p).
      expect(countAtoms(mol, 15), 2, 'expected 2 phosphates in modified');
      // No direct C-P bond (regression check from sp/Rsp fix).
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0,
        'expected zero direct C-P bonds');
    });
  });

  // Modified sugar — 2'-fluoro ribose. F must end up on a ring carbon of
  // a furanose (i.e., a sugar atom), not on an arbitrary aliphatic carbon.
  test('rna-modified-sugar', async () => {
    const {molfile: plain} = await helmRnaLinear(`RNA1{r(A)p.r(C)p}$$$$`);
    const {molfile: mod} = await helmRnaLinear(`RNA1{[fl2r](A)p.r(C)p}$$$$`);
    withMol(plain, (mol) => {
      expect(countAtoms(mol, 9), 0, 'plain RNA must contain no fluorine');
    });
    withMol(mod, (mol) => {
      expect(countAtoms(mol, 9), 1, 'fl2r contributes exactly one fluorine');
      // F is on a ring carbon of a furanose.
      expect(countSmarts(mol, SMARTS.FLUORO_ON_FURANOSE) >= 1, true,
        'fluorine must be on a furanose ring carbon (2\'-F)');
      // Furanose count unchanged (one ribose replaced by 2'-F ribose).
      expect(countSmarts(mol, SMARTS.FURANOSE) >= 2, true,
        'expected ≥ 2 furanose rings');
    });
  });

  // HELM omits the trailing phosphate (3'-OH terminus on the sugar). The
  // splitter must split the partial `r(C)` into [r, C], assembly must skip
  // the trailing P emit, and counts must agree.
  test('rna-no-trailing-phosphate', async () => {
    const {molfile: withTail} = await helmRnaLinear(`RNA1{r(A)p.r(C)p}$$$$`);
    const {molfile: noTail} = await helmRnaLinear(`RNA1{r(A)p.r(C)}$$$$`);
    const pCountWith = withMol(withTail, (mol) => countAtoms(mol, 15));
    const pCountNoTail = withMol(noTail, (mol) => countAtoms(mol, 15));
    expect(pCountWith, 2, 'with trailing P: 2 phosphates (1 linker + 1 trail)');
    expect(pCountNoTail, 1, 'no trailing P: 1 phosphate (the linker only)');
    withMol(noTail, (mol) => {
      // The remaining phosphate is still a proper diester (both bridging O
      // present, no direct C-P bond).
      expect(countSmarts(mol, SMARTS.PHOSPHODIESTER), 1,
        'inter-nucleotide diester must still be present');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0,
        'no direct C-P bond');
      // Both furanose rings still present.
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'both furanose rings present');
    });
  });

  // Missing trailing phosphate combined with modifications.
  test('rna-no-trailing-phosphate-with-modifications', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[fl2r]([m5C])[Rsp].r(A)}$$$$`);
    withMol(molfile, (mol) => {
      // 1 F (2'-F on the fl2r sugar), on a furanose carbon.
      expect(countAtoms(mol, 9), 1, 'expected exactly 1 fluorine');
      expect(countSmarts(mol, SMARTS.FLUORO_ON_FURANOSE), 1,
        '2\'-F must be on a furanose ring carbon');
      // 1 P, 1 S — single Rsp linker, no trailing P.
      expect(countAtoms(mol, 15), 1, 'expected exactly 1 phosphorus (Rsp)');
      expect(countAtoms(mol, 16), 1, 'expected exactly 1 sulfur (Rsp)');
      // Linker is a phosphorothioate diester (both bridging O present).
      expect(countSmarts(mol, SMARTS.PHOSPHOROTHIOATE_DIESTER), 1,
        'Rsp linker must remain a phosphorothioate diester');
      // m5C base present.
      expect(countSmarts(mol, SMARTS.METHYL_CYTOSINE), 1,
        'expected one 5-methylcytosine base');
    });
  });

  // All three modifications combined. End-to-end smoke test — every
  // modification's structural fingerprint must be detectable.
  test('rna-all-modifications', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[fl2r]([m5C])[Rsp].r(A)p}$$$$`);
    withMol(molfile, (mol) => {
      expect(countSmarts(mol, SMARTS.FLUORO_ON_FURANOSE), 1,
        'fl2r: 2\'-F on furanose');
      expect(countSmarts(mol, SMARTS.METHYL_CYTOSINE), 1,
        'm5C: 5-methylcytosine');
      expect(countSmarts(mol, SMARTS.PHOSPHOROTHIOATE_DIESTER), 1,
        'Rsp: phosphorothioate diester');
      expect(countAtoms(mol, 15), 2, 'two phosphates (Rsp + trailing p)');
      expect(countAtoms(mol, 16), 1, 'exactly one sulfur (from Rsp)');
    });
  });

  // 3'-end terminal modifier (GalNAc, R1 only). HELM puts it in the
  // "phosphate" slot of the last triple, but it's actually a chain end.
  // GalNAc carries an N-acetyl group — that's the structural fingerprint
  // the test should pin to (not "any nitrogen", which thymine satisfies).
  test('rna-helm-3p-terminal-galnac', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{r(T)[GalNAc]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      // No phosphate at all (GalNAc replaces the trailing P slot).
      expect(countAtoms(mol, 15), 0, 'GalNAc terminus: no P expected');
      // Acetamide group from GalNAc — must be present.
      expect(countSmarts(mol, SMARTS.N_ACETYL) >= 1, true,
        'expected N-acetyl group from GalNAc');
      // GalNAc is a hexopyranose (6-mem ring with one O). Plus thymine ring
      // and the ribose furanose, the molecule has more than one ring.
      // Pyranose: C-C-C-C-C-O 6-membered.
      expect(hasSmarts(mol, '[#6]1[#6][#6][#6][#6][O]1'), true,
        'expected a pyranose (6-membered) ring from GalNAc');
    });
  });

  // 5'-end terminal modifier (Chol, R2 only) at the start of the chain.
  // Cholesterol's structural fingerprint is the gonane: four fused rings
  // including a cyclopentane fused to a cyclohexane (D-C ring junction).
  test('rna-helm-5p-terminal-chol', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[Chol].r(T)}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(looksLikeSteroid(mol), true,
        'Chol terminus must produce the steroid (gonane) ring system');
      // Chol replaces the first sugar — only one furanose left (from r(T)).
      expect(countSmarts(mol, SMARTS.FURANOSE), 1,
        'expected exactly 1 furanose ring (from r(T))');
    });
  });

  // Chol at 5' with explicit trailing phosphate (the original failing case).
  // Chain: Chol → r(T) → P-OH. Steroid rings + ribose + 1 phosphate.
  test('rna-helm-5p-terminal-chol-with-trailing-phosphate', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[Chol].r(T)p}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(looksLikeSteroid(mol), true,
        'expected steroid ring system from Chol');
      expect(countAtoms(mol, 15), 1, 'expected exactly 1 phosphorus');
      expect(countSmarts(mol, SMARTS.FURANOSE), 1, 'expected 1 furanose');
    });
  });

  // Both terminals at once: Chol at 5', GalNAc at 3', single nucleotide
  // between. Both terminus markers must be present, no phosphate.
  test('rna-helm-both-terminals', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[Chol].r(T)[GalNAc]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 0, 'expected zero phosphates');
      expect(looksLikeSteroid(mol), true,
        'expected steroid (Chol) ring system');
      expect(countSmarts(mol, SMARTS.N_ACETYL) >= 1, true,
        'expected N-acetyl group from GalNAc');
      // r(T) brings exactly one furanose, GalNAc brings the pyranose.
      expect(hasSmarts(mol, '[#6]1[#6][#6][#6][#6][O]1'), true,
        'expected pyranose ring from GalNAc');
    });
  });

  // LNA (2',4'-BNA) regression. The structural marker is the bicyclic
  // sugar: every ring carbon of the LNA furanose is shared with a second
  // ring (the C2'-O-CH2-C4' bridge). Standard riboses produce zero such
  // R2-shared atoms — so this test is exclusive to LNA.
  //
  // Additionally, the depiction-level claim ("base above sugar") is
  // verified by reading molblock coordinates and confirming the base
  // atoms sit higher in Y than every sugar atom.
  test('rna-helm-lna-base-above-sugar', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[lna](A)p.[lna](T)}$$$$V2.0`);
    withMol(molfile, (mol) => {
      // Single connected fragment.
      expect(hasSmarts(mol, '[*]'), true, 'molecule must be non-empty');
      // LNA-specific bicyclic sugar: ring atoms shared between two rings.
      // Two LNA sugars × 3 bridgehead-class carbons each = ≥ 4.
      const r2 = countSmarts(mol, SMARTS.LNA_BRIDGEHEAD);
      expect(r2 >= 4, true,
        `expected ≥ 4 ring carbons in 2 rings (LNA bicyclic), got ${r2}`);
      // Inter-nucleotide phosphodiester present, no direct C-P.
      expect(countSmarts(mol, SMARTS.PHOSPHODIESTER) >= 1, true,
        'expected ≥ 1 phosphodiester linker');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0,
        'expected zero direct C-P bonds');
      // Adenine + thymine present.
      expect(countSmarts(mol, SMARTS.ADENINE_RING) >= 1, true,
        'expected adenine ring (purine)');
      expect(hasSmarts(mol, '[CH3][#6]1=[#6][#7]([!#1])[#6](=O)[#7][#6]1=O') ||
             hasSmarts(mol, 'Cc1cn([!#1])c(=O)[nH]c1=O'),
        true, 'expected thymine ring (5-methyluracil)');
    });
    // Depiction: base atoms above sugar atoms in Y.
    expectBaseAboveSugar(molfile);
  });

  // GalNAc oxygen-count regression. Previously the R1 placeholder atom
  // (substituted to 'O' from the "OH" cap) was left in the assembly,
  // adding a stray OH on the chain-attach carbon. lna(T)GalNAc has a known
  // expected oxygen count; an extra OH would push it to 11.
  test('rna-helm-3p-terminal-galnac-no-extra-oh', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[lna](T)[GalNAc]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      // Heavy oxygen atom count — RDKit doesn't double-count ring closures
      // or atoms inside brackets.
      expect(countAtoms(mol, 8), 10,
        'expected exactly 10 oxygen atoms in lna(T)GalNAc');
      // No phosphate (GalNAc replaces the trailing P slot).
      expect(countAtoms(mol, 15), 0, 'expected no phosphate');
      // GalNAc N-acetyl preserved.
      expect(hasSmarts(mol, SMARTS.N_ACETYL), true,
        'expected GalNAc N-acetyl group');
      // LNA still bicyclic.
      expect(countSmarts(mol, SMARTS.LNA_BRIDGEHEAD) >= 2, true,
        'expected LNA bicyclic bridgeheads');
    });
  });

  // sp (and similar phosphates with R-cap = H) used to disconnect the chain
  // because the H placeholder was removed by removeHydrogen, leaving
  // terminalNodes[0] pointing at the now-deleted atom. The result was a
  // molecule with two disconnected fragments. The fix promotes the H cap
  // to an O so the chain bond attaches at a real atom; the linker becomes
  // a true phosphorothioate diester.
  test('rna-helm-h-cap-phosphate-sp-connects', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{r(T)[sp].r(A)}$$$$V2.0`);
    // SMILES dot count is the canonical fragment-count test — keep it.
    expect(smiles.indexOf('.') === -1, true,
      `expected single connected fragment, got: ${smiles}`);
    withMol(molfile, (mol) => {
      // Sulfur is bonded to the phosphorus, not floating somewhere else.
      expect(countSmarts(mol, '[PX4][SX2,SX1H,SX1-]'), 1,
        'sp\'s sulfur must be on its phosphorus');
      expect(countAtoms(mol, 15), 1, 'one phosphorus from the sp linker');
      expect(countAtoms(mol, 16), 1, 'one sulfur from the sp linker');
      expect(countSmarts(mol, SMARTS.PHOSPHOROTHIOATE_DIESTER), 1,
        'sp linker must be a phosphorothioate diester (C-O-P-O-C)');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0,
        'no direct C-P bond');
    });
  });

  // Regression: H-cap phosphates (sp et al.) used to drop the bridging O
  // on the 3' side of the linkage. The previous sugar's 3'-O is removed
  // unconditionally during sugar processing on the assumption that the
  // following linker brings its own bridging oxygen via the R1 cap; with
  // an H cap that assumption breaks and the chain ended up as
  // C3'-P(=O)(SH)-O-C5' instead of the proper C3'-O-P(=O)(SH)-O-C5'.
  // The fix promotes the H cap to an O so the bridging atom always exists.
  // Use m(2'-OMe ribose) so we can also verify the methoxy group survives
  // the sp chain assembly.
  test('rna-helm-sp-bridging-o-preserved', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{m(A)[sp].r(A)[sp]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      // Element counts via RDKit (no SMILES regex).
      expect(countAtoms(mol, 15), 2, 'expected exactly 2 phosphorus atoms');
      expect(countAtoms(mol, 16), 2, 'expected exactly 2 sulfur atoms');
      // Each P carries its own sulfur (not floating somewhere else).
      expect(countSmarts(mol, '[PX4][SX2,SX1H,SX1-]'), 2,
        'both sulfurs must be bonded to a phosphorus atom');
      // Inter-nucleotide sp is a phosphorothioate diester (bridging O on
      // both sides). The trailing sp is a monoester (P-O-cap on the 3'
      // side), so we expect exactly ONE diester match.
      expect(countSmarts(mol, SMARTS.PHOSPHOROTHIOATE_DIESTER), 1,
        'inter-nucleotide sp must remain a phosphorothioate diester');
      // Bridging-O presence on the 5' side of every phosphorothioate.
      // The diester P has two C-O-P matches (5' and 3' bridges) and the
      // terminal monoester P has one — total 3 matches across both linkers.
      // The bug we guard against (lost 3'-O) would drop this to 1 or 2.
      expect(countSmarts(mol, '[CX4][OX2][PX4](=[OX1])[SX2,SX1H,SX1-]'), 3,
        'every C-O-P-P=O-S match must be present (3: 2 from diester, 1 from monoester)');
      // No direct C-P bond anywhere (the bug we're guarding against).
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0,
        'expected zero direct C-P bonds — bridging O must be present');
      // Methoxy group on the m sugar must survive — exactly one (m only at
      // position 0). 2'-OMe = OCH3 on a ring carbon. The 2nd nucleotide is
      // r(A), no methoxy.
      expect(countSmarts(mol, SMARTS.TWO_PRIME_OME), 1,
        'expected exactly one 2\'-OMe group on the m sugar');
    });
  });

  // R-group swap heuristic: a single-R-group terminal monomer can be placed
  // at either end of a HELM chain, even if its R-group label "should" only
  // belong at one end. The conversion swaps rNodes so the existing
  // TERMINAL_5P/3P role logic still works. Each test asserts the terminal
  // monomer's STRUCTURAL fingerprint as well as topology.
  //
  // Bio (R1 only) — naturally a 3'-terminal, but we accept it at 5' too.
  // Biotin's fingerprint is its bicyclic head: a thiophene (C-C-C-C-S 5-mem
  // ring) fused to an imidazolidone (N-C(=O)-N 5-mem ring with two NH).
  test('rna-helm-bio-terminal-at-end', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{r(T)[Bio]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 0, 'Bio terminus: no phosphate');
      // Biotin's cyclic urea (ureido) ring.
      expect(hasSmarts(mol, SMARTS.BIOTIN_UREIDO), true,
        'expected biotin ureido (cyclic urea) ring system');
      // Biotin's thiolane: a sulfur in a ring.
      expect(hasSmarts(mol, '[#16;R]'), true,
        'expected ring sulfur (biotin\'s thiolane)');
      // r(T) sugar still present.
      expect(countSmarts(mol, SMARTS.FURANOSE), 1, 'expected the r(T) furanose');
    });
  });

  test('rna-helm-bio-terminal-at-start', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[Bio].r(T)}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 0, 'no phosphates');
      // Biotin's ureido + ring-S marker (the thiolane).
      expect(hasSmarts(mol, SMARTS.BIOTIN_UREIDO), true,
        'expected biotin ureido ring system at the 5\' end');
      expect(hasSmarts(mol, '[#16;R]'), true,
        'expected biotin\'s thiolane ring sulfur');
      // r(T) sugar still present and connected (single fragment via R-swap).
      expect(countSmarts(mol, SMARTS.FURANOSE), 1, 'expected the r(T) furanose');
    });
  });

  // Chol (R2 only) — naturally a 5'-terminal, but we accept it at 3' too.
  // Chol's structural fingerprint is the steroid 4-ring core plus a
  // ring-fused junction, see `looksLikeSteroid()`.
  test('rna-helm-chol-terminal-at-start', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{[Chol].r(T)}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(looksLikeSteroid(mol), true,
        'expected steroid (gonane) ring system from Chol at 5\'');
      expect(countSmarts(mol, SMARTS.FURANOSE), 1, 'expected one r(T) furanose');
    });
  });

  test('rna-helm-chol-terminal-at-end', async () => {
    const {molfile} = await helmRnaLinear(`RNA1{r(T)[Chol]}$$$$V2.0`);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 0, 'no phosphate when Chol replaces trailing P');
      expect(looksLikeSteroid(mol), true,
        'expected steroid (gonane) ring system from Chol at 3\'');
      expect(countSmarts(mol, SMARTS.FURANOSE), 1, 'expected one r(T) furanose');
    });
  });

  // ===================================================================
  // Phosphate-position regressions.
  //
  // The linear HELM-RNA path used to assign sugar/base/phosphate roles by
  // POSITION (index mod 3), assuming a strict [sugar, base, phosphate]
  // repeat. Any phosphate that did not land where that triple index
  // expected — a 5'-leading phosphate, several phosphates in a row, a
  // linker dropped between nucleotides — shifted every role, so a phosphate
  // was laid out as a sugar (branch-angle math on a 2-R-group monomer →
  // NaN coordinates) and RDKit rejected the molfile on the OCL step.
  //
  // Roles are now derived from each monomer's library definition (monomer
  // type + R-groups), so a phosphate is a phosphate wherever it sits. These
  // tests pin: (1) the molfile has no NaN/Infinity coordinates, (2) it is a
  // single connected fragment, (3) phosphorus / sulfur / fluorine counts
  // are exactly right, and (4) the expected backbone chemistry (diesters,
  // P-O-P bridges, no spurious direct C-P bonds) is present.
  // ===================================================================

  // 5'-leading phosphate (single). Chain: HO-P-O-r(A)-p-r(C)-p-OH. The
  // leading p must become a real 5'-monophosphate, not be mistaken for a
  // sugar.
  test('rna-5p-leading-phosphate', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{p.r(A)p.r(C)p}$$$$`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // leading p + inter-nucleotide p + trailing p = 3 phosphorus.
      expect(countAtoms(mol, 15), 3, 'expected 3 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
      // The inter-nucleotide linker is a proper phosphodiester.
      expect(countSmarts(mol, SMARTS.PHOSPHODIESTER) >= 1, true,
        'expected ≥ 1 inter-nucleotide phosphodiester');
      // The 5'-leading phosphate is not part of a P-O-P bridge (only one P
      // at the 5' end).
      expect(countSmarts(mol, SMARTS.P_O_P), 0, 'no P-O-P bridge for a single 5\' phosphate');
    });
  });

  // The user's exact 5'-leading-phosphate report. 1 leading p + 18 r(X)p
  // + a trailing phosphate-less r(C) ⇒ 19 phosphorus, 19 furanoses.
  test('rna-5p-leading-phosphate-long', async () => {
    const helm = `RNA1{p.r(C)p.r(A)p.r(C)p.r(A)p.r(A)p.r(G)p.r(T)p.r(T)p.r(T)p.r(A)p.r(T)p.r(A)p.r(T)p.r(T)p.r(C)p.r(A)p.r(G)p.r(T)p.r(C)}$$$$`;
    const {molfile, smiles} = await helmRnaLinear(helm);
    expect(smiles.indexOf('.') === -1, true, 'expected single connected fragment');
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 19, 'expected 19 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 19, 'expected 19 furanose rings');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // Multiple phosphates at the 5' end → a 5'-triphosphate (two P-O-P
  // bridges). Chain: HO-P-O-P-O-P-O-r(C)-p-r(A)-p-OH.
  test('rna-5p-multiple-leading-phosphates', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{p.p.p.r(C)p.r(A)p}$$$$`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // 3 leading + 1 inter-nucleotide + 1 trailing = 5 phosphorus.
      expect(countAtoms(mol, 15), 5, 'expected 5 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      // The 5'-triphosphate run contributes two P-O-P bridges.
      expect(countSmarts(mol, SMARTS.P_O_P) >= 2, true,
        'expected ≥ 2 P-O-P bridges from the 5\' triphosphate run');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // The user's exact multiple-leading-phosphates report. 3 leading p + 18
  // r(X)p + trailing r(C) ⇒ 21 phosphorus, 19 furanoses.
  test('rna-5p-multiple-leading-phosphates-long', async () => {
    const helm = `RNA1{p.p.p.r(C)p.r(A)p.r(C)p.r(A)p.r(A)p.r(G)p.r(T)p.r(T)p.r(T)p.r(A)p.r(T)p.r(A)p.r(T)p.r(T)p.r(C)p.r(A)p.r(G)p.r(T)p.r(C)}$$$$`;
    const {molfile, smiles} = await helmRnaLinear(helm);
    expect(smiles.indexOf('.') === -1, true, 'expected single connected fragment');
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 21, 'expected 21 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 19, 'expected 19 furanose rings');
      expect(countSmarts(mol, SMARTS.P_O_P) >= 2, true,
        'expected ≥ 2 P-O-P bridges from the 5\' triphosphate run');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // Multiple thio-linkers in the middle of the chain. Chain:
  // r(A)-p-[sp]-[sp]-r(C)-p. The two consecutive sp linkers must chain
  // (P-O-P) without disconnecting the molecule or losing a bridging O.
  test('rna-middle-multiple-linkers', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{r(A)p.[sp].[sp].r(C)p}$$$$V2.0`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // p + sp + sp + p = 4 phosphorus; sp + sp = 2 sulfur.
      expect(countAtoms(mol, 15), 4, 'expected 4 phosphorus atoms');
      expect(countAtoms(mol, 16), 2, 'expected 2 sulfur atoms (from the two sp)');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      // A run of three phosphates (p-sp-sp) ⇒ at least two P-O-P bridges.
      expect(countSmarts(mol, SMARTS.P_O_P) >= 2, true,
        'expected ≥ 2 P-O-P bridges in the p-sp-sp run');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // The user's exact mid-chain linkers report. 18 r(X)p + trailing r(C) +
  // two [sp] linkers between positions ⇒ 20 phosphorus, 2 sulfur, 19
  // furanoses.
  test('rna-middle-multiple-linkers-long', async () => {
    const helm = `RNA1{r(C)p.r(A)p.r(C)p.r(A)p.r(A)p.r(G)p.r(T)p.r(T)p.r(T)p.r(A)p.[sp].[sp].r(T)p.r(A)p.r(T)p.r(T)p.r(C)p.r(A)p.r(G)p.r(T)p.r(C)}$$$$V2.0`;
    const {molfile, smiles} = await helmRnaLinear(helm);
    expect(smiles.indexOf('.') === -1, true, 'expected single connected fragment');
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 20, 'expected 20 phosphorus atoms');
      expect(countAtoms(mol, 16), 2, 'expected 2 sulfur atoms (from the two sp)');
      expect(countSmarts(mol, SMARTS.FURANOSE), 19, 'expected 19 furanose rings');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // Multiple phosphates at the 3' end → a 3'-triphosphate. Chain:
  // r(A)-p-r(C)-p-p-p, capped with OH. Trailing run of three P ⇒ two P-O-P.
  test('rna-3p-multiple-trailing-phosphates', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{r(A)p.r(C)p.p.p}$$$$`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // inter-nucleotide p + 3 trailing p = 4 phosphorus.
      expect(countAtoms(mol, 15), 4, 'expected 4 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      expect(countSmarts(mol, SMARTS.P_O_P) >= 2, true,
        'expected ≥ 2 P-O-P bridges from the 3\' triphosphate run');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // Deoxyribose (DNA) sugar with a 5'-leading phosphate, to confirm the
  // role classification is not RNA-ribose specific.
  test('rna-5p-leading-phosphate-deoxyribose', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{p.d(A)p.d(C)p}$$$$`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      expect(countAtoms(mol, 15), 3, 'expected 3 phosphorus atoms');
      // Deoxyribose still has a ring oxygen ⇒ matches the furanose pattern.
      // Each monocyclic (deoxy)furanose matches exactly once (RDKit uniquify),
      // so the count is exactly 2 — assert it exactly like the ribose tests.
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected exactly 2 (deoxy)furanose rings');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
    });
  });

  // Methylphosphonate linker (mp) in the middle. mp carries an intramonomer
  // C-P bond by design, so DIRECT_C_P is expected here — we instead assert
  // the molecule stays connected, the bridging oxygens survive on both
  // sides, and the phosphorus count is right.
  test('rna-middle-methylphosphonate', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{r(A)[mp].r(C)p}$$$$V2.0`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // mp linker + trailing p = 2 phosphorus.
      expect(countAtoms(mol, 15), 2, 'expected 2 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      // The mp linker is a methylphosphonate diester: P with two bridging
      // C-O-P arms (one to each sugar). Exactly one such diester here.
      expect(countSmarts(mol, '[#6][OX2][PX4](=[OX1])([CH3])[OX2][#6]'), 1,
        'expected the methylphosphonate diester linkage');
    });
  });

  // Broadest torture case: a 5'-leading phosphate, a run of three
  // non-canonical linkers (Rsp / s2p / sp) in the middle, non-canonical
  // sugars (fana 2'-F, ena bicyclic, moe 2'-MOE) and non-canonical bases
  // (m1A, N-acetyl ac4C, 2-thio s2C). If this assembles cleanly, the linear
  // path handles arbitrary backbone composition.
  test('rna-complex-noncanonical-mix', async () => {
    const helm = `RNA1{p.[fana]([m1A])p.[Rsp].[s2p].[ena]([ac4C])[sp].[moe]([s2C])p}$$$$V2.0`;
    const {molfile, smiles} = await helmRnaLinear(helm);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // Phosphorus: leading p + p + Rsp + s2p + sp + trailing p = 6.
      expect(countAtoms(mol, 15), 6, 'expected 6 phosphorus atoms');
      // Sulfur: Rsp(1) + s2p(2) + sp(1) + s2C(1) = 5.
      expect(countAtoms(mol, 16), 5, 'expected 5 sulfur atoms');
      // Fluorine: fana 2'-F = 1, on a furanose ring carbon.
      expect(countAtoms(mol, 9), 1, 'expected exactly 1 fluorine (fana 2\'-F)');
      expect(countSmarts(mol, SMARTS.FLUORO_ON_FURANOSE) >= 1, true,
        'fluorine must sit on a furanose ring carbon');
      // Furanoses (ring O): fana + ena + moe = 3 (≥ 3 to tolerate the
      // bicyclic ena ring-match multiplicity).
      expect(countSmarts(mol, SMARTS.FURANOSE) >= 3, true, 'expected ≥ 3 furanose rings');
      // ac4C contributes an N-acetyl group.
      expect(countSmarts(mol, SMARTS.N_ACETYL) >= 1, true,
        'expected N-acetyl group from ac4C');
      // No bridging oxygen lost anywhere along this dense backbone.
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
      // The non-canonical linker run forms exactly two P-O-P bridges at the
      // two consecutive-phosphate joints p–Rsp and Rsp–s2p.
      expect(countSmarts(mol, SMARTS.P_O_P) >= 2, true,
        'expected ≥ 2 P-O-P bridges across the linker run');
    });
  });

  // 5'-leading phosphate combined with a 3'-terminal modifier (GalNAc).
  // Chain: HO-P-O-r(A)-p-r(C)-GalNAc. This is the layout that exercises BOTH
  // a leading phosphate AND the has3pTerm branch of getResultingAtomBondCounts
  // (no trailing OH cap): the declared bond count must equal the emitted bond
  // lines so the produced molfile is well-formed before the OCL pass — not
  // only after it. The terminal modifier replaces the trailing phosphate.
  test('rna-5p-leading-phosphate-3p-galnac', async () => {
    const {molfile, smiles} = await helmRnaLinear(`RNA1{p.r(A)p.r(C)[GalNAc]}$$$$V2.0`);
    expect(smiles.indexOf('.') === -1, true, `expected single fragment, got: ${smiles}`);
    expectNoNaN(molfile);
    withMol(molfile, (mol) => {
      // leading p + inter-nucleotide p = 2 phosphorus (GalNAc replaces the
      // trailing phosphate, so no 3' P).
      expect(countAtoms(mol, 15), 2, 'expected 2 phosphorus atoms');
      expect(countSmarts(mol, SMARTS.FURANOSE), 2, 'expected 2 furanose rings');
      expect(countSmarts(mol, SMARTS.DIRECT_C_P), 0, 'expected zero direct C-P bonds');
      // GalNAc's N-acetyl group at the 3' terminus.
      expect(countSmarts(mol, SMARTS.N_ACETYL) >= 1, true,
        'expected N-acetyl group from the 3\' GalNAc');
      // GalNAc is a hexopyranose (6-membered ring with one O).
      expect(hasSmarts(mol, '[#6]1[#6][#6][#6][#6][O]1'), true,
        'expected a pyranose ring from GalNAc');
    });
  });
});


function polishMolfile(mol: string): string {
  return mol.replaceAll('\r\n', '\n')
    .replace(/\n$/, '')
    .split('\n').map((l) => l.trimEnd()).join('\n');
}
