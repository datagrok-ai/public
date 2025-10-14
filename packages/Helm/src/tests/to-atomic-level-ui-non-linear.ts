import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

type TestDataTargetType = { atomCount: number, bondCount: number };
type TestDataType = {
  src: { seq: string, units: NOTATION },
  tgt: TestDataTargetType,
};

category('toAtomicLevel-ui', () => {
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let monomerLib: IMonomerLib;
  let rdKitModule: RDModule;

  before(async () => {
    rdKitModule = await getRdKitModule();
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests(); // load default libraries

    monomerLib = monomerLibHelper.getMonomerLib();
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  const tests: { [testName: string]: TestDataType } = {
    'fasta': {
      src: {seq: 'MDYKETLLMPK', units: NOTATION.FASTA},
      tgt: {atomCount: 94, bondCount: 95},
    },
    'fasta-with-gap': {
      src: {seq: 'MD-YKETLLMPK', units: NOTATION.FASTA},
      tgt: {atomCount: 94, bondCount: 95},
    },
    'helm': {
      src: {seq: 'PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM},
      tgt: {atomCount: 68, bondCount: 68},
    },
    'helm-with-gap': {
      src: {seq: 'PEPTIDE1{meI.hHis.*.Aca.N.T.dK.Thr_PO3H2}$$$$', units: NOTATION.HELM},
      tgt: {atomCount: 68, bondCount: 68},
    },
  };

  const getDfAndSeqCol = async (testData: TestDataType): Promise<{df: DG.DataFrame, seqCol: DG.Column<string>}> => {
    const seq = testData.src.seq;
    const df = DG.DataFrame.fromColumns([DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'seq', [seq])]);
    await grok.data.detectSemanticTypes(df);
    return {df: df, seqCol: df.getCol('seq')};
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}-nonlinear`, async () => {
      const res = await getDfAndSeqCol(testData);
      await _testToAtomicLevelFunc(res.df, res.seqCol, true, testData.tgt);
    });
  }

  async function _testToAtomicLevelFunc(
    df: DG.DataFrame, seqCol: DG.Column<string>, nonlinear: boolean, tgt: TestDataTargetType,
  ): Promise<void> {
    await grok.functions.call('Bio:toAtomicLevel', {
      table: df, seqCol: seqCol, nonlinear: true, highlight: false
    });
    const molCol = df.col('molfile(seq)');
    expect(molCol?.semType, DG.SEMTYPE.MOLECULE);
    const resMolStr = molCol?.get(0)!;
    const resRdMol = rdKitModule.get_mol(resMolStr);
    expect(resRdMol != null, true, 'No molecule generated');
    try {
      const resAtomCount = resRdMol.get_num_atoms();
      const resBondCount = resRdMol.get_num_bonds();
      expect(resAtomCount, tgt.atomCount);
      expect(resBondCount, tgt.bondCount);
    } finally {
      resRdMol.delete();
    }
  }
});
