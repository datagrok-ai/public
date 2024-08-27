import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {after, before, category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {sequenceToMolfile} from '../utils/sequence-to-mol';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';

import {ConverterFunc} from './types';
import {_package} from '../package';

category('toAtomicLevel-ui', () => {

  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;
  let helmHelper: IHelmHelper;

  before(async () => {
    helmHelper = await getHelmHelper(); // init Helm package
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadMonomerLib(true); // load default libraries
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });
  const fastaCsv = `seq
MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW
MIEVFLFGIVLGLIPITLAGLFVTAYLQYRRGDQLDL
MMELVLKTIIGPIVVGVVLRIVDKWLNKDK
`;
  const helmCsv = `seq
PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.Aze.dV.E.N.dV.Phe_4Me}$$$$
PEPTIDE1{meI.hHis.Aca.N.T.dK.Thr_PO3H2.Aca.meM.D-Chg.dV.E.N.D-Orn.D-aThr.Phe_4Me}$$$$
PEPTIDE1{meI.Aca.N.T.dK.Thr_PO3H2.Aca.D-Tyr_Et.Tyr_ab-dehydroMe.dV.D-Cit.N.D-Orn.D-aThr.Phe_4Me}$$$$
`;

  test('toAtomicLevel-fasta-linear', async () => {
    const df = DG.DataFrame.fromCsv(fastaCsv);
    await grok.data.detectSemanticTypes(df);
    const seqCol = df.getCol('seq');
    await _testToAtomicLevelFunc(df, seqCol, false);
  });

  test('toAtomicLevel-fasta-nonlinear', async () => {
    const df = DG.DataFrame.fromCsv(fastaCsv);
    await grok.data.detectSemanticTypes(df);
    const seqCol = df.getCol('seq');
    await _testToAtomicLevelFunc(df, seqCol, true);
  });

  test('toAtomicLevel-helm', async () => {
    const df = DG.DataFrame.fromCsv(helmCsv);
    await grok.data.detectSemanticTypes(df);
    const seqCol = df.getCol('seq');
    await _testToAtomicLevelFunc(df, seqCol, true);
  });

  async function _testToAtomicLevelFunc(
    table: DG.DataFrame, seqCol: DG.Column<string>, nonlinear: boolean
  ): Promise<void> {
    const molCol = await sequenceToMolfile(table, seqCol, nonlinear, monomerLibHelper.getMonomerLib());
    expect(molCol!.semType, DG.SEMTYPE.MOLECULE);
  }
});
