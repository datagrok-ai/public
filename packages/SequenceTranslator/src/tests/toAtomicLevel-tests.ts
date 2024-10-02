import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray, testEvent, delay} from '@datagrok-libraries/utils/src/test';
import {Monomer, MonomerLibData} from '@datagrok-libraries/bio/src/types/index';

import {PolymerTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {getSeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {getRdKitModule} from '@datagrok-libraries/bio/src/chem/rdkit-module';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

import {_package} from '../package-test';
import {getNewMonomer} from '../polytool/pt-conversion';
import {getRules, RuleReaction} from '../polytool/pt-rules';

category('toAtomicLevel', () => {
  let userLibSettings: UserLibSettings;
  let monomerLibHelper: IMonomerLibHelper;
  let rdKitModule: RDModule;

  before(async () => {

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    rdKitModule = await getRdKitModule();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });

  test('override', async () => {
    const systemMonomerLib = monomerLibHelper.getMonomerLib();
    const rLibStr = await _package.files.readAsText('tests/polytool-reaction-lib.json');
    const rLib: Monomer[] = JSON.parse(rLibStr);
    const ggazM = rLib.find((m) => m.symbol === 'GGaz')!;
    expect(ggazM != null, true, `Monomer 'GGaz' not found.`);

    const overrideMonomerLibData: MonomerLibData = {[PolymerTypes.PEPTIDE]: {'GGaz': ggazM}};
    const overriddenMonomerLib = systemMonomerLib.override(overrideMonomerLibData);

    const seqHelper = await getSeqHelper();

    const helmCol = DG.Column.fromList(DG.COLUMN_TYPE.STRING, 'helm',
      ['PEPTIDE1{F.P.Y.[GGaz].H.A.A.G.G.A.C}|PEPTIDE2{A.A.A}$PEPTIDE1,PEPTIDE2,4:R4-1:R1|PEPTIDE1,PEPTIDE1,11:R2-4:R3$$$V2.0']);
    helmCol.semType = DG.SEMTYPE.MACROMOLECULE;
    helmCol.meta.units = NOTATION.HELM;
    const talRes = await seqHelper.helmToAtomicLevel(helmCol, false, false, overriddenMonomerLib);

    expect(talRes.molCol != null, true, 'Result molCol is null');
    const molfile = talRes.molCol!.get(0)!;
    expect(!!molfile, true, 'Molfile is empty');

    const mol = rdKitModule.get_mol(molfile);
    try {
      const molInchi = mol.get_inchi();
      const molInchiKey = rdKitModule.get_inchikey_for_inchi(molInchi);

      expect(mol.get_num_bonds(), 103);
      expect(mol.get_num_atoms(), 98);
      expect(molInchiKey, 'XPQUTLNSNIMERR-WKBTUBAPSA-N');
    } finally {
      mol.delete();
    }
  });

  test('getNewMonomer', async () => {
    const rdKitModule = await getRdKitModule();
    const systemMonomerLib = monomerLibHelper.getMonomerLib();

    const rules = await getRules(['rules_example.json']);
    const reactionRule = rules.reactionRules.find((r) => r.name == 'GGaz')!;

    const [newSymbol, newMonomer] = getNewMonomer(rdKitModule, systemMonomerLib, reactionRule);
    expect(newSymbol, reactionRule.name);

    const mol = rdKitModule.get_mol(newMonomer.molfile);
    try {
      const molInchi = mol.get_inchi();
      const molInchiKey = rdKitModule.get_inchikey_for_inchi(molInchi);
      expect(mol.get_num_bonds(), 18);
      expect(mol.get_num_atoms(), 18);
      // TODO: Check inchi key for the new monomer molfile
      // expect(molInchiKey, 'V2H10N2O3S-UHFFFAOYSA-N');
    } finally {
      mol.delete();
    }
  });
});
