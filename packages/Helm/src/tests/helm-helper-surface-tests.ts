import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmTypes} from '@datagrok-libraries/bio/src/helm/consts';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {getUserLibSettings, setUserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {HelmMol} from '@datagrok-libraries/bio/src/helm/types';

/** The laid-out molecule HelmInput hit-tests: the view-only editor's mol, atoms[].p on the editor surface. */
function getEditorMol(helmHelper: IHelmHelper, helm: string): HelmMol {
  const editor = helmHelper.createHelmWebEditor(ui.div()).editor;
  editor.setHelm(helm);
  return editor.m;
}

category('HelmHelper: getHoveredAtom', () => {
  let helmHelper: IHelmHelper;
  let libHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    helmHelper = await getHelmHelper();
    libHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await libHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await libHelper.loadMonomerLib(true);
  });

  const helm = 'PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$';

  test('at-atom', async () => {
    const mol = getEditorMol(helmHelper, helm);
    expect(mol.atoms.length, 9);
    for (const atom of [mol.atoms[0], mol.atoms[4], mol.atoms[8]])
      expect(helmHelper.getHoveredAtom(atom.p.x, atom.p.y, mol, 100) === atom, true);
  });

  test('far-away', async () => {
    const mol = getEditorMol(helmHelper, helm);
    expect(helmHelper.getHoveredAtom(-10000, -10000, mol, 100), null);
  });

  test('keeps-mol', async () => {
    const mol = getEditorMol(helmHelper, helm);
    const p0 = mol.atoms[0].p;
    for (const height of [40, 80, 200])
      helmHelper.getHoveredAtom(p0.x, p0.y, mol, height);
    expect(mol.atoms.length, 9);
    expect(mol.bonds.length, 8);
  });
});

category('HelmHelper: monomersFuncs', () => {
  let helmHelper: IHelmHelper;
  let libHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings;

  before(async () => {
    helmHelper = await getHelmHelper();
    libHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();
    await libHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await libHelper.loadMonomerLib(true);
  });

  test('build-from-lib', async () => {
    const lib = libHelper.getMonomerLib();
    const funcs = helmHelper.buildMonomersFuncsFromLib(lib);
    const a = funcs.getMonomer(HelmTypes.AA, 'A');
    expect(a?.id, 'A');
    expect(a?.id, lib.getWebEditorMonomer(HelmTypes.AA, 'A')?.id);
    const unknown = funcs.getMonomer(HelmTypes.AA, 'Xz_NotInLib');
    expect(unknown?.id, 'Xz_NotInLib');
    expect(unknown?.n, 'missing');
    const bracketed = funcs.getMonomer(HelmTypes.AA, '[meI]');
    expect(bracketed?.id, '[meI]');
    expect(bracketed?.n, 'missing');
    expect(funcs.getMonomerSet(HelmTypes.AA) != null, true);
  });
});
