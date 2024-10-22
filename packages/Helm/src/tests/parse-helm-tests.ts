import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {after, before, category, delay, expect, test, expectArray, timeout} from '@datagrok-libraries/utils/src/test';
import {HelmType, IHelmEditorOptions, Mol, OrgType} from '@datagrok-libraries/bio/src/helm/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {JSDraw2Module} from '../types';
import {initHelmMainPackage} from './utils';

declare const org: OrgType;
declare const JSDraw2: JSDraw2Module;

type TestTgtType = { atomCount: number, bondCount: number, helm: string };

category('parseHelm', () => {
  const testData: { [testName: string]: { src: string, tgt: TestTgtType } } = {
    'PT': {
      src: 'PEPTIDE1{[L-hArg(Et,Et)].E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0',
      tgt: {
        atomCount: 6, bondCount: 5,
        helm: 'PEPTIDE1{[L-hArg(Et,Et)].E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0'
      }
    },
    'PT-cyclic': {
      src: 'PEPTIDE1{[meI].C.[Aca].[Cys_SEt].T.C.[Thr_PO3H2].[Aca]}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0',
      tgt: {
        atomCount: 8, bondCount: 8,
        helm: 'PEPTIDE1{[meI].C.[Aca].[Cys_SEt].T.C.[Thr_PO3H2].[Aca]}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0'
      },
    },
    'RNA': {
      src: 'RNA1{[dhp]([m1A])p.r(T)p.r(T)p.r(G)p}$$$$V2.0',
      tgt: {
        atomCount: 12, bondCount: 11,
        helm: 'RNA1{[dhp]([m1A])p.r(T)p.r(T)p.r(G)p}$$$$V2.0',
      }
    },
    'PT-with-gap': {
      src: 'PEPTIDE1{[meY].*.C.R.N.P.C.T}$$$$V2.0',
      tgt: {
        atomCount: 8, bondCount: 7,
        helm: 'PEPTIDE1{[meY].*.C.R.N.P.C.T}$$$$V2.0',
      }
    },
    'PT-with-gap-at-connection': {
      src: 'PEPTIDE1{[meY].*.C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$V2.0',
      tgt: {
        atomCount: 8, bondCount: 8,
        helm: 'PEPTIDE1{[meY].*.C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$V2.0',
      }
    },
  };

  let helmHelper: IHelmHelper;
  let libHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  before(async () => {
    await initHelmMainPackage();

    helmHelper = await getHelmHelper();

    libHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // parseHelm is dependent on monomers RGroups available, test requires default monomer library
    await libHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await libHelper.loadMonomerLib(true);
  });

  for (const [testName, {src, tgt}] of Object.entries(testData)) {
    test(`Editor-${testName}`, async () => {
      _testParseHelmWithEditor(src, tgt);
    });
  }

  for (const [testName, {src, tgt}] of Object.entries(testData)) {
    test(`HelmHelper-${testName}`, async () => {
      _testParseHelmWithHelmHelper(src, tgt, helmHelper);
    });
  }

  for (const [testName, {src, tgt}] of Object.entries(testData)) {
    test(`woDOM-${testName}`, async () => {
      _testParseHelmWithoutDOM(src, tgt);
    });
  }
});

function _testParseHelmWithEditor(src: string, tgt: TestTgtType): void {
  const io = org.helm.webeditor.IO;

  const editorHost = ui.div();
  const editor = new JSDraw2.Editor<HelmType, IHelmEditorOptions>(editorHost, {viewonly: true});

  const plugin = new org.helm.webeditor.Plugin(editor);
  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, src, origin, undefined);

  const m: Mol<HelmType> = plugin.jsd.m;
  expect(m.atoms.length, tgt.atomCount);
  expect(m.bonds.length, tgt.bondCount);
}

function _testParseHelmWithHelmHelper(src: string, tgt: TestTgtType, helmHelper: IHelmHelper): void {
  const io = org.helm.webeditor.IO;

  const origin = new JSDraw2.Point(0, 0);
  const resMol = helmHelper.parse(src, origin);

  const resHelm = io.getHelm(resMol);
  expect(resMol.atoms.length, tgt.atomCount);
  expect(resMol.bonds.length, tgt.bondCount);
  expect(resHelm, tgt.helm);
  expect(resMol.atoms.every((a) => !a.elem.startsWith('[') && !a.elem.endsWith(']')), true,
    'Atoms should not contain square braces.');
}

function _testParseHelmWithoutDOM(src: string, tgt: TestTgtType): void {
  const io = org.helm.webeditor.IO;

  const molHandler = new JSDraw2.MolHandler<HelmType, IHelmEditorOptions>();

  const plugin = new org.helm.webeditor.Plugin(molHandler);
  const origin = new JSDraw2.Point(0, 0);
  io.parseHelm(plugin, src, origin, undefined);

  const m: Mol<HelmType> = plugin.jsd.m;
  const resHelm = io.getHelm(m);
  expect(m.atoms.length, tgt.atomCount);
  expect(m.bonds.length, tgt.bondCount);
  expect(resHelm, tgt.helm);
  expect(m.atoms.every((a) => !a.elem.startsWith('[') && !a.elem.endsWith(']')), true,
    'Atoms should not contain square braces.');
}
