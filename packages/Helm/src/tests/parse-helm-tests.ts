import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';


import {after, before, category, delay, expect, test, expectArray, timeout} from '@datagrok-libraries/test/src/test';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {IHelmHelper, getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {initHelmMainPackage} from './utils';

type TestTgtType = { atomCount: number, bondCount: number, helm: string };

// hwe migration (Phase 7): `helmHelper.parse` now delegates to the standalone
// `@datagrok-libraries/hwe` parser (immutable `Mol`, exposed as a legacy
// `HelmMol`-shaped view). The legacy variants that drove `org.helm.webeditor.IO`
// / `new JSDraw2.Editor(...)` directly were removed with the forked editor.
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
    test(`HelmHelper-${testName}`, async () => {
      _testParseHelmWithHelmHelper(src, tgt, helmHelper);
    });
  }
});

function _testParseHelmWithHelmHelper(src: string, tgt: TestTgtType, helmHelper: IHelmHelper): void {
  const resMol = helmHelper.parse(src);

  expect(resMol.atoms.length, tgt.atomCount);
  expect(resMol.bonds.length, tgt.bondCount);
  expect(resMol.atoms.every((a) => !a.elem.startsWith('[') && !a.elem.endsWith(']')), true,
    'Atoms should not contain square braces.');
}
