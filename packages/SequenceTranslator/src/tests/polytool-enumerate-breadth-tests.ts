import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes} from '../polytool/types';
import {doPolyToolEnumerateHelm} from '../polytool/pt-enumeration-helm';

import {_package} from '../package-test';

category('PolyTool: Enumerate', () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  let userLibSettings: UserLibSettings; //backup

  before(async () => {
    helmHelper = await getHelmHelper(); // initialize JSDraw2 and org

    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const tests: {
    [testName: string]: { src: string, params: PolyToolEnumeratorParams, tgt: { seq: string, name: string }[] }
  } = {
    'breadth1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        breadthPlaceholders: [
          {start: 2, end: 4, monomers: ['K']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.K.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-W3K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.K.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-G4K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
      ],
    },
    'breadth1-with-original': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        breadthPlaceholders: [
          {start: 2, end: 4, monomers: ['K']},
        ],
        keepOriginal: true,
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: ''},
        {seq: 'PEPTIDE1{[Ac(1)].F.K.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-W3K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.K.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-G4K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
      ],
    },
    'breadth2': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        breadthPlaceholders: [
          {start: 2, end: 4, monomers: ['K', 'Y']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.K.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-W3K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.K.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-G4K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.Y.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-W3Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.Y.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-G4Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.Y.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5Y'},
      ],
    },
    'breadth-double': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        breadthPlaceholders: [
          {start: 2, end: 4, monomers: ['K']},
          {start: 2, end: 4, monomers: ['Y']},
        ],
      },
      tgt: [
        {seq: "PEPTIDE1{[Ac(1)].F.Y.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-W3K-K3Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.K.Y.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-W3K-G4Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.K.G.Y.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-W3K-P5Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.Y.K.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-G4K-W3Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.W.Y.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-G4K-K4Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.W.K.Y.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-G4K-P5Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.Y.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-P5K-W3Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.W.Y.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-P5K-G4Y"},
        {seq: "PEPTIDE1{[Ac(1)].F.W.G.Y.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0", name: "-P5K-K5Y"}
      ]
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const res = doPolyToolEnumerateHelm(testData.src, '', testData.params);
      expectArray(res, testData.tgt.map((r) => [r.seq, r.name]));
    });
  }
});
