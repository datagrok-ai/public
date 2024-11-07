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
    // Clear settings to test default
    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true);
  });

  const tests: {
    [testName: string]: { src: string, params: PolyToolEnumeratorParams, tgt: { seq: string, name: string }[] }
  } = {
    'single1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        placeholders: [
          {position: 4, monomers: ['K', 'P', 'F4COO']},
          {position: 6, monomers: ['Y', 'T']},
        ],
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5P'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5[F4COO]'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7T'},
      ]
    },
    'single-with-original': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        placeholders: [
          {position: 4, monomers: ['K', 'P', 'F4COO']},
          {position: 6, monomers: ['Y', 'T']},
        ],
        keepOriginal: true,
      },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: ''},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5K'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5P'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0', name: '-P5[F4COO]'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-[Tic]7T'},
      ]
    },
    'matrix1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
      params:
        {
          type: PolyToolEnumeratorTypes.Matrix,
          placeholders: [
            {position: 1, monomers: ['D', 'L']},
            {position: 4, monomers: ['K', 'P', 'F4COO']},
            {position: 6, monomers: ['Y', 'T']},
          ]
        },
      tgt: [
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5K-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5K-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5P-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5P-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5[F4COO]-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2D-P5[F4COO]-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5K-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5K-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5P-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5P-[Tic]7T'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5[F4COO]-[Tic]7Y'},
        {seq: 'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0', name: '-F2L-P5[F4COO]-[Tic]7T'},
      ],
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      const res = doPolyToolEnumerateHelm(testData.src, '', testData.params);
      expectArray(res, testData.tgt.map((r) => [r.seq, r.name]));
    });
  }
});
