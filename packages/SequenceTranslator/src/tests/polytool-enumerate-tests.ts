import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect, test, expectArray} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes} from '../polytool/types';
import {polyToolEnumerateHelm} from '../package';
import {getPtEnumeratorHelm} from '../polytool/pt-enumeration-helm';
import {HelmMol} from '@datagrok-libraries/bio/src/helm/types';

category('PolyTool', () => {
  let helmHelper: IHelmHelper;

  before(async () => {
    helmHelper = await getHelmHelper(); // initialize JSDraw2 and org
  });

  after(async () => {

  });

  const tests: {
    [testName: string]: { src: string, params: PolyToolEnumeratorParams, tgt: string[] }
  } = {
    'single1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$',
      params: {
        type: PolyToolEnumeratorTypes.Single,
        placeholders: {
          [4]: ['K', 'P', 'F4COO'],
          [6]: ['Y', 'T'],
        }
      },
      tgt: [
        'PEPTIDE1{[Ac(1)].F.W.G.K.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].F.W.G.[F4COO].L.[Tic].[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].F.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].F.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0',
      ]
    },
    'matrix1': {
      src: 'PEPTIDE1{[Ac(1)].F.W.G.P.L.[Tic].[C(1)].G.[NH2]}$$$$',
      params:
        {
          type: PolyToolEnumeratorTypes.Matrix,
          placeholders: {
            [1]: ['D', 'L'],
            [4]: ['K', 'P', 'F4COO'],
            [6]: ['Y', 'T'],
          }
        },
      tgt: [
        'PEPTIDE1{[Ac(1)].D.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].D.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].D.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].D.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].D.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.K.L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.K.L.T.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.P.L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.P.L.T.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.Y.[C(1)].G.[NH2]}$$$$V2.0',
        'PEPTIDE1{[Ac(1)].L.W.G.[F4COO].L.T.[C(1)].G.[NH2]}$$$$V2.0',
      ],
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`enumerator-${testName}`, async () => {
      const res = getPtEnumeratorHelm(testData.src, testData.params);
      expectArray(res, testData.tgt);
    });
  }
});
