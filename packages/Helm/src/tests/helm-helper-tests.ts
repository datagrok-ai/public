import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, expectArray, test, testEvent} from '@datagrok-libraries/utils/src/test';
import {HelmNotSupportedError, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {_package} from '../package-test';

type TestSrcType = { helm: string };
type TestTgtType = { helm: string | null, map: [number, number][] };

category('HelmHelper: removeGaps', () => {
  let helmHelper: IHelmHelper;

  before(async () => {
    helmHelper = await getHelmHelper();
  });

  const tests: { [testName: string]: { src: TestSrcType, tgt: TestTgtType } } = {
    'single-linear-start': {
      src: {helm: 'PEPTIDE1{*.[meY].A.G.T}$$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].A.G.T}$$$$V2.0',
        map: [[1, 0], [2, 1], [3, 2], [4, 3]],
      }
    },
    'single-linear-end': {
      src: {helm: 'PEPTIDE1{[meY].A.G.T.*}$$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].A.G.T}$$$$V2.0',
        map: [[0, 0], [1, 1], [2, 2], [3, 3]],
      }
    },
    'single-linear-middle': {
      src: {helm: 'PEPTIDE1{[meY].A.*.G.T}$$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].A.G.T}$$$$V2.0',
        map: [[0, 0], [1, 1], [3, 2], [4, 3]],
      }
    },
    'single-cycle-start': {
      src: {helm: 'PEPTIDE1{*.[meY].C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,3:R3-7:R3$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-6:R3$$$V2.0',
        map: [[1, 0], [2, 1], [3, 2], [4, 3], [5, 4], [6, 5], [7, 6]],
      }
    },
    'single-cycle-end': {
      src: {helm: 'PEPTIDE1{[meY].C.R.N.P.C.T.*}$PEPTIDE1,PEPTIDE1,2:R3-6:R3$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-6:R3$$$V2.0',
        map: [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4], [5, 5], [6, 6]],
      }
    },
    'single-cycle-middle': {
      src: {helm: 'PEPTIDE1{[meY].C.R.N.*.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$V2.0'},
      tgt: {
        helm: 'PEPTIDE1{[meY].C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-6:R3$$$V2.0',
        map: [[0, 0], [1, 1], [2, 2], [3, 3], [5, 4], [6, 5], [7, 6]],
      }
    },
    'single-cycle-gap-at-connection': {
      src: {helm: 'PEPTIDE1{[meY].*.C.R.N.P.C.T}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$V2.0'},
      tgt: {helm: null, map: []}
    }
  };

  for (const [testName, testData] of Object.entries(tests)) {
    test(`${testName}`, async () => {
      let resErr: any = null;
      try {
        const res = helmHelper.removeGaps(testData.src.helm);
        expect(res.resHelm, testData.tgt.helm);
        expectArray(Array.from(res.monomerMap!.entries()), testData.tgt.map);
      } catch (err: any) {
        resErr = err;
      }

      if (testData.tgt.helm === null) { // err expected, for debug on GitHub CI
        expect(resErr != null, true, 'Error expected');
        const isErrorInstance = resErr instanceof HelmNotSupportedError;
        const errCtorName = resErr.constructor.name;
        const isErrorCtor = errCtorName === 'HelmNotSupportedError';
        _package.logger.debug(`Check error object. ` +
          `isErrorInstance: ${isErrorInstance}, isErrorCtor: ${isErrorCtor}, errCtorName: ${errCtorName}`);
      }

      expect((resErr instanceof HelmNotSupportedError) || resErr?.constructor.name === 'HelmNotSupportedError',
        testData.tgt.helm === null, 'HelmNotSupportedError thrown expected');
    }, testData.tgt.helm == null ? {skipReason: 'GitHub CI does not support testing exceptions'} : undefined);
  }
});
