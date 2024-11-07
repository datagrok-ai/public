import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {
  category, expect, test, awaitCheck, delay, testEvent, before, after
} from '@datagrok-libraries/utils/src/test';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {
  getHelmService, HelmServiceBase, HelmAux, HelmProps
} from '@datagrok-libraries/bio/src/viewers/helm-service';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';

import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {
  getUserLibSettings, setUserLibSettings,
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';

import {initHelmMainPackage} from './utils';

import {_package} from '../package-test';

category('HelmService', () => {
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;

  let helmSvc: HelmServiceBase;

  before(async () => {
    const logPrefix = `Helm: tests: HelmService: before`;
    _package.logger.debug(`${logPrefix}, start`);
    await initHelmMainPackage();
    try {
      monomerLibHelper = await getMonomerLibHelper();
      userLibSettings = await getUserLibSettings();

      await monomerLibHelper.loadMonomerLibForTests(); // load default libraries

      helmSvc = await getHelmService();
      await helmSvc.reset();

      _package.logger.debug(`${logPrefix}, end`);
    } finally {
      _package.logger.debug(`${logPrefix}, fina`);
    }
  });

  after(async () => {
    const logPrefix = `Helm: tests: HelmService: after`;
    _package.logger.debug(`${logPrefix}, start`);

    await helmSvc.reset();

    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
    _package.logger.debug(`${logPrefix}, end`);
  });

  test('helm', async () => {
    const logPrefix: string = `Helm: tests: HelmService: helm`;
    _package.logger.debug(`${logPrefix}, start`);
    const helmStr: string = 'PEPTIDE1{I.H.A.N.T.Thr_PO3H2.Aca.D-Tyr_Et}$$$$';
    const monomerLib = monomerLibHelper.getMonomerLib();
    let consumerId: number | null = null;
    let taskAux: HelmAux | null = null;
    const event = new Subject<void>();
    await testEvent(event,
      () => {
      },
      () => {
        const task: RenderTask<HelmProps, HelmAux> = {
          name: 'test1',
          props: {
            helm: helmStr, monomerLib: monomerLib,
            backColor: DG.Color.fromHtml('#FFFFFF'),
            width: 300, height: 100,
          },
          onAfterRender: (canvas: HTMLCanvasElement, aux: HelmAux) => {
            taskAux = aux;
            event.next();
          }
        };
        helmSvc.reset();
        consumerId = helmSvc.render(consumerId, task, 0);
      }, 7000);
    expect(taskAux !== null, true, 'Helm task aux is null');
    expect(taskAux!.mol.atoms.length > 0, true, 'Polymer mol atoms empty');
    expect(taskAux!.mol.atoms.length, 8);
    expect(taskAux!.mol.bonds.length, 7);
    expect(helmSvc.errorCount, 0, 'There was errors in NglGlService.');
    expect(consumerId !== null, true, 'consumerId not assigned');
    _package.logger.debug(`${logPrefix}, end`);
  }, {timeout: 10000});
});
