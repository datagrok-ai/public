import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';
import {getAutoDockService, GridSize, IAutoDockService} from '@datagrok-libraries/bio/src/pdb/auto-dock-service';
import {BiostructureData} from '@datagrok-libraries/bio/src/pdb/types';

import {_package} from '../package-test';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';


category('AutoDock', () => {
  let adSvc: IAutoDockService | null;

  before(async () => {
    try {
      adSvc = await getAutoDockService();
      if (!adSvc.ready) {
        _package.logger.warning('AutoDock docker container is not ready, trying to start.');
        await adSvc.startDockerContainer();
        _package.logger.warning('AutoDock docker container successfully started .');
        if (!adSvc.ready) {
          _package.logger.warning('AutoDock docker container can not start, skip tests');
          adSvc = null;
        }
      }
    } catch (err: any) {
      const [errMsg, errStack] = errInfo(err);
      _package.logger.error(errMsg, undefined, errStack);
    }
  });

  test('1bdq', async () => {
    if (!adSvc) return;

    const receptorPdb = await _package.files.readAsText('samples/1bdq-wo-ligands.pdb');
    const ligandPdb = await _package.files.readAsText('samples/1bdq-ligand-0.by-babel.pdb');
    if (receptorPdb === '' && ligandPdb === '')
      throw new Error('Empty test data');

    const ligandData: BiostructureData = {binary: false, ext: 'pdb', data: ligandPdb};

    const npts = new GridSize(20, 20, 20);
    const adRes = await adSvc.run(receptorPdb, ligandData, npts);
    expect(adRes.posesDf.rowCount, 30);
  }, {timeout: 60000});
});
