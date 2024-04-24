import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subject} from 'rxjs';

import {before, after, category, expect/*, expect*/, test, testEvent} from '@datagrok-libraries/utils/src/test';
import {
  NglGlServiceBase, getNglGlService, NglGlProps, NglGlAux
} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {RenderTask} from '@datagrok-libraries/bio/src/utils/cell-renderer-async-base';

import {_package} from '../package-test';

category('NglGlService', () => {
  let nglSvc: NglGlServiceBase;

  before(async () => {
    nglSvc = await getNglGlService();
    await nglSvc.reset();
  });

  after(async () => {
    await nglSvc.reset();
  });

  test('pdb', async () => {
    /** Tests rendering errors after NglGlService.reset() */
    _package.logger.debug('tests NglGlService/pdb, end');
    const pdbStr = await _package.files.readAsText('samples/1bdq.pdb');
    expect(pdbStr.length > 0, true, 'Empty test file.');
    const event = new Subject<void>();
    await testEvent(event,
      () => {
        _package.logger.debug('tests NglGlService/pdb, event handler');
      },
      () => {
        const task: RenderTask<NglGlProps, NglGlAux> = {
          name: 'test1',
          props: {
            pdb: pdbStr,
            backColor: DG.Color.fromHtml('#FFFFFF'),
            width: 300, height: 300,
          },
          onAfterRender: (canvas: HTMLCanvasElement) => {
            event.next();
          },
        };
        nglSvc.reset();
        nglSvc.render(task, 0);
      }, 15000);

    expect(nglSvc.errorCount, 0, 'There was errors in NglGlService.');
    _package.logger.debug('tests NglGlService/pdb, end');
  }, {timeout: 25000});
});
