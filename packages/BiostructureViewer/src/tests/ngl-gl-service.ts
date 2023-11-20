import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {before, after, category, expect/*, expect*/, test, testEvent} from '@datagrok-libraries/utils/src/test';
import {NglGlServiceBase, NglGlTask} from '@datagrok-libraries/bio/src/viewers/ngl-gl-service';
import {getNglGlService} from '../package';
import {Subject} from 'rxjs';
import {_package} from '../package-test';

category('NglGlService', () => {
  let svc: NglGlServiceBase;

  before(async () => {
    svc = await getNglGlService();
  });

  after(async () => {
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
        const task: NglGlTask = {
          name: 'test1',
          backColor: DG.Color.fromHtml('#FFFFFF'),
          props: {
            pdb: pdbStr,
            width: 300, height: 300,
          },
          onAfterRender: (canvas: HTMLCanvasElement) => {
            event.next();
          },
        };
        svc.reset();
        svc.render(task, 0);
      }, 10000);

    expect(svc.errorCount, 0, 'There was errors in NglGlService.');
    _package.logger.debug('tests NglGlService/pdb, end');
  }, {timeout: 20000});
});
