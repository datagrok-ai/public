import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category/*, expect*/, test} from '@datagrok-libraries/test/src/test';
import {NglViewerApp} from '../apps/ngl-viewer-app';

category('NglViewer', () => {
  test('NglViewerApp', async () => {
    const app = new NglViewerApp('NglViewerApp');
    await app.init();
  });
});
