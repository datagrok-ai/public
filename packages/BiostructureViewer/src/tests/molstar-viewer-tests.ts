import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category/*, expect*/, test} from '@datagrok-libraries/utils/src/test';
import { DockingApp } from '../apps/docking-app';
import { MolstarViewerApp } from '../apps/molstar-viewer-app';


category('Demo', () => {
  test('dockingDemo', async () => {
    const app = new DockingApp();
    app.init();
  });

  test('bioStructureDemo', async () => {
    const app = new MolstarViewerApp('molstarViewerApp');
    await app.init();
  });
});
