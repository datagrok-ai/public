import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test} from '@datagrok-libraries/utils/src/test';
import {testApp} from './gui-utils';

category('Apps', () => {
  let TM: DG.DockNode;

  before(async () => {
    TM = Array.from(grok.shell.dockManager.rootNode.children)[0].dart.a.cx;
    grok.shell.windows.showContextPanel = false;
  });

  const apps = DG.Func.find({tags: ['app']});
  for (const app of apps) {
    test(app.friendlyName, async () => {
      //@ts-ignore
      await testApp(app, TM);
    });
  }

  after(async () => {
    grok.shell.windows.showContextPanel = true;
    grok.shell.windows.showToolbox = true;
  });
});
