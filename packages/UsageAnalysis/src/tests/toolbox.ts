import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, expect} from '@datagrok-libraries/test/src/test';
import {ViewHandler} from '../view-handler';


category('Toolbox', () => {
  test('cancelled filter queries', async () => {
    const cancelled = ['PackagesCategories', 'EntitiesTags', 'ProjectsList'];
    const sub = grok.functions.onBeforeRunAction.subscribe((fc: DG.FuncCall) => {
      if (cancelled.includes(fc.func.name))
        setTimeout(() => fc.cancel(), 30);
    });
    const handler = new ViewHandler();
    try {
      grok.shell.addView(handler.view);
      await handler.init();
      expect(handler.view.tabs != null, true, 'tabs not initialized');
    }
    finally {
      sub.unsubscribe();
      handler.view.close();
    }
  });
}, {clear: false, timeout: 60000});
