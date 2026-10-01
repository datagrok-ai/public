import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {APP_NAME, MENU_PATH} from '../constants';
import {ForgeApp} from '../ui/forge-app';

category('UI', () => {
  test('Forge app opens', async () => {
    const view = await ForgeApp.create();
    grok.shell.addView(view);
    try {
      expect(view.name, APP_NAME);
      const text = view.root.textContent ?? '';
      expect(text.includes('XGBoost'), true, 'XGBoost is not listed');
      expect(text.includes('Models'), true, 'The Models section is missing');
      expect(view.models.col('name')?.meta.friendlyName, 'Name');
      expect(view.models.col('row_count')?.meta.friendlyName, 'Training rows');
      expect(view.models.col('created_on')?.meta.friendlyName, 'Created');
    } finally {
      view.close();
    }
  });

  test('app function is registered', async () => {
    const apps = DG.Func.find({package: 'Forge', tags: [DG.FUNC_TYPES.APP]});
    expect(apps.length, 1);
    expect(apps[0].friendlyName, APP_NAME);
  });

  test('menu function is registered', async () => {
    const topMenu = DG.Func.find({package: 'Forge', name: 'forgeModels'})[0]?.topMenu;
    expect(topMenu?.startsWith(MENU_PATH), true, `Unexpected top menu: ${topMenu}`);
  });
});
