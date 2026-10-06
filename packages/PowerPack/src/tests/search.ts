import * as DG from 'datagrok-api/dg';
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {initTemplates} from '../search/templates-search';
import {appSearch, exactAppFuncSearch} from '../search/entity-search';

category('Search', () => {
  test('InitTemplates', async () => await initTemplates());

  test('AppSearchIgnoresSpaces', async () => {
    const app = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}}).find((f) => f.friendlyName?.includes(' '));
    if (app == null)
      return;
    const found = await appSearch(app.friendlyName.replaceAll(' ', ''));
    expect(found.some((f) => f.id === app.id), true);
  });

  test('ExactAppSearchIgnoresSpaces', async () => {
    const app = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, returnType: 'view'})
      .find((f) => f.friendlyName?.includes(' '));
    if (app == null)
      return;
    expect(exactAppFuncSearch(app.friendlyName.replaceAll(' ', ''))?.id, app.id);
  });
});
