import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('ui: spaceSaveRow', () => {
  test('the provider is the spaceSaveRow function', async () => {
    expect(DG.Func.find({meta: {role: DG.FUNC_TYPES.SPACE_SAVE_ROW}}).length > 0, true,
      'a package registers the spaceSaveRow role (PowerPack)');
  });

  test('gives the Space row: a new project with no space stays in the user\'s projects', async () => {
    const row: any = await ui.spaceSaveRow(null, 'Kinase panel');
    expect(row instanceof DG.Widget, true);
    try {
      expect(row.root instanceof HTMLElement, true);
      expect(row.problem, null);
      expect(JSON.stringify(await row.resolveTarget()), JSON.stringify({spaceId: null, link: false, changed: false}));
    }
    finally {
      row?.detach();
    }
  });
});
