import {category, expect, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

category('ui: spaceSaveRow', () => {
  // The row comes from a package (PowerPack). A stand without one, such as the CI stand that has the test
  // packages only, has nothing to test here.
  const provided = () => DG.Func.find({meta: {role: DG.FUNC_TYPES.SPACE_SAVE_ROW}}).length > 0;
  const noProvider = 'skipped: no package on this stand registers the spaceSaveRow role';

  test('the provider is the spaceSaveRow function', async () => {
    if (!(await grok.dapi.packages.filter('PowerPack').list({pageSize: 5})).length) {
      console.log(noProvider);
      return;
    }
    expect(provided(), true, 'PowerPack registers the spaceSaveRow role');
  });

  test('gives the Space row: a new project with no space stays in the user\'s projects', async () => {
    if (!provided()) {
      console.log(noProvider);
      return;
    }
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
