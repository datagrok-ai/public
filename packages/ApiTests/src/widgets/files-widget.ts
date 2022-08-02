import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {before, category, expect, test} from '@datagrok-libraries/utils/src/test';


category('Widgets', () => {
  let testConnection: DG.DataConnection;
  let packageDataConnection: DG.DataConnection;

  before(async () => {
    testConnection = await grok.dapi.connections.filter('shortName = "Home"').first();
    packageDataConnection = await grok.dapi.connections.filter('shortName = "AppData"').first();
  });

  test('Files', async () => {
    const fw = DG.FilesWidget.create();
    expect(fw instanceof DG.FilesWidget, true);
    expect(fw.root instanceof HTMLElement, true);
    expect(fw.root.classList.contains('d4-tree-view-root'), true);
    expect(ui.fileBrowser() instanceof DG.FilesWidget, true)

  });
});
