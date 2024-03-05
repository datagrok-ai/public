import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import wu from 'wu';

import {category, expect, test, testViewer} from '@datagrok-libraries/utils/src/test';

import {awaitGrid} from './utils';
import {WebLogoViewer} from '../viewers/web-logo-viewer';

import {_package} from '../package-test';

category('WebLogo-layout', () => {
  test('fasta', async () => {
    const df = await _package.files.readCsv('tests/filter_FASTA.csv');
    const col = df.getCol('fasta');
    await grok.data.detectSemanticTypes(df);
    const view = grok.shell.addTableView(df);
    const wlViewer = await df.plot.fromType('WebLogo',
      {sequenceColumnName: col.name}) as unknown as WebLogoViewer;
    view.dockManager.dock(wlViewer);
    await wlViewer.awaitRendered();
    await awaitGrid(view.grid);

    const viewLayout = view.saveLayout();
    const viewLayoutJsonStr = viewLayout.toJson();
    view.loadLayout(viewLayout);
    await wlViewer.awaitRendered();
    await awaitGrid(view.grid);

    const viewersA = wu(view.viewers).toArray();
    expect(viewersA.length, 2 /* Grid and WebLogo */);
    expect(viewersA.filter((f) => f.type === 'Grid').length, 1);
    expect(viewersA.filter((f) => f.type === 'WebLogo').length, 1);
  });
});
