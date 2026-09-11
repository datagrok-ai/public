import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, assure, expect, awaitCheck, after} from '@datagrok-libraries/test/src/test';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {eda} from '../tracks/eda';

category('Tutorials App', () => {
  after(async () => grok.shell.dockManager.close(document.querySelector('.tutorials-root') as HTMLElement));

  test('Launch app', async () => {
    await grok.functions.call('Tutorials:trackOverview');
    assure.notNull(document.querySelector('div.panel-content > div.tutorials-root'));
  });

  test('Open tutorial by path', async () => {
    const tutorial = eda.tutorials.find((t) => t.name === 'Scatter Plot')!;
    try {
      await grok.functions.call('Tutorials:trackOverview', {path: '/eda/ScatterPlot'});
      await awaitCheck(() => grok.shell.v instanceof DG.TableView && grok.shell.v.path === tutorial.path,
        'Tutorial view path was not set', 5000);
      expect(tutorial.path, `${Tutorial.APP_PATH}/ExploratoryDataAnalysis/ScatterPlot`);
      expect(tutorial.url, window.location.origin + tutorial.path);
    } finally {
      tutorial.close();
    }
  });
});
