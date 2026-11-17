import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {category, test, assure, expect, waitUntil, after} from '@datagrok-libraries/test/src/test';
import {Tutorial} from '@datagrok-libraries/tutorials/src/tutorial';
import {eda} from '../tracks/eda';

category('Tutorials App', () => {
  after(async () => {
    for (const v of Array.from(grok.shell.views)) {
      if (v.root.classList.contains(Tutorial.VIEW_CLASS))
        v.close();
    }
  });

  test('Launch app', async () => {
    const view: DG.View = await grok.functions.call('Tutorials:trackOverview');
    expect(view.root.classList.contains(Tutorial.VIEW_CLASS), true);
    assure.notNull(view.root.querySelector('.tutorials-card'));
    expect(document.querySelector('.tutorials-root, .tutorials-bar') == null, true);
  });

  test('Open a track by path', async () => {
    const view: DG.View = await grok.functions.call('Tutorials:trackOverview', {path: '/eda'});
    expect(view.root.querySelectorAll('.tutorials-track[data-name]').length, 1);
    expect(view.root.querySelector('.tutorials-track')!.getAttribute('data-name'), eda.name);
    expect(view.path, `${Tutorial.APP_PATH}/ExploratoryDataAnalysis`);
    expect(document.querySelector('.tutorials-bar') == null, true);
  });

  test('Open tutorial by path', async () => {
    const tutorial = eda.tutorials.find((t) => t.name === 'Scatter Plot')!;
    try {
      await grok.functions.call('Tutorials:trackOverview', {path: '/eda/ScatterPlot'});
      await waitUntil(() => grok.shell.v instanceof DG.TableView && grok.shell.v.path === tutorial.path,
        'Tutorial view path was not set', 5000);
      expect(tutorial.path, `${Tutorial.APP_PATH}/ExploratoryDataAnalysis/ScatterPlot`);
      expect(tutorial.url, window.location.origin + tutorial.path);
    } finally {
      tutorial.close();
    }
  });

  test('Starting a tutorial closes the running one', async () => {
    const first = eda.tutorials.find((t) => t.name === 'Scatter Plot')!;
    const second = eda.tutorials.find((t) => t.name === 'Embedded Viewers')!;
    try {
      await grok.functions.call('Tutorials:trackOverview', {path: '/eda/ScatterPlot'});
      await waitUntil(() => grok.shell.v instanceof DG.TableView && grok.shell.v.path === first.path,
        'Tutorial did not start', 5000);
      await grok.functions.call('Tutorials:trackOverview', {path: '/eda/EmbeddedViewers'});
      await waitUntil(() => grok.shell.v instanceof DG.TableView && grok.shell.v.path === second.path,
        'The second tutorial did not start', 5000);
      expect(first.closed, true);
      expect(document.querySelectorAll('.tutorials-bar').length, 1);
    } finally {
      first.close();
      second.close();
    }
  });
});
