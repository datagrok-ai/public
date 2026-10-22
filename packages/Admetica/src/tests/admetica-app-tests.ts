import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, category, delay, expect, test} from '@datagrok-libraries/test/src/test';

import {AdmeticaViewApp} from '../utils/admetica-app';

/** The sketch mode debounces sketcher changes by a second before it rebuilds the table. */
const SKETCH_DEBOUNCE_MS = 1000;

async function openApp(): Promise<{app: AdmeticaViewApp, view: DG.TableView}> {
  const app = new AdmeticaViewApp(null);
  app.addPath = false;
  app.persistSketcherValue = false;
  await app.init();
  const view = grok.shell.addView(app.tableView!) as DG.TableView;
  await awaitCheck(() => app.formContainer.childElementCount > 0 &&
    app.formContainer.querySelector('.demo-splash-container') == null, 'Sketch form was not shown', 30000);
  return {app, view};
}

function clickMode(app: AdmeticaViewApp, mode: string): void {
  const headers = Array.from(app.modeContainer.querySelectorAll<HTMLElement>('.d4-tab-header'));
  headers.find((h) => h.textContent?.trim().startsWith(mode))!.click();
}

async function awaitFileTable(app: AdmeticaViewApp): Promise<DG.DataFrame> {
  await awaitCheck(() => app.tableView!.dataFrame.rowCount > 1, 'File table was not loaded', 30000);
  return app.tableView!.dataFrame;
}

category('Admetica: App', () => {
  test('File mode keeps its table when a sketch change was pending', async () => {
    const {app, view} = await openApp();
    try {
      app.sketcherInstance!.onChanged.next();
      clickMode(app, 'File');
      const fileTable = await awaitFileTable(app);
      await delay(SKETCH_DEBOUNCE_MS + 500);
      expect(app.tableView!.dataFrame.dart === fileTable.dart, true);
    } finally {
      view.close();
    }
  }, {timeout: 90000});

  test('File mode ignores sketcher changes', async () => {
    const {app, view} = await openApp();
    try {
      clickMode(app, 'File');
      const fileTable = await awaitFileTable(app);
      app.sketcherInstance!.onChanged.next();
      await delay(SKETCH_DEBOUNCE_MS + 500);
      expect(app.tableView!.dataFrame.dart === fileTable.dart, true);
    } finally {
      view.close();
    }
  }, {timeout: 90000});
});
