/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import {expect, type Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given, Then, When} from '../../src/registry.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';

declare const grok: any;
declare const DG: any;

/** Opening a table starts semantic-type detection in the background (package detectors, a few
 * hundred ms on a molecule table); the platform reports its end on the global event bus, and the
 * step is over only then — otherwise that work lands on whatever step comes next. */
async function openTable(page: Page, dataset: DatasetEntry, rows?: number, name?: string): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__bddDetected)
      return;
    w.__bddDetected = [];
    w.grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      w.__bddDetected = [...w.__bddDetected.slice(-19), a?.args?.dataFrame?.dart];
    });
  });
  if (rows === undefined) {
    await openTableFromFile(page, dataset.path);
  }
  else {
    // a subset is a clone: the file's rows are read the same way, the view gets the first N,
    // named as the file (readCsv names nothing) unless the step names it
    await page.evaluate(async ([p, r, n]) => {
      const src = await grok.dapi.files.readCsv(p);
      const df = r < src.rowCount ? src.clone(DG.BitSet.create(src.rowCount, (i: number) => i < r)) : src;
      df.name = n ?? p.replace(/^.*\//, '').replace(/\.[^.]+$/, '');
      grok.shell.addTableView(df);
    }, [dataset.path, rows, name ?? null] as [string, number, string | null]);
  }
  await page.locator('[name="viewer-Grid"]').first().waitFor();
  await page.waitForFunction(() => {
    const w = window as any;
    return w.__bddDetected.includes(w.grok.shell.tv?.dataFrame?.dart);
  }, undefined, {timeout: 15000}).catch(() => {
    throw new Error(`${dataset.name}: semantic types were not detected within 15 s (is auto-detection on?)`);
  });
  // a second after the grid is created the view makes row 0 current when no row is, and every
  // viewer repaints its marker mid-feature; done here, the view's timer skips it
  await page.evaluate(() => {
    const df = (window as any).grok.shell.tv?.dataFrame;
    if (df && df.currentRowIdx === -1 && df.rowCount > 0)
      df.currentRowIdx = 0;
  });
}

export const openDataset = Given('user opens {dataset} dataset', (page: Page, dataset: DatasetEntry) => openTable(page, dataset),
  {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI; done when semantic types are detected, with row 0 current as the view would make it a second later'});

export const openDatasetRows = Given('user opens {dataset} dataset keeping the first {int} rows', (page: Page, dataset: DatasetEntry, rows: number) =>
  openTable(page, dataset, rows), {tier: 'api', description: 'a clone with the first N rows, named as the file — for a feature whose commands cost a call per row'});

export const openDatasetRowsAs = Given('user opens {dataset} dataset keeping the first {int} rows as {string}', (page: Page, dataset: DatasetEntry, rows: number, name: string) =>
  openTable(page, dataset, rows, name), {tier: 'api', description: 'the same, with the table named — what "table {string} should …" and "switches to the {string} table view" then use'});

export const switchTableView = Given('user switches to (the ){string} table view', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const grok = (window as any).grok;
    const views = Array.from(grok.shell.tableViews) as any[];
    const view = views.find((x) => String(x.dataFrame?.name).toLowerCase() === n.toLowerCase());
    if (!view)
      throw new Error(`no table view for "${n}"; open: ${views.map((x) => x.dataFrame?.name).join(', ')}`);
    grok.shell.v = view;
  }, name);
  await page.waitForFunction((n) => (window as any).grok.shell.tv?.dataFrame?.name?.toLowerCase() === n.toLowerCase(), name);
}, {tier: 'api', description: 'by table name — the view of a dataset opened earlier in the scenario'});

export const switchView = Given('user switches to (the ){string} view', async (page: Page, name: string) => {
  await page.evaluate((n) => {
    const grok = (window as any).grok;
    const views = Array.from(grok.shell.views) as any[];
    const view = views.find((x) => String(x.name).toLowerCase() === n.toLowerCase());
    if (!view)
      throw new Error(`no "${n}" view; open: ${views.map((x) => x.name).join(', ')}`);
    grok.shell.v = view;
  }, name);
  await page.waitForFunction((n) => String((window as any).grok.shell.v?.name).toLowerCase() === n.toLowerCase(), name);
}, {tier: 'api', description: 'any open view by its name (an app view, Home) — the view tabs are hidden in the simple mode a bdd page runs in, so a click on one is not a step'});

export const closeAllViews = When('user closes all views', async (page: Page) => {
  await page.evaluate(() => { grok.shell.closeAll(); });
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok');
}, {tier: 'api', description: 'grok.shell.closeAll — tables, views and viewers gone, the Home view current'});

/** The project's table is uploaded and the view saved with its layout, as the ribbon's Save
 * dialog does it; a plain view info would drop the viewport. Deleted when the feature ends. */
export const saveAsProject = When('user saves the current view as project {string}', async (page: Page, name: string) => {
  const ids: {project: string; table: string} = await page.evaluate(async (n) => {
    const tv = grok.shell.tv;
    const project = DG.Project.create();
    project.name = n;
    const tableInfo = tv.dataFrame.getTableInfo();
    const viewInfo = tv.getInfo();
    viewInfo.viewState = tv.saveLayout({saveWithData: true}).viewState;
    project.addChild(tableInfo);
    project.addChild(viewInfo);
    await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
    await grok.dapi.tables.save(tableInfo);
    await grok.dapi.views.save(viewInfo);
    await grok.dapi.projects.save(project);
    const w = window as any;
    w.__bddProjects = {...(w.__bddProjects ?? {}), [n]: String(project.id)};
    return {project: String(project.id), table: String(tableInfo.id)};
  }, name);
  atFeatureEnd(page, () => page.evaluate(async (i) => {
    for (const [source, id] of [[grok.dapi.projects, i.project], [grok.dapi.tables, i.table]]) {
      const e = await source.find(id).catch(() => null);
      if (e)
        await source.delete(e);
    }
  }, ids));
}, {tier: 'api', description: 'the current table view as a project on the server, removed again when the feature ends'});

export const openProject = When('user opens the {string} project', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const id = ((window as any).__bddProjects ?? {})[n];
    const p = id ? await grok.dapi.projects.find(id) : await grok.dapi.projects.filter(`name = "${n}"`).first();
    if (!p)
      throw new Error(`no project "${n}" on the server`);
    await p.open();
  }, name);
  await page.waitForFunction(() => grok.shell.tv?.dataFrame != null);
}, {tier: 'api', description: 'the project this feature saved under that name, else the server\'s by name; done when a table view is current'});

export const viewIsCurrent = Then('the {string} view should be current', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => String((window as any).grok.shell.v?.name ?? '')), {message: 'the current view'}).toBe(name);
}, {description: 'grok.shell.v by name — what a command that opens a view leaves in front'});

export const closeCurrentView = When('user closes the current view', async (page: Page) => {
  await page.evaluate(() => { (window as any).grok.shell.v?.close(); });
}, {tier: 'api'});

export const openApp = Given('user opens the {string} app', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const DG = (window as any).DG;
    const grok = (window as any).grok;
    const apps = DG.Func.find({tags: ['app']});
    const app = apps.find((f: any) => f.friendlyName === n) ?? apps.find((f: any) => f.name === n);
    if (!app)
      throw new Error(`no "${n}" app is registered; the apps: ${apps.map((f: any) => f.friendlyName).sort().join(', ')}`);
    const view = await app.apply({});
    if (view?.root && !Array.from(grok.shell.views).some((v: any) => v.dart === view.dart))
      grok.shell.addView(view);
  }, name);
  await expect.poll(() => page.evaluate(() => String((window as any).grok.shell.v?.name ?? '')), {message: 'the current view'}).toBe(name);
}, {tier: 'api', description: 'runs the app function by name (the way the browse tree does) and shows the view it returns; done when that view is current'});

export const autostartsCompleted = Given('the package autostarts have completed', async (page: Page) => {
  await page.evaluate(async () => { await grok.shell.autostartsCompleted; });
}, {tier: 'api', description: 'grok.shell.autostartsCompleted — a viewer a package registers is not there before it'});
