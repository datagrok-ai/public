/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import type {Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given, When} from '../../src/registry.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';

declare const grok: any;
declare const DG: any;

/** Opening a table starts semantic-type detection in the background (package detectors, a few
 * hundred ms on a molecule table); the platform reports its end on the global event bus, and the
 * step is over only then — otherwise that work lands on whatever step comes next. */
export const openDataset = Given('user opens {dataset} dataset', async (page: Page, dataset: DatasetEntry) => {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__bddDetected)
      return;
    w.__bddDetected = [];
    w.grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      w.__bddDetected = [...w.__bddDetected.slice(-19), a?.args?.dataFrame?.dart];
    });
  });
  await openTableFromFile(page, dataset.path);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
  await page.waitForFunction(() => {
    const w = window as any;
    return w.__bddDetected.includes(w.grok.shell.tv?.dataFrame?.dart);
  }, undefined, {timeout: 15000}).catch(() => {
    throw new Error(`${dataset.name}: semantic types were not detected within 15 s (is auto-detection on?)`);
  });
}, {tier: 'api', description: 'OpenFile through the JS API — provenance as in the UI; done when semantic types are detected'});

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
