/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import {expect, type Page} from '@playwright/test';
import {openTableFromFile} from '@datagrok-libraries/test/src/playwright/openers.js';
import {DatasetEntry, Given, Then, When} from '../../src/registry.js';
import {el, type ElementRef} from '../../src/runtime/args.js';
import {editorOf} from '../../src/runtime/gestures.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';
import {exactText} from '../../src/runtime/locate.js';

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

/** A bare `input[type="checkbox"]` in a Dart dialog — the select-all of "Order or Hide Columns"
 * and its like, which carry no class and no label for the `checkbox` kind to match on. */
export const clickPlainCheckbox = When('user clicks the plain checkbox in the {string} dialog',
  async (page: Page, title: string) => {
    const dialog = page.locator('.d4-dialog').filter({has: page.locator('.d4-dialog-title', {hasText: exactText(title)})}).last();
    await dialog.waitFor({state: 'visible', timeout: 5000});
    const box = dialog.locator('input[type="checkbox"]').filter({visible: true}).first();
    await expect(box, `a checkbox in the "${title}" dialog`).toBeVisible({timeout: 5000});
    await box.click();
  }, {tier: 'ui', description: 'the only checkbox of that dialog the library\'s kinds cannot name'});

/* --- the browse panel -------------------------------------------------------------------------
   A bdd page runs in simple mode (set by `user is logged in`), where the browse panel is not built
   at all — so a feature about the Browse tree opens it first. Opening it is a shell setting, never
   a click on the Browse tab: that tab TOGGLES the panel, and clicking it on an open one closes it
   again. Simple mode is restored when the feature ends, so the next feature on the same page finds
   the shell as it expects it. */

async function showBrowsePanel(page: Page, open: boolean): Promise<void> {
  await page.evaluate((show) => {
    if (show)
      grok.shell.windows.simpleMode = false;
    grok.shell.windows.showBrowse = show;
  }, open);
}

export const browsePanelOpen = Given('the browse panel is open', async (page: Page) => {
  await showBrowsePanel(page, true);
  await expect(page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]').first(),
    'the browse tree').toBeVisible({timeout: 60000});
  atFeatureEnd(page, () => page.evaluate(() => { grok.shell.windows.simpleMode = true; }));
}, {tier: 'api', description: 'idempotent: leaves simple mode, shows the panel and waits for its tree; puts simple mode back at feature end'});

export const openBrowsePanel = When('user opens the browse panel', async (page: Page) => {
  await showBrowsePanel(page, true);
  await expect(page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]').first(),
    'the browse tree').toBeVisible({timeout: 60000});
}, {tier: 'api'});

export const closeBrowsePanel = When('user closes the browse panel', (page: Page) => showBrowsePanel(page, false), {tier: 'api'});

export const browsePanelShouldBeOpen = Then('the browse panel should be open', async (page: Page) => {
  await expect(page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]').first(),
    'the browse tree').toBeVisible();
});

export const browsePanelShouldBeClosed = Then('the browse panel should be closed', async (page: Page) => {
  await expect(page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]')
    .filter({visible: true}), 'the browse tree').toHaveCount(0);
});

/* --- the second account ------------------------------------------------------------------------
   A sharing feature needs a user other than the one running it. It is DATAGROK_SHARING_LOGIN — the
   same variable the hand-written suites read from playwright-tests/.env — and the platform's user
   typeahead offers it under a name with the punctuation stripped ("a+b@x" shows as "ab"), so the
   step types the local part and picks the row rather than trusting what it typed. */

export function sharingLogin(): string {
  const login = process.env.DATAGROK_SHARING_LOGIN;
  if (!login)
    throw new Error('no DATAGROK_SHARING_LOGIN in the environment: a sharing feature needs a second account');
  return login;
}

export const pickSharingUser = When('user picks the sharing user in {element}', async (page: Page, target: ElementRef) => {
  const login = sharingLogin();
  const editor = await editorOf(page, el(target.phrase));
  await editor.click();
  await editor.pressSequentially(login.split('@')[0]);
  const wanted = login.split('@')[0].replace(/[^a-z0-9]/gi, '').toLowerCase();
  const row = page.locator('.d4-user-selector-drop-down tr, .d4-tags-selector-drop-down tr')
    .filter({hasText: new RegExp(wanted, 'i')}).first();
  await expect(row, `the "${login}" row of the user typeahead`).toBeVisible({timeout: 15000});
  await row.click();
}, {tier: 'ui', description: 'types the login of DATAGROK_SHARING_LOGIN and takes it from the typeahead'});

/* The address bar is part of what the platform promises: a view that puts its state in the URL can
   be shared by copying it. */
export const urlShouldContain = Then('the page address should contain {string}', async (page: Page, part: string) => {
  await expect.poll(() => page.url(), {message: 'the page address'}).toContain(part);
});

export const openAddress = When('user opens the page address of the current view', async (page: Page) => {
  await page.goto(page.url(), {waitUntil: 'domcontentloaded', timeout: 180000});
}, {tier: 'ui', description: 'loads the address again from scratch — what pasting the copied link into a new tab does'});

/** How many viewers the current view holds — an analysis that is done is one that has put its
 * viewers on screen. */
export const viewHoldsViewers = Then('the current view should hold at least {int} viewer(s)',
  async (page: Page, count: number) => {
    await expect.poll(() => page.evaluate(() => Array.from(grok.shell.v?.viewers ?? []).length),
      {message: 'viewers of the current view', timeout: 300000}).toBeGreaterThanOrEqual(count);
  }, {tier: 'api', description: 'polls for up to five minutes: an analysis or a fit builds them when its run ends'});

/** A table in the workspace and nothing else: no view, so a form that offers the open tables in a
 * choice gains the option without losing the focus of the view it lives in. Named after the file,
 * which is the name such a choice shows. */
export const loadTable = Given('the {string} file is loaded as a table', async (page: Page, path: string) => {
  const name = await page.evaluate(async (p) => {
    const df = await grok.dapi.files.readCsv(p);
    df.name = p.replace(/^.*\//, '').replace(/\.[^.]+$/, '');
    grok.shell.addTable(df);
    return String(df.name);
  }, path);
  await expect.poll(() => page.evaluate((n) => (grok.shell.tables ?? []).some((t: any) => t.name === n), name),
    {message: `"${name}" among the open tables`}).toBe(true);
}, {tier: 'api', description: 'a file on the stand into the workspace, without a view of its own'});
