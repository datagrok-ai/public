/* Platform base steps: setup through the JS API (the openers of @datagrok-libraries/test keep the
   provenance tags the UI would set). Viewer steps live in the `viewers` tier. */
import {type Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {DatasetEntry, Given, Then, When} from '../../src/registry.js';
import {el, type ElementRef} from '../../src/runtime/args.js';
import {click, editorOf} from '../../src/runtime/gestures.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';
import {shellSimpleMode, silent} from '../../src/runtime/guide.js';
import {exactText, locate} from '../../src/runtime/locate.js';
import {armEvent} from '../../src/runtime/viewer-menus.js';
import {deleteChatsOf, serverRequests} from '../../src/runtime/server.js';

declare const grok: any;
declare const DG: any;

/** A file is read once per page and every feature gets a clone of it — the clone carries the
 * semantic types the first detection found, so the platform's detection on it skips the typed
 * columns — named as the file (readCsv names nothing) unless the step names it; a subset keeps
 * the first N rows. Opening a table starts semantic-type detection in the background (package
 * detectors, a few hundred ms on a molecule table); the platform reports its end on the global
 * event bus, and the step is over only then — otherwise that work lands on whatever step comes
 * next. */
async function openTable(page: Page, dataset: DatasetEntry, rows?: number, name?: string): Promise<void> {
  await page.evaluate(async ([p, r, n]) => {
    const w = window as any;
    if (!w.__bddDetected) {
      w.__bddDetected = [];
      w.grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
        w.__bddDetected = [...w.__bddDetected.slice(-19), a?.args?.dataFrame?.dart];
      });
    }
    w.__bddTables ??= {};
    // a file that is not a csv goes through the file handler its extension is registered for (sdf → Chem);
    // a .csv the comma parser reads as one tab-separated column goes through it too
    const read = async (path: string) => {
      if (!/\.(csv|tsv|txt)$/i.test(path))
        return grok.data.files.openTable(path);
      const df = await grok.dapi.files.readCsv(path);
      return df.columns.length === 1 && df.columns.names()[0].includes('\t') ? grok.data.files.openTable(path) : df;
    };
    const src = (w.__bddTables[p] ??= await read(p));
    const df = r !== null && r < src.rowCount ? src.clone(DG.BitSet.create(src.rowCount, (i: number) => i < r)) : src.clone();
    df.name = n ?? p.replace(/^.*\//, '').replace(/\.[^.]+$/, '');
    grok.shell.addTableView(df);
  }, [dataset.path, rows ?? null, name ?? null] as [string, number | null, string | null]);
  await page.locator('[name="viewer-Grid"]').first().waitFor();
  // resolved the moment the detection event lands, not on the next poll
  await page.evaluate(([timeout, what, handled]) => new Promise<void>((resolve, reject) => {
    const w = window as any;
    const df = w.grok.shell.tv?.dataFrame;
    const dart = df?.dart;
    // a file handler that types its own column (an sdf, a mol) fires no detection event, so such a
    // table that already carries a semantic type is as ready as one detection has run over
    if (w.__bddDetected.includes(dart) || (handled && (df?.columns?.toList() ?? []).some((c: any) => c.semType)))
      return resolve();
    const timer = setTimeout(() => { sub.unsubscribe(); reject(new Error(`${what}: semantic types were not detected (is auto-detection on?)`)); }, timeout);
    const sub = w.grok.events.onEvent('ddt-semantic-type-detected').subscribe((a: any) => {
      if (a?.args?.dataFrame?.dart !== dart)
        return;
      clearTimeout(timer);
      sub.unsubscribe();
      resolve();
    });
  }), [pollMs(60000), dataset.name, !/\.(csv|tsv|txt)$/i.test(dataset.path)] as [number, string, boolean]);
  // a second after the grid is created the view makes row 0 of its first column current when no row
  // is, and every viewer repaints its marker mid-feature; done here, the same cell, the view's timer
  // skips it — and a "current column" claim does not depend on which of the two got there first
  await page.evaluate(() => {
    const tv = (window as any).grok.shell.tv;
    const df = tv?.dataFrame;
    if (!df || df.currentRowIdx !== -1 || df.rowCount === 0)
      return;
    const first = tv.grid?.columns.byIndex(1)?.column ?? df.columns.byIndex(0);
    df.currentCell = df.cell(0, first.name);
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

/** A project of this name, or of this family and older than an hour, is what a run that never
 * reached its feature end left behind: it goes, with the table and the view it holds. */
async function deleteLeftoverProjects(page: Page, name: string): Promise<void> {
  const families = fixtureFamilies([name]);
  const leftovers = (await serverEntities(page, 'projects'))
    .filter((project) => [project.name, project.friendlyName].includes(name) || isStaleFixture(project, families));
  for (const leftover of leftovers)
    await page.evaluate(async (id) => {
      // find types the children (TableInfo, ViewInfo); a listing's include('children') leaves them plain entities
      const project = await grok.dapi.projects.find(id);
      if (!project)
        return;
      for (const child of project.children) {
        const source = child instanceof DG.TableInfo ? grok.dapi.tables : child instanceof DG.ViewInfo ? grok.dapi.views : null;
        if (source)
          await source.delete(child);
      }
      await grok.dapi.projects.delete(project);
    }, leftover.id);
}

/** The project's tables are uploaded and every view saved with its layout, as the ribbon's Save
 * dialog does it; a plain view info would drop the viewport. Deleted when the feature ends, and
 * whatever an earlier run left under the name goes first. */
async function saveProject(page: Page, name: string, everyView: boolean): Promise<void> {
  await deleteLeftoverProjects(page, name);
  const ids: {project: string; tables: string[]; views: string[]} = await page.evaluate(async ([n, every]) => {
    const project = DG.Project.create();
    project.name = n;
    const tables: string[] = [];
    const views: string[] = [];
    for (const tv of every ? Array.from(grok.shell.tableViews) as any[] : [grok.shell.tv]) {
      const tableInfo = tv.dataFrame.getTableInfo();
      const viewInfo = DG.ViewInfo.fromJson(tv.saveLayout({saveWithData: true}).toJson());
      project.addChild(tableInfo);
      project.addChild(viewInfo);
      await grok.dapi.tables.uploadDataFrame(tv.dataFrame);
      await grok.dapi.tables.save(tableInfo);
      await grok.dapi.views.save(viewInfo);
      tables.push(String(tableInfo.id));
      views.push(String(viewInfo.id));
    }
    await grok.dapi.projects.save(project);
    const w = window as any;
    w.__bddProjects = {...(w.__bddProjects ?? {}), [n]: String(project.id)};
    return {project: String(project.id), tables, views};
  }, [name, everyView] as [string, boolean]);
  atFeatureEnd(page, () => page.evaluate(async (i) => {
    for (const [source, entityIds] of [[grok.dapi.projects, [i.project]], [grok.dapi.views, i.views], [grok.dapi.tables, i.tables]] as [any, string[]][]) {
      for (const id of entityIds) {
        const e = await source.find(id).catch(() => null);
        if (e)
          await source.delete(e);
      }
    }
  }, ids));
}

export const noProjectOnServer = Given('no project named {string} is on the server', async (page: Page, name: string) => {
  silent(page);
  const cleanup = () => deleteLeftoverProjects(page, name);
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'deletes the project an earlier run left under that name (with its table and view), and again when the feature ends — for a save made through the Save dialog'});

export const projectsOnServer = Then('{int} project(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'projects', 'projects', name, count),
{tier: 'api', description: 'what the server holds, not what the dialog said'});

export const saveAsProject = When('user saves the current view as project {string}', (page: Page, name: string) =>
  saveProject(page, name, false),
{tier: 'api', description: 'the current table view as a project on the server, removed again when the feature ends'});

export const saveAllAsProject = When('user saves all open table views as project {string}', (page: Page, name: string) =>
  saveProject(page, name, true),
{tier: 'api', description: 'every open table view in one project, as the ribbon\'s Save does — for a project that has to hold more than the current table'});

export const openProject = When('user opens the {string} project', async (page: Page, name: string) => {
  await page.evaluate(async (n) => {
    const id = ((window as any).__bddProjects ?? {})[n];
    const p = id ? await grok.dapi.projects.find(id) :
      await grok.dapi.projects.filter(`friendlyName = "${n}" or name = "${n}"`).first();
    if (!p)
      throw new Error(`no project "${n}" on the server`);
    await p.open();
  }, name);
  await page.waitForFunction(() => grok.shell.tv?.dataFrame != null);
}, {tier: 'api', description: 'the project this feature saved under that name, else the server\'s by name or friendly name; done when a table view is current'});

/** The same, awaited to the rows: a project whose table the save uploaded opens its view first and
 * fills it afterwards, which on a loaded stand is well past the budget of the step above. */
export const openProjectWithTable = When('user opens the {string} project and waits for its table', async (page: Page, name: string) => {
  await openProject(page, name);
  await expect.poll(() => page.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? -1),
    {message: `the rows of the table the "${name}" project opened`, timeout: pollMs(60000)}).toBeGreaterThan(0);
}, {tier: 'api', description: 'opens the project and is done when its table view holds rows'});

/** How the server keeps a table of a project: "sync: <creation script>" when opening the project
 * re-runs the script (the Save dialog's Data sync), "snapshot" when it loads the uploaded data. */
async function savedTableMode(page: Page, project: string, table: string): Promise<string> {
  return page.evaluate(async ([p, t]) => {
    const listed = await grok.dapi.projects.filter(`friendlyName = "${p}" or name = "${p}"`).first();
    if (!listed)
      return `no project "${p}" on the server`;
    const found = await grok.dapi.projects.find(listed.id);
    const infos = found.children.filter((c: any) => c instanceof DG.TableInfo);
    const info = infos.find((c: any) => c.friendlyName === t || c.name === t);
    if (!info)
      return `the project holds no "${t}" table; it holds: ${infos.map((c: any) => c.friendlyName).join(', ') || 'none'}`;
    const tags = (await grok.dapi.tables.find(info.id)).tags;
    return tags['.data-sync'] === 'sync' ? `sync: ${tags['.script'] ?? ''}` : 'snapshot';
  }, [project, table]);
}

export const savedWithDataSync = Then('the {string} table of the {string} project should be saved with data sync', async (page: Page, table: string, project: string) => {
  await expect.poll(() => savedTableMode(page, project, table),
    {message: `how the server keeps the "${table}" table of the "${project}" project`, timeout: pollMs(30000)}).toMatch(/^sync: \S/);
}, {tier: 'api', description: 'the table carries the data-sync flag and a creation script on the server, so opening the project re-runs it'});

export const savedAsSnapshot = Then('the {string} table of the {string} project should be saved as a snapshot', async (page: Page, table: string, project: string) => {
  await expect.poll(() => savedTableMode(page, project, table),
    {message: `how the server keeps the "${table}" table of the "${project}" project`, timeout: pollMs(30000)}).toBe('snapshot');
}, {tier: 'api', description: 'the table has no data-sync flag on the server, so opening the project loads the uploaded data'});

/** `TableInfo.execDataSync` marks the frame it rebuilt from the creation script; a frame loaded from
 * the uploaded data carries no mark. */
export const reloadedByDataSync = Then('the table should have been reloaded by data sync', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => grok.shell.tv?.dataFrame?.getTag('.data-sync') ?? 'no mark'),
    {message: 'the data-sync mark of the current table', timeout: pollMs(30000)}).toBe('success');
});

export const loadedAsSnapshot = Then('the table should have been loaded as a snapshot', async (page: Page) => {
  expect(await page.evaluate(() => grok.shell.tv?.dataFrame?.getTag('.data-sync') ?? 'no mark'),
    'the data-sync mark of the current table').toBe('no mark');
});

/** The kind of view in front, when its name does not tell them apart: a query editor
 * (DataQueryView) and the table view its Run leaves behind carry the same name. */
export const currentViewType = Then('the current view should be a {word} view', async (page: Page, type: string) => {
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.type ?? '')),
    {message: 'the type of the current view'}).toBe(type);
}, {description: 'grok.shell.v.type'});

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

export const browsePanelOpen = Given('the browse panel is open', async (page: Page) => {
  const shown = await page.evaluate(() => {
    const was = grok.shell.windows.showBrowse;
    grok.shell.windows.simpleMode = false;
    grok.shell.windows.showBrowse = true;
    return was;
  });
  await expect(page.locator('.grok-view-browse [role="tree"], .layout-browse [role="tree"]').first(), 'the browse tree').toBeVisible({timeout: 60000});
  // showBrowse is a user setting: left on, it opens the panel in every later page of the account
  atFeatureEnd(page, () => page.evaluate(([simple, browse]) => {
    grok.shell.windows.showBrowse = browse;
    grok.shell.windows.simpleMode = simple;
  }, [shellSimpleMode(), shown] as const));
}, {tier: 'api', description: 'idempotent: leaves simple mode, shows the panel and waits for its tree; puts the panel and simple mode back at feature end'});

export const toolboxPaneShown = Given('the toolbox pane is shown', async (page: Page) => {
  await page.evaluate(() => {
    grok.shell.windows.simpleMode = false;
    // a toolbox restored open by the user's layout and hidden at startup while empty still reads
    // as shown, and the setter ignores a value it already has: off, then on, docks it again
    grok.shell.windows.showToolbox = false;
    grok.shell.windows.showToolbox = true;
  });
  await expect(page.locator('.d4-toolbox[caption]').first(), 'the toolbox pane').toBeVisible({timeout: 15000});
  atFeatureEnd(page, () => page.evaluate((simple) => { grok.shell.windows.showToolbox = false; grok.shell.windows.simpleMode = simple; }, shellSimpleMode()));
}, {tier: 'api', description: 'idempotent: leaves simple mode and docks the toolbox pane afresh (off by default for a user, hidden at startup while empty); puts both back at feature end'});

/* Which sketcher a molecule input, a filter card or a dialog opens is the account's choice, kept on
   the server: a feature that draws or types a molecule names the one it was written against, so an
   account that picked another one elsewhere does not change what the feature sees. */
export const sketcherIs = Given('the molecule sketcher is {string}', async (page: Page, name: string) => {
  silent(page);
  const was = await page.evaluate((n) => {
    const known = DG.Func.find({meta: {role: 'moleculeSketcher'}}).map((f: any) => f.friendlyName);
    if (!known.includes(n))
      throw new Error(`no molecule sketcher "${n}"; the stand has: ${known.join(', ')}`);
    const before = grok.userSettings.getValue(DG.chem.STORAGE_NAME, DG.chem.KEY) ?? null;
    grok.userSettings.add(DG.chem.STORAGE_NAME, DG.chem.KEY, n);
    DG.chem.currentSketcherType = n;
    return before;
  }, name);
  atFeatureEnd(page, () => page.evaluate((b) => {
    if (b === null)
      grok.userSettings.delete(DG.chem.STORAGE_NAME, DG.chem.KEY);
    else
      grok.userSettings.add(DG.chem.STORAGE_NAME, DG.chem.KEY, b);
    DG.chem.currentSketcherType = b ?? DG.DEFAULT_SKETCHER;
  }, was));
}, {tier: 'api', description: 'the sketcher every molecule editor opens from then on (OpenChemLib is the platform\'s default); the account\'s own choice comes back at feature end; not in the video'});

/** Every guide's second step (the compiler insists): the shell as a person has it, view tabs and
 * menu bar included, in a plain run as much as in a filmed one. Silent, like the login. */
export const simpleModeOff = Given('simple mode is off', async (page: Page) => {
  silent(page);
  await page.evaluate(() => { grok.shell.windows.simpleMode = false; });
  await expect(page.locator('.d4-view-handle, [name^="view-handle: "]').first(), 'a view tab').toBeVisible({timeout: 15000});
  atFeatureEnd(page, () => page.evaluate((simple) => { grok.shell.windows.simpleMode = simple; }, shellSimpleMode()));
}, {tier: 'api', description: 'the full shell — view tabs, menu bar, panels — for the feature; not in the video; simple mode is back at feature end'});

/* --- the context panel -------------------------------------------------------------------------
   The panel renders the current object (`grok.shell.o`) and nothing else: a click that did not
   change it — the setter drops a change to the object already current, one within 2 s of a
   property edit, one while the object is frozen — leaves the panel as it was. So a claim about
   the panel names the object first, and reads the panel only once that is the current one. */

/* The setter is ignored while the shell is rearranging the side panels (a view opening, the toolbox
   taking the panel), and the panel then stays hidden with `showContextPanel` already true — so the
   flag is set again until the panel is really up. */
export const contextPanelOpen = Given('the context panel is open', async (page: Page) => {
  await expect.poll(async () => page.evaluate(() => {
    const shown = document.querySelector('.grok-prop-panel') as HTMLElement | null;
    if (shown?.offsetParent != null)
      return 'shown';
    grok.shell.windows.showContextPanel = true;
    return shown ? 'in the DOM, not shown' : 'not in the DOM';
  }), {message: 'the context panel', timeout: pollMs(30000)}).toBe('shown');
}, {tier: 'api', description: 'idempotent: the shell setting, re-applied until the panel is shown — not a click on its toggle'});

export const contextPanelShows = Then('the context panel should show {string}', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => {
    const o = grok.shell.o;
    return o == null ? 'nothing' : `${o.constructor?.name ?? typeof o} "${o.friendlyName ?? o.name ?? ''}"`;
  }), {message: `the current object (grok.shell.o), which the context panel renders`}).toMatch(new RegExp(`"${name.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}"$`));
  await expect(page.locator('.grok-prop-panel'), 'the context panel').toContainText(name);
}, {description: 'the current object (grok.shell.o) is the entity of that name, and the panel shows it'});

/* --- the second account ------------------------------------------------------------------------
   A sharing feature needs a user other than the one running it: DATAGROK_SHARING_LOGIN — the same
   variable the hand-written suites read from playwright-tests/.env — or, when it is unset, the
   "bddsecond" user the global setup creates on the stand with the dev key. The platform's user
   typeahead offers a user under a name with the punctuation stripped ("a+b@x" shows as "ab"), so
   the step types the local part and picks the row rather than trusting what it typed. */

export function sharingLogin(): string {
  const login = process.env.DATAGROK_SHARING_LOGIN;
  if (!login)
    throw new Error('no second account: set DATAGROK_SHARING_LOGIN, or run with a dev key so the setup can create one');
  return login;
}

/** The second account as the platform shows it: the local part of the login, punctuation stripped. */
export function sharingShownName(): string {
  return sharingLogin().split('@')[0].replace(/[^a-z0-9]/gi, '');
}

export const pickSharingUser = When('user picks the sharing user in {element}', async (page: Page, target: ElementRef) => {
  const login = sharingLogin();
  const editor = await editorOf(page, el(target.phrase));
  await editor.click();
  await editor.pressSequentially(login.split('@')[0]);
  const row = page.locator('.d4-user-selector-drop-down tr, .d4-tags-selector-drop-down tr')
    .filter({hasText: new RegExp(sharingShownName(), 'i')}).first();
  await expect(row, `the "${login}" row of the user typeahead`).toBeVisible({timeout: pollMs(15000)});
  await row.click();
}, {tier: 'ui', description: 'types the login of DATAGROK_SHARING_LOGIN and takes it from the typeahead'});

/* The grant row of the second account in a Share dialog; its Remove button shows while the row is
   hovered. The row goes at once, the grant only when the dialog is confirmed. */
export const removeSharingUser = When('user removes the sharing user from {element}', async (page: Page, target: ElementRef) => {
  const shown = sharingShownName();
  const row = (await locate(page, target)).locator('[name^="div-permissions-row-"]').filter({hasText: new RegExp(shown, 'i')}).first();
  await expect(row, `the grant row of "${shown}"`).toBeVisible({timeout: pollMs(15000)});
  await row.hover();
  await row.locator('[name="button-Remove"]').click();
  await expect(row, `the grant row of "${shown}" after Remove`).toBeHidden();
}, {tier: 'ui', description: 'hovers the grant row of DATAGROK_SHARING_LOGIN and clicks its Remove button'});

/* Whom an entity is shared with, read where the platform shows it. grok.dapi.permissions.get answers
   with the edit and view buckets only, and a share made through the dialog lands in neither — the
   Sharing pane calls it "has special permissions" — so the API cannot see it and the pane is the
   claim. The pane fills its grants in after it is built and appends its Share... button last
   (db_entity_meta.dart renderSharingSection), so a claim reads it only once the button is there: a
   "does not list" read off the loading pane would pass whatever the grants are. Two expressions
   rather than one with "(not )": an optional literal is not a parameter, so a single step would
   always take the positive branch. */
async function sharingPane(page: Page): Promise<ReturnType<Page['locator']>> {
  const header = page.locator('.grok-prop-panel [name="div-section--Sharing"]').first();
  await expect(header, 'the Sharing pane of the context panel').toBeVisible({timeout: pollMs(30000)});
  if (await header.getAttribute('aria-expanded') !== 'true')
    await header.click();
  const pane = page.locator('.grok-prop-panel .d4-pane-sharing').first();
  await expect(pane.getByRole('button', {name: /^share\.\.\.$/i}), 'the Sharing pane, loaded (its Share... button)')
    .toBeVisible({timeout: pollMs(30000)});
  return pane;
}

export const sharingPaneLists = Then('the sharing pane should list the sharing user', async (page: Page) => {
  await expect(await sharingPane(page)).toContainText(new RegExp(sharingShownName(), 'i'));
}, {tier: 'ui', description: 'the Sharing section of the context panel, opened if it is closed, read once it has loaded'});

export const sharingPaneListsNot = Then('the sharing pane should not list the sharing user', async (page: Page) => {
  await expect(await sharingPane(page)).not.toContainText(new RegExp(sharingShownName(), 'i'));
}, {tier: 'ui', description: 'read once the pane has loaded'});

// Space and group name filters miss existing entities, so names are matched after reading every page.
type NamedSource = 'spaces' | 'models' | 'groups' | 'queries' | 'scripts' | 'connections';
type CleanupSource = NamedSource | 'projects' | 'tables';
type ServerEntity = {id: string; name: string; friendlyName: string; createdOn: number; children?: string[]};
type CleanupStage = {source: CleanupSource; ids: string[]};

async function serverEntities(page: Page, source: CleanupSource, filter = ''): Promise<ServerEntity[]> {
  return page.evaluate(async ([src, query]) => {
    try {
      // The tables gallery hides system tables, including training artifacts.
      let data = (src === 'tables' ? grok.dapi.entities : grok.dapi[src]).order('id');
      const filters = [src === 'tables' ? 'entityType.name = "TableInfo"' : '', query].filter(Boolean);
      if (filters.length > 0)
        data = data.filter(filters.map((filter) => `(${filter})`).join(' and '));
      if (src === 'projects')
        data = data.include('children');
      const result: ServerEntity[] = [];
      for (let pageNumber = 1; ; pageNumber++) {
        const entities = await data.list({pageSize: 1000, pageNumber});
        for (const entity of entities)
          result.push({id: entity.id, name: entity.name, friendlyName: entity.friendlyName,
            createdOn: entity.createdOn?.valueOf() ?? 0,
            // a project can hold a child whose entity is gone: one of those must not fail the listing
            children: src === 'projects' ? entity.children.filter(Boolean).map((child: any) => child.id) : undefined});
        if (entities.length < 1000)
          return result;
      }
    }
    catch (error) {
      // Dart ApiException loses its message when Playwright serializes it directly.
      throw new Error(`${src} list: ${(error as any)?.message ?? String(error)}`);
    }
  }, [source, filter] as [CleanupSource, string]);
}

/* groups.delete refuses a group holding a global permission, and an entity delete orphans the grant. */
async function deleteGlobalGrantsOf(page: Page, entity: ServerEntity): Promise<void> {
  const api = await serverRequests(page);
  for (const grant of await api.get<{id: string; userGroup?: {id: string}}[]>(`/privileges/permissions/?groupId=${entity.id}&global=true`))
    if (grant.userGroup?.id === entity.id)
      await api.remove(`/privileges/permissions/${grant.id}`);
}

/* A fixture name ends in its run's {run} or {time}. A run that was killed never reached its
   feature-end cleanup, so the fixtures of the same family that are older than any live feature go too. */
const RUN_SUFFIX = /-(\d{13,}|[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})$/;
const STALE_AFTER_MS = 60 * 60 * 1000;

export const fixtureFamilies = (names: string[]): string[] =>
  names.filter((name) => RUN_SUFFIX.test(name)).map((name) => name.replace(RUN_SUFFIX, ''));

export function isStaleFixture(entity: {name: string; friendlyName: string; createdOn: number}, families: string[],
  now = Date.now()): boolean {
  if (entity.createdOn === 0 || now - entity.createdOn < STALE_AFTER_MS)
    return false;
  return [entity.friendlyName, entity.name]
    .some((name) => RUN_SUFFIX.test(name) && families.includes(name.replace(RUN_SUFFIX, '')));
}

function namedCleanup(page: Page, source: NamedSource, what: string, names: string[]): () => Promise<void> {
  const families = fixtureFamilies(names);
  let createdAfter: number | undefined;
  let authorId: string;
  const pending = new Map<string, CleanupStage[]>();
  return async () => {
    await expect.poll(async () => {
      try {
        if (source === 'models' && createdAfter === undefined) {
          const owner = await page.evaluate(async () => {
            try {
              return {root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
                authorId: (await grok.dapi.users.current()).id};
            }
            catch (error) {
              throw new Error(`cleanup owner lookup: ${(error as any)?.message ?? String(error)}`);
            }
          });
          // No public dapi clock getter: use the server's HTTP Date and preserve its whole second.
          const response = await page.request.get(`${owner.root}/info/server`);
          if (!response.ok())
            throw new Error(`server clock request failed: HTTP ${response.status()}`);
          const serverTime = Date.parse(response.headers()['date'] ?? '');
          if (!Number.isFinite(serverTime))
            throw new Error('server clock response has no valid Date header');
          createdAfter = serverTime + 1000;
          authorId = owner.authorId;
        }
        const named = (await serverEntities(page, source))
          .filter((entity) => names.includes(entity.friendlyName) || names.includes(entity.name) ||
            isStaleFixture(entity, families));
        for (const entity of named) {
          if (pending.has(entity.id))
            continue;
          const stages: CleanupStage[] = [{source, ids: [entity.id]}];
          if (source === 'models') {
            const training: {id: string; owned: boolean} | null = await page.evaluate(async ([id, ownerId, cutoff]) => {
              try {
                const model = await grok.dapi.models.include('trainedOn').find(id);
                const table = await grok.functions.call('Get', {object: model, propertyName: 'trainedOn'});
                if (!table)
                  return null;
                const saved = await grok.dapi.tables.find(table.id);
                return {id: table.id, owned: saved != null && model.author.id === ownerId && saved.author.id === ownerId &&
                  model.createdOn.valueOf() > cutoff && saved.createdOn.valueOf() > cutoff};
              }
              catch (error) {
                throw new Error(`model ${id} trainedOn lookup: ${(error as any)?.message ?? String(error)}`);
              }
            }, [entity.id, authorId, createdAfter!] as [string, string, number]);
            const wrappers = await serverEntities(page, 'projects',
              `isEntity = true and isPackage = false and relations.entity.id = "${entity.id}"`);
            for (const wrapper of wrappers)
              if (wrapper.children?.length !== 1 || wrapper.children[0] !== entity.id)
                throw new Error(`project ${wrapper.id} has children outside model ${entity.id}; it was retained`);
            stages.unshift({source: 'projects', ids: wrappers.map((wrapper) => wrapper.id)});
            // trainedOn is a used table: require a new artifact by this run's user and no other model.
            if (training?.owned) {
              const users = await serverEntities(page, 'models', `trainedOn.id = "${training.id}"`);
              if (!users.some((model) => model.id === entity.id))
                throw new Error(`the training-table lookup did not return model ${entity.id}; ownership is unverified`);
              if (users.length === 1)
                stages.push({source: 'tables', ids: [training.id]});
            }
          }
          pending.set(entity.id, stages);
        }
        const left: string[] = [];
        for (const [modelId, stages] of pending) {
          while (stages.length > 0) {
            const stage = stages[0];
            // Only projects and tables answer an ID filter reliably; every other source is verified
            // against the complete listing, as spaces and groups had to be.
            const filter = stage.source === 'projects' || stage.source === 'tables' ?
              stage.ids.map((id) => `id = "${id}"`).join(' or ') : '';
            const entities = stage.ids.length === 0 ? [] : (await serverEntities(page, stage.source, filter))
              .filter((entity) => stage.ids.includes(entity.id));
            if (entities.length === 0) {
              stages.shift();
              continue;
            }
            for (const entity of entities) {
              try {
                if (stage.source === 'tables' &&
                  (await serverEntities(page, 'models', `trainedOn.id = "${entity.id}"`)).length > 0)
                  throw new Error(`training table ${entity.id} is still used by another model; it was retained`);
                // a chat outlives the entity it is about and then throws in every profile's chat
                // listing: the sources whose features post one delete it first
                if (['groups', 'queries', 'scripts', 'connections'].includes(stage.source))
                  await deleteChatsOf(page, entity.id);
                if (stage.source === 'groups')
                  await deleteGlobalGrantsOf(page, entity);
                await page.evaluate(async ([src, id, ownerId]) => {
                  let operation = 'lookup';
                  try {
                    // Project find expands every child; the listed wrapper is sufficient to delete it.
                    const entity = src === 'projects' ?
                      await grok.dapi.projects.filter(`id = "${id}"`).include('children').first() :
                      await grok.dapi[src].find(id);
                    if (!entity)
                      return;
                    if (entity.id !== id)
                      throw new Error(`the lookup returned ${entity.id} instead of ${id}; it was retained`);
                    if (src === 'projects' && (entity.children.length !== 1 || entity.children[0].id !== ownerId))
                      throw new Error(`project ${id} has children outside model ${ownerId}; it was retained`);
                    operation = 'delete';
                    await grok.dapi[src].delete(entity);
                  }
                  catch (error) {
                    throw new Error(`${src} ${id} ${operation}: ${(error as any)?.message ?? String(error)}`);
                  }
                }, [stage.source, entity.id, modelId] as [CleanupSource, string, string]);
                left.push(`${stage.source} ${entity.id}`);
              }
              catch (error) {
                left.push(`${stage.source} ${entity.id}: ${String(error)}`);
              }
            }
            break;
          }
        }
        return left;
      }
      catch (error) {
        return [`the listing or artifact lookup failed: ${String(error)}`];
      }
    }, {message: `${what} still on the server under ${names.join(', ')}`, timeout: pollMs(60000)}).toEqual([]);
  };
}

async function expectNamedCount(page: Page, source: CleanupSource, what: string, name: string, count: number): Promise<void> {
  await expect.poll(async () => {
    try {
      return (await serverEntities(page, source)).filter((entity) => entity.friendlyName === name || entity.name === name).length;
    }
    catch (error) {
      return `the listing failed: ${String(error)}`;
    }
  }, {message: `${what} the server holds under "${name}"`, timeout: pollMs(60000)}).toBe(count);
}

const namesOf = (list: string): string[] => list.split(',').map((n) => n.trim()).filter(Boolean);

/** The open Browse tree caches its nodes, so what the API added or deleted shows after a Refresh; the
 * tree rebuilds after it, and a node right-clicked mid-rebuild opens no menu, or another node's. */
async function refreshBrowseTree(page: Page): Promise<void> {
  if (await (await locate(page, el('browse panel'))).filter({visible: true}).count() === 0)
    return;
  const refreshed = await armEvent(page, 'onBrowseTreeRefreshed', pollMs(15000));
  await click(page, el('"Refresh" icon inside browse toolbar'));
  await refreshed();
}

export const noSpaceOnServer = Given('no space named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'spaces', 'spaces', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
  await refreshBrowseTree(page);
}, {tier: 'api', description: 'deletes earlier fixtures by name (comma-separated), refreshes the open Browse tree and waits for it to rebuild, and deletes them again at feature end'});

/* A space is listed once its save returns, and the save of a ROOT space is slow: 4.8 s alone and
   18 s with four features creating at once on a local stand (2026-09-10); the claim right after OK
   owns the budget the dialog-close claim does. */
export const spacesOnServer = Then('{int} space(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'spaces', 'spaces', name, count),
{tier: 'api', description: 'what the server holds, not what the tree draws — the refusal of a duplicate is a space that was never created'});

/* --- queries, scripts and connections -----------------------------------------------------------
   The entities a Queries, Scripts or Connections feature saves. Each is cleaned like a space: the
   complete listing (a name filter misses rows, and a grok name differs from the friendly one —
   "BDD-Q-x" is saved as "BDDQX"), the chats posted on it first, and the {run}/{time} family of a
   run that was killed swept with it. */

export const noQueryOnServer = Given('no query named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'queries', 'queries', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'deletes what an earlier run left under those names (comma-separated) with their chats, and again at feature end — verified gone'});

export const queriesOnServer = Then('{int} query/queries named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'queries', 'queries', name, count),
{tier: 'api', description: 'what the server holds, by friendly or grok name — not what the tree draws'});

export const noScriptOnServer = Given('no script named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'scripts', 'scripts', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'deletes it (and its chats) now and again at feature end, whatever a scenario saved under the name'});

export const scriptsOnServer = Then('{int} script(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'scripts', 'scripts', name, count), {tier: 'api'});

export const noConnectionOnServer = Given('no connection named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'connections', 'connections', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'deletes what an earlier run left under those names (comma-separated) with their chats, and again at feature end — verified gone'});

export const connectionsOnServer = Then('{int} connection(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'connections', 'connections', name, count),
{tier: 'api', description: 'what the server holds, not what the tree draws'});

/** A script saved through the JS API: the doc string is the script, with the `name:` header set to
 * the name the feature uses (its comment character follows the language of the body). */
export const scriptOnServer = Given('a script {string} is on the server:', async (page: Page, name: string, body: string) => {
  const cleanup = namedCleanup(page, 'scripts', 'scripts', [name]);
  atFeatureEnd(page, cleanup);
  await cleanup();
  const lines = body.split(/\r?\n/).filter((line) => !/^\s*(#|\/\/)name:/.test(line));
  const comment = /^\s*\/\//.test(lines.find((line) => /^\s*(#|\/\/)language:/.test(line)) ?? '') ? '//' : '#';
  await page.evaluate(async (text) => { await grok.dapi.scripts.save(DG.Script.create(text)); },
    [`${comment}name: ${name}`, ...lines].join('\n'));
  await expectNamedCount(page, 'scripts', 'scripts', name, 1);
}, {tier: 'api', description: 'saved through the JS API under that name; deleted with its chats at feature end'});

/** The coordinates of a data source without its credentials: Postgres points at the Northwind of
 * the stand's test server, any other source copies the parameters of its Samples Northwind. A
 * connection saved without `connString` and `ssl` answers every test with a bare HTTP 500, so both
 * are set the way the dialog saves them. */
export const connectionOnServer = Given('a {string} connection named {string} is on the server', async (page: Page, dataSource: string, name: string) => {
  const cleanup = namedCleanup(page, 'connections', 'connections', [name]);
  atFeatureEnd(page, cleanup);
  await cleanup();
  const failed = await page.evaluate(async ([ds, n]) => {
    let params: Record<string, unknown> = {server: 'db.datagrok.ai', port: 54322, db: 'northwind', connString: '', ssl: false};
    if (ds !== 'Postgres') {
      const sample = (await grok.dapi.connections.list({pageSize: 5000}))
        .find((c: any) => c.dataSource === ds && c.friendlyName === 'Northwind' && String(c.nqName).startsWith('Samples:'));
      if (!sample)
        return `no Samples Northwind connection of ${ds} to copy the coordinates from`;
      params = {...(await grok.dapi.connections.find(sample.id)).parameters};
    }
    try {
      await grok.dapi.connections.save(DG.DataConnection.create(n, {dataSource: ds, ...params}));
      return '';
    }
    catch (error: any) {
      // a Dart ApiException loses its message when Playwright serializes it directly
      return `saving the ${ds} connection "${n}": ${error?.message ?? String(error)}`;
    }
  }, [dataSource, name] as [string, string]);
  if (failed)
    throw new Error(failed);
  await expectNamedCount(page, 'connections', 'connections', name, 1);
  await refreshBrowseTree(page);
}, {tier: 'api', description: 'a connection of that data source without credentials, deleted at feature end'});

export const connectionDataSource = Then('the {string} connection on the server should have the data source {string}', async (page: Page, name: string, source: string) => {
  await expect.poll(async () => {
    const list = await page.evaluate(async (n) => (await grok.dapi.connections.list({pageSize: 5000}))
      .filter((c: any) => c.friendlyName === n || c.name === n).map((c: any) => String(c.dataSource)), name);
    return list.length === 1 ? list[0] : `${list.length} connections named "${name}"`;
  }, {message: `the data source of the connection "${name}" on the server`, timeout: pollMs(30000)}).toBe(source);
}, {tier: 'api', description: 'the provider the connection was saved under, not the tree branch it is shown in'});

export const noModelOnServer = Given('no predictive model named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'models', 'predictive models', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'deletes what an earlier run left under those names (comma-separated), and deletes them again when the feature ends'});

export const modelsOnServer = Then('{int} predictive model(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'models', 'predictive models', name, count),
{tier: 'api', description: 'what the server holds, not what the gallery draws'});

/** A login typed as the text of an input — a field that lists users by login, not a typeahead. */
async function enterLogin(page: Page, target: ElementRef, login: string): Promise<void> {
  const editor = await editorOf(page, el(target.phrase));
  await editor.fill(login);
  await editor.press('Enter');
  await expect(editor, `${target.phrase} after typing a login`).toHaveValue(login);
}

export const enterOwnLogin = When('user enters the current user\'s login into {element}', async (page: Page, target: ElementRef) => {
  const login: string = await page.evaluate(() => String(grok.shell.user.login));
  await enterLogin(page, target, login);
}, {tier: 'ui', description: 'the login of the account the run is signed in with, typed and committed with Enter'});

export const enterSharingLogin = When('user enters the sharing user\'s login into {element}', (page: Page, target: ElementRef) =>
  enterLogin(page, target, sharingLogin()),
{tier: 'ui', description: 'the login of the second account (DATAGROK_SHARING_LOGIN, else the setup\'s "bddsecond"), typed and committed with Enter'});
/* --- users, groups and roles ----------------------------------------------------------------------
   A role is a group on the server (the Roles view lists the groups flagged as roles), so the group
   steps serve both, under either word. A user can never be deleted: a feature that needs one shares
   a fixture user made once per stand, and only users-create, which tests the creation, adds any. */

export const noGroupOnServer = Given('no group/role named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'groups', 'groups or roles', namesOf(name));
  atFeatureEnd(page, cleanup);
  await cleanup();
}, {tier: 'api', description: 'a role is a group on the server: deletes earlier fixtures by name (comma-separated), and deletes them again at feature end'});

export const groupsOnServer = Then('{int} group(s)/role(s) named {string} should be on the server', (page: Page, count: number, name: string) =>
  expectNamedCount(page, 'groups', 'groups or roles', name, count), {tier: 'api', description: 'a role is a group on the server'});

/** The friendly name is set with the name: the server derives it by splitting camel case otherwise,
 * and "BDD-Group" would be listed as "BD D-Group". */
export const groupOnServer = Given('a group named {string} is on the server', async (page: Page, name: string) => {
  const cleanup = namedCleanup(page, 'groups', 'groups', [name]);
  atFeatureEnd(page, cleanup);
  await cleanup();
  await page.evaluate(async (n) => {
    const group = DG.Group.create(n);
    group.friendlyName = n;
    await grok.dapi.groups.save(group);
  }, name);
  await expectNamedCount(page, 'groups', 'groups', name, 1);
}, {tier: 'api', description: 'a new group under that name, deleted at feature end'});

const serverUsers = (page: Page, login: string): Promise<{id: string; status: string}[]> =>
  page.evaluate(async (l) => (await grok.dapi.users.filter(`login = "${l}"`).list()).map((u: any) => ({id: u.id, status: u.status})), login);

/** Like the global setup's bddsecond: the login is the name the gallery shows, and `<login>@datagrok.ai`
 * the mail. The favorites are the running account's; a run that died between Add and Remove left one. */
export const userOnServer = Given('a user {string} is on the server', async (page: Page, login: string) => {
  const status = await page.evaluate(async (l) => {
    const user = await grok.dapi.users.filter(`login = "${l}"`).first();
    if (!user) {
      const made = DG.User.create();
      made.login = l;
      made.email = `${l}@datagrok.ai`;
      made.firstName = l;
      made.lastName = '';
      made.status = 'active';
      await grok.dapi.users.save(made);
      return 'active';
    }
    for (const favorite of await grok.dapi.entities.getFavorites())
      if (favorite.id === user.id)
        await (window as any).grok_Favorites_Remove(favorite.dart, null);
    return user.status;
  }, login);
  if (status !== 'active') {
    const api = await serverRequests(page);
    await api.post('/public/v1/users/unblock', await api.get(`/public/v1/users/${login}`));
  }
  await expect.poll(() => serverUsers(page, login).then((users) => users.map((u) => u.status).join(', ') || 'no such user'),
    {message: `the user "${login}"`, timeout: pollMs(30000)}).toBe('active');
}, {tier: 'api', description: 'a fixture user made once per stand, since users cannot be deleted: found by login or created; put back to active and out of the favorites'});

export const usersOnServer = Then('{int} user(s) with login {string} should be on the server', async (page: Page, count: number, login: string) => {
  await expect.poll(() => serverUsers(page, login).then((users) => users.length), {message: `users with login "${login}"`, timeout: pollMs(60000)})
    .toBe(count);
}, {tier: 'api', description: 'what the server holds, not what the gallery draws'});

export const userStatusOnServer = Then('the user {string} should be {word} on the server', async (page: Page, login: string, status: string) => {
  const want = status === 'disabled' ? 'blocked' : status;
  await expect.poll(() => serverUsers(page, login).then((users) => users.map((u) => u.status).join(', ') || 'no such user'),
    {message: `the status of "${login}"`, timeout: pollMs(30000)}).toBe(want);
}, {tier: 'api', description: '"active", or "disabled" (the server says blocked)'});

export const personalGroupOnServer = Then('the user {string} should have a personal group on the server', async (page: Page, login: string) => {
  await expect.poll(() => page.evaluate(async (l) => {
    const user = await grok.dapi.users.include('group').filter(`login = "${l}"`).first();
    if (!user?.group)
      return user ? 'the user has no group yet' : 'no such user';
    const group = await grok.dapi.groups.find(user.group.id);
    return group?.personal ? `personal group "${group.friendlyName}"` : `group ${user.group.id} is not personal`;
  }, login), {message: `the personal group of "${login}"`, timeout: pollMs(30000)}).toBe(`personal group "${login}"`);
}, {tier: 'api', description: 'the security group every user has, named by the login and flagged personal'});

/** Who belongs to a group or holds a role: the member is a user by login or a group by name, and an
 * admin member is a group's Admin or a role's Can assign. */
async function membership(page: Page, member: string, group: string): Promise<string> {
  return page.evaluate(async ([m, g]) => {
    try {
      const all = async (source: any): Promise<any[]> => {
        const out: any[] = [];
        for (let pageNumber = 1; ; pageNumber++) {
          const items = await source.order('id').list({pageSize: 1000, pageNumber});
          out.push(...items);
          if (items.length < 1000)
            return out;
        }
      };
      const groups = await all(grok.dapi.groups);
      const target = groups.find((x) => x.friendlyName === g || x.name === g);
      if (!target)
        return `no group or role "${g}"`;
      const ids = new Set<string>(groups.filter((x) => x.friendlyName === m || x.name === m).map((x) => x.id));
      for (const user of await grok.dapi.users.include('group').filter(`login = "${m}"`).list())
        if (user.group)
          ids.add(user.group.id);
      if (ids.size === 0)
        return `no user or group "${m}"`;
      const full = await grok.dapi.groups.find(target.id);
      return full.adminMembers.some((x: any) => ids.has(x.id)) ? 'admin member' :
        full.members.some((x: any) => ids.has(x.id)) ? 'member' : 'not a member';
    }
    catch (error) {
      return `the lookup failed: ${(error as any)?.message ?? String(error)}`;
    }
  }, [member, group]);
}

const expectMembership = (page: Page, member: string, group: string, want: string[]) =>
  expect.poll(() => membership(page, member, group), {message: `"${member}" in "${group}"`, timeout: pollMs(30000)})
    .toMatch(new RegExp(`^(${want.join('|')})$`));

export const memberOnServer = Then('{string} should be a member of {string} on the server', (page: Page, member: string, group: string) =>
  expectMembership(page, member, group, ['member', 'admin member']), {tier: 'api', description: 'a group\'s members or a role\'s assignees, admins included'});

export const plainMemberOnServer = Then('{string} should be a plain member of {string} on the server', (page: Page, member: string, group: string) =>
  expectMembership(page, member, group, ['member']), {tier: 'api', description: 'a member without Admin (a role: without Can assign)'});

export const adminMemberOnServer = Then('{string} should be an admin member of {string} on the server', (page: Page, member: string, group: string) =>
  expectMembership(page, member, group, ['admin member']), {tier: 'api', description: 'Admin of a group, Can assign of a role'});

export const notMemberOnServer = Then('{string} should not be a member of {string} on the server', (page: Page, member: string, group: string) =>
  expectMembership(page, member, group, ['not a member']), {tier: 'api'});

/* --- the gallery ------------------------------------------------------------------------------------ */

/* The counter reads "N", "N of M" (M is the list) or "shown / total", "..." before it knows. Other
   features may add items meanwhile, so a comparison with the remembered count is one-sided. */
const rememberedCounts = new WeakMap<Page, number>();

async function galleryCount(page: Page): Promise<number | string> {
  const text = ((await (await locate(page, el('gallery counter'))).filter({visible: true}).first().textContent()
    .catch(() => null)) ?? '').trim();
  const m = /^(?:\d+\s+of\s+)?(\d+)(?:\s*\/\s*\d+)?$/.exec(text);
  return m ? Number(m[1]) : `not a count: "${text}"`;
}

export const rememberGalleryCount = When('user remembers the gallery counter', async (page: Page) => {
  let count: number | string = '';
  await expect.poll(async () => String(count = await galleryCount(page)), {message: 'the gallery counter'}).toMatch(/^[0-9]+$/);
  if (!rememberedCounts.has(page))
    atFeatureEnd(page, async () => { rememberedCounts.delete(page); });
  rememberedCounts.set(page, Number(count));
}, {tier: 'ui', description: 'the number the counter shows once the gallery has loaded, until the feature ends'});

type CountRelation = 'lower' | 'not lower' | 'higher';

async function expectCountVersusRemembered(page: Page, relation: CountRelation): Promise<void> {
  const remembered = rememberedCounts.get(page);
  if (remembered === undefined)
    throw new Error('no gallery counter remembered: "user remembers the gallery counter" first');
  await expect.poll(async () => {
    const count = await galleryCount(page);
    if (typeof count !== 'number')
      return count;
    const holds = relation === 'lower' ? count < remembered : relation === 'higher' ? count > remembered : count >= remembered;
    return `${holds ? '' : 'not '}${relation} (${count} vs ${remembered})`;
  }, {message: 'the gallery counter against the remembered one'}).toMatch(new RegExp(`^${relation} \\(`));
}

export const galleryCountLower = Then('the gallery counter should be lower than remembered', (page: Page) =>
  expectCountVersusRemembered(page, 'lower'));

export const galleryCountNotLower = Then('the gallery counter should not be lower than remembered', (page: Page) =>
  expectCountVersusRemembered(page, 'not lower'), {description: 'back to at least the remembered count — another feature may have added items meanwhile'});

export const galleryCountHigher = Then('the gallery counter should be higher than remembered', (page: Page) =>
  expectCountVersusRemembered(page, 'higher'), {description: 'a cleared search against the count it showed: a jump no item or two of another feature can fake'});

/* The first item says a reordering happened; which order the rest is in no reading exposes. */
const rememberedFirstItems = new WeakMap<Page, string>();

async function firstGalleryItem(page: Page): Promise<string> {
  if (typeof await galleryCount(page) !== 'number')
    return '';
  const gallery = await locate(page, el('gallery'));
  return ((await gallery.filter({visible: true}).first().locator('.d4-link-label').first().textContent({timeout: 1000})
    .catch(() => null)) ?? '').trim();
}

export const rememberFirstGalleryItem = When('user remembers the first item in gallery', async (page: Page) => {
  let first = '';
  await expect.poll(async () => first = await firstGalleryItem(page), {message: 'the first item in the gallery'}).not.toBe('');
  if (!rememberedFirstItems.has(page))
    atFeatureEnd(page, async () => { rememberedFirstItems.delete(page); });
  rememberedFirstItems.set(page, first);
}, {tier: 'ui', description: 'the name of the first item once the gallery has loaded, until the feature ends'});

async function expectFirstVersusRemembered(page: Page, same: boolean): Promise<void> {
  const remembered = rememberedFirstItems.get(page);
  if (remembered === undefined)
    throw new Error('no first item remembered: "user remembers the first item in gallery" first');
  await expect.poll(async () => {
    const first = await firstGalleryItem(page);
    return first === '' ? 'still loading' : first === remembered ? 'the remembered one' : `another one ("${first}")`;
  }, {message: `the first item in the gallery against "${remembered}"`}).toMatch(same ? /^the remembered one$/ : /^another one/);
}

export const firstGalleryItemSame = Then('the first item in gallery should be the remembered one', (page: Page) =>
  expectFirstVersusRemembered(page, true));

export const firstGalleryItemOther = Then('the first item in gallery should not be the remembered one', (page: Page) =>
  expectFirstVersusRemembered(page, false), {description: 'a reordering: another item leads the list'});

export const galleryMode = Then('the gallery should be in {word} mode', async (page: Page, mode: string) => {
  const gallery = await locate(page, el('gallery'));
  await expect.poll(() => gallery.filter({visible: true}).first().getAttribute('mode')
    .then((m) => (m ?? 'none').toLowerCase()).catch(() => 'no gallery'), {message: 'the render mode of the gallery'}).toBe(mode.toLowerCase());
}, {description: 'brief, card or grid — what the gallery renders its items as, not which toggle is lit'});

export const urlShouldContain = Then('the page address should contain {string}', async (page: Page, part: string) => {
  await expect.poll(() => page.url(), {message: 'the page address'}).toContain(part);
});

export const urlShouldNotContain = Then('the page address should not contain {string}', async (page: Page, part: string) => {
  await expect.poll(() => page.url(), {message: 'the page address'}).not.toContain(part);
});

/** How many viewers the current view holds — an analysis that is done is one that has put its
 * viewers on screen. */
export const viewHoldsViewers = Then('the current view should hold at least {int} viewer(s)',
  async (page: Page, count: number) => {
    await expect.poll(() => page.evaluate(() => Array.from(grok.shell.v?.viewers ?? []).length),
      {message: 'viewers of the current view', timeout: pollMs(120000)}).toBeGreaterThanOrEqual(count);
  }, {tier: 'api', description: 'a claim with the budget of the run that builds them: an analysis puts its viewers up when it ends'});

/** A table in the workspace and nothing else: no view, so a form that offers the open tables in a
 * choice gains the option without losing the focus of the view it lives in. Named after the file,
 * which is the name such a choice shows. */
export const openTableOf = Given('user opens a table {string} with:', async (page: Page, name: string, rows: string[][]) => {
  const [header, ...body] = rows;
  if (!header || body.length === 0)
    throw new Error('the table needs a header row and at least one row of values');
  await page.evaluate(([n, csv]) => {
    const df = (window as any).DG.DataFrame.fromCsv(csv);
    df.name = n;
    grok.shell.addTableView(df);
  }, [name, [header, ...body].map((r) => r.map((c) => /[",\n]/.test(c) ? `"${c.replace(/"/g, '""')}"` : c).join(',')).join('\n')] as [string, string]);
  await expect.poll(() => page.evaluate((n) => grok.shell.tv?.dataFrame?.name === n && !!grok.shell.tv.grid, name),
    {message: `a table view of "${name}"`}).toBe(true);
}, {tier: 'api', description: 'a small table written in the feature — the header row names the columns, types are detected as from a CSV — in a table view of its own'});

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

/** A dialog that commits to the server stays up until the server answers: creating a space took
 * 6-18 s on dev, and a confirmation dialog on a loaded stand the same, which straddles the shared
 * 15 s budget — so "should be hidden" passed or failed by luck. This claim owns its budget. */
export const dialogCloses = Then('the {string} dialog should close', async (page: Page, title: string) => {
  const dialog = await locate(page, el(`${JSON.stringify(title)} dialog`));
  await expect(dialog.filter({visible: true}), `the "${title}" dialog`).toHaveCount(0, {timeout: pollMs(60000)});
}, {tier: 'ui', description: 'the platform closes it when the work it started is done'});

/* --- what the server holds for a script or a query ----------------------------------------------
   The saved entity, not what the editor shows: the claim after a save is about what reached the
   server. A name matches the friendly name or the grok name, as the cleanup does. */

async function serverEntityNamed(page: Page, source: 'scripts' | 'queries', name: string): Promise<ServerEntity> {
  const found = (await serverEntities(page, source)).filter((e) => e.friendlyName === name || e.name === name);
  if (found.length !== 1)
    throw new Error(`${found.length} ${source} named "${name}" on the server, expected one`);
  return found[0];
}

export const scriptContains = Then('the script {string} on the server should contain {string}', async (page: Page, name: string, text: string) => {
  await expect.poll(async () => {
    try {
      const script = await serverEntityNamed(page, 'scripts', name);
      return page.evaluate(async (id) => String((await grok.dapi.scripts.find(id))?.script ?? ''), script.id);
    }
    catch (error) {
      return String(error);
    }
  }, {message: `the body of the script "${name}" on the server`, timeout: pollMs(30000)}).toContain(text);
}, {tier: 'api', description: 'the script body the save sent, not the text in the editor'});

export const scriptHasParam = Then('the script {string} on the server should have (an ){word} {string} of type {string}',
  async (page: Page, name: string, direction: string, param: string, type: string) => {
    if (!['input', 'output'].includes(direction))
      throw new Error(`a parameter is an input or an output, not "${direction}"`);
    const script = await serverEntityNamed(page, 'scripts', name);
    const params: string[] = await page.evaluate(async (id) => {
      const full = await grok.dapi.scripts.find(id);
      return [...full.inputs.map((p: any) => `input ${p.name}: ${p.propertyType}`),
        ...full.outputs.map((p: any) => `output ${p.name}: ${p.propertyType}`)];
    }, script.id);
    expect(params, `the parameters the server parsed from "${name}"`).toContain(`${direction} ${param}: ${type}`);
  }, {tier: 'api', description: 'the parameters the server parsed from the script header, not the header text'});

/* The transformation script a saved query carries, read from the server, not from the editor: the
   view runs a query from the editor's own copy (data_query_view.dart), so without this a save that
   never reached the server stays invisible. The JS API's DataQuery does not carry the script, so it
   is read from the query's REST entity. */
async function savedTransformations(page: Page, name: string): Promise<string> {
  const query = await serverEntityNamed(page, 'queries', name);
  const saved = await (await serverRequests(page)).get<{script?: string}>(`/connectors/queries/${query.id}`);
  return String(saved?.script ?? '');
}

export const queryTransformations = Then('the query {string} on the server should have transformations containing {string}', async (page: Page, name: string, text: string) => {
  await expect.poll(() => savedTransformations(page, name),
    {message: `the transformations of the query "${name}" on the server`, timeout: pollMs(30000)}).toContain(text);
}, {tier: 'api', description: 'query.script — what the Transformations tab saved, read back from the server'});

export const queryNoTransformations = Then('the query {string} on the server should not have transformations containing {string}', async (page: Page, name: string, text: string) => {
  await expect.poll(() => savedTransformations(page, name),
    {message: `the transformations of the query "${name}" on the server`, timeout: pollMs(30000)}).not.toContain(text);
}, {tier: 'api', description: 'polled like its positive twin: a step the save has not written yet is not an absent step'});

export const queryPostProcess = Then('the query {string} on the server should have a post-process containing {string}', async (page: Page, name: string, text: string) => {
  await expect.poll(async () => {
    try {
      const query = await serverEntityNamed(page, 'queries', name);
      return page.evaluate(async (id) => String((await grok.dapi.queries.find(id))?.postProcessScript ?? ''), query.id);
    }
    catch (error) {
      return String(error);
    }
  }, {message: `the post-process of the query "${name}" on the server`, timeout: pollMs(30000)}).toContain(text);
}, {tier: 'api', description: 'what the save sent to the server, not what the editor shows'});

export const entityHasLayout = Then('the {word} {string} on the server should have a layout', async (page: Page, kind: string, name: string) => {
  const source = kind === 'query' ? 'queries' : kind === 'script' ? 'scripts' : null;
  if (source == null)
    throw new Error(`a layout is saved with a query or a script, not with a ${kind}`);
  await expect.poll(async () => {
    try {
      const entity = await serverEntityNamed(page, source, name);
      // the view writes the layout id into the options of the output table parameter, and shares
      // the layout itself: a reference to a layout the server no longer holds is not a layout
      return page.evaluate(async ([id, src]) => {
        const full = await (src === 'queries' ? grok.dapi.queries : grok.dapi.scripts).find(id);
        const ids = full.outputs.map((p: any) => p.options?.['layout']).filter((x: any) => x);
        if (ids.length === 0)
          return 'no layout';
        const layout = await grok.dapi.layouts.find(ids[0]).catch(() => null);
        return layout ? 'a layout' : `a layout id the server does not hold (${ids[0]})`;
      }, [entity.id, source]);
    }
    catch (error) {
      return String(error);
    }
  }, {message: `the layout of the ${kind} "${name}" on the server`, timeout: pollMs(30000)}).toBe('a layout');
}, {tier: 'api', description: 'the layout id the save wrote on the output parameter, and the layout it points at'});

/* --- the script view ------------------------------------------------------------------------- */

export const scriptsView = Given('user opens the Scripts view', async (page: Page) => {
  await page.evaluate(() => {
    // the route opens another Scripts view each time; one already open is made current instead
    const open = Array.from(grok.shell.views as Iterable<any>).find((v) => v.type === 'scripts');
    if (open)
      grok.shell.v = open;
    else
      grok.shell.route('/scripts');
  });
  await page.waitForFunction(() => grok.shell.v?.type === 'scripts', null, {timeout: pollMs(30000)});
  // the search text is a user setting: a run that died mid-search leaves the next visit filtered
  const search = page.locator('.grok-gallery-search-bar input').filter({visible: true}).first();
  if (await search.inputValue() !== '')
    await search.fill('');
  await expect.poll(() => page.locator('.grok-items-view-counts').filter({visible: true}).first().textContent()
    .then((t) => (t ?? '').trim()).catch(() => ''), {message: 'the Scripts gallery counter once it has loaded every script',
    timeout: pollMs(60000)}).toMatch(/^\d+ of \d+$/);
}, {tier: 'api', description: 'Browse > Platform > Functions > Scripts, by its route (or made current when open), with its search empty and every script loaded'});

/* --- the console ---------------------------------------------------------------------------- */

const consoleText = (page: Page): Promise<string> =>
  page.evaluate(() => (document.querySelector('.d4-console-body') as HTMLElement | null)?.innerText ?? '');
const consoleMark = new WeakMap<Page, number>();

/** The console logs the calls only while it is open. */
async function openConsole(page: Page): Promise<void> {
  const input = page.locator('.d4-console-wrapper input.ui-input-editor').filter({visible: true});
  if (await input.count() === 0) {
    // a focused text field would take the backquote as a character
    await page.evaluate(() => (document.activeElement as HTMLElement | null)?.blur());
    await page.keyboard.press('Backquote');
  }
  await expect(input, 'the console input').toBeVisible();
}

export const noteConsole = When('user notes the console output', async (page: Page) => {
  await openConsole(page);
  consoleMark.set(page, (await consoleText(page)).length);
}, {tier: 'ui', description: 'opens the console (Backquote) when it is closed — it logs calls only while open — and marks what it shows now: a later "the console should show" reads only what came after'});

export const consoleShows = Then('the console should show {string}', async (page: Page, text: string) => {
  await expect.poll(async () => (await consoleText(page)).slice(consoleMark.get(page) ?? 0),
    {message: 'the console output since it was noted', timeout: pollMs(30000)}).toContain(text);
}, {description: 'the console logs every function call the UI makes and its outputs; read from the point "user notes the console output" marked'});

export const consoleShowsTimes = Then('the console should show {string} {int} time(s)', async (page: Page, text: string, times: number) => {
  const count = async () => (await consoleText(page)).slice(consoleMark.get(page) ?? 0).split(text).length - 1;
  await expect.poll(count, {message: `occurrences of "${text}" in the console output since it was noted`, timeout: pollMs(30000)})
    .toBeGreaterThanOrEqual(times);
  expect(await count(), `occurrences of "${text}" in the console output since it was noted`).toBe(times);
}, {description: 'exactly that many — a count that reaches the number and does not pass it'});

/** A call typed into the console as a person would, under the qualified name the platform gave the
 * saved script (its namespace is the owner's login, capitalised). */
export const consoleCall = When('user calls the script {string} from the console with {string}', async (page: Page, name: string, args: string) => {
  const script = await serverEntityNamed(page, 'scripts', name);
  const nqName: string = await page.evaluate(async (id) => String((await grok.dapi.scripts.find(id)).nqName), script.id);
  await openConsole(page);
  const input = page.locator('.d4-console-wrapper input.ui-input-editor').filter({visible: true});
  // put in whole: typed key by key, the console's completion rewrites the arguments as they come
  await input.fill(`${nqName}(${args})`);
  await input.press('Enter');
}, {tier: 'ui', description: 'opens the console (Backquote) when it is closed and enters the qualified call'});

/* --- the context panel's counting panes ------------------------------------------------------ */

const paneCountOf = (page: Page, pane: string): Promise<string> =>
  page.evaluate((p) => document.querySelector(`.grok-prop-panel [name="pane-${p}"]`)?.getAttribute('d4-info') ?? '', pane);

export const paneCountAtLeast = Then('the {string} pane of the context panel should count at least {int}', async (page: Page, pane: string, count: number) => {
  await expect.poll(async () => Number(await paneCountOf(page, pane) || '0'),
    {message: `the count the "${pane}" pane of the context panel shows`, timeout: pollMs(60000)}).toBeGreaterThanOrEqual(count);
}, {description: 'a counting pane with at least that many entries — the platform records an entity made through the JS API with a delay of its own'});

/* --- the side panel ------------------------------------------------------------------------------ */

/** Browse and the toolbox share the side panel: once the toolbox has taken it, "the browse panel is
 * open" finds `showBrowse` already true and the tree stays hidden until the toolbox lets it go. */
export const toolboxPaneHidden = Given('the toolbox pane is hidden', async (page: Page) => {
  await page.evaluate(() => { grok.shell.windows.showToolbox = false; });
  await expect(page.locator('.d4-toolbox[caption]').filter({visible: true}), 'the toolbox pane').toHaveCount(0, {timeout: pollMs(15000)});
}, {tier: 'api', description: 'the side panel goes back to what it showed before the toolbox (the browse tree when it is open)'});
