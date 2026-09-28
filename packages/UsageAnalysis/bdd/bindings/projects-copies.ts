/* Projects: copies, links, uploading. What the platform vocabulary has no phrase for: a project
   card's own thumbnail, the second browser tab a copied link is opened in, whose table a copy
   holds, a space of demo files, presentation mode, and the personal view customizations the Save
   dialog keeps in the account's settings. */
import {type Page} from '@playwright/test';
import {element, Given, Then, When} from '@datagrok-libraries/bdd';
import {atFeatureEnd, expect, gestures, locate, pollMs} from '@datagrok-libraries/bdd/runtime';
import type {ElementRef} from '@datagrok-libraries/bdd/runtime';
import {noSpaceOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';

declare const grok: any;
declare const DG: any;

/* The Save project dialog's description: a bare textarea with a placeholder, which no text-area
   phrase reaches (the Name field has an aria-label, the description has none). */
element('project description field', {selector: '[name="dialog-Save-project"] textarea#description'});

/* The Save button the Dashboards panel shows on the New Dashboard (scratchpad) row. */
element('new dashboard save button', {selector: '.grok-view-browse [name="tree-New-Dashboard"] [name="button-Save"]'});

/* The Link Tables dialog (Data > Link Tables...) names its two table choices and the key column
   selectors of its first key row, but under one "Tables" and one "Key Columns" host each. */
const LINK = '[name="dialog-Link-Tables"]';
element('link tables first table', {selector: `${LINK} select[name="input-selectTableLeft"]`});
element('link tables second table', {selector: `${LINK} select[name="input-selectTableRight"]`});
element('link tables first key', {selector: `${LINK} [name="div-selectKeyCol0Row1"]`,
  description: 'the key column of the first table, first key row: a Dart column selector'});
element('link tables second key', {selector: `${LINK} [name="div-selectKeyCol1Row1"]`,
  description: 'the key column of the second table, first key row'});

/* --- a project card's thumbnail ------------------------------------------------------------------
   A card shows the picture the server keeps for the project (`/api/entities/picture/<id>`), which
   the ribbon's Save dialog renders from the view; a project saved without one shows the gallery's
   placeholder (`/images/datasets/no_picture.png`). The claim is on the image the card points at:
   the entity's own, and one the browser decodes. */
export const cardShowsPicture = Then('{element} should show the project\'s own picture', async (page: Page, target: ElementRef) => {
  const card = (await locate(page, target)).first();
  await expect(card, target.phrase).toBeVisible();
  await expect.poll(() => card.evaluate(async (e) => {
    const thumb = e.querySelector('.grok-gallery-grid-item-thumbnail') as HTMLElement | null;
    if (!thumb)
      return 'the card has no thumbnail';
    const url = getComputedStyle(thumb).backgroundImage.replace(/^url\("?|"?\)$/g, '');
    if (!url.includes('/api/entities/picture/'))
      return `the thumbnail is not the project's picture: ${url || 'none'}`;
    return new Promise<string>((resolve) => {
      const img = new Image();
      img.onload = () => resolve(img.naturalWidth > 0 && img.naturalHeight > 0 ? 'decoded' : 'empty image');
      img.onerror = () => resolve(`the picture does not load: ${url}`);
      img.src = url;
    });
  }), {message: `the thumbnail of ${target.phrase}`, timeout: pollMs(15000)}).toBe('decoded');
}, {description: 'the card\'s thumbnail is the entity picture the server keeps for the project (not the gallery placeholder), and the browser decodes it to an image'});

/* The link opened in a second tab of the same browser context — the session cookie and the stored
   token are the context's, as a person's new tab shares them. The tab boots the shell from the
   link; the step is done when the table the project opens has rows. Balloons are recorded from the
   tab's first script, since a project that fails to open says so in one. The tab closes when the
   scenario closes it, or when the feature ends. */
const tabs = new WeakMap<Page, Page>();

function newTab(page: Page): Page {
  const tab = tabs.get(page);
  if (!tab || tab.isClosed())
    throw new Error('no second tab is open: "user opens the copied link in a new tab" first');
  return tab;
}

export const openCopiedLink = When('user opens the copied link in a new tab', async (page: Page) => {
  const url = await gestures.readClipboard(page);
  expect(url, 'the clipboard text to open').toMatch(/^https?:\/\//);
  const tab = await page.context().newPage();
  tabs.set(page, tab);
  atFeatureEnd(page, async () => { await tab.close().catch(() => undefined); });
  // the balloons are read off the DOM as they are added: a `d4-balloon-shown` subscription made
  // while the shell boots is dropped (it recorded nothing, not even the boot's own balloons)
  await tab.addInitScript(() => {
    const w = window as any;
    w.__bddTabBalloons = [];
    new MutationObserver((records) => {
      for (const r of records)
        for (const n of Array.from(r.addedNodes))
          if (n instanceof HTMLElement && n.classList.contains('d4-balloon')) {
            const type = ['error', 'warning', 'info'].find((t) => n.classList.contains(t)) ?? '';
            w.__bddTabBalloons.push(`${type}: ${n.querySelector('.d4-balloon-content')?.textContent?.trim() ?? ''}`);
          }
    }).observe(document, {childList: true, subtree: true});
  });
  await tab.goto(url, {waitUntil: 'domcontentloaded', timeout: 180000});
  await expect.poll(() => tab.evaluate(() => (window as any).grok?.shell?.tv?.dataFrame?.rowCount ?? 0).catch(() => 0),
    {message: `rows of the table the link ${url} opened in the new tab`, timeout: pollMs(180000)}).toBeGreaterThan(0);
}, {tier: 'ui', description: 'a second page of the same browser context goes to the URL on the clipboard; done when the project it opens shows a table with rows'});

type TabState = {project: string; viewers: string[]; passing: number; sort: string; hidden: string[]; balloons: string[]};

function tabState(tab: Page): Promise<TabState> {
  return tab.evaluate(() => {
    const g = (window as any).grok;
    const tv = g.shell.tv;
    const grid = tv?.grid;
    const sortCols = (grid?.sortByColumns ?? []).map((c: any) => c.name);
    const sortTypes = grid?.sortTypes ?? [];
    return {
      project: String(g.shell.project?.friendlyName ?? g.shell.project?.name ?? ''),
      viewers: tv ? Array.from(tv.viewers).map((v: any) => String(v.type)).filter((t: string) => t !== 'Filters').sort() : [],
      passing: tv ? tv.dataFrame.filter.trueCount : -1,
      sort: sortCols.map((c: string, i: number) => `${c} ${sortTypes[i] ? 'ascending' : 'descending'}`).join(', '),
      hidden: grid ? tv.dataFrame.columns.names().filter((n: string) => grid.columns.byName(n)?.visible === false) : [],
      balloons: (window as any).__bddTabBalloons ?? [],
    };
  });
}

export const tabShowsProject = Then('the new tab should show the {string} project with the viewers {string}', async (page: Page, name: string, viewers: string) => {
  const tab = newTab(page);
  const want = viewers.split(',').map((s) => s.trim()).sort();
  await expect.poll(async () => { const s = await tabState(tab); return `${s.project}: ${s.viewers.join(', ')}`; },
    {message: 'the project the new tab opened and the viewers of its table view', timeout: pollMs(30000)})
    .toBe(`${name}: ${want.join(', ')}`);
}, {description: 'the project open in the second tab, and exactly these viewers in its table view (in any order); the filter panel is not counted — it docks late and is claimed by the rows it passes'});

export const tabPassing = Then('the table in the new tab should have {int} rows passing the filter', async (page: Page, rows: number) => {
  await expect.poll(async () => (await tabState(newTab(page))).passing, {message: 'the rows passing the filter of the table in the new tab'}).toBe(rows);
});

export const tabGridSorted = Then('the grid in the new tab should be sorted by {string}', async (page: Page, sort: string) => {
  await expect.poll(async () => (await tabState(newTab(page))).sort, {message: 'the sort of the grid in the new tab'}).toBe(sort);
}, {description: '"<column> ascending|descending", from the grid\'s sortByColumns and sortTypes'});

export const tabGridHides = Then('the grid in the new tab should hide the {string} column', async (page: Page, column: string) => {
  await expect.poll(async () => (await tabState(newTab(page))).hidden, {message: 'the hidden columns of the grid in the new tab'}).toContain(column);
});

/* Where a viewer sits in the second tab: its panel's box against the area the table view's docked
   panels cover, the same reading as "… should be docked along the left edge of the view" of the
   main page (6 px slack), read in the tab's own freshly loaded layout. */
export const tabViewerDockedLeft = Then('the {string} viewer in the new tab should be docked along the left edge of the view', async (page: Page, type: string) => {
  const tab = newTab(page);
  await expect.poll(() => tab.evaluate((t) => {
    const tv = (window as any).grok.shell.tv;
    const viewer = tv ? Array.from(tv.viewers as any[]).find((v: any) => v.type === t) : null;
    const panel = viewer?.root?.closest('.panel-base') as HTMLElement | null;
    if (!panel)
      return `no docked "${t}" viewer in the table view`;
    const box = (e: Element) => e.getBoundingClientRect();
    const root = box(tv.root);
    const panels = Array.from(tv.root.querySelectorAll('.panel-base') as NodeListOf<HTMLElement>).map(box).filter((r) => r.width > 0);
    const left = Math.min(root.x, ...panels.map((r) => r.x));
    const top = Math.min(...panels.map((r) => r.y));
    const bottom = Math.max(root.y + root.height, ...panels.map((r) => r.y + r.height));
    const a = box(panel);
    const ok = Math.abs(a.x - left) <= 6 && a.y <= top + 6 && a.y + a.height >= bottom - 6;
    return ok ? 'along the left edge' : `spans ${Math.round(a.x)},${Math.round(a.y)}..${Math.round(a.x + a.width)},${Math.round(a.y + a.height)}; ` +
      `the docked viewers span ${Math.round(left)},${Math.round(top)}..bottom ${Math.round(bottom)}`;
  }, type), {message: `where the ${type} viewer of the new tab is docked`, timeout: pollMs(15000)}).toBe('along the left edge');
}, {description: 'the viewer\'s panel touches the left edge of the docked area and runs its whole height (6 px slack), read in the second tab'});

export const tabNoErrorBalloon = Then('the new tab should have shown no error balloon', async (page: Page) => {
  const balloons = (await tabState(newTab(page))).balloons.filter((b) => b.startsWith('error'));
  expect(balloons, 'error balloons the new tab showed since it opened').toEqual([]);
});

export const closeTab = When('user closes the new tab', async (page: Page) => {
  await newTab(page).close();
  tabs.delete(page);
  await page.bringToFront();
}, {tier: 'ui'});

/* --- the grid of the current table view --------------------------------------------------------------
   Whether a column is shown is the grid column's own visibility (Order or Hide Columns clears it),
   not whether its header is drawn: a header scrolled out of a narrow grid, or not laid out yet, is
   no reading of a hidden column. */
function gridColumnShown(page: Page, column: string): Promise<string> {
  return page.evaluate((c) => {
    const grid = grok.shell.tv?.grid;
    const col = grid?.columns.byName(c);
    return !grid ? 'no table view is current' : !col ? `the grid has no "${c}" column` : col.visible ? 'shown' : 'hidden';
  }, column);
}

export const gridHidesColumn = Then('the grid should hide the {string} column', async (page: Page, column: string) => {
  await expect.poll(() => gridColumnShown(page, column), {message: `the "${column}" column of the grid`}).toBe('hidden');
}, {description: 'the grid column of the current table view exists and is not visible'});

export const gridShowsColumn = Then('the grid should show the {string} column', async (page: Page, column: string) => {
  await expect.poll(() => gridColumnShown(page, column), {message: `the "${column}" column of the grid`}).toBe('shown');
}, {description: 'the grid column of the current table view exists and is visible'});

/* The cell renderer the grid resolved for a column: the renderer a package registers for the
   column's semantic type (Chem's "Molecule"), or the grid's default when none is found. */
export const gridColumnRenderer = Then('the grid should draw the {string} column with the {string} renderer', async (page: Page, column: string, renderer: string) => {
  await expect.poll(() => page.evaluate((c) => {
    const col = grok.shell.tv?.grid?.col(c);
    return col ? `${col.cellType} / ${col.renderer?.name ?? 'no renderer'}` : `the grid has no "${c}" column`;
  }, column), {message: `the cell type and renderer of the "${column}" grid column`, timeout: pollMs(30000)}).toBe(`${renderer} / ${renderer}`);
}, {description: 'the grid column\'s cell type and the name of the renderer it resolved are both that name'});

/* --- where a copy's table lives on the server -------------------------------------------------------
   A copy with link keeps the original's table (the same TableInfo, among its links); a copy with
   clone uploads a table of its own (a new TableInfo, among its children). Read from the projects'
   children and links. */
async function projectTableId(page: Page, project: string, table: string): Promise<string> {
  return page.evaluate(async ([p, t]) => {
    const listed = await grok.dapi.projects.filter(`friendlyName = "${p}" or name = "${p}"`).first();
    if (!listed)
      return `no project "${p}" on the server`;
    const found = await grok.dapi.projects.find(listed.id);
    const infos = [...found.children, ...found.links].filter((c: any) => c instanceof DG.TableInfo);
    const info = infos.find((c: any) => c.friendlyName === t || c.name === t);
    return info ? `table ${info.id}` : `the project "${p}" holds no "${t}" table; it holds: ${infos.map((c: any) => c.friendlyName).join(', ') || 'none'}`;
  }, [project, table]);
}

async function tableIdsOf(page: Page, table: string, project: string, other: string): Promise<[string, string]> {
  const ids = [await projectTableId(page, project, table), await projectTableId(page, other, table)] as [string, string];
  for (const id of ids)
    expect(id, 'a table of the project on the server').toMatch(/^table /);
  return ids;
}

export const sameTableAs = Then('the {string} table of the {string} project should be the one of the {string} project',
  async (page: Page, table: string, project: string, other: string) => {
    const [a, b] = await tableIdsOf(page, table, project, other);
    expect(a, `the "${table}" table of "${project}" against the one of "${other}"`).toBe(b);
  }, {description: 'both projects hold the same TableInfo (the same id) among their children or links on the server'});

export const ownTableNotOf = Then('the {string} table of the {string} project should be its own, not the one of the {string} project',
  async (page: Page, table: string, project: string, other: string) => {
    const [a, b] = await tableIdsOf(page, table, project, other);
    expect(a, `the "${table}" table of "${project}" against the one of "${other}"`).not.toBe(b);
  }, {description: 'both projects hold a table of that name on the server, and they are different TableInfos (different ids)'});

/* --- a space holding demo files -----------------------------------------------------------------
   The fixture of the uploading matrix: a root space whose Files hold copies of Browse > Files > Demo
   > northwind files, made through the JS API (the space's own file client), so the scenario's
   subject — opening them from the space and saving them — starts from a known space. The space
   (with its files) goes through the platform's "no space named …" before it is made and when the
   feature ends; the Create Space dialog and the drag with Copy are claimed in
   projects-lifecycle-spaces. */
export const spaceWithNorthwindFiles = Given('a space {string} holds copies of the northwind files {string}', async (page: Page, name: string, files: string) => {
  await noSpaceOnServer(page, name);
  const held = await page.evaluate(async ([n, list]) => {
    const space = await grok.dapi.spaces.createRootSpace(n);
    const client = grok.dapi.spaces.id(space.id).files;
    for (const f of list) {
      await client.writeString(f, await grok.dapi.files.readAsText(`System:DemoFiles/northwind/${f}`));
      if (!await client.exists(f))
        throw new Error(`${f} was not written into the space "${n}"`);
    }
    return (await grok.dapi.files.list(`${space.name}:Files/`, false)).map((f: any) => f.fileName).sort();
  }, [name, files.split(',').map((f) => f.trim())] as [string, string[]]);
  expect(held, `the files of the space "${name}"`).toEqual(files.split(',').map((f) => f.trim()).sort());
}, {tier: 'api', description: 'a root space made through the JS API, with the named files of System:DemoFiles/northwind copied into its Files; removed (with its files) before and at feature end'});

/* --- presentation mode ----------------------------------------------------------------------------
   A project saved with the Save dialog's Presentation mode switch opens in presentation mode, and
   the save itself puts the shell into it: the shell's root gets the "presentation" class and a
   "back to design mode" link, which is how a person leaves it. The mode is the page's, not the
   account's, so a feature that turns it on puts it back when it ends. */
element('presentation back link', {selector: '.presentation-back-button',
  description: 'the "back to design mode" link the shell shows while in presentation mode'});

export const presentationModeOff = Given('presentation mode is off, now and when the feature ends', async (page: Page) => {
  const off = () => page.evaluate(() => { grok.shell.windows.presentationMode = false; });
  atFeatureEnd(page, off);
  await off();
}, {tier: 'api', description: 'grok.shell.windows.presentationMode = false, now and at feature end'});

/* --- personal view customizations ----------------------------------------------------------------
   "Save personal view customizations" keeps the changed views in the account's settings, under
   `projectCustomViews` and the project's id — not on the project, and not removed with it. A run
   removes the entries whose project is no longer on the server, when it starts and when it ends
   (after the projects are deleted). */
async function sweepCustomViews(page: Page): Promise<void> {
  const removed = await page.evaluate(async () => {
    const gone: string[] = [];
    for (const id of Object.keys(grok.userSettings.get('projectCustomViews') ?? {})) {
      if (!await grok.dapi.projects.find(id).catch(() => null)) {
        grok.userSettings.delete('projectCustomViews', id);
        gone.push(id);
      }
    }
    return gone;
  });
  if (removed.length === 0)
    return;
  // the settings reach the server's user data storage after the call returns
  await expect.poll(() => page.evaluate(async (ids) => {
    const stored = await grok.dapi.userDataStorage.get('projectCustomViews', true) ?? {};
    return ids.filter((id) => id in stored).join(', ') || 'gone';
  }, removed), {message: 'the personal view customizations of deleted projects on the server', timeout: pollMs(15000)}).toBe('gone');
}

/* The Save dialog puts the customizations into the page's copy of the settings only
   (`UserSettingsStorage.add`); the page sends it to the server on its 10 s sync timer. Until then
   another session of the account — a reload, a second tab — opens the project without them. */
export const customViewsOnServer = Then('the personal view customizations of the {string} project should be on the server', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(async (n) => {
    const project = await grok.dapi.projects.filter(`friendlyName = "${n}" or name = "${n}"`).first();
    if (!project)
      return `no project "${n}" on the server`;
    const stored = await grok.dapi.userDataStorage.get('projectCustomViews', true) ?? {};
    if (!(project.id in stored))
      return `the account's projectCustomViews on the server has no entry for ${project.id}; it has: ${Object.keys(stored).join(', ') || 'none'}`;
    const views = JSON.parse(stored[project.id]);
    return Array.isArray(views) && views.length > 0 ? 'stored' : `the entry for ${project.id} holds no view: ${stored[project.id]}`;
  }, name), {message: `the personal view customizations of "${name}" in the account's settings on the server`, timeout: pollMs(30000)}).toBe('stored');
}, {description: 'the account\'s projectCustomViews user data on the server holds a non-empty view list under the project\'s id — what another session of the account opens the project with'});

export const tabBalloonShown = Then('the new tab should have shown a warning balloon containing {string}', async (page: Page, text: string) => {
  const warnings = (await tabState(newTab(page))).balloons.filter((b) => b.startsWith('warning'));
  expect(warnings.some((b) => b.includes(text)), `warning balloons the new tab showed since it opened: ${warnings.join(' | ') || 'none'}`).toBe(true);
}, {description: 'read from the balloons the second tab recorded since its first script, after the project it opened has a table with rows'});

export const noOrphanCustomViews = Given('no personal view customizations of deleted projects are kept', async (page: Page) => {
  atFeatureEnd(page, () => sweepCustomViews(page));
  await sweepCustomViews(page);
}, {tier: 'api', description: 'removes the account\'s personal view customizations of projects no longer on the server, now and when the feature ends (after the projects are deleted)'});

