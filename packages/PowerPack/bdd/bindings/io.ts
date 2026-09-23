/* Files of the bdd project put where a user keeps files — the "My files" share of the account, a
   space's storage — and a file dropped onto the platform window the way the browser delivers a file
   dragged in from the desktop. Each put file is deleted when the feature ends. */
import {readFileSync} from 'node:fs';
import {resolve} from 'node:path';
import {Page} from '@playwright/test';
import {Given, Then, When, element} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, el, expect, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

element('file drop overlay', {selector: 'xpath=//div[./*[local-name()="svg"]//*[local-name()="text" and normalize-space()="Incoming!"]]',
  description: 'the "Incoming! Drop your CSV files to open them locally" layer the platform lays over the window while files are dragged over it'});

const bytesOf = (file: string): string => readFileSync(resolve(process.env.BDD_ROOT ?? process.cwd(), file)).toString('base64');

export const fixtureInHome = Given('the {string} file of the project is in the home folder as {string}', async (page: Page, file: string, name: string) => {
  const path: string = await page.evaluate(async ([b64, n]) => {
    const project = grok.shell.user.project.name;
    const home = async () => (await grok.dapi.connections.list()).find((c: any) => c.dataSource === 'Files' && c.nqName === `${project}:Home`);
    // an account made before the stand kept personal folders has none until the account is saved again,
    // which is what creates the folder for a new account
    if (!await home())
      await grok.dapi.users.save(await grok.dapi.users.find(grok.shell.user.id));
    const share = await home();
    if (!share)
      throw new Error(`no home folder ${project}:Home on this stand`);
    const target = `${share.nqName}/${n}`;
    await grok.dapi.files.write(target, Array.from(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))));
    return target;
  }, [bytesOf(file), name] as [string, string]);
  atFeatureEnd(page, () => page.evaluate(async (p) => {
    if (await grok.dapi.files.exists(p))
      await grok.dapi.files.delete(p);
  }, path));
}, {tier: 'api', description: 'a file of the bdd project (a path under its root) written into the "My files" share of the current user, and deleted when the feature ends'});

export const fixtureInSpace = Given('the {string} file of the project is in the space {string} as {string}', async (page: Page, file: string, space: string, name: string) => {
  await page.evaluate(async ([b64, s, n]) => {
    const found = (await grok.dapi.spaces.list({pageSize: 1000})).find((p: any) => p.friendlyName === s || p.name === s);
    const target = found ?? await grok.dapi.spaces.createRootSpace(s);
    await grok.dapi.spaces.id(target.id).files.write(n, Array.from(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))));
  }, [bytesOf(file), space, name] as [string, string, string]);
}, {tier: 'api', description: 'the space made (a root space) when there is none of that name, and the file written into its storage; "no space named" before it takes the space away at feature end'});

/* A DataTransfer made in the page carries no file-system entry (`webkitGetAsEntry()` is null), and
   the platform reads a drop through those entries: the drag is made by the browser itself, through
   the DevTools protocol, with the file on disk, as the operating system hands a dragged file over. */
export const dropFile = When('user drops the {string} file of the project onto {element}', async (page: Page, file: string, target: ElementRef) => {
  const box = await (await locate(page, target)).first().boundingBox();
  if (!box)
    throw new Error(`${target.phrase} has no box to drop onto`);
  const cdp = await page.context().newCDPSession(page);
  const data = {items: [], files: [resolve(process.env.BDD_ROOT ?? process.cwd(), file)], dragOperationsMask: 1};
  const at = {x: box.x + box.width / 2, y: box.y + box.height / 2};
  try {
    await cdp.send('Input.dispatchDragEvent', {type: 'dragEnter', ...at, data});
    const overlay = await locate(page, el('file drop overlay'));
    await expect(overlay.first(), 'the drop layer the platform shows while files are dragged over it').toBeVisible();
    await cdp.send('Input.dispatchDragEvent', {type: 'dragOver', ...at, data});
    await cdp.send('Input.dispatchDragEvent', {type: 'drop', ...at, data});
    await expect(overlay, 'the drop layer, gone once the file is dropped').toHaveCount(0);
  }
  finally {
    await cdp.detach().catch(() => undefined);
  }
}, {tier: 'ui', description: 'the browser\'s own drag of the file from disk: it enters over the element, the platform lays its drop layer over the window, and the file is dragged over it and dropped there'});

export const tableViewsOpen = Then('the table views {string} should be open', async (page: Page, list: string) => {
  const wanted = list.split(',').map((n) => n.trim()).filter(Boolean);
  await expect.poll(() => page.evaluate(() => [...grok.shell.tableViews].map((v: any) => String(v.dataFrame?.name))),
    {message: 'the table views open in the workspace', timeout: pollMs(30000)}).toEqual(wanted);
}, {description: 'exactly these table views, in the order they were opened — no other table view, none twice'});

export const cellOfTable = Then('the value of {string} column in row {int} of table {string} should be {string}', async (page: Page, column: string, row: number, table: string, value: string) => {
  await expect.poll(() => page.evaluate(([c, r, t]) => {
    const df = grok.shell.tables.find((x: any) => x.name === t);
    if (!df)
      return `no table ${t}`;
    const col = df.col(c);
    return col ? String(col.getString(r - 1)) : `no column ${c} in ${t}`;
  }, [column, row, table] as [string, number, string]), {message: `${column} of row ${row} of ${table}`}).toBe(value);
}, {description: 'the value as the column formats it; rows are 1-based'});

export const columnTypeOfTable = Then('{string} column of table {string} should have type {string}', async (page: Page, column: string, table: string, type: string) => {
  await expect.poll(() => page.evaluate(([c, t]) => {
    const df = grok.shell.tables.find((x: any) => x.name === t);
    return df?.col(c)?.type ?? `no column ${c} in ${t}`;
  }, [column, table] as [string, string]), {message: `the type of ${column} in ${table}`}).toBe(type);
});


export const noTableViewOpen = Given('no table view is open', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => [...grok.shell.tableViews].map((v: any) => String(v.dataFrame?.name))),
    {message: 'the table views open before the file is opened'}).toEqual([]);
}, {description: 'a claim, not a cleanup: what opens afterwards is all the file brought'});

export const menuFileChooser = When('user picks {string} from the open menu and chooses the {string} file of the project', async (page: Page, path: string, file: string) => {
  const chooser = page.waitForEvent('filechooser', {timeout: pollMs(15000)});
  await viewers.pickMenuPath(page, path);
  await (await chooser).setFiles(resolve(process.env.BDD_ROOT ?? process.cwd(), file));
}, {tier: 'ui', description: 'the menu command opens the browser\'s file chooser, which is answered with a file of the bdd project'});

/* The Browse tree rebuilds itself after it is shown (the nodes an account had expanded come back),
   and a node double-clicked while that happens is replaced under the pointer. Refresh rebuilds it
   once more and says when it is done (onBrowseTreeRefreshed). */
export const refreshBrowseTree = When('user refreshes the browse tree', async (page: Page) => {
  const refreshed = await viewers.armEvent(page, 'onBrowseTreeRefreshed', pollMs(15000));
  await (await locate(page, el('"Refresh" icon inside browse toolbar'))).first().click();
  if (!await refreshed())
    throw new Error('the Browse tree did not report the end of its refresh (onBrowseTreeRefreshed)');
}, {tier: 'ui', description: 'the Refresh icon of the Browse toolbar; done when the tree reports it has been rebuilt'});
