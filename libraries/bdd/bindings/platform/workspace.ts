/* The shell's workspace and what surrounds it: the tables and table views that are open, a project
   reached by its direct link, and files on the server (the user's own files, a space's storage, a
   file dropped onto the window or chosen from the open menu). */
import {readFileSync} from 'node:fs';
import {resolve} from 'node:path';
import {type Page} from '@playwright/test';
import {expect, pollMs} from '../../src/runtime/patience.js';
import {Given, Then, When} from '../../src/registry.js';
import {el, type ElementRef} from '../../src/runtime/args.js';
import {atFeatureEnd} from '../../src/runtime/harness.js';
import {locate} from '../../src/runtime/locate.js';
import {fixtureFamilies, isStaleFixture} from '../../src/runtime/server.js';
import {installViewerRuntime, pickMenuPath, settleAll} from '../../src/runtime/viewers.js';
import {currentSimpleMode} from '../common/session.js';

declare const grok: any;

const listOf = (text: string): string[] => text.split(',').map((s) => s.trim()).filter(Boolean);

/* --- the tables and table views that are open ------------------------------------------------------ */

export const noTableLeft = Then('no table should be left in the workspace', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => (grok.shell.tables as any[]).map((t) => t.name).join(' | ')),
    {message: 'the tables still open in the workspace'}).toBe('');
}, {description: 'grok.shell.tables is empty, polled — the claim after a Close All made through the UI, before anything is reopened'});

const tableViewTables = (page: Page): Promise<string[]> =>
  page.evaluate(() => [...grok.shell.tableViews].map((v: any) => String(v.dataFrame?.name)));

export const openTableViewsExactly = Then('the open table views should be exactly {string}', async (page: Page, names: string) => {
  await expect.poll(() => tableViewTables(page), {message: 'the tables of the open table views', timeout: pollMs(30000)})
    .toEqual(listOf(names));
}, {description: 'the tables of the open table views (comma-separated), in the order the views opened: no other, none twice'});

export const noTableViewOpen = Given('no table view is open', async (page: Page) => {
  await expect.poll(() => tableViewTables(page), {message: 'the table views open before the file is opened'}).toEqual([]);
}, {description: 'a claim, not a cleanup: what opens afterwards is all the gesture brought'});

/** A file opened from the Files tree arrives as a table view once it is read and parsed, which for a
 * large one (SPGI) is well past a claim's default budget. */
export const tableViewOpened = Then('the {string} table view should open with {int} rows', async (page: Page, name: string, rows: number) => {
  await expect.poll(() => page.evaluate((n) => {
    const tv = (Array.from(grok.shell.tableViews) as any[]).find((v) => v.dataFrame?.name === n);
    return tv ? tv.dataFrame.rowCount : -1;
  }, name), {timeout: pollMs(120000), message: `the rows of the "${name}" table view (-1: not open)`}).toBe(rows);
  await settleAll(page);
}, {description: 'polls for a table view of that table for as long as a large file takes to load (up to two minutes), then its row count'});

/* --- a project reached by its address ---------------------------------------------------------------
   The page is loaded afresh on the project's direct link (`/p/<namespace>.<name>`, the address the
   platform gives a project), as pasting the link does. */
async function projectPath(page: Page, name: string): Promise<string> {
  const path = await page.evaluate(async (n) => {
    const filter = `name = ${JSON.stringify(n)} or friendlyName = ${JSON.stringify(n)}`;
    const found = (await grok.dapi.projects.filter(filter).list())
      .filter((p: any) => [p.friendlyName, p.name].some((x) => String(x).toLowerCase() === n.toLowerCase()));
    return found.length === 1 ? String(found[0].path) : `${found.length} projects named "${n}"`;
  }, name);
  if (!path.startsWith('/p/'))
    throw new Error(`no direct link for the project: ${path}`);
  return path;
}

export const openDirectLink = When('user loads the direct link of project {string}', async (page: Page, name: string) => {
  const path = await projectPath(page, name);
  const simple = await currentSimpleMode(page);
  await page.goto(new URL(path, page.url()).href, {waitUntil: 'domcontentloaded', timeout: 180000});
  await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
  await expect.poll(() => page.evaluate(() => grok.shell.tv?.dataFrame?.name ?? 'no table view'),
    {message: `the table view the direct link ${path} opens`, timeout: pollMs(60000)}).not.toBe('no table view');
  await page.evaluate((s) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = s;
  }, simple);
  await installViewerRuntime(page);
}, {tier: 'ui', description: 'the page loaded anew on /p/<namespace>.<name>, the address the platform gives the project; done when a table view is current'});

const openProjectName = (page: Page): Promise<string> => page.evaluate(() => {
  const p = grok.shell.project;
  return p ? String(p.friendlyName ?? p.name).toLowerCase() : 'no project';
});

export const projectIsOpen = Then('the project {string} should be open', async (page: Page, name: string) => {
  await expect.poll(() => openProjectName(page),
    {message: 'the project the shell has open (its name as the platform capitalizes it aside)'}).toBe(name.toLowerCase());
}, {description: 'grok.shell.project: the project the views belong to'});

export const projectNotOpen = Then('the project {string} should not be open', async (page: Page, name: string) => {
  await expect.poll(() => openProjectName(page), {message: 'the project the shell has open'}).not.toBe(name.toLowerCase());
}, {description: 'grok.shell.project is another project (or none)'});

export const noLoader = Then('no loading indicator should be visible', async (page: Page) => {
  await expect(page.locator('#grok-preloader, .grok-preloader, .grok-loader').filter({visible: true}),
    'the start-up splash and the loaders of the page').toHaveCount(0);
}, {description: 'the start-up splash of the platform and every loader spinner are gone'});

/* --- files on the server ----------------------------------------------------------------------------
   Browse > Files > My files is the account's Home storage, "<namespace>:Home" (the namespace is the
   user's own project). A file a feature puts there is written through the JS API — the file is the
   scene, not the subject — deleted when the feature ends, and read back gone; one of the same {run} or
   {time} family over an hour old is a killed run's, and goes first. A file of the bdd project is a path
   under its root. */
const bytesOf = (file: string): string => readFileSync(resolve(process.env.BDD_ROOT ?? process.cwd(), file)).toString('base64');

async function putIntoUsersFiles(page: Page, name: string, write: (path: string) => Promise<void>): Promise<void> {
  const home: string = await page.evaluate(async () => {
    const project = grok.shell.user.project.name;
    const share = (await grok.dapi.connections.list())
      .find((c: any) => c.dataSource === 'Files' && c.nqName === `${project}:Home`);
    if (!share)
      throw new Error(`the account has no home folder ${project}:Home on this stand`);
    return String(share.nqName);
  });
  const extension = /\.[^.]+$/.exec(name)?.[0] ?? '';
  const base = (file: string): string => file.endsWith(extension) ? file.slice(0, file.length - extension.length) : file;
  const families = fixtureFamilies([base(name)]);
  if (families.length > 0) {
    const files: {name: string; changed: number}[] = await page.evaluate(async (dir) =>
      (await grok.dapi.files.list(`${dir}/`, false)).filter((f: any) => !f.isDirectory)
        .map((f: any) => ({name: String(f.name), changed: f.updatedOn?.valueOf() ?? 0})), home);
    const stale = files.filter((f) => f.name.endsWith(extension) &&
      isStaleFixture({name: base(f.name), friendlyName: base(f.name), createdOn: f.changed}, families));
    await page.evaluate(async ([dir, names]) => {
      for (const n of names)
        await grok.dapi.files.delete(`${dir}/${n}`);
    }, [home, stale.map((f) => f.name)] as [string, string[]]);
  }
  const path = `${home}/${name}`;
  atFeatureEnd(page, async () => {
    await page.evaluate(async (p) => {
      if (await grok.dapi.files.exists(p))
        await grok.dapi.files.delete(p);
    }, path);
    await expect.poll(() => page.evaluate((p) => grok.dapi.files.exists(p), path),
      {message: `${path} still in My files`, timeout: pollMs(15000)}).toBe(false);
  });
  await write(path);
}

export const copyInUsersFiles = Given('a copy of the {string} file is in the user\'s files as {string}', (page: Page, source: string, name: string) =>
  putIntoUsersFiles(page, name, (path) => page.evaluate(async ([src, p]) => {
    await grok.dapi.files.write(p, await grok.dapi.files.readAsBytes(src));
  }, [source, path])),
{tier: 'api', description: 'a server file (System:DemoFiles/…) copied into My files of the signed-in account; deleted when the feature ends'});

export const fixtureInUsersFiles = Given('the {string} file of the project is in the user\'s files as {string}', (page: Page, file: string, name: string) =>
  putIntoUsersFiles(page, name, (path) => page.evaluate(async ([b64, p]) => {
    await grok.dapi.files.write(p, Array.from(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))));
  }, [bytesOf(file), path] as [string, string])),
{tier: 'api', description: 'a file of the bdd project written into My files of the signed-in account; deleted when the feature ends'});

export const fixtureInSpace = Given('the {string} file of the project is in the space {string} as {string}', async (page: Page, file: string, space: string, name: string) => {
  await page.evaluate(async ([b64, s, n]) => {
    const found = (await grok.dapi.spaces.list({pageSize: 1000})).find((p: any) => p.friendlyName === s || p.name === s);
    const target = found ?? await grok.dapi.spaces.createRootSpace(s);
    await grok.dapi.spaces.id(target.id).files.write(n, Array.from(Uint8Array.from(atob(b64), (c) => c.charCodeAt(0))));
  }, [bytesOf(file), space, name] as [string, string, string]);
  // the file goes on its own, so a space the feature does not take away ("no space named") keeps nothing of it
  atFeatureEnd(page, async () => {
    const left = await page.evaluate(async ([s, n]) => {
      const found = (await grok.dapi.spaces.list({pageSize: 1000})).find((p: any) => p.friendlyName === s || p.name === s);
      if (!found)
        return false;
      const files = grok.dapi.spaces.id(found.id).files;
      if (await files.exists(n))
        await files.delete(n);
      return files.exists(n);
    }, [space, name] as [string, string]);
    expect(left, `the file ${name} in the space ${space}`).toBe(false);
  });
}, {tier: 'api', description: 'the space made (a root space) when there is none of that name, and the file written into its storage; the file deleted at feature end and read back gone ("no space named" before it takes a space the feature made away)'});


/* A DataTransfer made in the page carries no file-system entry (`webkitGetAsEntry()` is null), and the
   platform reads a drop through those entries: the drag is made by the browser itself, through the
   DevTools protocol, with the file on disk, as the operating system hands a dragged file over. */
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

export const menuFileChooser = When('user picks {string} from the open menu and chooses the {string} file of the project', async (page: Page, path: string, file: string) => {
  const chooser = page.waitForEvent('filechooser', {timeout: pollMs(15000)});
  await pickMenuPath(page, path);
  await (await chooser).setFiles(resolve(process.env.BDD_ROOT ?? process.cwd(), file));
}, {tier: 'ui', description: 'the menu command opens the browser\'s file chooser, which is answered with a file of the bdd project'});
