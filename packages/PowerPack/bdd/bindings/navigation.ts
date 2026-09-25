/* A project reached by its address: the page is loaded afresh on the project's direct link
   (`/p/<namespace>.<name>`, the address the platform gives a project), as pasting the link does. */
import {Page} from '@playwright/test';
import {Then, When} from '@datagrok-libraries/bdd';
import {expect, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;

async function projectPath(page: Page, name: string): Promise<string> {
  const path = await page.evaluate(async (n) => {
    const found = (await grok.dapi.projects.list({pageSize: 1000}))
      .filter((p: any) => [p.friendlyName, p.name].some((x) => String(x).toLowerCase() === n.toLowerCase()));
    return found.length === 1 ? String(found[0].path) : `${found.length} projects named "${n}"`;
  }, name);
  if (!path.startsWith('/p/'))
    throw new Error(`no direct link for the project: ${path}`);
  return path;
}

export const openDirectLink = When('user loads the direct link of project {string}', async (page: Page, name: string) => {
  const path = await projectPath(page, name);
  const simple = await page.evaluate(() => grok.shell.windows.simpleMode);
  await page.goto(new URL(path, page.url()).href, {waitUntil: 'domcontentloaded', timeout: 180000});
  await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
  await expect.poll(() => page.evaluate(() => grok.shell.tv?.dataFrame?.name ?? 'no table view'),
    {message: `the table view the direct link ${path} opens`, timeout: pollMs(60000)}).not.toBe('no table view');
  await page.evaluate((s) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = s;
  }, simple);
  await viewers.installViewerRuntime(page);
}, {tier: 'ui', description: 'the page loaded anew on /p/<namespace>.<name>, the address the platform gives the project; done when a table view is current'});

export const projectIsOpen = Then('the project {string} should be open', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => {
    const p = grok.shell.project;
    return p ? String(p.friendlyName ?? p.name).toLowerCase() : 'no project';
  }), {message: 'the project the shell has open (its name as the platform capitalizes it aside)'}).toBe(name.toLowerCase());
}, {description: 'grok.shell.project: the project the views belong to'});

export const projectNotOpen = Then('the project {string} should not be open', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => {
    const p = grok.shell.project;
    return p ? String(p.friendlyName ?? p.name).toLowerCase() : 'no project';
  }), {message: 'the project the shell has open'}).not.toBe(name.toLowerCase());
}, {description: 'grok.shell.project is another project (or none)'});

export const viewNotCurrent = Then('the {string} view should not be current', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')), {message: 'the current view'}).not.toBe(name);
});

export const noLoader = Then('no loading indicator should be visible', async (page: Page) => {
  await expect(page.locator('#grok-preloader, .grok-preloader, .grok-loader').filter({visible: true}),
    'the start-up splash and the loaders of the page').toHaveCount(0);
}, {description: 'the start-up splash of the platform and every loader spinner are gone'});
