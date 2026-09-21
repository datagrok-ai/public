import {Page, expect} from '@playwright/test';
import {loginToDatagrok} from '../spec-login';

export const BASE_URL = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

export const projectsTestOptions = {
  baseURL: BASE_URL,
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 120_000,
};

export async function evalJs<T = any>(page: Page, script: string): Promise<T> {
  return page.evaluate(script) as Promise<T>;
}

export async function gotoApp(page: Page) {

  await loginToDatagrok(page);
}

export async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

export async function setupSession(page: Page) {
  await evalJs(page, `(() => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  })()`);
}

export async function openCsv(page: Page, path: string) {
  await evalJs(page, `(async () => {
    const df = await grok.dapi.files.readCsv('${path}');
    grok.shell.addTableView(df);
    return df.rowCount;
  })()`);
  await page.waitForTimeout(1500);
}

export async function saveProjectViaDialog(page: Page, name: string) {

  const saveBtn = page.locator('button[name="button-Save"]:visible').first();
  await saveBtn.waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500);  

  try {
    await saveBtn.click({timeout: 5_000});
  } catch (_) {
    await page.evaluate(() => {
      const candidates = Array.from(document.querySelectorAll('button[name="button-Save"]'));
      const visible = candidates.find(b => (b as HTMLElement).offsetParent !== null);
      if (visible) (visible as HTMLElement).click();
    });
  }
  const dialog = page.locator(
    '.d4-dialog[name="dialog-Save-project"], .d4-dialog:has-text("Save project")',
  ).first();
  await dialog.waitFor({timeout: 15_000});
  const nameInput = dialog.locator(
    'input[name="input-Name"], input[type="text"].ui-input-editor',
  ).first();
  await nameInput.fill(name);
  await dialog.locator('button.ui-btn-ok, [name="button-OK"]').first().click();

  const shareDialog = page.locator(
    '.d4-dialog[name^="dialog-Share-"], .d4-dialog:has-text("Share ")',
  ).first();
  if (await shareDialog.isVisible({timeout: 10_000}).catch(() => false)) {
    const cancel = shareDialog.locator(
      '[name="button-CANCEL"], button.ui-btn-cancel, button:has-text("Cancel")',
    ).first();
    if (await cancel.isVisible({timeout: 2_000}).catch(() => false))
      await cancel.click();
    else
      await page.keyboard.press('Escape');
    await expect(shareDialog).toBeHidden({timeout: 10_000});
  }
}

export async function saveProjectViaApi(page: Page, name: string): Promise<string> {
  return await evalJs<string>(page, `(async () => {
    const tv = grok.shell.tv;
    const layout = tv.saveLayout();
    const proj = grok.shell.project;
    proj.name = '${name}';
    const saved = await grok.dapi.projects.save(proj);
    return saved.id;
  })()`);
}

export async function projectExists(page: Page, name: string): Promise<boolean> {
  return await evalJs<boolean>(page,
    `(async () => (await grok.dapi.projects.filter('name = "${name}"').first()) != null)()`);
}

export async function deleteProjectByName(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`).catch(() => {});
}

export async function pollUntilProjectExists(page: Page, name: string, timeout = 60_000) {
  await expect.poll(async () => projectExists(page, name),
    {timeout, intervals: [500, 1000, 2000, 5000]}).toBe(true);
}

export async function navigateToDashboards(page: Page) {
  await page.goto(BASE_URL + '/projects');
  await page.waitForSelector('.grok-gallery-grid', {timeout: 30_000});
}

export interface SpaceFixture {
  spaceId: string;
  spaceName: string;            
  fileName: string;             
  filePath: string;             
}

export interface ProvisionSpaceOptions {
  namePrefix: string;          
  sourceDirectory?: string;    
  fileName: string;            
  asName?: string;             
}

export interface ProvisionSpaceResult {
  blocked: boolean;
  reason: string;

  fixture: SpaceFixture | null;
}

export async function provisionSpaceFixture(
  page: Page, opts: ProvisionSpaceOptions): Promise<ProvisionSpaceResult> {
  const payload = {
    namePrefix: opts.namePrefix,
    sourceDirectory: opts.sourceDirectory ?? 'System:DemoFiles/',
    fileName: opts.fileName,
    asName: opts.asName ?? null,
  };
  return await evalJs<ProvisionSpaceResult>(page, `(async () => {
    const opts = ${JSON.stringify(payload)};
    try {
      if (typeof grok.dapi.spaces?.createRootSpace !== 'function')
        return {blocked: true, reason: 'spaces.createRootSpace not implemented on this env', fixture: null};
      const space = await grok.dapi.spaces.createRootSpace(opts.namePrefix + '-' + Date.now());
      const client = grok.dapi.spaces.id(space.id);
      try {
        const list = (await grok.dapi.files.list(opts.sourceDirectory)) || [];
        const src = list.find(x => x.name === opts.fileName);
        if (!src) throw new Error('source file not found in ' + opts.sourceDirectory + ': ' + opts.fileName);
        const target = opts.asName || opts.fileName;
        const bytes = await src.readAsBytes();
        await client.files.write(target, Array.from(bytes));
        if (!(await client.files.exists(target))) throw new Error('files.exists returned false after write: ' + target);
        // Resolve the global path via the space's children — the FileInfo
        // child carries the canonical fullPath ("<spaceName>:Files/<file>").
        // Fall back to the synthesized form if children.list misses it.
        let filePath = space.name + ':Files/' + target;
        try {
          const children = await client.children.list();
          const match = (children || []).find(c => c?.name === target && (c?.fullPath || '').includes('/' + target));
          if (match?.fullPath) filePath = match.fullPath;
        } catch (_) { /* fall through */ }
        return {
          blocked: false,
          reason: 'ok (' + filePath + ')',
          fixture: {
            spaceId: space.id,
            spaceName: space.name,
            fileName: target,
            filePath,
          },
        };
      } catch (e) {
        // Best-effort cleanup of the empty Space before returning a blocker.
        await grok.dapi.spaces.delete(space).catch(() => {});
        return {blocked: true, reason: String((e && e.message) || e).slice(0, 250), fixture: null};
      }
    } catch (e) {
      return {blocked: true, reason: 'Spaces fixture provisioning threw: ' + String((e && e.message) || e).slice(0, 250), fixture: null};
    }
  })()`);
}

export async function releaseSpaceFixture(page: Page, fixture: SpaceFixture | null) {
  if (!fixture) return;
  await evalJs(page, `(async () => {
    try {
      const sp = await grok.dapi.spaces.find(${JSON.stringify(fixture.spaceId)});
      if (sp) await grok.dapi.spaces.delete(sp);
    } catch (_) { /* best-effort cleanup */ }
  })()`).catch(() => {});
}
