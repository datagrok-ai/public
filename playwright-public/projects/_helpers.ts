// Shared helpers for Projects regression specs.
// Auth comes from playwright.config.ts (`storageState: 'e2e/.auth.json'`,
// written by ../e2e/global-setup.ts at suite start). Do NOT override
// storageState here — a `path.resolve(__dirname, ...)` override breaks
// under the CI runner's esbuild transpile path and the test lands on the
// login page (build #45 page snapshot: "Login or Email / Token is empty").
// Helpers transcribed from .claude/skills/grok-browser/references/projects.md
// and .claude/skills/grok-browser/references/widgets/dialog.md.
import {Page, expect} from '@playwright/test';

export const BASE_URL = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';

export const projectsTestOptions = {
  viewport: {width: 1920, height: 1080},
  launchOptions: {args: ['--window-size=1920,1080', '--window-position=0,0']},
  actionTimeout: 15_000,
  navigationTimeout: 120_000,
};

export async function evalJs<T = any>(page: Page, script: string): Promise<T> {
  return page.evaluate(script) as Promise<T>;
}

export async function gotoApp(page: Page) {
  await page.goto(BASE_URL);
  await page.locator('[name="Browse"]').first().waitFor({state: 'visible', timeout: 60_000});
  // Cold CI Datlas: the Browse sidebar attaches before `#grok-preloader`
  // detaches; while the preloader is up it covers the rootDiv and intercepts
  // every click. Mirrors what connections/helpers.ts goHome() does.
  await page.waitForFunction(
    () => document.querySelector('#grok-preloader, .grok-preloader') == null,
    null, {timeout: 90_000},
  ).catch(() => { /* best-effort — some flows tolerate a lingering preloader */ });
  await page.waitForTimeout(500);
  await page.addStyleTag({content: `
    .d4-tooltip { display: none !important; }
    #grok-preloader, .grok-preloader { pointer-events: none !important; }
  `}).catch(() => {});
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

// Drive the Save Project dialog. Per the platform's own UITests pattern
// (scripts-layout.test.ts:415) the toolbar Save button is reachable as
// `button[name="button-Save"]:visible` — explicit `button` element +
// `:visible` filter + `.first()`. Bare `[name="button-Save"]` returns
// undefined in Playwright (DOM has hidden duplicates that confuse the
// resolver).
//
// Match the dialog in two stages:
//   * `[name="dialog-Save-project"]` — name attribute (preferred when present)
//   * `.d4-dialog:has-text("Save project")` — text-based fallback
export async function saveProjectViaDialog(page: Page, name: string) {
  // Wait for the toolbar to settle — under Playwright the toolbar render
  // can lag the dataframe render by several seconds. Also guard against
  // overlay banners that intercept pointer events on the toolbar.
  const saveBtn = page.locator('button[name="button-Save"]:visible').first();
  await saveBtn.waitFor({timeout: 30_000, state: 'visible'});
  await page.waitForTimeout(500);  // settle render
  // Try regular click first; fall back to JS DOM click if intercepted.
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
  // Save dialog auto-opens a Share dialog after OK. Match by name-attribute
  // prefix `dialog-Share-*` (preferred) or by title-text prefix `^Share /`
  // — the user-typed `name` is server-normalized (PascalCase, no dashes) in
  // the share dialog title, so exact `Share ${name}` text match doesn't fire.
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

// Save via JS API (no UI). Use for setup state where the scenario doesn't
// own the Save dialog.
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

// Spaces fixture provisioning. Creates a root Space and physically uploads a
// file into the space's own storage via `client.files.write(name, bytes)`,
// making it openable via `<spaceName>:Files/<name>` — the path-based open
// flow scenarios use.
//
// Path format note: a Space exposes a child DataConnection literally named
// `Files` (full name `<spaceName>:Files`); the file written via
// `client.files.write(name, ...)` ends up at `<spaceName>:Files/<name>`.
// `Spaces:<spaceName>/<name>` does NOT resolve — the global `Spaces:` namespace
// is not the storage root. Verified empirically 2026-05-07 against dev:
// `client.children.list()` returns the FileInfo with
// `fullPath = "<spaceName>:Files/<fileName>"`, and that path opens cleanly
// through the canonical `OpenFile` opener.
//
// Linked-entity provisioning (`addEntity(uuid, link=true)`) is intentionally
// omitted: on dev/sandbox the test account doesn't own System:DemoFiles
// FileInfos, so the server rejects the call with `Permission denied to change
// entity`. Physical copy bypasses this — the bytes are uploaded fresh, the
// space owns the resulting file outright.
export interface SpaceFixture {
  spaceId: string;
  spaceName: string;            // server-normalized (PascalCase, no dashes)
  fileName: string;             // filename inside the space (e.g. 'cars.csv')
  filePath: string;             // `<spaceName>:Files/<fileName>` — usable with OpenFile
}

export interface ProvisionSpaceOptions {
  namePrefix: string;          // suffixed with `-${Date.now()}` to avoid collisions
  sourceDirectory?: string;    // default 'System:DemoFiles/'
  fileName: string;            // file in sourceDirectory to copy into the space
  asName?: string;             // override stored name in the space (default = fileName)
}

export interface ProvisionSpaceResult {
  blocked: boolean;
  reason: string;
  // Populated when the Space itself was created, even if the copy failed.
  // Always pass to releaseSpaceFixture for cleanup, regardless of `blocked`.
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
