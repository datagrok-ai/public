import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, fn: string): Promise<any> {
  return page.evaluate(fn);
}

async function projectExists(page: Page, name: string): Promise<boolean> {
  return await evalJs(page,
    `(async () => (await grok.dapi.projects.filter('name = "${name}"').first()) !== null)()`);
}

async function waitProjectGone(page: Page, name: string, timeoutMs: number): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    if (!(await projectExists(page, name))) return true;
    await page.waitForTimeout(500);
  }
  return false;
}

test('Projects / Deleting', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const t = Date.now();
  const projectNames = [`AutoTestDelete1${t}`, `AutoTestDelete2${t}`];

  await loginToDatagrok(page);

  await evalJs(page, `(async () => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
  })()`);

  try {
    await softStep('Setup: Create test projects via API', async () => {
      for (const name of projectNames) {
        await evalJs(page, `(async () => {
          const project = DG.Project.create();
          project.name = ${JSON.stringify(name)};
          project.friendlyName = ${JSON.stringify(name)};
          const table = grok.data.demo.demog(20);
          table.name = 'data_' + ${JSON.stringify(name)};
          const tableInfo = table.getTableInfo();
          project.addChild(tableInfo);
          await grok.dapi.tables.uploadDataFrame(table);
          await grok.dapi.tables.save(tableInfo);
          await grok.dapi.projects.save(project);
        })()`);
        await page.waitForTimeout(1500);
      }
      await evalJs(page, `grok.shell.closeAll()`);
    });

    await softStep('Step 1: Find the projects from previous setup', async () => {
      // Navigate directly to /projects to render the Dashboards gallery view
      await page.goto(`${process.env.DATAGROK_URL ?? 'http://localhost:8888'}/projects`);
      await page.waitForFunction(() =>
        document.querySelector('.grok-card-view') !== null ||
        document.querySelector('.grok-gallery-grid-item') !== null, null, {timeout: 30_000});

      // Gallery may not have the latest API-created projects until we click the
      // refresh icon (top of the Dashboards toolbar) — the Dart view caches.
      for (const name of projectNames) {
        const tile = page.locator(`[name="div-${name}"]`);
        let visible = await tile.count() > 0;
        for (let i = 0; i < 8 && !visible; i++) {
          await page.evaluate(() => {
            const refreshIcons = Array.from(document.querySelectorAll<HTMLElement>(
              'i.fa-sync, i.fal.fa-sync, i.far.fa-sync, [name="icon-sync"]'));
            // Pick the visible refresh icon inside the gallery view (skip navbar-level ones)
            const visIcon = refreshIcons.find((el) => el.offsetParent !== null);
            if (visIcon) visIcon.click();
          });
          await page.waitForTimeout(1500);
          visible = await tile.count() > 0;
        }
        await tile.waitFor({timeout: 15_000});
        expect(await projectExists(page, name)).toBe(true);
      }
    });

    for (const name of projectNames) {
      await softStep(`Steps 2-3: Right-click ${name} > Delete Project > DELETE`, async () => {
        // Make sure no leftover dialog from previous iteration
        await page.evaluate(() => {
          document.querySelectorAll('.d4-dialog').forEach((d) => {
            const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
            if (cancel) cancel.click();
          });
        });
        await page.waitForTimeout(300);

        // Right-click the tile via dispatched contextmenu event
        await page.evaluate((tileName) => {
          const tile = document.querySelector(`[name="div-${tileName}"]`) as HTMLElement | null;
          if (!tile) throw new Error('tile not found: ' + tileName);
          tile.scrollIntoView({block: 'center'});
          const r = tile.getBoundingClientRect();
          tile.dispatchEvent(new MouseEvent('contextmenu', {
            bubbles: true, cancelable: true, button: 2,
            clientX: Math.round(r.left + r.width / 2),
            clientY: Math.round(r.top + r.height / 2),
            view: window,
          }));
        }, name);

        // Click the Delete Project menu item
        await page.locator('[name="div-Delete-Project"]').first().waitFor({timeout: 5000});
        await page.evaluate(() => {
          (document.querySelector('[name="div-Delete-Project"]') as HTMLElement).click();
        });

        // Wait for the most recently opened confirm dialog (the one for THIS project)
        await page.waitForFunction((expected) => {
          const dialogs = Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog'));
          return dialogs.some((d) => d.textContent?.includes(`Delete project "${expected}"`));
        }, name, {timeout: 5000});

        // Verify dialog title/body
        const dialogInfo = await page.evaluate((expected) => {
          const dialogs = Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog'));
          const d = dialogs.find((x) => x.textContent?.includes(`Delete project "${expected}"`));
          if (!d) return null;
          return {
            title: d.querySelector('.d4-dialog-title')?.textContent?.trim() ?? null,
            body: d.querySelector('.d4-dialog-contents')?.textContent?.trim() ?? null,
          };
        }, name);
        expect(dialogInfo?.title).toBe('Are you sure?');
        expect(dialogInfo?.body ?? '').toContain(`Delete project "${name}"`);

        // Click the DELETE button inside this specific dialog using Playwright's
        // native click (not synthetic events) — Dart command-bar buttons require
        // a real pointer click to fire the OK action.
        const dialogLocator = page.locator('.d4-dialog').filter({hasText: `Delete project "${name}"`});
        const deleteBtn = dialogLocator.locator('[name="button-DELETE"]');
        await deleteBtn.waitFor({state: 'visible', timeout: 5000});
        await expect(deleteBtn).toBeEnabled({timeout: 5000});
        await deleteBtn.click();

        // Wait for THIS project's confirm dialog to close (not all dialogs)
        await page.waitForFunction((expected) => {
          const dialogs = Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog'));
          return !dialogs.some((d) => d.textContent?.includes(`Delete project "${expected}"`));
        }, name, {timeout: 30_000});

        // Wait for the API to actually have removed it (tile may linger briefly)
        const gone = await waitProjectGone(page, name, 30_000);
        expect(gone).toBe(true);
      });
    }

    await softStep('Step 4: Verify projects are deleted (UI tiles gone, API filter empty)', async () => {
      // Refresh dashboards view by clicking the tree node again (forces gallery re-render)
      await page.evaluate(() => {
        const labels = Array.from(document.querySelectorAll<HTMLElement>(
          '.d4-tree-view-group-label, .d4-tree-view-item-label'));
        const dash = labels.find((l) => l.textContent?.trim() === 'Dashboards');
        if (dash) dash.click();
      });
      await page.waitForTimeout(1500);

      for (const name of projectNames) {
        // API reports gone (authoritative)
        expect(await projectExists(page, name)).toBe(false);
        // Tile is also gone from gallery DOM
        const tileCount = await page.locator(`[name="div-${name}"]`).count();
        expect(tileCount).toBe(0);
      }
    });
  } finally {
    // Idempotent cleanup: best-effort delete in case any soft step left a project behind
    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`).catch(() => {});
    }
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
