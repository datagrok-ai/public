/* ---
sub_features_covered: [projects.upload, projects.api.save]
--- */
// Selector sources (grok-browser/references):
//   widgets/dialog.md:22 — `.d4-dialog` root selector
//   widgets/dialog.md:29,61 — `[name="button-OK"]` (preferred over .ui-btn-ok)
//   widgets/dialog.md:69 — `[name="button-CANCEL"]`
//   widgets/dialog.md:74-92 — Dart inputs require focus + keyboard.type, NOT .fill()
//   projects.md:15,48 — File → Save Project / Ctrl+S opens the Save dialog;
//     in Playwright the SAVE ribbon button is the reliable trigger (Ctrl+S is
//     intercepted by the browser's native save). For an unsaved project the
//     ribbon button opens the dialog (per projects.md:232 it only skips the
//     dialog on re-save of an already-saved project).
//   projects.md:28 — Preferred default is grok.dapi for verification, not UI scraping
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

async function deleteProject(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`);
}

test('Projects / Upload Project', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await closeAll(page);

  const projectName = 'AutoTest-UploadProject-' + Date.now();

  await softStep('Save demog with 3 viewers, Data Sync ON, dismiss Share dialog', async () => {
    try {
      // Setup steps 1-4: open demog.csv, add scatter/bar/line viewers
      await evalJs(page, `(async () => {
        const df = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        const tv = grok.shell.addTableView(df);
        tv.addViewer('Scatter plot');
        tv.addViewer('Bar chart');
        tv.addViewer('Line chart');
      })()`);

      // technical: wait for viewers to mount before driving the dialog
      await page.waitForTimeout(2000);
      const viewerCount = await evalJs(page,
        `Array.from(grok.shell.tv.viewers).length`);
      expect(viewerCount).toBeGreaterThanOrEqual(4); // grid + 3 added

      // Step 1: Trigger Save Project via the SAVE ribbon button (Ctrl+S is
      // intercepted by the browser's native save in Playwright; ribbon button
      // is the reliable trigger for an unsaved project — see header notes).
      await page.click('button:has-text("SAVE"), .ui-btn:has-text("SAVE")');
      const saveDialog = page.locator('.d4-dialog').first();
      await saveDialog.waitFor({timeout: 15000});

      // Step 2: leave Data Sync toggled ON (default per projects.md:61, no action needed)

      // Step 3: fill project name + click OK
      // Dart input — focus + Ctrl+A + type, NOT .fill() (dialog.md:74-92)
      const nameInput = saveDialog.locator('input[type="text"]').first();
      await nameInput.focus();
      await page.keyboard.press('Control+a');
      await page.keyboard.type(projectName);

      const okBtn = saveDialog.locator('[name="button-OK"]');
      await okBtn.click();

      // Verify project exists on server (canonical verification per projects.md:28)
      await expect.poll(async () => evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p !== null;
        })()`,
      ), {timeout: 30000}).toBe(true);

      // Step 4: Share dialog auto-opens; click Cancel; verify it closes
      const shareDialog = page.locator('.d4-dialog').filter({hasText: `Share ${projectName}`});
      await shareDialog.waitFor({timeout: 30000});
      await shareDialog.locator('[name="button-CANCEL"]').click();
      await expect(shareDialog).toBeHidden({timeout: 10000});

      // Step 5: Close All — handled in finally{}
    } finally {
      await deleteProject(page, projectName).catch(() => {});
      await closeAll(page);
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
