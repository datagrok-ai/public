/* ---
sub_features_covered: [projects.api.search, projects.api.delete]
--- */
// Selector sources (grok-browser/references):
//   widgets/dialog.md:22,29,61,69,74-92 — d4-dialog, button-OK, button-CANCEL, Dart input pattern
//   projects.md:171-184 — Delete via JS API (grok.dapi.projects.delete)
//   projects.md:232 — SAVE ribbon button
//
// SCOPE_REDUCTION rationale:
//   Original scenario uses right-click → "Delete project" + DELETE confirmation
//   button (no name= attribute per projects.md:178). Spec uses JS API delete
//   per projects.md:181-184 — this is the documented preferred default; the
//   right-click path adds no projects-side coverage beyond the delete API.
//
//   This is the chain's terminal cleanup scenario (must_run_last per chain
//   rev 2). For a standalone Validator run, it creates its own dummy projects
//   to delete (so it doesn't depend on upstream chain state).
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function evalJs(page: Page, script: string): Promise<any> {
  return page.evaluate(script);
}

async function closeAll(page: Page) {
  await evalJs(page, 'grok.shell.closeAll()');
}

async function saveProject(page: Page, name: string) {
  await page.click('button:has-text("SAVE"), .ui-btn:has-text("SAVE")');
  const dialog = page.locator('.d4-dialog').first();
  await dialog.waitFor({timeout: 15000});
  const nameInput = dialog.locator('input[type="text"]').first();
  await nameInput.focus();
  await page.keyboard.press('Control+a');
  await page.keyboard.type(name);
  await dialog.locator('[name="button-OK"]').click();
  const shareDialog = page.locator('.d4-dialog').filter({hasText: `Share ${name}`});
  if (await shareDialog.isVisible({timeout: 30000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }
}

test('Projects / Deleting', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectNames = [
    `AutoTest-Delete-1-${stamp}`,
    `AutoTest-Delete-2-${stamp}`,
  ];

  await loginToDatagrok(page);

  try {
    await softStep('Setup: Create test projects via SAVE-button flow', async () => {
      for (const name of projectNames) {
        await closeAll(page);
        await evalJs(page, `(async () => {
          grok.shell.addTableView(grok.data.demo.demog());
        })()`);
        await page.waitForTimeout(2000);
        await saveProject(page, name);
        await expect.poll(async () => evalJs(page,
          `(async () => {
            const p = await grok.dapi.projects.filter('name = "${name}"').first();
            return p != null;
          })()`,
        ), {timeout: 30000}).toBe(true);
        await closeAll(page);
      }
    });

    await softStep('Case 1: Find test projects on server', async () => {
      for (const name of projectNames) {
        const exists = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          return p != null;
        })()`);
        expect(exists).toBe(true);
      }
    });

    await softStep('Case 2-3: Delete projects via JS API', async () => {
      for (const name of projectNames) {
        await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          if (p) await grok.dapi.projects.delete(p);
        })()`);
      }
      await page.waitForTimeout(2000);
    });

    await softStep('Case 4: Verify projects are deleted', async () => {
      for (const name of projectNames) {
        const exists = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${name}"').first();
          return p != null;
        })()`);
        expect(exists).toBe(false);
      }
    });
  } finally {
    // Best-effort cleanup of anything still lingering
    for (const name of projectNames) {
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${name}"').first();
        if (p) await grok.dapi.projects.delete(p);
      })()`).catch(() => {});
    }
    await closeAll(page);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
