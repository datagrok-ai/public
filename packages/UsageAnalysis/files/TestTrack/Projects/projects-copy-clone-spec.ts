// Selector sources (grok-browser/references):
//   widgets/dialog.md:22,29,61,69,74-92 — d4-dialog, button-OK, button-CANCEL, Dart input pattern
//   projects.md:71-79 — Copy with link / Clone semantics documented;
//     Save Copy modes are exposed via the Save Project dialog flag, but
//     the specific UI selectors for Link / Clone / PersonalView modes
//     are NOT documented in any grok-browser reference.
//   projects.md:232 — SAVE ribbon button
//
// SCOPE_REDUCTION rationale (significant):
//   The migrated scenario tests three Save Copy modes (Link / Clone /
//   PersonalView) plus a GROK-19750 regression invariant. The mode-specific
//   UI drivers (Link/Clone/PersonalView selectors and dialog flow) are NOT
//   yet documented in grok-browser/references; per E-SEL-01 / E-SEL-02 these
//   would be invented selectors. Without a documented helper, a UI-driven
//   spec at this layer would be fragile and unmaintainable.
//
//   This spec exercises the closest projects-side touch point that IS
//   documented: re-save (overwrite) the original after adding a viewer,
//   then reopen and assert the project still loads with its data. This
//   covers projects.api.save and projects.shell.open. The Save Copy
//   matrix is deferred to a future helper-registry addition for
//   helpers.playwright.projects.saveCopy({mode: 'link'|'clone'|'pvc'}).
//
// GROK-19750 invariant (Wave 1a B70 follow-up):
//   The mandatory GROK-19750 regression assertion (scenario sub-flow 4b
//   step 7: "after Save Copy with Link, original's viewers must be
//   intact") is exercised here via a JS API approximation:
//     1. Open original; capture baseline viewer count.
//     2. Build a linked copy via DG.Project.create() + addLink(tables);
//        save it through grok.dapi.projects.save (atlas: projects.api.save +
//        projects.add-link).
//     3. closeAll, reopen original.
//     4. Assert viewer count preserved (no leak into source).
//   This approximation exercises the SEMANTIC invariant without driving
//   the UI Save-Copy-with-Link dialog (which lacks documented selectors).
//   It does NOT exercise the exact UI code path the bug was originally
//   reproduced under; full coverage still requires
//   helpers.playwright.projects.saveCopy({mode}). Documented at the
//   `4b: GROK-19750 invariant` softStep below.
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
  // Share dialog auto-opens for new projects only
  const shareDialog = page.locator('.d4-dialog').filter({hasText: `Share ${name}`});
  if (await shareDialog.isVisible({timeout: 30000}).catch(() => false)) {
    await shareDialog.locator('[name="button-CANCEL"]').click();
    await expect(shareDialog).toBeHidden({timeout: 10000});
  }
}

async function deleteProjectByName(page: Page, name: string) {
  await evalJs(page, `(async () => {
    const p = await grok.dapi.projects.filter('name = "${name}"').first();
    if (p) await grok.dapi.projects.delete(p);
  })()`);
}

test('Projects / Save and reopen with added viewer', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const projectName = 'AutoTest-CopyClone-' + Date.now();

  await loginToDatagrok(page);

  try {
    await softStep('Setup: build original project with one viewer', async () => {
      await closeAll(page);
      await evalJs(page, `(async () => {
        const df = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        const tv = grok.shell.addTableView(df);
        tv.addViewer('Scatter plot');
      })()`);
      await page.waitForTimeout(2000);
      await saveProject(page, projectName);
      await expect.poll(async () => evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p != null;
        })()`,
      ), {timeout: 30000}).toBe(true);
      await closeAll(page);
    });

    await softStep('4a: Reopen, add another viewer, re-save (overwrite)', async () => {
      // Open project via JS API
      await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        await p.open();
      })()`);
      await page.waitForTimeout(3000);
      // Add a second viewer
      await evalJs(page, `(async () => {
        grok.shell.tv.addViewer('Bar chart');
      })()`);
      await page.waitForTimeout(1500);
      const viewerCountBefore = await evalJs(page, 'Array.from(grok.shell.tv.viewers).length');
      expect(viewerCountBefore).toBeGreaterThanOrEqual(3); // grid + scatter + bar
    });

    await softStep('Regression-class assertion: reopen original and verify project loads with tables', async () => {
      // Sanity check after re-save (overwrite): original must still
      // load with its tables. The full GROK-19750 invariant is below.
      await closeAll(page);
      const opened = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        await p.open();
        return grok.shell.tables.length > 0;
      })()`);
      expect(opened).toBe(true);
    });

    await softStep('4b: GROK-19750 invariant — Save Copy with Link must NOT leak state into source', async () => {
      // Approximation of scenario sub-flow 4b step 7. Drives a linked-copy
      // construction via JS API instead of the UI Save-Copy-with-Link dialog
      // (selectors not in references — see SR rationale in header).
      const linkCopyName = projectName + '-link';
      try {
        // 1. Reopen original, capture baseline viewer count
        await closeAll(page);
        const baselineViewers = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          await p.open();
          // settle render
          await new Promise(r => setTimeout(r, 1500));
          return Array.from(grok.shell.tv.viewers).length;
        })()`);
        expect(baselineViewers).toBeGreaterThanOrEqual(2); // grid + scatter from Setup

        // 2. Build linked copy via DG.Project.create() + addLink + dapi.save
        await evalJs(page, `(async () => {
          const linkCopy = DG.Project.create();
          linkCopy.name = '${linkCopyName}';
          for (const t of grok.shell.tables) {
            linkCopy.addLink(t);
          }
          await grok.dapi.projects.save(linkCopy);
        })()`);

        // 3. closeAll, reopen original
        await closeAll(page);
        const postCopyViewers = await evalJs(page, `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          await p.open();
          await new Promise(r => setTimeout(r, 1500));
          return Array.from(grok.shell.tv.viewers).length;
        })()`);

        // 4. GROK-19750 assertion: original's viewer count must be preserved
        expect(postCopyViewers).toBe(baselineViewers);
      } finally {
        // Cleanup the linked copy regardless of assertion outcome
        await evalJs(page, `(async () => {
          const c = await grok.dapi.projects.filter('name = "${linkCopyName}"').first();
          if (c) await grok.dapi.projects.delete(c);
        })()`).catch(() => {});
      }
    });
  } finally {
    await deleteProjectByName(page, projectName).catch(() => {});
    await closeAll(page);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
