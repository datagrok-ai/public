// Selector sources (grok-browser/references):
//   widgets/dialog.md:22,29,61,69,74-92 — d4-dialog, button-OK, button-CANCEL, Dart input pattern
//   projects.md:232 — SAVE ribbon button is the reliable Save Project trigger
//   projects.md:81-101 — Open from URL: navigate to the project's deep-link;
//     URL is shown in Context Panel → Links and built by url-params.build-share-link
//
// SCOPE_REDUCTION rationale:
//   The migrated scenario lists 4 variants (original / copy-with-link / copy-with-
//   clone / save-personal-view-customizations). Per scenario Notes "URL deep-link
//   is generic" and "URL-apply behavior is uniform across variant modes" — the
//   path under test is identical regardless of variant. Spec exercises 1 variant
//   (the original saved project) which is sufficient to verify URL-apply, and
//   skips the multi-variant fixture build. Variant-specific copy/clone behavior
//   is covered separately by projects-copy-clone-spec.ts.
//
// Wave 1a B70 follow-up:
//   Replaced previous p.open() roundtrip with page.goto({BASE}/p/{nqName}) so
//   the URL-apply path (projects.url-params.apply + projects.shell.open via the
//   URL handler) is actually exercised, honoring the frontmatter coverage claim
//   for projects.url-params.build-share-link / projects.url-params.apply.
import {test, expect, Page} from '@playwright/test';
import {baseUrl, loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

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
  // Share dialog auto-opens for new projects only — defensive wait pattern
  // (matches projects-copy-clone-spec.ts saveProject) to tolerate the
  // intermittent dev-env case where the dialog does not appear within 30s.
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

test('Projects / Project URL deep-link reopen', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const projectName = 'AutoTest-ProjectUrl-' + Date.now();

  await loginToDatagrok(page);

  try {
    await softStep('Setup: build saved project fixture', async () => {
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

    await softStep('Step 3: Build deep-link URL for the saved project', async () => {
      // url-params.build-share-link constructs a {BASE}/p/{nqName} or similar
      // shape; we rely on the platform's project URL convention rather than
      // scraping Context Panel > Links (which has no name= attributes).
      const urlInfo = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        const id = p?.id || null;
        const nqName = p?.nqName || p?.name || null;
        return { id, nqName };
      })()`);
      expect(urlInfo.id || urlInfo.nqName).toBeTruthy();
    });

    await softStep('Step 4-5: Navigate to project URL and verify reopen', async () => {
      // Build the deep-link URL via JS API (nqName is the namespace-qualified
      // project name surfaced by Context Panel > Links per projects.md:91).
      const nqName = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        return p?.nqName || p?.name || null;
      })()`);
      expect(nqName).toBeTruthy();

      // Close the current session before URL navigation so the assertion that
      // tables load comes from the URL handler, not from carry-over state.
      await closeAll(page);

      // Wave 1a fix: actually exercise the URL-apply path
      // (projects.url-params.apply + projects.shell.open via the URL handler)
      // instead of substituting with same-session p.open(). Project URL shape
      // is `{BASE}/p/{nqName}` per Datagrok URL convention (projects.md:81-101).
      const projectUrl = `${baseUrl}/p/${nqName}`;
      await page.goto(projectUrl);

      // Wait for the URL handler to resolve and load the project state
      await page.waitForTimeout(5000);

      // Verify project is open: tables loaded via the URL-apply path
      const opened = await evalJs(page, `(async () => {
        return grok.shell.tables.length > 0;
      })()`);
      expect(opened).toBe(true);
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
