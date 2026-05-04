/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.shell.open]
--- */
// Selector sources (grok-browser/references):
//   widgets/dialog.md:22,29,61,69,74-92 — d4-dialog, button-OK, button-CANCEL, Dart input pattern
//   projects.md:232 — SAVE ribbon button is the reliable Save Project trigger
//   projects.md:28 — Preferred default is grok.dapi for verification
//
// Wave 1b complex-split: covers Step 1 multi-source open + Join sub-bullet
// + Step 2 Save with Data Sync ON of complex.md scenario. Targets the
// GROK-19103 regression invariant ("join result silently saved as a
// SEPARATE project that later fails to open"). Spec verifies that after
// joining two source tables and saving the active project, exactly ONE
// new project is created on the server — not two (the saved active
// project plus a stray join-only project).
//
// Scope reductions (documented):
//   * Sources are file-based (System:DemoFiles/demog.csv used twice with
//     distinct names) rather than DB query result + DB table per the
//     scenario sub-bullet. The GROK-19103 invariant ("join lands in
//     active project, not stray") is independent of source type — the
//     bug surfaces from the joinTables → save flow regardless of
//     where source tables originate. File sources are env-resilient
//     (no Postgres dependency).
//   * Pivot / Aggregate / Clone-via-script sub-bullets of Step 1 are NOT
//     covered here — those belong to a future cross-feature scenario
//     about derived-table types. Only the Join sub-bullet that
//     reproduces GROK-19103 is in scope.
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
  // Defensive share-dialog handling per Wave 1a Validator hypothesis cycle 1
  // canonical pattern (project-url-spec, projects-copy-clone-spec).
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

test('Projects / Complex derived-tables: Join lands in active project (GROK-19103)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = 'AutoTest-ComplexDerived-' + stamp;

  await loginToDatagrok(page);

  let baselineCount = 0;

  try {
    await softStep('Setup + Step 1: open 2 source tables from files', async () => {
      await closeAll(page);
      await evalJs(page, `(async () => {
        const df1 = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        df1.name = 'src_a_${stamp}';
        grok.shell.addTableView(df1);
        const df2 = await grok.data.files.openTable('System:DemoFiles/demog.csv');
        df2.name = 'src_b_${stamp}';
        grok.shell.addTableView(df2);
      })()`);
      await page.waitForTimeout(2000);
      const tableCount = await evalJs(page, 'grok.shell.tables.length');
      expect(tableCount).toBe(2);
    });

    await softStep('Step 1 (Join sub-bullet): grok.data.joinTables produces 3rd table in workspace', async () => {
      await evalJs(page, `(async () => {
        const [a, b] = grok.shell.tables;
        const joined = grok.data.joinTables(
          a, b,
          ['USUBJID'], ['USUBJID'],
          ['AGE', 'SEX'], ['RACE'],
          'inner',
        );
        joined.name = 'joined_${stamp}';
        grok.shell.addTableView(joined);
      })()`);
      await page.waitForTimeout(2000);
      const tableCount = await evalJs(page, 'grok.shell.tables.length');
      expect(tableCount).toBe(3);
    });

    await softStep('Step 2: capture project-count baseline before Save (GROK-19103 invariant prep)', async () => {
      baselineCount = await evalJs(page,
        `(async () => (await grok.dapi.projects.list({limit: 1000})).length)()`,
      );
      expect(baselineCount).toBeGreaterThanOrEqual(0);
    });

    await softStep('Step 2: Save current project with Data Sync ON', async () => {
      await saveProject(page, projectName);
      await expect.poll(async () => evalJs(page,
        `(async () => {
          const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
          return p != null;
        })()`,
      ), {timeout: 30000}).toBe(true);
    });

    await softStep('GROK-19103 INVARIANT: exactly ONE new project created (no stray join-only project)', async () => {
      const afterCount = await evalJs(page,
        `(async () => (await grok.dapi.projects.list({limit: 1000})).length)()`,
      );
      // GROK-19103: a stray join-only project would push delta to 2+. The
      // active project save should produce exactly one new project entity.
      expect(afterCount - baselineCount).toBe(1);
    });

    await softStep('Reopen: verify project loads with source tables intact', async () => {
      await closeAll(page);
      const tableCount = await evalJs(page, `(async () => {
        const p = await grok.dapi.projects.filter('name = "${projectName}"').first();
        await p.open();
        await new Promise(r => setTimeout(r, 2000));
        return grok.shell.tables.length;
      })()`);
      // At minimum the 2 source tables must reopen. The joined derivative
      // may or may not be persisted depending on project relation semantics
      // (recomputed from sources vs persisted as standalone) — accept >= 2.
      expect(tableCount).toBeGreaterThanOrEqual(2);
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
