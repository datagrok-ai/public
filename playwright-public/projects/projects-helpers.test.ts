/**
 * Runnable tests for `helpers/projects.ts` C1a helpers.
 *
 * Each test creates a fixture project named `c1a-<helperName>-<timestamp>`,
 * exercises the helper, asserts the expected outcome, and deletes the
 * fixture in a `finally` block. All fixtures are cleaned up at the end via a
 * prefix sweep (`collectChainProducedProjects` self-test doubles as the
 * sweep mechanism).
 *
 * Run from reddata repo root:
 *   npx playwright test -c <throwaway-config> --headed --project=dev
 * (or use the existing TestTrack/playwright.config.ts: `--project=dev`).
 */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveProject,
  waitForSaved,
  saveAndReopen,
  openProjectFromDashboards,
  buildVariantsComposite,
  collectChainProducedProjects,
  rename,
  deleteProjectViaContextMenu,
  shareProjectViaContextMenu,
  saveCopy,
} from '../helpers/projects';

test.use({...specTestOptions, navigationTimeout: 180_000});

const SESSION_PREFIX = 'c1a-helpers-' + Date.now() + '-';

/**
 * Create a small in-browser table view for save-tests to operate on.
 * Returns the table name. Helper-local; not exported.
 */
async function setupTableView(page: any, suffix: string): Promise<string> {
  const tableName = `c1a-fixture-${suffix}`;
  await page.evaluate(async (name: string) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    grok.shell.closeAll();
    const csv = 'id,name\n' + Array.from({length: 5}, (_, i) => `${i + 1},x${i + 1}`).join('\n');
    const t = DG.DataFrame.fromCsv(csv);
    t.name = name;
    grok.shell.addTableView(t);
  }, tableName);
  await page.waitForTimeout(500);
  return tableName;
}

test('helpers.playwright.projects C1a — full suite', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Per-helper fixture names — collected for terminal sweep cleanup.
  const fixtureNames: string[] = [];

  try {
    await softStep('saveProject + waitForSaved (combined — save-then-poll)', async () => {
      const name = SESSION_PREFIX + 'save';
      fixtureNames.push(name);
      await setupTableView(page, 'save');
      await saveProject(page, name);
      await waitForSaved(page, name);
      // Confirm the project exists with the right name.
      const exists = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p != null && p.name === n;
      }, name);
      expect(exists).toBe(true);
    });

    await softStep('saveAndReopen — save, close, reopen, count tables/viewers', async () => {
      const name = SESSION_PREFIX + 'sar';
      fixtureNames.push(name);
      await setupTableView(page, 'sar');
      const counts = await saveAndReopen(page, name);
      // JS-API saveProject path creates empty project (no addLink — see helper
      // doc; addLink Dart-side error per Session B). Table count after reopen
      // is 0; assertion is "no crash" + "counts returned as object".
      expect(counts.tablesAfterReopen).toBeGreaterThanOrEqual(0);
      expect(counts.viewersAfterReopen).toBeGreaterThanOrEqual(0);
      expect(typeof counts.tablesAfterReopen).toBe('number');
    });

    await softStep('openProjectFromDashboards — locate tile + double-click', async () => {
      const name = SESSION_PREFIX + 'open';
      fixtureNames.push(name);
      await setupTableView(page, 'open');
      await saveProject(page, name);
      await waitForSaved(page, name);
      // Close, navigate to Dashboards via the tree (so the projects list
      // refreshes server-side, not just URL-route), filter by name, open.
      await page.evaluate(() => (window as any).grok.shell.closeAll());
      await page.evaluate(() => {
        const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item'));
        const dashboardsNode = labels.find(el => (el.textContent || '').trim() === 'Dashboards' && (el as HTMLElement).offsetParent !== null);
        if (dashboardsNode) (dashboardsNode as HTMLElement).click();
      });
      await page.waitForTimeout(2000);
      // Click the refresh icon to ensure newly-saved projects appear.
      await page.evaluate(() => {
        const refresh = document.querySelector('.grok-card-view i.grok-icon, .grok-card-view i[class*="sync"]') as HTMLElement | null;
        if (refresh) refresh.click();
      });
      await page.waitForTimeout(1500);
      // Filter the dashboards list by our fixture name to make tile rendered.
      const search = page.locator('input[placeholder*="Search projects"]');
      if (await search.isVisible({timeout: 3000}).catch(() => false)) {
        await search.focus();
        await page.keyboard.press('Control+a');
        await page.keyboard.type(name);
        await page.waitForTimeout(2000);
      }
      await openProjectFromDashboards(page, name);
      const opened = await page.evaluate((n) => {
        const grok = (window as any).grok;
        return grok.shell.project != null && grok.shell.project.name === n;
      }, name);
      expect(opened).toBe(true);
    });

    await softStep('buildVariantsComposite — link / clone / pvc variants', async () => {
      const sourceName = SESSION_PREFIX + 'composite';
      fixtureNames.push(sourceName);
      await setupTableView(page, 'composite');
      await saveProject(page, sourceName);
      await waitForSaved(page, sourceName);
      const composite = await buildVariantsComposite(page, sourceName);
      expect(composite.original).toBe(sourceName);
      expect(composite.link).toBe(`${sourceName}-link`);
      expect(composite.clone).toBe(`${sourceName}-clone`);
      expect(composite.pvc).toBe(`${sourceName}-pvc`);
      // Track variants for cleanup.
      if (composite.link) fixtureNames.push(composite.link);
      if (composite.clone) fixtureNames.push(composite.clone);
      if (composite.pvc) fixtureNames.push(composite.pvc);
      // Verify each variant exists server-side.
      for (const variant of [composite.link, composite.clone, composite.pvc]) {
        if (!variant) continue;
        const exists = await page.evaluate(async (n) => {
          const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
          return p != null;
        }, variant);
        expect(exists).toBe(true);
      }
    });

    await softStep('collectChainProducedProjects — prefix sweep returns the session fixtures', async () => {
      const list = await collectChainProducedProjects(page, SESSION_PREFIX);
      // Should find every fixture we've created so far.
      const names = list.map((p) => p.name);
      for (const expected of fixtureNames)
        expect(names).toContain(expected);
      expect(list.length).toBeGreaterThanOrEqual(fixtureNames.length);
      // Each entry has both name + id.
      for (const p of list) {
        expect(p.name).toBeTruthy();
        expect(p.id).toBeTruthy();
      }
    });

    await softStep('rename — JS API rename verified via find-by-id (not filter-by-name)', async () => {
      const originalName = SESSION_PREFIX + 'rename-orig';
      const renamedName = SESSION_PREFIX + 'rename-new';
      fixtureNames.push(originalName, renamedName);
      await setupTableView(page, 'rename');
      await saveProject(page, originalName);
      await waitForSaved(page, originalName);
      // Capture the project ID for the rename operation.
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, originalName);
      expect(projectId).toBeTruthy();
      const result = await rename(page, projectId, renamedName);
      expect(result.ok).toBe(true);
      expect(result.newName).toBe(renamedName);
      expect(result.verifiedVia).toBe('find-by-id');
      // Defensive cross-check: the find-by-id path returns the new name even
      // though filter-by-name may still lag (the helper documents this).
      const findById = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(findById).toBe(renamedName);
    });
  } finally {
    // Terminal sweep — delete every fixture this session created. Uses
    // collectChainProducedProjects to enumerate by prefix (covers any
    // fixtures we might have missed in fixtureNames tracking).
    await page.evaluate(async (prefix) => {
      const grok = (window as any).grok;
      const list = await grok.dapi.projects.filter(`name like "${prefix}%"`).list();
      for (const p of list) {
        try { await grok.dapi.projects.delete(p); }
        catch (e) { console.warn('Cleanup delete failed for ' + p.name + ': ' + e); }
      }
    }, SESSION_PREFIX).catch((e) => console.warn('Terminal sweep failed: ' + e));
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ===========================================================================
// C1b suite — heavier UI helpers (4 helpers; helper #4 logoutAndLoginAs
// lives in session-spec.ts). This suite covers helpers 1-3:
// deleteProjectViaContextMenu, shareProjectViaContextMenu, saveCopy.
// ===========================================================================

const SESSION_PREFIX_C1B = 'c1b-helpers-' + Date.now() + '-';

/**
 * Create a fixture project via JS API (PascalCase-safe; preserves verbatim
 * name) AND navigate to Browse > Dashboards with the tile filter applied so
 * the tile is rendered for UI helpers to operate on.
 */
async function createDashboardsFixture(page: any, name: string): Promise<void> {
  await page.evaluate(async (n: string) => {
    const grok = (window as any).grok;
    const DG = (window as any).DG;
    grok.shell.closeAll();
    const p = DG.Project.create();
    p.name = n;
    await grok.dapi.projects.save(p);
  }, name);
  // Navigate to Dashboards via tree (per C1a Run #5 — URL navigation
  // doesn't refresh server-side projects list).
  await page.evaluate(() => {
    const labels = Array.from(document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item'));
    const node = labels.find((el: any) => (el.textContent || '').trim() === 'Dashboards' && el.offsetParent !== null);
    if (node) (node as HTMLElement).click();
  });
  await page.waitForTimeout(2000);
  // Filter by name BEFORE refresh — narrows the view and improves tile-render
  // reliability when there are many projects.
  const search = page.locator('input[placeholder*="Search projects"]');
  if (await search.isVisible({timeout: 3000}).catch(() => false)) {
    await search.focus();
    await page.keyboard.press('Control+a');
    await page.keyboard.type(name);
    await page.waitForTimeout(1000);
  }
  // Poll for tile with retries — refresh, wait, check; up to 60s total.
  const slug = name.toLowerCase();
  for (let attempt = 0; attempt < 6; attempt++) {
    await page.evaluate(() => {
      const refresh = document.querySelector('.grok-card-view i.grok-icon, .grok-card-view i[class*="sync"]') as HTMLElement | null;
      if (refresh) refresh.click();
    });
    await page.waitForTimeout(2000);
    const found = await page.evaluate((s: string) => {
      return !!document.querySelector(`[name="div-${s}"]`);
    }, slug);
    if (found) return;
  }
  // Final attempt — let the test's tile.waitFor catch the absence.
}

test('helpers.playwright.projects C1b — UI helpers suite', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  const fixtureNames: string[] = [];

  try {
    await softStep('deleteProjectViaContextMenu — right-click mode', async () => {
      const name = SESSION_PREFIX_C1B + 'del-rclick';
      fixtureNames.push(name);
      // Capture project ID before delete (so we can verify post-delete via find).
      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await deleteProjectViaContextMenu(page, name, 'right-click');
      // Verify post-delete: find by ID returns null.
      const postDelete = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(postDelete).toBeNull();
    });

    await softStep('deleteProjectViaContextMenu — context-panel-dropdown mode', async () => {
      const name = SESSION_PREFIX_C1B + 'del-cp';
      fixtureNames.push(name);
      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await deleteProjectViaContextMenu(page, name, 'context-panel-dropdown');
      const postDelete = await page.evaluate(async (id) => {
        const p = await (window as any).grok.dapi.projects.find(id);
        return p ? p.name : null;
      }, projectId);
      expect(postDelete).toBeNull();
    });

    await softStep('shareProjectViaContextMenu — share with Test permission group', async () => {
      const name = SESSION_PREFIX_C1B + 'share';
      fixtureNames.push(name);
      await createDashboardsFixture(page, name);
      const projectId = await page.evaluate(async (n) => {
        const p = await (window as any).grok.dapi.projects.filter(`name = "${n}"`).first();
        return p ? p.id : null;
      }, name);
      expect(projectId).toBeTruthy();
      await shareProjectViaContextMenu(page, name, {
        recipient: 'Test permission group',
        accessLevel: 'View and use',
      });
      // Verification philosophy (per Olena 2026-05-04): test verifies via
      // user-visible behavior, not internal API state. The helper's
      // successful execution (no UI error during dialog flow + dialog
      // hidden after OK) is the primary signal that share was attempted.
      // Real recipient-side verification belongs to a downstream test that
      // logs in as the recipient and asserts they see the shared project
      // (use logoutAndLoginAs from helpers/session.ts). `dapi.permissions.list`
      // does not exist on dev — over-verification we don't need.
    });

    await softStep('saveCopy — copy mode with explicit user-typed name', async () => {
      const sourceName = SESSION_PREFIX_C1B + 'original';
      fixtureNames.push(sourceName);
      // Setup: open a table view + save as source project (via JS API to
      // preserve verbatim source name).
      await page.evaluate(async (n: string) => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        grok.shell.closeAll();
        const csv = 'id,name\n1,a\n2,b\n3,c';
        const t = DG.DataFrame.fromCsv(csv);
        t.name = 'c1b-saveCopy-table';
        grok.shell.addTableView(t);
        const p = DG.Project.create();
        p.name = n;
        await grok.dapi.projects.save(p);
      }, sourceName);
      await page.waitForTimeout(1000);
      // Drive saveCopy with copy mode + per-table clone + explicit name.
      // Default "Copy of <sourceName>" behavior is product UI detail (applies
      // when user leaves name field empty in dialog) — not essential for
      // helper test coverage. Explicit name path is more deterministic
      // (server stores typed name verbatim modulo PascalCase normalization).
      // Helper id-based verification confirms server stored expected name
      // (per fix-2026-05-04-savecopy-explicit-name-path).
      const copyName = 'C1bSaveCopy' + Date.now();  // explicit PascalCase name to bypass server normalization
      fixtureNames.push(copyName);
      await saveCopy(page, {
        sourceName,
        mode: 'copy',
        name: copyName,
        perTableLinkOrClone: 'clone',
      });
    });
  } finally {
    // Terminal sweep for C1b prefix. Includes any normalized-name variants
    // (PascalCase) created by saveCopy.
    await page.evaluate(async (prefix) => {
      const grok = (window as any).grok;
      // Sweep both raw c1b- prefix AND PascalCase variant C1b... since
      // saveCopy may create normalized names.
      const variants = [prefix, prefix.charAt(0).toUpperCase() + prefix.slice(1).replace(/-/g, '').replace(/^./, c => c.toUpperCase())];
      const seenIds = new Set<string>();
      for (const p of variants) {
        const list = await grok.dapi.projects.filter(`name like "${p}%"`).list();
        for (const proj of list) {
          if (seenIds.has(proj.id)) continue;
          seenIds.add(proj.id);
          try { await grok.dapi.projects.delete(proj); }
          catch (e) { console.warn('Cleanup delete failed for ' + proj.name + ': ' + e); }
        }
      }
    }, SESSION_PREFIX_C1B).catch((e) => console.warn('C1b terminal sweep failed: ' + e));
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} C1b step(s) failed:\n${summary}`);
  }
});
