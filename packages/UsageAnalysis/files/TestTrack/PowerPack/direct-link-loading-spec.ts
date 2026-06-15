/* ---
sub_features_covered: [powerpack.dashboards, powerpack.lifecycle.init, powerpack.welcome.view]
--- */
import {test, expect, BrowserContext, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, baseUrl} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
async function injectTokenInNewContext(ctx: BrowserContext): Promise<Page> {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test` so the wrapper exchanges the dev key from ~/.grok/config.yaml.');
  const u = new URL(baseUrl);
  await ctx.addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  const page = await ctx.newPage();
  await page.goto(baseUrl + '/oauth/');
  await page.evaluate((t) => window.localStorage.setItem('auth', t), token);
  return page;
}
async function getBoxOf(page: Page, selector: string): Promise<{found: boolean; w: number; h: number}> {
  return await page.evaluate((sel) => {
    const el = document.querySelector(sel) as HTMLElement | null;
    if (!el) return {found: false, w: 0, h: 0};
    const r = el.getBoundingClientRect();
    return {found: true, w: r.width, h: r.height};
  }, selector);
}
async function waitForPreloaderGone(page: Page, timeout = 120_000): Promise<void> {
  await page.waitForFunction(() => document.querySelector('.grok-preloader') == null, null, {timeout});
}
test('PowerPack: Direct-link entry renders loading window fully (GROK-18721 regression)', async ({browser, page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  const stamp = Date.now();
  const projectName = `DirectLinkLoading${stamp}`;
  let projectId: string | null = null;
  let tableInfoId: string | null = null;
  let ownerLogin: string | null = null;
  let secondaryContext: BrowserContext | null = null;
  try {
    await softStep('Setup: create a project with a known direct-link URL (Setup step 2 of scenario)', async () => {
      await page.evaluate(async () => {
        const grok = (window as any).grok;
        document.body.classList.add('selenium');
        grok.shell.settings.showFiltersIconsConstantly = true;
        grok.shell.windows.simpleMode = true;
        try { grok.shell.closeAll(); } catch (_) {}
        const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
        grok.shell.addTableView(df);
        await new Promise<void>((resolve) => {
          const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
          setTimeout(resolve, 3000);
        });
      });
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1000);
      const saved = await page.evaluate(async (n) => {
        const grok = (window as any).grok;
        const DG = (window as any).DG;
        const df = grok.shell.t;
        const project = DG.Project.create();
        project.name = n;
        const ti = df.getTableInfo();
        project.addChild(ti);
        await grok.dapi.tables.uploadDataFrame(df);
        await grok.dapi.tables.save(ti);
        const tv = grok.shell.tv;
        const layout = tv?.saveLayout?.();
        if (layout) {
          project.addChild(layout);
          await grok.dapi.layouts.save(layout);
        }
        await grok.dapi.projects.save(project);
        return {
          projectId: project.id,
          tableInfoId: ti.id,
          ownerLogin: grok.shell.user?.login ?? null,
        };
      }, projectName);
      projectId = saved.projectId;
      tableInfoId = saved.tableInfoId;
      ownerLogin = saved.ownerLogin;
      expect(projectId).toBeTruthy();
      expect(ownerLogin).toBeTruthy();
    });
    const directLinkPath = `/p/${ownerLogin}.${projectName}`; // direct-link URL form: /p/<owner>.<project>

    // Scenario 1: Direct-link entry (fresh browser context).
    await softStep('Scenario 1 Step 1+2: open fresh context and navigate to direct-link URL', async () => {
      // Fresh context (no shared cookies/localStorage) so PowerPack's powerPackInit runs from scratch.
      secondaryContext = await browser.newContext({
        viewport: specTestOptions.viewport,
      });
      const freshPage = await injectTokenInNewContext(secondaryContext);
      await freshPage.goto(baseUrl + directLinkPath);
      (secondaryContext as any)._freshPage = freshPage;
    });
    await softStep('Scenario 1 Step 3: observe PowerPack loading window during page load (no zero-dimension cropping)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;
      // GROK-18721 invariant: while loading, the welcome view must not render as a degenerate box (one dim
      // zero, the other non-zero — the pre-fix "cropped" presentation). Snapshot rects across the load window.
      const snapshots: Array<{t: number; box: {found: boolean; w: number; h: number}}> = [];
      for (let i = 0; i < 30; i++) {
        const box = await getBoxOf(freshPage, '.power-pack-welcome-view');
        snapshots.push({t: i * 200, box});
        const gridMounted = await freshPage.evaluate(() =>
          !!document.querySelector('[name="viewer-Grid"] canvas'));
        if (gridMounted) break;
        await freshPage.waitForTimeout(200);
      }
      const degenerate = snapshots.filter((s) =>
        s.box.found &&
        ((s.box.w === 0 && s.box.h > 0) || (s.box.h === 0 && s.box.w > 0)));
      expect(degenerate, `GROK-18721 invariant: no cropped loading window with one-zero-one-nonzero dimensions. Offending snapshots: ${JSON.stringify(degenerate)}`).toEqual([]);
      // Also: a mounted welcome view must clear ≥100px in both dimensions (no 1-px sliver crop).
      const tinyBoxes = snapshots.filter((s) =>
        s.box.found && s.box.w > 0 && s.box.h > 0 &&
        (s.box.w < 100 || s.box.h < 100));
      expect(tinyBoxes, `GROK-18721 invariant: welcome view must not render at sub-100px in either dimension. Offending snapshots: ${JSON.stringify(tinyBoxes)}`).toEqual([]);
    });
    await softStep('Scenario 1 Step 4: wait for the load to complete (preloader gone + grid mounted)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;
      await waitForPreloaderGone(freshPage);
      await freshPage.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await freshPage.waitForTimeout(1500);
    });
    await softStep('Scenario 1 Step 5: verify post-load rendering (grid has non-zero dimensions, no zombie welcome fragments)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;
      // Project view rendered correctly: grid present with non-zero box.
      const gridBox = await getBoxOf(freshPage, '[name="viewer-Grid"]');
      expect(gridBox.found).toBe(true);
      expect(gridBox.w).toBeGreaterThan(100);
      expect(gridBox.h).toBeGreaterThan(100);
      // The dataframe under the grid is the demog table (the direct-link
      // target). Confirm by reading the table info via JS API.
      const tableMeta = await freshPage.evaluate(() => {
        const grok = (window as any).grok;
        const df = grok.shell.tv?.dataFrame;
        return df ? {name: df.name, rowCount: df.rowCount, colCount: df.columns.length} : null;
      });
      expect(tableMeta).not.toBeNull();
      expect(tableMeta!.rowCount).toBeGreaterThan(0);
      expect(tableMeta!.colCount).toBeGreaterThan(0);
      // Welcome view must have yielded — removed, hidden, or zero-box (not occupying the project-view area).
      const welcomeStillActive = await freshPage.evaluate(() => {
        const w = document.querySelector('.power-pack-welcome-view') as HTMLElement | null;
        if (!w) return false;
        const cs = getComputedStyle(w);
        if (cs.display === 'none' || cs.visibility === 'hidden') return false;
        const r = w.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      });
      expect(welcomeStillActive, 'Welcome view should have yielded to the project view after load').toBe(false);
    });
    // Scenario 2: in-app navigation control (warm session). URL goto is equivalent to a Recent Projects click.
    await softStep('Scenario 2 Step 1+2: from inside Datagrok, navigate to the same project via direct-link URL (warm session)', async () => {
      await page.evaluate(async () => {
        const grok = (window as any).grok;
        try { grok.shell.closeAll(); } catch (_) {}
        await new Promise((r) => setTimeout(r, 500));
      });
      await page.goto(baseUrl + directLinkPath);
    });
    await softStep('Scenario 2 Step 3: observe loading window during in-app open (control case — no cropping)', async () => {
      // Same invariant against the warm-session page (control: pre-fix this path never reproduced the glitch).
      const snapshots: Array<{t: number; box: {found: boolean; w: number; h: number}}> = [];
      for (let i = 0; i < 30; i++) {
        const box = await getBoxOf(page, '.power-pack-welcome-view');
        snapshots.push({t: i * 200, box});
        const gridMounted = await page.evaluate(() =>
          !!document.querySelector('[name="viewer-Grid"] canvas'));
        if (gridMounted) break;
        await page.waitForTimeout(200);
      }
      const degenerate = snapshots.filter((s) =>
        s.box.found &&
        ((s.box.w === 0 && s.box.h > 0) || (s.box.h === 0 && s.box.w > 0)));
      expect(degenerate, `Scenario 2 control: in-app navigation must not produce cropped loading window. Offending snapshots: ${JSON.stringify(degenerate)}`).toEqual([]);
    });
    await softStep('Scenario 2 Step 4: wait for in-app open to complete', async () => {
      await waitForPreloaderGone(page);
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.waitForTimeout(1500);
    });
    await softStep('Scenario 2 Step 5: compare against Scenario 1 outcome (same visual quality)', async () => {
      const gridBox = await getBoxOf(page, '[name="viewer-Grid"]');
      expect(gridBox.found).toBe(true);
      expect(gridBox.w).toBeGreaterThan(100);
      expect(gridBox.h).toBeGreaterThan(100);
      const tableMeta = await page.evaluate(() => {
        const grok = (window as any).grok;
        const df = grok.shell.tv?.dataFrame;
        return df ? {name: df.name, rowCount: df.rowCount, colCount: df.columns.length} : null;
      });
      expect(tableMeta).not.toBeNull();
      expect(tableMeta!.rowCount).toBeGreaterThan(0);
      expect(tableMeta!.colCount).toBeGreaterThan(0);
    });
  } finally {
    // Cleanup: delete the project + table info, close both contexts.
    try {
      if (projectId || tableInfoId) {
        await page.evaluate(async (ids) => {
          const grok = (window as any).grok;
          if (ids.projectId) {
            try {
              const p = await grok.dapi.projects.find(ids.projectId);
              if (p) await grok.dapi.projects.delete(p);
            } catch (_) { /* best effort */ }
          }
          if (ids.tableInfoId) {
            try {
              const ti = await grok.dapi.tables.find(ids.tableInfoId);
              if (ti) await grok.dapi.tables.delete(ti);
            } catch (_) { /* best effort */ }
          }
        }, {projectId: projectId ?? undefined, tableInfoId: tableInfoId ?? undefined});
      }
    } catch (_) { /* best-effort cleanup, do not mask test outcome */ }
    try { await page.evaluate(() => (window as any).grok?.shell?.closeAll?.()); } catch (_) {}
    if (secondaryContext) {
      try { await secondaryContext.close(); } catch (_) {}
    }
  }
  finishSpec();
});
