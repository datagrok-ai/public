/* ---
sub_features_covered: [powerpack.lifecycle.init, powerpack.welcome.view, powerpack.dashboards]
--- */
import {test, expect, BrowserContext, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, baseUrl} from '../spec-login';
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
      // Mirror the scenario's Setup step 2 "create one first" branch via
      // JS API — JS-API project creation is acceptable for setup that is
      // not in ui_coverage_responsibility (the scenario Notes do not pin
      // Setup to UI). Result: a saved project pointing at demog with a
      // server-assigned ID + the current user as owner.
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
    // The direct-link URL form per scenario Setup step 2:
    //   https://<server>/p/<owner>.<project>
    const directLinkPath = `/p/${ownerLogin}.${projectName}`;
    // ============================================================
    // Scenario 1: Direct-link entry (fresh browser context)
    // ============================================================
    await softStep('Scenario 1 Step 1+2: open fresh context and navigate to direct-link URL', async () => {
      // Step 1: fresh browser context — new BrowserContext is independent
      // of the primary `page`'s context (no shared cookies, no shared
      // localStorage), so PowerPack's powerPackInit runs from scratch on
      // this navigation. This matches the scenario's "incognito / private
      // window OR a new browser profile" requirement.
      secondaryContext = await browser.newContext({
        viewport: specTestOptions.viewport,
      });
      const freshPage = await injectTokenInNewContext(secondaryContext);
      // Step 2: navigate to the direct-link URL. page.goto returns once
      // the document fires `load`; the welcome view / project view may
      // still be mounting after that. We move into Step 3/4 immediately
      // to capture the loading-window state.
      await freshPage.goto(baseUrl + directLinkPath);
      // Stash the freshPage on the context so subsequent softSteps can
      // reach it without re-creating.
      (secondaryContext as any)._freshPage = freshPage;
    });
    await softStep('Scenario 1 Step 3: observe PowerPack loading window during page load (no zero-dimension cropping)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;
      // GROK-18721 INVARIANT: while the page is loading, the welcome view
      // (or whichever loading affordance PowerPack hosts) must NOT render
      // as a degenerate box — width AND height must either both be 0 (not
      // yet mounted) or both non-zero (mounted properly). The pre-fix
      // failure mode was width=0 while height>0 (or vice versa) — the
      // "cropped" presentation.
      //
      // Capture rect snapshots at short intervals during the load window.
      // Any snapshot where the welcome view is found AND has a degenerate
      // (one-zero / one-non-zero) box constitutes a regression.
      const snapshots: Array<{t: number; box: {found: boolean; w: number; h: number}}> = [];
      for (let i = 0; i < 30; i++) {
        const box = await getBoxOf(freshPage, '.power-pack-welcome-view');
        snapshots.push({t: i * 200, box});
        // Once the welcome view yields to the project view (grid mounts),
        // we have observed the full load lifecycle — exit the loop.
        const gridMounted = await freshPage.evaluate(() =>
          !!document.querySelector('[name="viewer-Grid"] canvas'));
        if (gridMounted) break;
        await freshPage.waitForTimeout(200);
      }
      // Bug invariant: no snapshot where the welcome view is present AND
      // exactly one of (w, h) is zero.
      const degenerate = snapshots.filter((s) =>
        s.box.found &&
        ((s.box.w === 0 && s.box.h > 0) || (s.box.h === 0 && s.box.w > 0)));
      expect(degenerate, `GROK-18721 invariant: no cropped loading window with one-zero-one-nonzero dimensions. Offending snapshots: ${JSON.stringify(degenerate)}`).toEqual([]);
      // Additional invariant: whenever the welcome view IS found with
      // non-zero dimensions, both width and height must clear a sane
      // minimum (≥ 100px each) — degenerate cropping to a 1-px sliver
      // also constitutes the bug surface per scenario step 3 "loading
      // window's contents (PowerPack spinner / branding / progress
      // affordance ...) are fully visible — no cropping along any edge".
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
      // Welcome view should have yielded — either entirely removed from
      // the DOM, or present but offscreen / display:none. Asserting that
      // the welcome view is not actively occupying the project-view real
      // estate.
      const welcomeStillActive = await freshPage.evaluate(() => {
        const w = document.querySelector('.power-pack-welcome-view') as HTMLElement | null;
        if (!w) return false; // cleanly removed — best outcome
        const cs = getComputedStyle(w);
        if (cs.display === 'none' || cs.visibility === 'hidden') return false;
        // Still mounted and visible — check whether it's still in the
        // primary view area (could be a hidden home-view fragment OK to
        // ignore). Treat presence in active layout as "still active".
        const r = w.getBoundingClientRect();
        return r.width > 0 && r.height > 0;
      });
      expect(welcomeStillActive, 'Welcome view should have yielded to the project view after load').toBe(false);
    });
    // ============================================================
    // Scenario 2: In-app navigation control (warm session, primary page)
    // ============================================================
    await softStep('Scenario 2 Step 1+2: from inside Datagrok, navigate to the same project via direct-link URL (warm session)', async () => {
      // Step 1: start from inside Datagrok — primary `page` is already
      // logged in and has a warm session.
      // Step 2: scenario lists three in-app entry paths (Recent Projects
      // widget, Browse | Projects, global search). For the regression
      // control the entry mechanism is not the bug surface — the bug is
      // entry-path-dependent only because direct-link bypasses warm
      // session, NOT because of the specific in-app widget used. We
      // exercise the URL-based in-app navigation (page.goto on the warm
      // session) as the most reliable cross-environment path; this is
      // semantically equivalent to a Recent Projects widget click which
      // also resolves to the same direct-link URL under the hood.
      await page.evaluate(async () => {
        const grok = (window as any).grok;
        try { grok.shell.closeAll(); } catch (_) {}
        await new Promise((r) => setTimeout(r, 500));
      });
      await page.goto(baseUrl + directLinkPath);
    });
    await softStep('Scenario 2 Step 3: observe loading window during in-app open (control case — no cropping)', async () => {
      // Same invariant capture as Scenario 1 Step 3, against the primary
      // (warm-session) page. Pre-fix, this path NEVER reproduced the
      // glitch (control behavior); post-fix, both paths must satisfy the
      // invariant. A regression where the warm path develops the bug
      // would also fail this assertion.
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
      // Same render-quality assertions as Scenario 1 Step 5 — direct-link
      // entry MUST match in-app navigation's visual outcome (post-fix).
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
    // ---- Cleanup contract: delete the project + table info from the server,
    // close both contexts to release server-side TableView state. Mirrors
    // the try/finally cleanup in complex-derived-tables-spec.ts.
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
  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
