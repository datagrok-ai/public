import {expect, BrowserContext, Page} from '@playwright/test';
import {test} from '../shared-page';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, baseUrl} from '../spec-login';
import {finishSpec, openTable} from '../helpers/viewers';
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
interface BoxSample { t: number; box: {found: boolean; w: number; h: number}; }

/** Samples the loading window's box every 200ms until the grid canvas mounts, entirely
 *  in the page — the Node-side version paid two round trips per sample. */
async function sampleWelcomeBox(page: Page): Promise<BoxSample[]> {
  return await page.evaluate(async () => {
    const out: BoxSample[] = [];
    for (let i = 0; i < 30; i++) {
      const el = document.querySelector('.power-pack-welcome-view') as HTMLElement | null;
      const r = el?.getBoundingClientRect();
      out.push({t: i * 200, box: r ? {found: true, w: r.width, h: r.height} : {found: false, w: 0, h: 0}});
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((res) => setTimeout(res, 200));
    }
    return out;
  });
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
      await openTable(page, {path: 'System:DemoFiles/demog.csv'});
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
    const directLinkPath = `/p/${ownerLogin}.${projectName}`; 

    await softStep('Scenario 1 Step 1+2: open fresh context and navigate to direct-link URL', async () => {

      secondaryContext = await browser.newContext({
        viewport: specTestOptions.viewport,
      });
      const freshPage = await injectTokenInNewContext(secondaryContext);
      await freshPage.goto(baseUrl + directLinkPath);
      (secondaryContext as any)._freshPage = freshPage;
    });
    await softStep('Scenario 1 Step 3: observe PowerPack loading window during page load (no zero-dimension cropping)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;

      const snapshots = await sampleWelcomeBox(freshPage);
      const degenerate = snapshots.filter((s) =>
        s.box.found &&
        ((s.box.w === 0 && s.box.h > 0) || (s.box.h === 0 && s.box.w > 0)));
      expect(degenerate, `GROK-18721 invariant: no cropped loading window with one-zero-one-nonzero dimensions. Offending snapshots: ${JSON.stringify(degenerate)}`).toEqual([]);

      const tinyBoxes = snapshots.filter((s) =>
        s.box.found && s.box.w > 0 && s.box.h > 0 &&
        (s.box.w < 100 || s.box.h < 100));
      expect(tinyBoxes, `GROK-18721 invariant: welcome view must not render at sub-100px in either dimension. Offending snapshots: ${JSON.stringify(tinyBoxes)}`).toEqual([]);
    });
    await softStep('Scenario 1 Step 4: wait for the load to complete (preloader gone + grid mounted)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;
      await waitForPreloaderGone(freshPage);
      await freshPage.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await freshPage.locator('[name="viewer-Grid"] canvas').first()
        .waitFor({timeout: 1500}).catch(() => {});
    });
    await softStep('Scenario 1 Step 5: verify post-load rendering (grid has non-zero dimensions, no zombie welcome fragments)', async () => {
      const freshPage: Page = (secondaryContext as any)._freshPage;

      const gridBox = await getBoxOf(freshPage, '[name="viewer-Grid"]');
      expect(gridBox.found).toBe(true);
      expect(gridBox.w).toBeGreaterThan(100);
      expect(gridBox.h).toBeGreaterThan(100);

      const tableMeta = await freshPage.evaluate(() => {
        const grok = (window as any).grok;
        const df = grok.shell.tv?.dataFrame;
        return df ? {name: df.name, rowCount: df.rowCount, colCount: df.columns.length} : null;
      });
      expect(tableMeta).not.toBeNull();
      expect(tableMeta!.rowCount).toBeGreaterThan(0);
      expect(tableMeta!.colCount).toBeGreaterThan(0);

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

    await softStep('Scenario 2 Step 1+2: from inside Datagrok, navigate to the same project via direct-link URL (warm session)', async () => {
      await page.evaluate(async () => {
        const grok = (window as any).grok;
        try { grok.shell.closeAll(); } catch (_) {}
        await new Promise((r) => setTimeout(r, 500));
      });
      await page.goto(baseUrl + directLinkPath);
    });
    await softStep('Scenario 2 Step 3: observe loading window during in-app open (control case — no cropping)', async () => {

      const snapshots = await sampleWelcomeBox(page);
      const degenerate = snapshots.filter((s) =>
        s.box.found &&
        ((s.box.w === 0 && s.box.h > 0) || (s.box.h === 0 && s.box.w > 0)));
      expect(degenerate, `Scenario 2 control: in-app navigation must not produce cropped loading window. Offending snapshots: ${JSON.stringify(degenerate)}`).toEqual([]);
    });
    await softStep('Scenario 2 Step 4: wait for in-app open to complete', async () => {
      await waitForPreloaderGone(page);
      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
      await page.locator('[name="viewer-Grid"] canvas').first()
        .waitFor({timeout: 1500}).catch(() => {});
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

    await page.evaluate((ids) => {
      const w = window as any;
      const drop = async (dapi: any, id: string | null | undefined) => {
        if (!id) return;
        try { const e = await dapi.find(id); if (e) await dapi.delete(e); } catch (_) {  }
      };
      w.__pendingDeletes = w.__pendingDeletes ?? [];
      w.__pendingDeletes.push(Promise.all([
        drop(w.grok.dapi.projects, ids.projectId), drop(w.grok.dapi.tables, ids.tableInfoId),
      ]));
    }, {projectId, tableInfoId}).catch(() => {});
    try { await page.evaluate(() => (window as any).grok?.shell?.closeAll?.()); } catch (_) {}
    if (secondaryContext) {
      try { await secondaryContext.close(); } catch (_) {}
    }
  }
  finishSpec();
});
