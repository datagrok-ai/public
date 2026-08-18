import {test, expect, Page} from '@playwright/test';
import {specTestOptions, softStep, stepErrors} from '../spec-login';

test.use({...specTestOptions, viewport: {width: 1280, height: 2400}});

const BASE = process.env.DATAGROK_URL ?? 'http://localhost:8888';

const DWELL_MS = process.env.PREVIEW_DWELL_MS ? Number(process.env.PREVIEW_DWELL_MS) : 7000;

async function expandToAllFormats(page: Page): Promise<string[]> {
  return page.evaluate(async () => {
    const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
    const findChildGroupByLabel = (parent: ParentNode, label: string): Element | null => {
      for (const c of Array.from(parent.querySelectorAll('.d4-tree-view-group-label'))
        .filter((el) => el.textContent?.trim() === label)) {
        const g = c.closest('.d4-tree-view-group');
        if (g && (parent === document || (parent as Element).contains(g))) return g;
      }
      return null;
    };
    const expandGroup = (group: Element) => {
      const tri = group.querySelector(':scope > .d4-tree-view-node > .d4-tree-view-tri') as HTMLElement | null;
      if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
    };
    let scope: ParentNode = document;
    for (const [label, next] of [['Files', 'My files'], ['My files', 'all_formats']] as [string, string][]) {
      const grp = findChildGroupByLabel(scope, label);
      if (!grp) throw new Error(`Tree node not found: ${label}`);
      expandGroup(grp);
      const dl = Date.now() + 20_000;
      while (Date.now() < dl && !findChildGroupByLabel(grp, next)) await wait(200);
      if (!findChildGroupByLabel(grp, next)) throw new Error(`Expand of ${label} did not reveal ${next}`);
      scope = grp;
    }
    const allFormats = findChildGroupByLabel(scope, 'all_formats');
    if (!allFormats) throw new Error('all_formats not found under My files');
    expandGroup(allFormats);

    const collect = () => (Array.from(allFormats.querySelectorAll('.d4-tree-view-item-label')) as HTMLElement[])
      .map((e) => e.textContent?.trim() ?? '').filter((t) => /\.\w{2,}$/.test(t));
    const dl = Date.now() + 20_000;
    let prev = -1; let stableSince = Date.now();
    while (Date.now() < dl) {
      await wait(300);
      const n = collect().length;
      if (n !== prev) { prev = n; stableSince = Date.now(); }
      else if (n > 0 && Date.now() - stableSince > 1000) break;
    }
    return collect();
  });
}

async function openBrowseDirect(page: Page): Promise<void> {
  const token = process.env.DATAGROK_AUTH_TOKEN;
  if (!token || token.length === 0)
    throw new Error('DATAGROK_AUTH_TOKEN is not set. Run via `grok test` or set it from the dev key.');
  const u = new URL(BASE);
  await page.context().addCookies([{name: 'auth', value: token, domain: u.hostname, path: '/'}]);
  await page.addInitScript((t) => { try { window.localStorage.setItem('auth', t); } catch (e) {  } }, token);
  await page.goto(`${BASE}/browse`, {waitUntil: 'domcontentloaded'});
  await page.waitForFunction(() => !document.querySelector('.grok-preloader'), null, {timeout: 120_000});
  await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
}

async function clickFileNode(page: Page, name: string): Promise<boolean> {
  for (let i = 0; i < 50; i++) {
    const r = await page.evaluate((t: string) => {
      const items = Array.from(document.querySelectorAll('.d4-tree-view-item-label')) as HTMLElement[];
      const el = items.find((e) => e.textContent?.trim() === t);
      if (el) { el.scrollIntoView({block: 'center', behavior: 'instant'}); el.click(); return 'clicked'; }
      items[items.length - 1]?.scrollIntoView({block: 'end', behavior: 'instant'}); 
      return 'scrolled';
    }, name);
    if (r === 'clicked') return true;
    await page.waitForTimeout(200);
  }
  return false;
}

test('File formats: preview every file by clicking through the Browse tree', async ({page}) => {

  test.setTimeout(360_000);
  stepErrors.length = 0;

  await openBrowseDirect(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
    grok.shell.windows.showBrowse = true;
  });
  await page.waitForFunction(() => document.querySelectorAll('.d4-tree-view-group-label').length > 5,
    undefined, {timeout: 30_000});

  await page.waitForTimeout(2000);

  let files = await expandToAllFormats(page);
  expect(files.length > 0, 'No files found under My files / all_formats in the Browse tree').toBe(true);

  if (process.env.FILES_SUBSET) {
    const want = process.env.FILES_SUBSET.split(',').map((s) => s.trim()).filter(Boolean);
    files = files.filter((f) => want.some((w) => f.includes(w)));
  }
  console.log(`Found ${files.length} files in the tree: ${files.join(', ')}`);

  const closeToHome = () => page.evaluate(() => {
    for (const v of Array.from((window as any).grok.shell.views) as any[])
      if (v?.type !== 'datagrok' && (v?.name ?? '') !== 'Home') { try { v.close(); } catch (e) {  } }
  });

  for (const name of files) {
    await softStep(`Preview ${name}: open its preview from the tree`, async () => {
      await closeToHome();

      await page.waitForFunction(() => Array.from((window as any).grok.shell.views)
        .every((v: any) => v?.type === 'datagrok' || (v?.name ?? '') === 'Home'),
      undefined, {timeout: 10_000}).catch(() => {});
      const clicked = await clickFileNode(page, name);
      expect(clicked, `Tree node not found for ${name}`).toBe(true);

      const previewed = await page.waitForFunction(() => {
        const g = (window as any).grok;
        const newView = Array.from(g.shell.views).some((v: any) => v?.type !== 'datagrok' && (v?.name ?? '') !== 'Home');
        const curNotHome = !!g.shell.v && (g.shell.v.name ?? '') !== 'Home';
        const inSidebar = (el: Element) => !!el.closest('.d4-tree-view, .grok-sidebar, [name="Browse"]');
        const bigMedia = (Array.from(document.querySelectorAll('canvas, img, iframe, .d4-grid')) as HTMLElement[])
          .some((e) => !inSidebar(e) && e.clientWidth > 150 && e.clientHeight > 80);
        return newView || curNotHome || bigMedia ? true : null;
      }, undefined, {timeout: 20_000}).then(() => true).catch(() => false);

      await page.waitForTimeout(DWELL_MS);

      const diag = await page.evaluate(() => {
        const g = (window as any).grok;
        const err = !!document.querySelector('.d4-balloon-error') ||
          Array.from(document.querySelectorAll('.d4-balloon, .d4-error')).some((e) =>
            /error|failed|cannot|unable/i.test((e as HTMLElement).textContent ?? ''));
        return {curView: `${g.shell.v?.name ?? ''}|${g.shell.v?.type ?? ''}`, err};
      });

      console.log(`[PREVIEW ${name}] previewOpened=${previewed} curView=${diag.curView} error=${diag.err}`);
      expect(previewed, `No preview opened for ${name}`).toBe(true);
      expect(diag.err, `Error balloon appeared while previewing ${name}`).toBe(false);
    });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} file(s) failed to preview:\n${summary}`);
  }
});
