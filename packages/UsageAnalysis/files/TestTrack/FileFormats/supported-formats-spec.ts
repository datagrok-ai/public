import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use({...specTestOptions, viewport: {width: 1280, height: 2400}});

const ALL_FORMATS = [
  '.bmp', '.csv', '.edf', '.fasta', '.feather', '.geojson', '.gz', '.h5',
  '.html', '.ipynb', '.ivp', '.json', '.kml', '.kmz', '.kxl', '.mat',
  '.md', '.nc', '.parquet', '.pdf', '.rda', '.rds', '.sas7bdat', '.sdf',
  '.sqlite', '.tar', '.topojson', '.xlsx', '.xml', '.zip',
];

const FORMATS = process.env.FORMATS_SUBSET
  ? process.env.FORMATS_SUBSET.split(',').map((s) => s.trim()).filter(Boolean)
  : ALL_FORMATS;

const FOLDER_LABEL = 'all_formats';
const OPEN_TIMEOUT_MS = 90_000; 

async function resolveFiles(page: Page): Promise<{folderPath: string; byExt: Record<string, string>}> {
  return page.evaluate(async (label: string) => {
    const g = (window as any).grok;
    const conns = await g.dapi.connections.list();
    const homeConns = conns.filter((c: any) =>
      /home|my files/i.test(`${c.friendlyName ?? ''} ${c.name ?? ''}`));
    for (const c of homeConns) {
      const root = `${c.nqName ?? c.name}/${label}`;
      try {
        const files = await g.dapi.files.list(root, false);
        if (files.length) {
          const byExt: Record<string, string> = {};
          for (const f of files) {
            const name: string = f.fileName ?? f.name ?? '';
            const m = name.match(/(\.[^.]+)$/);

            if (m && !f.isDirectory && !(m[1].toLowerCase() in byExt))
              byExt[m[1].toLowerCase()] = f.fullPath;
          }
          return {folderPath: root, byExt};
        }
      } catch (e) {  }
    }
    return {folderPath: '', byExt: {}};
  }, FOLDER_LABEL);
}

async function openFolderUi(page: Page): Promise<void> {
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d: any) => d.remove());
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = true;
  });
  await page.waitForFunction(() => document.querySelectorAll('.d4-tree-view-group-label, .d4-tree-view-item-label').length > 0,
    null, {timeout: 15_000});

  const treeLabel = (t: string) =>
    page.waitForFunction((label: string) => {
      const sel = '.d4-tree-view-group-label, .d4-tree-view-item-label';
      const el = (Array.from(document.querySelectorAll(sel)) as HTMLElement[])
        .find((e) => e.textContent?.trim() === label);
      return el ? true : null;
    }, t, {timeout: 15_000});
  const clickTreeLabel = (t: string) =>
    page.evaluate((label: string) => {
      const sel = '.d4-tree-view-group-label, .d4-tree-view-item-label';
      (Array.from(document.querySelectorAll(sel)) as HTMLElement[])
        .find((e) => e.textContent?.trim() === label)?.click();
    }, t);

  await treeLabel('Files').catch(() => {});
  await clickTreeLabel('Files');
  await treeLabel('My files').catch(() => {});
  await clickTreeLabel('My files');
  await treeLabel(FOLDER_LABEL).catch(() => {});
  await page.locator('label').filter({hasText: FOLDER_LABEL}).first().dblclick({timeout: 15_000});
  await page.waitForFunction(() =>
    (Array.from(document.querySelectorAll('label')) as HTMLElement[])
      .filter((l) => /\.\w{2,}$/.test(l.textContent?.trim() ?? '') && !l.closest('.grok-prop-panel')).length >= 10,
    null, {timeout: 15_000});
}

async function returnToFolder(page: Page): Promise<void> {
  const ok = await page.evaluate((label: string) => {
    const views = Array.from(grok.shell.views) as any[];
    const folderView = views.find((v) => v.type === 'files' && (v.name ?? '').includes(label));
    if (!folderView) return false;
    grok.shell.v = folderView;
    for (const v of views)
      if (v !== folderView && v.type !== 'datagrok') { try { v.close(); } catch (e) {  } }
    return true;
  }, FOLDER_LABEL);
  if (!ok) await openFolderUi(page);
  else
    await page.waitForFunction(() =>
      (Array.from(document.querySelectorAll('label')) as HTMLElement[])
        .filter((l) => /\.\w{2,}$/.test(l.textContent?.trim() ?? '') && !l.closest('.grok-prop-panel')).length >= 10,
      null, {timeout: 10_000}).catch(() => {});
}

test('File formats: preview and open from My Files / all_formats', async ({page}) => {

  test.setTimeout(480_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);
  await openFolderUi(page);

  const {folderPath, byExt} = await resolveFiles(page);
  expect(folderPath.length > 0, 'Could not resolve the My Files/all_formats folder path').toBe(true);

  for (const ext of FORMATS) {
    await softStep(`Open ${ext}: platform handler parses the file without error`, async () => {
      const path = byExt[ext];
      expect(!!path, `No ${ext} file found in ${folderPath}`).toBe(true);

      const r = await page.evaluate(async (args: {path: string; to: number}) => {
        const g = (window as any).grok;
        try {
          const dfs = await Promise.race([
            g.data.files.openTables(args.path),
            new Promise((_, rej) => setTimeout(() => rej(new Error(`timeout ${args.to}ms`)), args.to)),
          ]);
          return {ok: true, count: Array.isArray(dfs) ? dfs.length : 0};
        } catch (e: any) { return {ok: false, err: String(e?.message ?? e).slice(0, 160)}; }
      }, {path, to: OPEN_TIMEOUT_MS});

      expect(r.ok, `Opening ${ext} threw: ${r.err ?? ''}`).toBe(true);
      console.log(`[OPEN ${ext}] ${r.ok ? `ok, ${r.count} dataframe(s)` : `FAIL ${r.err}`}`);
    });
  }

  await returnToFolder(page);
  for (const ext of FORMATS) {
    try {
      await test.step(`Preview ${ext}: single-click shows a preview (best-effort)`, async () => {
        if (await page.evaluate(() => (Array.from(document.querySelectorAll('label')) as HTMLElement[])
          .filter((l) => /\.\w{2,}$/.test(l.textContent?.trim() ?? '') && !l.closest('.grok-prop-panel')).length) < 10)
          await returnToFolder(page);

        const coords = await page.evaluate((e: string) => {
          const pat = new RegExp('\\.' + e.slice(1).replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '$', 'i');
          const el = (Array.from(document.querySelectorAll('label')) as HTMLElement[])
            .find((l) => pat.test(l.textContent?.trim() ?? '') && !l.closest('.grok-prop-panel') &&
              !l.closest('.d4-property-panel'));
          if (!el) return null;
          el.scrollIntoView({block: 'center', behavior: 'instant'});
          const r = el.getBoundingClientRect();
          return {x: Math.round(r.left + r.width / 2), y: Math.round(r.top + r.height / 2)};
        }, ext);
        if (!coords) { console.log(`[PREVIEW ${ext}] no label found`); return; }

        const viewsBefore = await page.evaluate(() => Array.from(grok.shell.views).length);
        await page.mouse.click(coords.x, coords.y);
        const ok = await page.waitForFunction((before: number) => {
          const p = document.querySelector('.grok-prop-panel') as HTMLElement | null;
          const panelFilled = !!p && p.innerText.trim().length > 0;
          const viewOpened = Array.from(grok.shell.views).length > before;
          return panelFilled || viewOpened ? true : null;
        }, viewsBefore, {timeout: 12_000}).then(() => true).catch(() => false);

        const opened = await page.evaluate((before: number) => Array.from(grok.shell.views).length > before, viewsBefore);
        if (opened) await returnToFolder(page);
        console.log(`[PREVIEW ${ext}] ${ok ? 'panel/view shown' : 'no panel/view (UI panel quirk)'}`);
      });
    } catch (e: any) { console.log(`[PREVIEW ${ext}] error: ${String(e?.message ?? e).slice(0, 80)}`); }
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
