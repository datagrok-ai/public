import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import { canvasHash } from './helpers/canvas-hash';
import {
  openDiffStudio, openModelFromLibrary, listTabs, clickTab,
  setInputValue, selectChoice, inputEditor, inputHost,
} from './helpers/diff-studio';

test('DiffStudio Open Model — Bioreactor: Multiaxis/Facet tabs, Switch at, Process mode cascade', async ({ page }) => {
  test.setTimeout(300_000);
  const { softStep, assertAllPassed } = createSoftStepCollector();
  const monitor = attachErrorMonitor(page);

  await softStep('Step 1: Open Diff Studio (Apps > Diff Studio)', async () => {
    await openDiffStudio(page);
    await expect(page.locator('.d4-ribbon').first()).toBeVisible();
  });

  await softStep('Step 2: Load Bioreactor from Library', async () => {
    await openModelFromLibrary(page, 'Bioreactor');
    await expect(page.locator(inputHost('Process-mode'))).toBeVisible({ timeout: 15_000 });
  });

  await softStep('Step 3: Multiaxis and Facet tabs are visible under the line chart', async () => {
    const tabs = await listTabs(page);
    expect(tabs).toEqual(expect.arrayContaining(['Multiaxis', 'Facet']));
    await clickTab(page, 'Multiaxis');
    await page.waitForTimeout(800);
    await clickTab(page, 'Facet');
    await page.waitForTimeout(800);

    const canvasCount = await page.locator('.d4-viewer canvas').count();
    expect(canvasCount).toBeGreaterThan(2);
  });

  await softStep('Step 4: Facet plot curves are not the same color (pixel histogram)', async () => {

    await clickTab(page, 'Facet');
    await page.waitForTimeout(1500);
    const distinctColours = await page.evaluate(() => {
      const canvases = Array.from(document.querySelectorAll('.d4-viewer canvas')) as HTMLCanvasElement[];
      const buckets = new Set<string>();
      for (const c of canvases) {
        const ctx = c.getContext('2d');
        if (!ctx) continue;
        let img: ImageData;
        try { img = ctx.getImageData(0, 0, c.width, c.height); }
        catch { continue; }
        const d = img.data;
        for (let i = 0; i < d.length; i += 4) {
          const r = d[i], g = d[i + 1], b = d[i + 2], a = d[i + 3];
          if (a < 128) continue;
          if (r > 240 && g > 240 && b > 240) continue;   
          if (r < 30 && g < 30 && b < 30) continue;       
          if (Math.abs(r - g) < 15 && Math.abs(g - b) < 15) continue;  
          buckets.add(`${r >> 4}|${g >> 4}|${b >> 4}`);
        }
      }
      return buckets.size;
    });

    expect(distinctColours).toBeGreaterThanOrEqual(10);
  });

  await softStep('Step 5: Adjust Switch at input; URL + chart redraw live', async () => {
    const ed = page.locator(inputEditor('switch-at'));
    await ed.waitFor({ timeout: 10_000 });

    await clickTab(page, 'Multiaxis');
    await page.waitForTimeout(800);
    const chartBefore = await canvasHash(page, '.d4-viewer');
    await setInputValue(page, 'switch-at', '150');
    await expect(ed).toHaveValue('150');
    expect(page.url()).toContain('switchat=150');

    await page.waitForTimeout(2500);
    const chartAfter = await canvasHash(page, '.d4-viewer');
    expect(chartAfter).not.toBe(chartBefore);
    expect(chartAfter).not.toBe('<missing>');
  });

  await softStep('Step 6: Modify Process mode; FFox/KKox/others cascade + charts redraw', async () => {
    const readAll = async () => ({
      ffox: await page.locator(inputEditor('FFox')).inputValue().catch(() => ''),
      kkox: await page.locator(inputEditor('KKox')).inputValue().catch(() => ''),
      ffred: await page.locator(inputEditor('FFred')).inputValue().catch(() => ''),
      mea: await page.locator(inputEditor('MEAthiol')).inputValue().catch(() => ''),
      temp: await page.locator(inputEditor('temperature')).inputValue().catch(() => ''),
      gas: await page.locator(inputEditor('Gas')).inputValue().catch(() => ''),
    });
    const before = await readAll();
    await clickTab(page, 'Multiaxis');
    await page.waitForTimeout(800);
    const chartBefore = await canvasHash(page, '.d4-viewer');
    await selectChoice(page, 'Process-mode', 'Mode 1');
    await page.waitForTimeout(2500);
    const after = await readAll();
    const changed = Object.keys(before).filter(k => (before as any)[k] !== (after as any)[k]);
    expect(changed.length).toBeGreaterThanOrEqual(4);
    const chartAfter = await canvasHash(page, '.d4-viewer');
    expect(chartAfter).not.toBe(chartBefore);
  });

  assertAllPassed();
  monitor.assertNone();
});
