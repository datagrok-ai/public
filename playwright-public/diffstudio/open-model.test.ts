import { test, expect } from '@playwright/test';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import { canvasHash } from './helpers/canvas-hash';
import {
  openDiffStudio, openModelFromLibrary, listTabs, clickTab,
  setInputValue, selectChoice, inputEditor, inputHost,
} from './helpers/diff-studio';

/**
 * Test Track scenario: DiffStudio/open-model.md
 * 1. Open Diff Studio.
 * 2. Load Bioreactor from Library.
 * 3. Check Multiaxis and Facet tabs (under linechart).
 * 4. Check that curves in Facet plot are not the same color.
 * 5. Adjust Switch at input; table and line chart update on the fly.
 * 6. Modify Process mode; FFox/KKox (and other) inputs cascade-update.
 */
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
    // At least the Multiaxis + Facet line-chart canvases mount
    const canvasCount = await page.locator('.d4-viewer canvas').count();
    expect(canvasCount).toBeGreaterThan(2);
  });

  await softStep('Step 4: Facet plot curves are not the same color (pixel histogram)', async () => {
    // Read pixels from every canvas under `.d4-viewer`, drop near-white/black/grey, bucket each
    // remaining channel to 4 bits to collapse antialiasing noise, and count distinct buckets.
    // A Facet plot drawing 12 differently coloured curves should produce at least ~10 distinct
    // buckets; a single-colour plot would produce only 1–2.
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
          if (r > 240 && g > 240 && b > 240) continue;   // near-white
          if (r < 30 && g < 30 && b < 30) continue;       // near-black
          if (Math.abs(r - g) < 15 && Math.abs(g - b) < 15) continue;  // greys
          buckets.add(`${r >> 4}|${g >> 4}|${b >> 4}`);
        }
      }
      return buckets.size;
    });
    // Threshold: 12 curves with antialiasing comfortably produce 50+ buckets in practice.
    // We assert >= 10 to give headroom against future palette changes.
    expect(distinctColours).toBeGreaterThanOrEqual(10);
  });

  await softStep('Step 5: Adjust Switch at input; URL + chart redraw live', async () => {
    const ed = page.locator(inputEditor('switch-at'));
    await ed.waitFor({ timeout: 10_000 });
    // Hash the line-chart canvas BEFORE the change, to prove the chart actually redrew.
    await clickTab(page, 'Multiaxis');
    await page.waitForTimeout(800);
    const chartBefore = await canvasHash(page, '.d4-viewer');
    await setInputValue(page, 'switch-at', '150');
    await expect(ed).toHaveValue('150');
    expect(page.url()).toContain('switchat=150');
    // The solver re-runs asynchronously; give it room to redraw, then hash again.
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
