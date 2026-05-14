import { test, expect } from '@playwright/test';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import { canvasHash } from './helpers/canvas-hash';
import { readIvpTooltips } from './helpers/ivp-tooltips';
import {
  openDiffStudio, openModelFromLibrary, listTabs, clickTab,
  setInputValue, readInputTooltip, inputEditor, inputHost,
} from './helpers/diff-studio';

/**
 * Test Track scenario: DiffStudio/stages.md
 * 1. Open Diff Studio + Acid Production from Library.
 *    NOTE: the library menu item is titled 'Acid production' (TITLE.ACID).
 * 2. Multiaxis and Facet tabs are updated.
 * 3. Modify an input value; the solution updates in real time and the URL reflects the new value.
 * 4. Hover inputs — tooltips appear with relevant info.
 */
test('DiffStudio Stages — Acid Production: Multiaxis/Facet tabs, input modify live update, tooltips',
  async ({ page }) => {
    test.setTimeout(300_000);
    const { softStep, assertAllPassed } = createSoftStepCollector();
    const monitor = attachErrorMonitor(page);

    await softStep('Step 1: Open Diff Studio + Acid Production from Library', async () => {
      await openDiffStudio(page);
      await openModelFromLibrary(page, 'Acid production');
      const inputs = await page.locator('[name^="input-host-"]').count();
      expect(inputs).toBeGreaterThan(5);
    });

    await softStep('Step 2: Multiaxis and Facet tabs are visible', async () => {
      const tabs = await listTabs(page);
      expect(tabs).toEqual(expect.arrayContaining(['Multiaxis', 'Facet']));
    });

    await softStep('Step 3: Modify "1-st stage" input; value, URL AND chart redraw live', async () => {
      const ed = page.locator(inputEditor('1-st-stage'));
      await ed.waitFor({ timeout: 10_000 });
      const before = await ed.inputValue();
      await clickTab(page, 'Multiaxis');
      await page.waitForTimeout(800);
      const chartBefore = await canvasHash(page, '.d4-viewer');
      await setInputValue(page, '1-st-stage', '50');
      await page.waitForTimeout(2000);
      const after = await ed.inputValue();
      expect(after).toBe('50');
      expect(after).not.toBe(before);
      expect(page.url()).toContain('1-ststage=50');
      const chartAfter = await canvasHash(page, '.d4-viewer');
      expect(chartAfter).not.toBe(chartBefore);
      expect(chartAfter).not.toBe('<missing>');
    });

    await softStep('Step 4: Tooltips on 1-st stage, biomass, glucose MATCH expected text', async () => {
      // Try parsing the IVP source for the canonical tooltip text. If the file isn't
      // reachable on this environment, fall back to the known expected strings.
      let tooltips = await readIvpTooltips(page, 'System:AppData/DiffStudio/library/ga-production.ivp');
      if (tooltips.size === 0) {
        // Expected values per the IVP file checked into the repo
        // (public/packages/DiffStudio/files/library/ga-production.ivp).
        tooltips = new Map([
          ['1-st stage', 'Duration of the 1-st stage'],
          ['biomass', 'Aspergillus niger biomass'],
          ['glucose', 'Glucose'],
        ]);
      }
      const e1 = tooltips.get('1-st stage') ?? '';
      const e2 = tooltips.get('biomass') ?? '';
      const e3 = tooltips.get('glucose') ?? '';
      expect(e1.length).toBeGreaterThan(0);
      expect(e2.length).toBeGreaterThan(0);
      expect(e3.length).toBeGreaterThan(0);

      const t1 = await readInputTooltip(page, '1-st-stage');
      const t2 = await readInputTooltip(page, 'biomass');
      const t3 = await readInputTooltip(page, 'glucose');
      expect(t1).toContain(e1);
      expect(t2).toContain(e2);
      expect(t3).toContain(e3);
    });

    assertAllPassed();
    monitor.assertNone();
  });
