import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import { canvasHash } from './helpers/canvas-hash';
import { readIvpTooltips } from './helpers/ivp-tooltips';
import {
  openDiffStudio, openModelFromLibrary, listTabs, clickTab, clickerIncrement,
  readInputTooltip, inputEditor, inputHost,
} from './helpers/diff-studio';

test('DiffStudio Cyclic Models — PK-PD: Multiaxis/Facet tabs, Count clickers live update, input tooltips',
  async ({ page }) => {
    test.setTimeout(300_000);
    const { softStep, assertAllPassed } = createSoftStepCollector();
    const monitor = attachErrorMonitor(page);

    await softStep('Step 1: Open Diff Studio + PK-PD from Library', async () => {
      await openDiffStudio(page);
      await openModelFromLibrary(page, 'PK-PD');
      await expect(page.locator(inputHost('count'))).toBeVisible({ timeout: 15_000 });
    });

    await softStep('Step 2: Multiaxis and Facet tabs are visible', async () => {
      const tabs = await listTabs(page);
      expect(tabs).toEqual(expect.arrayContaining(['Multiaxis', 'Facet']));
      const canvases = await page.locator('.d4-viewer canvas').count();
      expect(canvases).toBeGreaterThan(2);
    });

    await softStep('Step 3: Modify Count via clickers; URL, value AND chart redraw live', async () => {
      const before = await page.locator(inputEditor('count')).inputValue();

      await clickTab(page, 'Multiaxis');
      await page.waitForTimeout(800);
      const chartBefore = await canvasHash(page, '.d4-viewer');
      await clickerIncrement(page, 'count', 3);
      await page.waitForTimeout(2000);
      const after = await page.locator(inputEditor('count')).inputValue();
      expect(parseInt(after)).toBeGreaterThan(parseInt(before));
      expect(page.url()).toContain('count=');
      const chartAfter = await canvasHash(page, '.d4-viewer');
      expect(chartAfter).not.toBe(chartBefore);
      expect(chartAfter).not.toBe('<missing>');
    });

    await softStep('Step 4: Tooltips on Begin, End, Step inputs MATCH expected text', async () => {

      let tooltips = await readIvpTooltips(page, 'System:AppData/DiffStudio/library/pk-pd.ivp');
      if (tooltips.size === 0) {
        tooltips = new Map([
          ['begin', 'Begin of dosing interval'],
          ['end', 'End of dosing interval'],
          ['step', 'Time step of simulation'],
        ]);
      }
      const expectBegin = tooltips.get('begin') ?? '';
      const expectEnd = tooltips.get('end') ?? '';
      const expectStep = tooltips.get('step') ?? '';
      expect(expectBegin.length).toBeGreaterThan(0);
      expect(expectEnd.length).toBeGreaterThan(0);
      expect(expectStep.length).toBeGreaterThan(0);

      const beginTip = await readInputTooltip(page, 'begin');
      const endTip = await readInputTooltip(page, 'end');
      const stepTip = await readInputTooltip(page, 'step');
      expect(beginTip).toContain(expectBegin);
      expect(endTip).toContain(expectEnd);
      expect(stepTip).toContain(expectStep);
    });

    assertAllPassed();
    monitor.assertNone();
  });
