import { test, expect } from '@playwright/test';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import { canvasHash } from './helpers/canvas-hash';
import {
  BASE, setInputValue, clickerIncrement, clickerDecrement, listTabs,
  inputEditor, inputHost,
} from './helpers/diff-studio';

/**
 * Test Track scenario: DiffStudio/files-and-sharing.md
 * 1. File preview: Browse > Files > App Data > Diff Studio > Library > pk.ivp.
 * 2. Modify inputs: Step = 0.1 via the input editor, Count = 4 via clickers.
 * 3. Sharing: open the resulting URL in a new tab — the model loads with the same inputs.
 * 4. REMARK: with 1–3 curves there is just a line chart (no Multiaxis & Facet); with ≥4 they appear.
 */
test('DiffStudio Files & Sharing — pk.ivp preview, modify Step/Count, URL share, curve count remark',
  async ({ page, context }) => {
    test.setTimeout(300_000);
    const { softStep, assertAllPassed } = createSoftStepCollector();
    const monitor = attachErrorMonitor(page);

    await softStep('Step 1: Open pk.ivp preview via Browse > Files > App Data > Diff Studio > library', async () => {
      // Navigate to the directory listing via address bar, then click the file row to trigger
      // the .ivp previewer (`DiffStudio:previewIvp`). Browsing to a directory URL is UI-level
      // (address bar); selecting a file via mouse click is UI-level (mouse).
      // The filesystem-side folder is `DiffStudio` (capital D) — DG does not lowercase it.
      await page.goto(`${BASE}/files/system.appdata/DiffStudio/library`);
      await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
      await page.waitForTimeout(2500);

      // The directory grid renders each file as a row containing the filename as text.
      // Find a clickable element whose text is exactly "pk.ivp" and click it.
      const clicked = await page.evaluate(() => {
        const all = Array.from(document.querySelectorAll('*')) as HTMLElement[];
        // Visible elements whose own text is exactly "pk.ivp" (no nested children)
        const candidates = all.filter(el =>
          el.children.length === 0 && el.textContent?.trim() === 'pk.ivp');
        if (candidates.length === 0) return false;
        const target = candidates[0];
        target.scrollIntoView({ block: 'center' });
        // Walk up to find the clickable row (file rows are usually `.grok-gallery-item` or list-item)
        let clickTarget: HTMLElement | null = target;
        for (let i = 0; clickTarget && i < 6; i++) {
          (clickTarget as any).click?.();
          clickTarget = clickTarget.parentElement;
        }
        return true;
      });
      expect(clicked).toBe(true);

      // Wait for the IVP previewer to render its input form (Initial/Final/Step/Count/...).
      await page.locator('[name^="input-host-"]').first().waitFor({ timeout: 30_000 });
      await page.waitForTimeout(2000);
      const inputs = await page.locator('[name^="input-host-"]').count();
      expect(inputs).toBeGreaterThan(5);
    });

    await softStep('Step 2: Step via slider drag, Count via clicker; URL + chart redraw live', async () => {
      // The first canvas in the pk.ivp preview is the Multiaxis line chart.
      const chartBefore = await canvasHash(page, '.d4-viewer');

      // Real slider drag — `<input type="range">` lives inside the Step input host alongside
      // the text editor. Per MD: "Set Step to 0.1 using the slider". pk.ivp Step has min=0.01,
      // max=0.1 — drag the thumb to the right-most position.
      const stepSlider = page.locator(`${inputHost('step')} input[type="range"]`).first();
      if (await stepSlider.count() > 0) {
        const box = await stepSlider.boundingBox();
        if (box) {
          await page.mouse.move(box.x + 4, box.y + box.height / 2);
          await page.mouse.down();
          // Move in small steps for realism — DG slider fires `input` events on each move.
          for (let i = 0; i <= 10; i++) {
            const x = box.x + (box.width - 8) * (i / 10) + 4;
            await page.mouse.move(x, box.y + box.height / 2);
          }
          await page.mouse.up();
          await page.waitForTimeout(800);
        }
      } else {
        // Slider not present in this build — fall back to typed value
        await setInputValue(page, 'step', '0.1');
      }

      // pk.ivp default Count is 1 — three clicks of "+" → 4
      await clickerIncrement(page, 'count', 3);
      await page.waitForTimeout(2000);

      const stepVal = await page.locator(inputEditor('step')).inputValue();
      const countVal = await page.locator(inputEditor('count')).inputValue();
      expect(parseFloat(stepVal)).toBeCloseTo(0.1, 1);  // 0.1 ± 0.1 — slider granularity
      expect(countVal).toBe('4');

      const url = page.url();
      expect(url).toMatch(/step=0\.1/);
      expect(url).toContain('count=4');

      const chartAfter = await canvasHash(page, '.d4-viewer');
      expect(chartAfter).not.toBe(chartBefore);
      expect(chartAfter).not.toBe('<missing>');
    });

    await softStep('Step 3: Open the URL in a new tab — model loads with the same inputs', async () => {
      const url = page.url();
      const newTab = await context.newPage();
      await newTab.goto(url);
      await newTab.waitForSelector('.d4-ribbon', { timeout: 60_000 });
      await newTab.locator(`${inputHost('step')} input.ui-input-editor`).waitFor({ timeout: 30_000 });
      await newTab.waitForTimeout(2000);
      const stepVal = await newTab.locator(`${inputHost('step')} input.ui-input-editor`).inputValue();
      const countVal = await newTab.locator(`${inputHost('count')} input.ui-input-editor`).inputValue();
      expect(parseFloat(stepVal)).toBeCloseTo(0.1);
      expect(countVal).toBe('4');
      await newTab.close();
    });

    await softStep('Step 4: REMARK — with 1–3 curves there is no Multiaxis/Facet', async () => {
      test.info().annotations.push({
        type: 'remark',
        description: 'REMARK from MD: with 1, 2 or 3 curves there is just the linechart — Multiaxis & Facet do not appear.',
      });
      // Already at count=4 → Multiaxis/Facet expected. Drop back to 1 via the "-" clicker
      // (the count input has `.ui-input-options` overlaying the editor, blocking typing).
      await clickerDecrement(page, 'count', 3);
      await page.waitForTimeout(1500);
      const tabsAt1 = await listTabs(page);
      expect(tabsAt1).not.toContain('Multiaxis');
      expect(tabsAt1).not.toContain('Facet');
    });

    assertAllPassed();
    monitor.assertNone();
  });
