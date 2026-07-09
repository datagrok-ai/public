import { test, expect } from './helpers/diff-studio';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  BASE, openDiffStudio, openModelFromLibrary, toggleRibbonSwitch,
  ribbonSwitchOn, setInputValue, inputEditor, inputHost,
} from './helpers/diff-studio';

/**
 * Test Track scenario: DiffStudio/scripting.md
 * 1. Open Diff Studio. Turn on Edit toggle, click </> to open JS script view.
 * 2. Run the script and move the slider for Final at — table and line chart update live.
 * 3. Add "//tags: model" to JS body and save the script.
 * 4. Access the saved model in Model Hub.
 * 5. Interact with the saved model (move Final-at slider).
 *
 * Implementation note on ordering: in the platform ScriptView, clicking Run can shift
 * focus to a Compute2 RichFunctionView (the script declares `//editor: ...RichFunctionViewEditor`),
 * after which the CodeMirror editor is no longer reachable from the active view. To keep all
 * steps UI-only while still covering both observable behaviours (script runs, tag persists),
 * we add `//tags: model` *first* (while the editor is fresh and focused) and then Run.
 */
test('DiffStudio Scripting — Edit toggle, </> JS view, Run, Save with //tags: model, Model Hub', async ({ page }) => {
  test.setTimeout(300_000);
  const { softStep, assertAllPassed } = createSoftStepCollector();
  const monitor = attachErrorMonitor(page);

  await softStep('Setup: Open Diff Studio + Bioreactor', async () => {
    await openDiffStudio(page);
    await openModelFromLibrary(page, 'Bioreactor');
  });

  await softStep('Step 1: Enable Edit toggle, then click </> to open the JS script view', async () => {
    const wasOn = await ribbonSwitchOn(page, 'Edit');
    if (!wasOn) await toggleRibbonSwitch(page, 'Edit');
    const editorVisible = await page.locator('.diff-studio-eqs-editor').count();
    expect(editorVisible).toBeGreaterThan(0);

    // Click the </> ribbon item by its label text
    await page.locator('.d4-ribbon-name', { hasText: '</>' }).first().click();
    // After the click, the platform switches to the Script editor — CodeMirror appears
    await page.locator('.CodeMirror').first().waitFor({ timeout: 30_000 });
    await page.waitForTimeout(1500);
  });

  await softStep('Step 2: Add "//tags: model" to the script body and Save', async () => {
    // Diff Studio's `</>` opens a platform ScriptView, which uses CodeMirror 5
    // (`/js/common/codemirror/codemirror.js`). CM5 renders editor as `.CodeMirror`.
    // The generated script always opens with this two-line header:
    //   line 1: //name: <model name>
    //   line 2: //language: javascript
    // Focus the editor, navigate to the end of line 2, insert a new line, type the tag.
    const cm = page.locator('.CodeMirror').first();
    await cm.waitFor({ timeout: 10_000 });
    await cm.click();
    await page.keyboard.press('Control+Home'); // top of file (line 1)
    await page.keyboard.press('ArrowDown');    // line 2 (//language: javascript)
    await page.keyboard.press('End');
    await page.keyboard.press('Enter');
    await page.keyboard.type('//tags: model');
    await page.waitForTimeout(500);

    // Verify the tag is now visible somewhere in the editor's rendered DOM
    const tagLine = page.locator('.CodeMirror-line', { hasText: '//tags: model' });
    await expect(tagLine.first()).toBeVisible({ timeout: 5_000 });

    // Click Save — the script editor's ribbon Save button is `[name="button-Save"]`.
    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(3000);
  });

  await softStep('Step 3: Run the script; move Final-at slider — live update', async () => {
    // The script editor exposes a play button as the "icon-play" name
    await page.locator('[name="icon-play"]').first().click();
    // Final input host appears once the script has run
    await page.locator(inputHost('Final')).waitFor({ timeout: 30_000 });
    await page.waitForTimeout(1500);

    const inp = page.locator(inputEditor('Final'));
    const before = await inp.inputValue();
    await setInputValue(page, 'Final', '500');
    await page.waitForTimeout(2500);
    const after = await inp.inputValue();
    expect(after).toBe('500');
    expect(after).not.toBe(before);
    // NOTE: chart redraw verification is manual (see ui-only.md M-1.5). The ScriptView
    // renders the chart inside a Compute2 RichFunctionView Vue component that exposes
    // neither a stable `.d4-viewer` container nor an updated `grok.shell.t` dataframe.

    // REMARK from MD: this UI exposes neither Process-mode nor Facet (only Multiaxis)
    await expect(page.locator(inputHost('Process-mode'))).toHaveCount(0);
  });

  await softStep('Step 4: Access the saved model in Model Hub (Apps > Compute > Model Hub)', async () => {
    // Navigate home, then walk the Browse tree: Apps > Compute > Model Hub.
    await page.goto(BASE);
    await page.waitForSelector('.d4-ribbon', { timeout: 30_000 });
    await page.waitForTimeout(1500);

    // Helper: scroll an element into view and click it. The Browse tree is sometimes virtualised
    // and labels can sit outside the viewport, so we drive the click via DOM scrollIntoView + click.
    const clickTreeLabel = async (label: string): Promise<boolean> => {
      return await page.evaluate((text) => {
        const candidates = Array.from(document.querySelectorAll(
          '.d4-tree-view-group-label, .d4-tree-view-node-label, .d4-tree-view-item-label')) as HTMLElement[];
        const el = candidates.find(e => e.textContent?.trim() === text);
        if (!el) return false;
        el.scrollIntoView({ behavior: 'instant', block: 'center' });
        el.click();
        return true;
      }, label);
    };

    expect(await clickTreeLabel('Apps')).toBe(true);
    await page.waitForTimeout(800);
    expect(await clickTreeLabel('Compute')).toBe(true);
    await page.waitForTimeout(800);
    expect(await clickTreeLabel('Model Hub')).toBe(true);
    await page.waitForTimeout(3000);

    // Model Hub renders saved models as cards (Compute2 Vue components) in the central grid.
    // Plain text match — the model name is the card's visible label.
    await expect(page.getByText('Bioreactor', { exact: true }).first())
      .toBeVisible({ timeout: 20_000 });
  });

  await softStep('Step 5: Interact with the saved model; adjust Final', async () => {
    // Pick the last "Bioreactor" entry — the newly saved one
    const items = page.getByText('Bioreactor', { exact: true });
    const count = await items.count();
    expect(count).toBeGreaterThan(0);
    await items.nth(count - 1).dblclick();

    await page.locator(inputHost('Final')).waitFor({ timeout: 30_000 });
    await setInputValue(page, 'Final', '800');
    await page.waitForTimeout(2500);
    const after = await page.locator(inputEditor('Final')).inputValue();
    expect(after).toBe('800');
    // NOTE: chart redraw verification is manual (see ui-only.md M-1.6) — same Vue isolation
    // as Step 3 above.
  });

  assertAllPassed();
  monitor.assertNone();
});
