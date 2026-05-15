import { test, expect } from '@playwright/test';
import { createSoftStepCollector } from './helpers/soft-step';
import { attachErrorMonitor } from './helpers/error-monitor';
import {
  openDiffStudio, openModelFromLibrary, openModelHub, openModelHubCard,
  modelHubCardCount, waitForModelScript,
  toggleRibbonSwitch, ribbonSwitchOn, setInputValue, inputEditor, inputHost,
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
    // The platform Save button fires `grok.dapi.scripts.save(...)` and may not
    // be settled by the time Step 4 navigates with `page.goto(BASE)`. Block here
    // until the server confirms the script exists with tag `model`, so the
    // subsequent goto does not abort the still-in-flight save POST.
    await waitForModelScript(page, 'Bioreactor');
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
    // Use the JS-API helper — Browse tree navigation is unreliable on cold CI Datlas
    // because the 'Apps' / 'Compute' / 'Model Hub' labels may not be mounted when the
    // test reaches this step.
    await openModelHub(page);
    // Model Hub renders saved models as `.d4-list-item` entries with the model name
    // as the visible label.
    expect(await modelHubCardCount(page, 'Bioreactor')).toBeGreaterThan(0);
  });

  await softStep('Step 5: Interact with the saved model; adjust Final', async () => {
    // Pick the last "Bioreactor" card (the newly saved one) and open it; cards require
    // both `click()` and a `dispatchEvent('dblclick')` to navigate.
    expect(await modelHubCardCount(page, 'Bioreactor')).toBeGreaterThan(0);
    await openModelHubCard(page, 'Bioreactor');

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
