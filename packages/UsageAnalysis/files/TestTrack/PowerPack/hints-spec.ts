/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
test('PowerPack: Add new column — hover inserted function name surfaces signature tooltip (demog)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await softStep('Step 1: open demog dataset; verify grid renders', async () => {
    await page.evaluate(async () => {
      const grok = (window as any).grok;
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      grok.shell.windows.simpleMode = true;
      try { grok.shell.closeAll(); } catch (_) {  }
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    await page.waitForTimeout(1000);
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? df.columns.names() : [];
    });
    expect(cols.length).toBeGreaterThan(0);
  });
  await softStep('Step 2: open Add New Column dialog via toolbar icon', async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    await expect(dlg).toBeVisible();
    const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cm.waitFor({timeout: 15_000, state: 'visible'});
  });
  const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await softStep('Step 3: type "a", accept Abs from autocomplete tooltip; verify "Abs(num)" inserted', async () => {
    await cm.click();
    await page.waitForTimeout(150);
    await page.keyboard.press('Control+A');
    await page.keyboard.press('Delete');
    await page.waitForTimeout(100);
    await page.keyboard.type('a', {delay: 60});
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    const clickedAbs = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return false;
      const items = Array.from(tip.querySelectorAll('li'));
      const abs = items.find((li) => /\bAbs\b/.test(li.textContent || ''));
      if (!abs) return false;
      abs.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      abs.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      (abs as HTMLElement).click();
      return true;
    });
    if (!clickedAbs)
      await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
    const doc = await page.evaluate(() => {
      const cmDiv = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return '';
      const view = (cmDiv as any).cmView?.view ?? (cmDiv.parentElement as any)?.cmView?.view ?? null;
      return view ? view.state.doc.toString() : (cmDiv.innerText || '');
    });
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(/^Abs\([^)]*\)/.test(normalized)).toBe(true);
    await page.keyboard.press('Escape');
    await page.waitForTimeout(150);
  });
  await softStep('Step 4: hover over inserted function name; verify signature tooltip surfaces', async () => {
    const rect = await page.evaluate(() => {
      const cmDiv = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
      if (!cmDiv) return null;
      const walker = document.createTreeWalker(cmDiv, NodeFilter.SHOW_TEXT);
      let node: Node | null;
      while ((node = walker.nextNode())) {
        const t = (node as Text).data;
        const idx = t.indexOf('Abs');
        if (idx >= 0) {
          const range = document.createRange();
          range.setStart(node, idx);
          range.setEnd(node, idx + 3);
          const r = range.getBoundingClientRect();
          if (r.width > 0 && r.height > 0)
            return {x: r.left + r.width / 2, y: r.top + r.height / 2};
        }
      }
      return null;
    });
    expect(rect).not.toBeNull();
    await page.mouse.move(10, 10);
    await page.waitForTimeout(120);
    await page.mouse.move(rect!.x, rect!.y, {steps: 10});
    const hover = page.locator('.cm-tooltip.cm-tooltip-hover').first();
    await hover.waitFor({timeout: 5_000, state: 'visible'});
    await expect(hover).toBeVisible();
    const tooltipText = (await hover.textContent() ?? '').trim();
    expect(tooltipText.length).toBeGreaterThan(0);
    expect(/\bAbs\b/.test(tooltipText)).toBe(true);
    expect(/Abs\s*\([^)]*\)/.test(tooltipText)).toBe(true);
  });
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {  });
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) {  }
  }).catch(() => {  });
  finishSpec();
});
