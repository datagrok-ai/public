/* ---
sub_features_covered: [powerpack.dialogs.add-new-column, powerpack.dialogs.add-new-column-func]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
test.use(specTestOptions);
test('PowerPack: Add new column autocomplete (demog — type, Ctrl+Space, $)', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  await loginToDatagrok(page);
  await page.evaluate(async () => {
    const grok = (window as any).grok;
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    try { grok.shell.closeAll(); } catch (_) {}
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
  expect(cols).toContain('HEIGHT');
  expect(cols).toContain('WEIGHT');
  expect(cols).toContain('AGE');
  const openDialog = async () => {
    const icon = page.locator('[name="icon-add-new-column"]').first();
    await icon.waitFor({timeout: 30_000, state: 'visible'});
    await icon.click({timeout: 10_000});
    const dlg = page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first();
    await dlg.waitFor({timeout: 30_000});
    return dlg;
  };
  const dlg = await openDialog();
  await expect(dlg).toBeVisible();
  const cm = dlg.locator('.add-new-column-dialog-cm-div .cm-content').first();
  await cm.waitFor({timeout: 15_000, state: 'visible'});
  const readDoc = async () => page.evaluate(() => {
    const cmContent = document.querySelector('.d4-dialog .add-new-column-dialog-cm-div .cm-content') as HTMLElement | null;
    if (!cmContent) return '';
    const lines = Array.from(cmContent.querySelectorAll('.cm-line')) as HTMLElement[];
    const text = lines.length > 0
      ? lines.map((l) => l.textContent || '').join('\n')
      : (cmContent.textContent || '');
    return text.replace(/​/g, '');
  });
  const closeDialog = async () => {
    await page.evaluate(() => {
      const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
      if (cancel) { cancel.click(); return; }
      const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
      if (anyCancel) { anyCancel.click(); return; }
    }).catch(() => {});
    await page.keyboard.press('Escape').catch(() => {});
    await page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).first()
      .waitFor({state: 'detached', timeout: 10_000})
      .catch(async () => {
        for (let i = 0; i < 5; i++) {
          if (await page.locator('.d4-dialog').filter({hasText: 'Add New Column'}).count() === 0) break;
          await page.keyboard.press('Escape');
          await page.waitForTimeout(150);
        }
      });
    await page.waitForTimeout(200);
  };
  let dlgRef = dlg;
  let cmRef = cm;
  const resetDialog = async () => {
    await closeDialog();
    dlgRef = await openDialog();
    await expect(dlgRef).toBeVisible();
    cmRef = dlgRef.locator('.add-new-column-dialog-cm-div .cm-content').first();
    await cmRef.waitFor({timeout: 15_000, state: 'visible'});
    await cmRef.click();
    await page.waitForTimeout(150);
    const doc = await readDoc();
    if (doc.length > 0) {
      await page.keyboard.press('Control+A');
      await page.keyboard.press('Delete');
      await page.waitForTimeout(120);
    }
  };
  const dismissTooltip = async () => {
    await page.keyboard.press('Escape');
    await page.waitForTimeout(150);
  };
  await softStep('Scenario 1 Step 1-2: focus editor and type "a"', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    await page.keyboard.type('a', {delay: 60});
  });
  await softStep('Scenario 1 Step 3: verify autocomplete tooltip appears listing functions starting with "a"', async () => {
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    const entries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      const items = Array.from(tip.querySelectorAll('li'));
      return items.map((li) => (li.textContent || '').trim()).filter((t) => t.length > 0);
    });
    expect(entries.length).toBeGreaterThan(0);
    const hasFnStartingWithA = entries.some((t) => /^[Aa]/.test(t));
    expect(hasFnStartingWithA).toBe(true);
  });
  await softStep('Scenario 1 Step 4a: select function from list via Enter on highlighted entry', async () => {
    await page.keyboard.press('Enter');
    await page.waitForTimeout(300);
  });
  await softStep('Scenario 1 Step 5 (via Enter): verify function inserted as Name(params)', async () => {
    const doc = await readDoc();
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(normalized.length).toBeGreaterThan(0);
    expect(/^[Aa][A-Za-z0-9_]*\([^)]*\)/.test(normalized)).toBe(true);
  });
  await softStep('Scenario 1 Step 6 (post-Enter): remove inserted function from editor', async () => {
    await dismissTooltip();
    await resetDialog();
    const doc = await readDoc();
    expect(doc.length).toBe(0);
  });
  await softStep('Scenario 1 Step 4b: re-trigger by typing "a", then select via mouse click on entry', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    await page.keyboard.type('a', {delay: 60});
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    const clicked = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return false;
      const items = Array.from(tip.querySelectorAll('li'))
        .filter((li) => ((li as HTMLElement).offsetParent !== null) && (li.textContent || '').trim().length > 0);
      if (items.length === 0) return false;
      const first = items[0] as HTMLElement;
      first.dispatchEvent(new MouseEvent('mousedown', {bubbles: true}));
      first.dispatchEvent(new MouseEvent('mouseup', {bubbles: true}));
      first.click();
      return true;
    });
    expect(clicked).toBe(true);
    await page.waitForTimeout(300);
  });
  await softStep('Scenario 1 Step 5 (via mouse click): verify function inserted as Name(params)', async () => {
    const doc = await readDoc();
    const normalized = doc.replace(/\s+/g, ' ').trim();
    expect(normalized.length).toBeGreaterThan(0);
    expect(/^[Aa][A-Za-z0-9_]*\([^)]*\)/.test(normalized)).toBe(true);
  });
  await softStep('Scenario 1 Step 6 (post-click): remove inserted function from editor', async () => {
    await dismissTooltip();
    await resetDialog();
    const doc = await readDoc();
    expect(doc.length).toBe(0);
  });
  await softStep('Scenario 2 Step 1: with editor focused and empty, press Ctrl+Space', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    const beforeDoc = await readDoc();
    expect(beforeDoc.length).toBe(0);
    await page.keyboard.press('Control+Space');
  });
  await softStep('Scenario 2 Step 2: verify autocomplete tooltip appears with full function list', async () => {
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    const entryCount = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return 0;
      return Array.from(tip.querySelectorAll('li'))
        .filter((li) => (li.textContent || '').trim().length > 0)
        .length;
    });
    expect(entryCount).toBeGreaterThan(0);
    const sampleEntries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      return Array.from(tip.querySelectorAll('li'))
        .map((li) => (li.textContent || '').trim())
        .filter((t) => t.length > 0)
        .slice(0, 20);
    });
    expect(sampleEntries.length).toBeGreaterThan(0);
  });
  await softStep('Scenario 3 Step 1 (pre): dismiss any open tooltip with Escape', async () => {
    await dismissTooltip();
    await resetDialog();
  });
  await softStep('Scenario 3 Step 1: with editor focused, type "$"', async () => {
    await cmRef.click();
    await page.waitForTimeout(120);
    await page.keyboard.type('$', {delay: 60});
  });
  await softStep('Scenario 3 Step 2: verify tooltip lists dataset columns (HEIGHT, WEIGHT, AGE) — distinct from function list', async () => {
    const tooltip = page.locator('.cm-tooltip-autocomplete').first();
    await tooltip.waitFor({timeout: 5_000, state: 'visible'});
    await expect(tooltip).toBeVisible();
    const entries = await page.evaluate(() => {
      const tip = document.querySelector('.cm-tooltip-autocomplete');
      if (!tip) return [];
      return Array.from(tip.querySelectorAll('li'))
        .map((li) => (li.textContent || '').trim())
        .filter((t) => t.length > 0);
    });
    expect(entries.length).toBeGreaterThan(0);
    const hasHeight = entries.some((e) => /\bHEIGHT\b/.test(e));
    const hasWeight = entries.some((e) => /\bWEIGHT\b/.test(e));
    const hasAge = entries.some((e) => /\bAGE\b/.test(e));
    expect(hasHeight).toBe(true);
    expect(hasWeight).toBe(true);
    expect(hasAge).toBe(true);
    const looksLikeColumnList = hasHeight && hasWeight && hasAge;
    expect(looksLikeColumnList).toBe(true);
  });
  await page.evaluate(() => {
    const cancel = document.querySelector('.d4-dialog [name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    if (cancel) cancel.click();
    const anyCancel = document.querySelector('.d4-dialog [name="button-CANCEL"]') as HTMLElement | null;
    if (anyCancel) anyCancel.click();
  }).catch(() => {});
  await page.evaluate(() => {
    try { (window as any).grok.shell.closeAll(); } catch (_) {  }
  }).catch(() => {});
  finishSpec();
});
