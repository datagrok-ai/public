import {test, expect, chromium} from '@playwright/test';

const baseUrl = process.env.DATAGROK_URL ?? 'https://dev.datagrok.ai';
const datasetPath = 'System:DemoFiles/demog.csv';

const stepErrors: {step: string; error: string}[] = [];

async function softStep(name: string, fn: () => Promise<void>) {
  try {
    await test.step(name, fn);
  } catch (e: any) {
    stepErrors.push({step: name, error: e.message ?? String(e)});
    console.error(`[STEP FAILED] ${name}: ${e.message ?? e}`);
  }
}

test('PowerPack: Add New Columns', async () => {
  // Connect to existing Chrome via CDP and reuse the Datagrok page
  const browser = await chromium.connectOverCDP('http://localhost:9222');
  const context = browser.contexts()[0];
  let page = context.pages().find(p => p.url().includes('datagrok'));
  if (!page) {
    page = await context.newPage();
    await page.goto(baseUrl, {waitUntil: 'networkidle', timeout: 60000});
    await page.waitForFunction(() => {
      try { return typeof grok !== 'undefined' && typeof grok.shell.closeAll === 'function'; }
      catch { return false; }
    }, {timeout: 45000});
  }

  // Close any open dialogs, then open dataset
  await page.evaluate(async (path) => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click();
    });
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    const tv = grok.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some(c => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 5000));
    }
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  // Step 1: Verify dataset opened correctly
  await softStep('Open the Demog Dataset', async () => {
    const info = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const cols = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name);
      return { rows: df.rowCount, cols };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toContain('HEIGHT');
    expect(info.cols).toContain('WEIGHT');
  });

  // Step 2: Press "Add new column" icon
  await softStep('Press Add New Column icon, dialog opens', async () => {
    await page!.evaluate(() => {
      document.querySelector('[name="icon-add-new-column"]')!.click();
    });
    await page!.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 5000});
  });

  // Step 3: UI Check
  await softStep('UI Check: no overlapping, resize', async () => {
    const dialog = page!.locator('[name="dialog-Add-New-Column"]');
    await expect(dialog).toBeVisible();

    // Check no overflow in dialog contents
    const hasOverflow = await page!.evaluate(() => {
      const d = document.querySelector('[name="dialog-Add-New-Column"]');
      const contents = d!.querySelector('.d4-dialog-contents');
      return contents ? (contents.scrollHeight > contents.clientHeight + 5) : false;
    });
    expect(hasOverflow).toBe(false);

    // Check CodeMirror editor is visible
    await expect(page!.locator('[name="dialog-Add-New-Column"] .cm-content')).toBeVisible();
  });

  // Step 4: Add column with formula
  await softStep('Add column "New" with Round(HEIGHT+WEIGHT)', async () => {
    // Set column name
    const nameInput = page!.locator('[name="input-Add-New-Column---Name"]');
    await nameInput.click();
    await page!.keyboard.press('Control+A');
    await page!.keyboard.type('New');

    // Click CodeMirror editor and type formula
    const cmEditor = page!.locator('[name="dialog-Add-New-Column"] .cm-content');
    await cmEditor.click();
    await page!.keyboard.type('Round(${HEIGHT} + ${WEIGHT})');
    await page!.waitForTimeout(500);

    // Click OK
    await page!.evaluate(() => {
      document.querySelector('[name="button-Add-New-Column---OK"]')!.click();
    });
    await page!.waitForTimeout(1000);

    // Verify column was created with correct values
    const result = await page!.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const hasNew = df.columns.contains('New');
      if (!hasNew) return {hasNew: false};
      const h = df.col('HEIGHT');
      const w = df.col('WEIGHT');
      const n = df.col('New');
      return {
        hasNew: true,
        row0: {height: h.get(0), weight: w.get(0), newVal: n.get(0), expected: Math.round(h.get(0) + w.get(0))},
      };
    });
    expect(result.hasNew).toBe(true);
    expect(result.row0!.newVal).toBe(result.row0!.expected);
  });

  // Step 5: Recent Activity
  await softStep('Recent Activities: reopen dialog, history autofill', async () => {
    // Reopen the Add New Column dialog
    await page!.evaluate(() => {
      document.querySelector('[name="icon-add-new-column"]')!.click();
    });
    await page!.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 5000});

    // Click history icon
    await page!.evaluate(() => {
      document.querySelector('[name="dialog-Add-New-Column"] [name="icon-history"]')!.click();
    });
    await page!.waitForTimeout(500);

    // Click the most recent history menu item
    const historyItem = page!.locator('.d4-menu-popup .d4-menu-item-vert').first();
    await historyItem.click();
    await page!.waitForTimeout(500);

    // Verify autofill
    const autofilled = await page!.evaluate(() => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]');
      const nameInput = document.querySelector('[name="input-Add-New-Column---Name"]') as HTMLInputElement;
      const cmContent = dialog!.querySelector('.cm-content');
      return {
        name: nameInput?.value,
        formula: cmContent?.textContent
      };
    });
    expect(autofilled.name).toBe('New');
    expect(autofilled.formula).toContain('Round');
    expect(autofilled.formula).toContain('HEIGHT');
    expect(autofilled.formula).toContain('WEIGHT');

    // Close the dialog
    await page!.keyboard.press('Escape');
  });

  // Summary check
  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
