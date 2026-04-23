import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('AddNewColumn/formula-refreshing — dependency propagation across calculated columns', async ({page}) => {
  test.setTimeout(300_000);

  await page.bringToFront();
  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    const tv = g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    tv.getFiltersGroup();
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30_000});

  async function openAddNewColumnDialog() {
    await page.evaluate(() => {
      const icon = document.querySelector('[name="icon-add-new-column"]') as HTMLElement;
      const r = icon.getBoundingClientRect();
      ['mousedown', 'mouseup', 'click'].forEach(t =>
        icon.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 0, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2})),
      );
    });
    await page.locator('[name="dialog-Add-New-Column"]').last().waitFor({timeout: 15_000});
    await page.waitForFunction(() => {
      const dlgs = document.querySelectorAll('[name="dialog-Add-New-Column"]');
      const dlg = dlgs[dlgs.length - 1];
      return !!dlg?.querySelector('.cm-content') && !!dlg?.querySelector('[name="button-Add-New-Column---OK"]');
    }, null, {timeout: 15_000});
  }

  async function dispatchClickInLastDialog(selectorInsideDialog: string) {
    await page.evaluate((sel) => {
      const dlgs = document.querySelectorAll('[name="dialog-Add-New-Column"]');
      const dlg = dlgs[dlgs.length - 1];
      const el = dlg?.querySelector(sel) as HTMLElement | null;
      if (!el) throw new Error('Not found: ' + sel);
      const r = el.getBoundingClientRect();
      const opts = {bubbles: true, cancelable: true, button: 0, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2};
      el.dispatchEvent(new MouseEvent('mousedown', opts));
      el.dispatchEvent(new MouseEvent('mouseup', opts));
      el.dispatchEvent(new MouseEvent('click', opts));
      el.focus?.();
    }, selectorInsideDialog);
  }

  async function addCalcColumn(colName: string, formula: string) {
    await openAddNewColumnDialog();

    // CodeMirror focus — direct JS dispatch bypasses any overlay
    await dispatchClickInLastDialog('.cm-content');
    await page.waitForTimeout(200);
    await page.keyboard.type(formula);

    await dispatchClickInLastDialog('[name="input-Add-New-Column---Name"]');
    await page.keyboard.press('Control+A');
    await page.keyboard.type(colName);

    await page.waitForFunction((f) => {
      const dlgs = document.querySelectorAll('[name="dialog-Add-New-Column"]');
      const cm = dlgs[dlgs.length - 1]?.querySelector('.cm-content');
      return cm?.textContent?.trim() === f;
    }, formula, {timeout: 5_000});

    await dispatchClickInLastDialog('[name="button-Add-New-Column---OK"]');
    await page.waitForFunction((n) => {
      const g: any = (window as any).grok;
      const col = g.shell.tv?.dataFrame?.col?.(n);
      return !!col && document.querySelectorAll('[name="dialog-Add-New-Column"]').length === 0;
    }, colName, {timeout: 30_000});

    const info = await page.evaluate((n) => {
      const g: any = (window as any).grok;
      const col = g.shell.tv.dataFrame.col(n);
      return {exists: !!col, formula: col?.getTag('formula')};
    }, colName);
    expect(info.exists).toBe(true);
    expect(info.formula).toBe(formula);
  }

  async function verifySample(
    targetCol: string,
    srcCol: string,
    compute: (src: number) => number,
    tolerance = 1e-3,
  ) {
    const rows = await page.evaluate(([t, s]) => {
      const g: any = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      const tc = df.col(t);
      const sc = df.col(s);
      return Array.from({length: 5}, (_, i) => ({src: sc.get(i), target: tc.get(i)}));
    }, [targetCol, srcCol]);
    for (const r of rows)
      expect(Math.abs(r.target - compute(r.src))).toBeLessThan(tolerance);
  }

  async function editFormulaViaContextPanel(colName: string, newFormula: string) {
    // Capture the target column's CURRENT formula before selecting — the pane refresh is async,
    // so we must wait for CM content to match this value before typing (otherwise we edit stale text).
    const currentFormula = await page.evaluate((n) => {
      const g: any = (window as any).grok;
      return g.shell.tv.dataFrame.col(n)?.getTag('formula') ?? '';
    }, colName);

    await page.evaluate((n) => {
      const g: any = (window as any).grok;
      g.shell.o = g.shell.tv.dataFrame.col(n);
    }, colName);

    // Expand the Formula pane if collapsed, then wait for its CodeMirror to render
    await page.waitForFunction(() => !!document.querySelector('.d4-pane-formula'), null, {timeout: 15_000});
    await page.evaluate(() => {
      const pane = document.querySelector('.d4-pane-formula') as HTMLElement | null;
      if (pane && !pane.classList.contains('expanded')) {
        const header = pane.parentElement?.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        header?.click();
      }
    });
    // Wait until CM content reflects this column's current formula (not the previously-selected column's)
    await page.waitForFunction((f) => {
      const cm = document.querySelector('.d4-pane-formula .cm-content');
      return cm?.textContent?.trim() === f;
    }, currentFormula, {timeout: 30_000});

    const dispatchInFormulaPane = async (sel: string) => {
      await page.evaluate((s) => {
        const el = document.querySelector('.d4-pane-formula ' + s) as HTMLElement | null;
        if (!el) throw new Error('Not found in Formula pane: ' + s);
        const r = el.getBoundingClientRect();
        const opts = {bubbles: true, cancelable: true, button: 0, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2};
        el.dispatchEvent(new MouseEvent('mousedown', opts));
        el.dispatchEvent(new MouseEvent('mouseup', opts));
        el.dispatchEvent(new MouseEvent('click', opts));
        el.focus?.();
      }, sel);
    };

    await dispatchInFormulaPane('.cm-content');
    await page.keyboard.press('Control+A');
    await page.keyboard.type(newFormula);

    await page.waitForFunction((f) => {
      const cm = document.querySelector('.d4-pane-formula .cm-content');
      return cm?.textContent?.trim() === f;
    }, newFormula, {timeout: 5_000});

    await dispatchInFormulaPane('[name="button-Apply"]');

    // Wait until the formula tag on the column reflects the new formula (platform has applied it)
    await page.waitForFunction(([n, f]) => {
      const g: any = (window as any).grok;
      return g.shell.tv.dataFrame.col(n)?.getTag('formula') === f;
    }, [colName, newFormula], {timeout: 15_000});

    // In headed mode, the Formula pane's CM retains stale content until the next user interaction.
    // Dismiss it actively: move selection off the column, so the next edit starts with a fresh pane.
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      g.shell.o = g.shell.tv.dataFrame;
    });
    await page.waitForFunction(() => !document.querySelector('.d4-pane-formula'), null, {timeout: 5_000})
      .catch(() => {/* pane may persist briefly; not fatal */});

    const applied = await page.evaluate((n) => {
      const g: any = (window as any).grok;
      return g.shell.tv.dataFrame.col(n)?.getTag('formula');
    }, colName);
    expect(applied).toBe(newFormula);
  }

  await softStep('Step 1: demog dataset open with WEIGHT column', async () => {
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      const df = g.shell.tv.dataFrame;
      return {rows: df.rowCount, hasWeight: !!df.col('WEIGHT')};
    });
    expect(info.hasWeight).toBe(true);
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('Step 2: create Weight2 = ${WEIGHT} + 100', async () => {
    await addCalcColumn('Weight2', '${WEIGHT} + 100');
    await verifySample('Weight2', 'WEIGHT', w => w + 100);
  });

  await softStep('Step 3: create Weight3 = ${Weight2} + 100', async () => {
    await addCalcColumn('Weight3', '${Weight2} + 100');
    await verifySample('Weight3', 'Weight2', w2 => w2 + 100);
  });

  await softStep('Step 4: create Weight4 = Log10(${Weight3}) - 0.2', async () => {
    await addCalcColumn('Weight4', 'Log10(${Weight3}) - 0.2');
    await verifySample('Weight4', 'Weight3', w3 => Math.log10(w3) - 0.2);
  });

  await softStep('Step 5a: editing Weight2 formula propagates to Weight3 and Weight4', async () => {
    await editFormulaViaContextPanel('Weight2', '${WEIGHT} + 200');
    await verifySample('Weight2', 'WEIGHT', w => w + 200);
    await verifySample('Weight3', 'WEIGHT', w => (w + 200) + 100);
    await verifySample('Weight4', 'WEIGHT', w => Math.log10((w + 200) + 100) - 0.2, 1e-2);
  });

  await softStep('Step 5b: editing Weight3 formula propagates to Weight4', async () => {
    await editFormulaViaContextPanel('Weight3', '${Weight2} * 2');
    await verifySample('Weight3', 'Weight2', w2 => w2 * 2);
    await verifySample('Weight4', 'Weight2', w2 => Math.log10(w2 * 2) - 0.2, 1e-2);
  });

  await softStep('Step 5c: editing Weight4 formula applies cleanly', async () => {
    await editFormulaViaContextPanel('Weight4', 'Log10(${Weight3}) + 1');
    await verifySample('Weight4', 'Weight3', w3 => Math.log10(w3) + 1);
  });

  await page.evaluate(() => {
    const g: any = (window as any).grok;
    g.shell.closeAll();
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n'));
});
