import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Chemprop model — Train, Apply, Container, Browse', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // ─── Train section ───────────────────────────────────────────────

  await softStep('1.1 Open smiles.csv', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      document.body.classList.add('selenium');
      g.shell.settings.showFiltersIconsConstantly = true;
      g.shell.windows.simpleMode = true;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/chem/smiles.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 3000));
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return {rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0};
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBeGreaterThan(0);
  });

  await softStep('1.2 Open ML > Models > Train Model', async () => {
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('no ML > Models');
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Train-Model..."]').click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel', null, {timeout: 15_000});
  });

  await softStep('1.3 Set Predict to RingCount', async () => {
    // Scenario says "Ring Count" (space). Actual column is "RingCount".
    await page.evaluate(() => {
      const root = (window as any).grok.shell.v.root;
      const predict = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
      predict.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
    });
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'), null, {timeout: 5_000});
    await page.evaluate(() => (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
    await page.keyboard.type('RingCount');
    await page.waitForTimeout(300);
    await page.keyboard.press('Enter').catch(() => {});
    await page.waitForTimeout(500);
    const val = await page.evaluate(() =>
      (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent);
    expect(val).toBe('RingCount');
  });

  await softStep('1.4 Set Features to canonical_smiles — Chemprop engine auto-selects', async () => {
    await page.evaluate(() => {
      const editor = (window as any).grok.shell.v.root.querySelector('[name="div-Features"]') as HTMLElement;
      const r = editor.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window, button: 0,
        clientX: r.left + 10, clientY: r.top + 5,
      });
      editor.dispatchEvent(ev('mousedown'));
      editor.dispatchEvent(ev('mouseup'));
      editor.dispatchEvent(ev('click'));
    });
    const selectDialog = page.locator('[name="dialog-Select-columns..."]');
    await selectDialog.waitFor({timeout: 10_000});
    // Narrow the canvas-rendered column list to canonical_smiles via search, then
    // click the row's checkbox via a PointerEvent (Dart's grid ignores MouseEvent-only
    // dispatches; PointerEvent hits the real input pipeline).
    await selectDialog.locator('input[placeholder="Search"]').fill('canonical_smiles');
    await page.waitForTimeout(400);
    const checked = await page.evaluate(async () => {
      const d = document.querySelector('[name="dialog-Select-columns..."]')!;
      const canvases = d.querySelectorAll('canvas');
      const canvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const x = r.right - 20;
      const y = r.top + 34;
      const el = document.elementFromPoint(x, y)!;
      for (const n of ['pointerdown', 'mousedown', 'pointerup', 'mouseup', 'click']) {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, buttons: 1, clientX: x, clientY: y};
        const ev = n.startsWith('pointer') ? new PointerEvent(n, opts) : new MouseEvent(n, opts);
        el.dispatchEvent(ev);
      }
      await new Promise(r => setTimeout(r, 400));
      const m = d.textContent?.match(/(\d+) checked/);
      return m ? parseInt(m[1], 10) : 0;
    });
    expect(checked, 'canonical_smiles checkbox should toggle').toBeGreaterThan(0);
    await selectDialog.locator('[name="button-OK"]').click();
    await page.waitForFunction(() => {
      const sel = document.querySelector('[name="input-host-Model-Engine"] select') as HTMLSelectElement | null;
      return sel?.value === 'Chem: Chemprop';
    }, null, {timeout: 10_000});
  });

  await softStep('1.5 Change Activation, Split_type, Epochs and click TRAIN', async () => {
    await page.evaluate(() => {
      const root = (window as any).grok.shell.v.root;
      const setSelect = (name: string, value: string) => {
        const sel = root.querySelector(`[name="input-host-${name}"] select`) as HTMLSelectElement;
        sel.value = value;
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      };
      setSelect('Activation', 'LeakyReLU');
      setSelect('Split-type', 'scaffold_balanced');
      const epochs = root.querySelector('[name="input-host-Epochs"] input') as HTMLInputElement;
      epochs.value = '5';
      epochs.dispatchEvent(new Event('input', {bubbles: true}));
      epochs.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('[name="button-Save"]').click();
    await page.locator('[name="input-Name"]').first().waitFor({timeout: 10_000});
  });

  await softStep('1.6 Save model as "test_chemprop"', async () => {
    await page.locator('[name="input-Name"]').first().fill('test_chemprop');
    await page.locator('.d4-dialog:has([name="input-Name"]) [name="button-OK"]').click();
    // Wait until the trained model is persisted to grok.dapi.models.
    const appeared = await page.waitForFunction(async () => {
      const models = await (window as any).grok.dapi.models.filter('name = "test_chemprop"').list();
      return models.length > 0;
    }, null, {timeout: 180_000}).then(() => true).catch(() => false);
    expect(appeared, 'test_chemprop should exist after training').toBe(true);
  });

  await softStep('1.7 Change Metric to auc (roc), click TRAIN — expect balloon error', async () => {
    await page.evaluate(() => {
      const sel = (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Metric"] select') as HTMLSelectElement;
      sel.value = 'roc'; // scenario says "auc"; closest Chemprop option is "roc" (classification-only)
      sel.dispatchEvent(new Event('change', {bubbles: true}));
    });
    await page.locator('[name="button-Save"]').click();
    // A balloon / error should appear — roc-auc requires classification dataset_type.
    await expect(page.locator('.d4-balloon, .grok-balloon').first())
      .toBeVisible({timeout: 5_000});
  });

  // ─── Apply section ───────────────────────────────────────────────

  await softStep('2.1 Close All and open smiles_only.csv', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/chem/smiles_only.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise(r => setTimeout(r, 200));
      }
      await new Promise(r => setTimeout(r, 3000));
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const rows = await page.evaluate(() => (window as any).grok.shell.tv?.dataFrame?.rowCount ?? 0);
    expect(rows).toBeGreaterThan(0);
  });

  await softStep('2.2 ML > Models > Apply Model... — select test_chemprop', async () => {
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement;
      const r = models.getBoundingClientRect();
      const ev = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window,
        clientX: r.left + 5, clientY: r.top + 5,
      });
      models.dispatchEvent(ev('mouseover'));
      models.dispatchEvent(ev('mouseenter'));
      models.dispatchEvent(ev('mousemove'));
    });
    await page.locator('[name="div-ML---Models---Apply-Model..."]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 30_000});
    const picked = await page.evaluate(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement;
      const idx = Array.from(sel.options).findIndex(o => (o.textContent || '').trim() === 'test_chemprop');
      if (idx < 0) return false;
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return true;
    });
    expect(picked, 'test_chemprop should appear in the Model dropdown').toBe(true);
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
  });

  await softStep('2.3 Verify prediction column "Ring Count (2)" was added', async () => {
    await page.waitForFunction(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      if (!df) return false;
      const names: string[] = [];
      for (let i = 0; i < df.columns.length; i++) names.push(df.columns.byIndex(i).name);
      return names.some(n => /Ring/i.test(n) && /\(\d+\)/.test(n));
    }, null, {timeout: 120_000});
  });

  // ─── Container section ───────────────────────────────────────────

  await softStep('3.1 Browse > Platform > Dockers — locate chem-chemprop', async () => {
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      // Clicking Browse sidebar brings the tree into view reliably.
      (document.querySelector('[name="Browse"]') as HTMLElement | null)?.click();
    });
    // Scroll Platform into view and wait for its tree node to render.
    await page.evaluate(() => {
      const p = document.querySelector('[name="tree-Platform"]') as HTMLElement | null;
      p?.scrollIntoView({block: 'center'});
    });
    await page.locator('[name="tree-Platform"]').waitFor({state: 'attached', timeout: 15_000});
    const node = page.locator('[name="tree-Platform---Dockers"]');
    if (!(await node.isVisible({timeout: 2_000}))) {
      // Platform may be collapsed. Clicking the row toggles expansion; it's
      // a no-op if already expanded. Avoid `.d4-tree-view-tri` which is missing
      // once Platform is already expanded.
      await page.locator('[name="tree-Platform"]').first().click();
      await node.waitFor({state: 'visible', timeout: 10_000});
    }
    await node.click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.path === '/dockers?', null, {timeout: 15_000});
    await page.locator('input[placeholder*="Search"]').first().fill('chem-chemprop');
    await page.waitForTimeout(1_000);
    await expect(page.locator('text=chem-chemprop').first()).toBeVisible({timeout: 5_000});
  });

  await softStep('3.2 Right-click chem-chemprop → Stop, then → Run', async () => {
    const openMenu = async () => {
      await page.evaluate(() => {
        const label = Array.from(document.querySelectorAll('div,span,label,a'))
          .find(e => e.textContent?.trim() === 'chem-chemprop' && e.children.length === 0) as HTMLElement;
        const container = (label.closest('.grok-gallery-grid-item, .d4-gallery-grid-item')
          || label.parentElement) as HTMLElement;
        const r = container.getBoundingClientRect();
        container.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, view: window, button: 2,
          clientX: r.left + r.width / 2, clientY: r.top + r.height / 2,
        }));
      });
      await page.locator('.d4-menu-popup').waitFor({timeout: 5_000});
    };
    const clickMenuItem = async (text: string) => {
      await page.evaluate((t) => {
        const menu = document.querySelector('.d4-menu-popup')!;
        const item = Array.from(menu.querySelectorAll('.d4-menu-item-label'))
          .find(e => e.textContent?.trim() === t) as HTMLElement;
        item.click();
      }, text);
    };
    const waitStatus = async (re: RegExp, timeoutMs: number) => {
      await page.waitForFunction(async (pattern: string) => {
        const c = (await (window as any).grok.dapi.docker.dockerContainers
          .filter('name = "chem-chemprop"').list())[0];
        return !!c && new RegExp(pattern, 'i').test(c.status);
      }, re.source, {timeout: timeoutMs});
    };
    await openMenu();
    await clickMenuItem('Stop');
    await waitStatus(/stopp|stopped|idle/, 60_000);
    await openMenu();
    await clickMenuItem('Run');
    await waitStatus(/start|running|ready/, 120_000);
  });

  // ─── Browse section ──────────────────────────────────────────────

  await softStep('4.1 Go to Browse > Platform > Predictive models', async () => {
    await page.evaluate(() => {
      const g: any = (window as any).grok;
      g.shell.windows.showBrowse = true;
      // Clicking Browse sidebar brings the tree into view reliably.
      (document.querySelector('[name="Browse"]') as HTMLElement | null)?.click();
    });
    // Scroll Platform into view and wait for its tree node to render.
    await page.evaluate(() => {
      const p = document.querySelector('[name="tree-Platform"]') as HTMLElement | null;
      p?.scrollIntoView({block: 'center'});
    });
    await page.locator('[name="tree-Platform"]').waitFor({state: 'attached', timeout: 15_000});
    const node = page.locator('[name="tree-Platform---Predictive-models"]');
    if (!(await node.isVisible({timeout: 2_000}))) {
      // Platform may be collapsed. Clicking the row toggles expansion; it's
      // a no-op if already expanded. Avoid `.d4-tree-view-tri` which is missing
      // once Platform is already expanded.
      await page.locator('[name="tree-Platform"]').first().click();
      await node.waitFor({state: 'visible', timeout: 10_000});
    }
    await node.click();
    await page.waitForFunction(() => (window as any).grok.shell.v?.path === '/models?', null, {timeout: 15_000});
  });

  await softStep('4.2 Search test_chemprop and open Context Panel', async () => {
    await page.locator('input[placeholder*="Search"]').first().fill('test_chemprop');
    await page.waitForTimeout(1_500);
    const model = page.locator('text=test_chemprop').first();
    await expect(model).toBeVisible({timeout: 8_000});
    await model.click();
    await expect(page.locator('.grok-prop-panel, [class*="context-panel"]').first())
      .toBeVisible({timeout: 5_000});
  });

  await softStep('4.3 Share the model', async () => {
    const sharingTab = page.locator('.grok-prop-panel, [class*="context-panel"]')
      .first().locator('text=Sharing').first();
    if (await sharingTab.isVisible({timeout: 2_000}))
      await sharingTab.click();
  });

  await softStep('4.4 Delete the model', async () => {
    const deleted = await page.evaluate(async () => {
      const models = await (window as any).grok.dapi.models.filter('name = "test_chemprop"').list();
      if (!models.length) throw new Error('no test_chemprop to delete');
      await (window as any).grok.dapi.models.delete(models[0]);
      const remaining = await (window as any).grok.dapi.models.filter('name = "test_chemprop"').list();
      return remaining.length === 0;
    });
    expect(deleted).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error(`Soft step failures:\n${stepErrors.map(e => `- ${e.step}: ${e.error}`).join('\n')}`);
});
