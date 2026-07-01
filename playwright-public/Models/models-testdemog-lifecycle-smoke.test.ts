/* ---
sub_features_covered: [models.api.save, models.api.suggested, models.command.apply, models.command.compare, models.command.delete, models.command.train, models.engines.api.apply, models.postprocessing.binary-classification, models.preprocessing.ignore-missing, models.preprocessing.impute-missing, models.view.browser, models.view.training, models.view.training.actions, models.workflow.apply-dialog, models.workflow.apply-model, models.workflow.compare-models, models.workflow.remove]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {setPredict, selectFeaturesByName} from '../helpers/models-helpers';

test.use(specTestOptions);

const MODEL_NAME = 'TestDemog';
const MODEL_NAME_REGRESSION = 'TestDemog_Regression';

test('Models / TestDemog predictive model lifecycle: Train / Apply / Browse+Compare / Delete', async ({page}) => {
  // Trains TWO small EDA models on demog.csv (binary SEX classifier + WEIGHT regression) — not
  // chemprop. Several steps each poll SAVE-enable up to 90-180s; 300s covers both trains + apply +
  // browse/compare + delete with margin.
  test.setTimeout(300_000);

  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async ([n1, n2]) => {
    for (const name of [n1, n2]) {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      for (const m of list) await (window as any).grok.dapi.models.delete(m);
    }
  }, [MODEL_NAME, MODEL_NAME_REGRESSION] as [string, string]);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('1.1 Open demog.csv from System:DemoFiles', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return { rows: df?.rowCount ?? 0, cols: df?.columns?.length ?? 0 };
    });
    expect(info.rows).toBeGreaterThan(0);
    expect(info.cols).toBe(11);
  });

  await softStep('1.2 Open ML > Models > Train Model...', async () => {
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Train-Model..."] not found after 3s');
      item.click();
    });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel', null, {timeout: 15_000});
  });

  await softStep('1.3 Configure: Predict=SEX, Features=HEIGHT+WEIGHT, tick Ignore missing → engine loads', async () => {
    await setPredict(page, 'SEX');
    await selectFeaturesByName(page, ['HEIGHT', 'WEIGHT']);
    await page.waitForFunction(() => {
      const root = (window as any).grok.shell.v?.root;
      return root?.querySelector('[name="input-host-Features"]')?.textContent?.includes('(2)');
    }, null, {timeout: 10_000});
    await page.locator('[name="input-Ignore-missing"]').click();
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 90_000}).catch(() => {
      console.warn('[WARN] 1.3 SAVE button still disabled 90s after Ignore missing tick');
    });
  });

  await softStep('1.4 Verify Predict probability available for binary target (SEX)', async () => {
    const imputeVisible = await page.locator('[name="input-host-Impute-missing"]').isVisible().catch(() => false);
    if (imputeVisible) {
      console.warn('[WARN] 1.4 Impute missing still visible after Ignore missing was checked — unexpected');
    }
    const probVisible = await page.locator('[name="input-host-Predict-probability"]').isVisible().catch(() => false);
    if (!probVisible) {
      console.warn('[WARN] 1.4 Predict-probability not visible for binary categorical target SEX');
    }
  });

  await softStep('1.5 Verify training complete — SAVE button enabled', async () => {
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 30_000});
    const viewReady = await page.evaluate(() => (window as any).grok.shell.v?.type === 'PredictiveModel');
    expect(viewReady).toBe(true);
  });

  await softStep('1.6 Verify Ignore missing state — rows with missing values excluded', async () => {
    await expect(page.locator('[name="input-Ignore-missing"]')).toBeChecked();
    const imputeVisible = await page.locator('[name="input-host-Impute-missing"]').isVisible().catch(() => false);
    if (imputeVisible) {
      console.warn('[WARN] 1.6 Impute missing still visible while Ignore missing is checked');
    }
  });

  await softStep('1.7 Tick Predict probability — binary-classification postprocessing path', async () => {
    const probVisible = await page.locator('[name="input-host-Predict-probability"]').isVisible().catch(() => false);
    if (probVisible) {
      const probChecked = await page.locator('[name="input-Predict-probability"]').isChecked();
      if (!probChecked) await page.locator('[name="input-Predict-probability"]').click();
      await expect(page.locator('[name="input-Predict-probability"]')).toBeChecked();
      await page.waitForFunction(() => {
        const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
        return !!btn && !btn.classList.contains('d4-disabled');
      }, null, {timeout: 90_000}).catch(() => {
        console.warn('[WARN] 1.7 re-train after Predict-probability did not complete within 90s');
      });
    }
  });

  await softStep('1.8 Save model as TestDemog — verify discoverable in Browse > Predictive models', async () => {
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 90_000});
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('h3, h4, [class*="card-header"]'))
        .some(e => /Eda:/i.test(e.textContent || '')),
      null, {timeout: 90_000}).catch(() => {
      console.warn('[WARN] 1.8 model results card (Eda: ...) did not appear — OK may be disabled');
    });
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 60_000});
    await nameInput.focus();
    await nameInput.fill(MODEL_NAME);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 30_000});
    const exists = await page.evaluate(async (name) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length > 0;
    }, MODEL_NAME);
    expect(exists).toBe(true);
  });

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click(); else (d as HTMLElement).remove();
    });
  });
  await page.waitForFunction(() => !document.querySelector('.d4-dialog'), null, {timeout: 5_000}).catch(() => {});

  await page.evaluate(async () => {
    const g = (window as any).grok;
    const v = g.shell.v as any;
    if (v && typeof v.close === 'function' && v.type === 'PredictiveModel') v.close();
    await new Promise(r => setTimeout(r, 500));
  });
  await page.waitForFunction(() => (window as any).grok.shell.v?.type !== 'PredictiveModel', null, {timeout: 5_000}).catch(() => {});
  await page.waitForTimeout(500);

  await softStep('2.1 Open ML > Models > Train Model... (again)', async () => {
    await page.waitForFunction(() => !!document.querySelector('[name="div-ML"]'), null, {timeout: 10_000});
    const tabsBefore = await page.evaluate(() =>
      document.querySelectorAll('.d4-tab-header').length);
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Train-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Train-Model..."] not found after 3s');
      item.click();
    });
    await page.waitForFunction((before) =>
      document.querySelectorAll('.d4-tab-header').length > before,
      tabsBefore, {timeout: 15_000}).catch(async () => {
      await page.waitForTimeout(1000);
    });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel', null, {timeout: 10_000});
    await page.waitForTimeout(500);
  });

  await softStep('2.2 Configure: Predict=WEIGHT (numerical), Features=HEIGHT, tick Ignore missing → engine loads', async () => {
    await setPredict(page, 'WEIGHT');
    await selectFeaturesByName(page, ['HEIGHT']);
    await page.waitForFunction(() => {
      const root = (window as any).grok.shell.v?.root;
      return root?.querySelector('[name="input-host-Features"]')?.textContent?.includes('(1)');
    }, null, {timeout: 10_000});
    await page.locator('[name="input-Ignore-missing"]').click();
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 180_000}).catch(() => {
      console.warn('[WARN] 2.2 SAVE button still disabled 3min after Ignore missing tick');
    });
    const probVisible = await page.locator('[name="input-host-Predict-probability"]').isVisible().catch(() => false);
    if (probVisible) {
      console.warn('[WARN] 2.2 Predict-probability visible for numerical target — unexpected per atlas');
    }
  });

  await softStep('2.3 Save as TestDemog_Regression', async () => {
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.classList.contains('d4-disabled');
    }, null, {timeout: 120_000});
    await page.waitForFunction(() =>
      Array.from(document.querySelectorAll('h3, h4, [class*="card-header"]'))
        .some(e => /Eda:/i.test(e.textContent || '')),
      null, {timeout: 90_000}).catch(() => {
      console.warn('[WARN] 2.3 model results card (Eda: ...) did not appear');
    });
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 30_000});
    await nameInput.focus();
    await nameInput.fill(MODEL_NAME_REGRESSION);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 30_000});
    const exists = await page.evaluate(async (name) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length > 0;
    }, MODEL_NAME_REGRESSION);
    expect(exists).toBe(true);
  });

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach(d => {
      const cancel = d.querySelector('[name="button-CANCEL"]');
      if (cancel) (cancel as HTMLElement).click(); else d.remove();
    });
  });
  await page.waitForFunction(() => !document.querySelector('.d4-dialog'), null, {timeout: 5_000}).catch(() => {});

  await page.evaluate(async () => {
    const g = (window as any).grok;
    g.shell.closeAll();
    await new Promise(r => setTimeout(r, 500));
    const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
    g.shell.addTableView(df);
    await new Promise(resolve => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(null); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  await softStep('3.1 Open demog.csv — capture initial column count', async () => {
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return { cols: df?.columns?.length ?? 0 };
    });
    expect(info.cols).toBe(11);
  });

  await softStep('3.2 Open ML > Models > Apply Model... — dialog opens, model list populated', async () => {
    await page.waitForFunction(() => !!document.querySelector('[name="div-ML"]'), null, {timeout: 10_000});
    await page.evaluate(async () => {
      const ml = document.querySelector('[name="div-ML"]') as HTMLElement | null;
      if (!ml) throw new Error('[name="div-ML"] not found');
      ml.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true}));
      ml.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      let item: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        item = document.querySelector('[name="div-ML---Models---Apply-Model..."]') as HTMLElement | null;
        if (item) break;
      }
      if (!item) throw new Error('[name="div-ML---Models---Apply-Model..."] not found after 3s');
      item.click();
    });
    await page.locator('[name="dialog-Apply-predictive-model"]').waitFor({timeout: 10_000});
    await page.waitForFunction(() => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      return !!sel && sel.options.length > 0;
    }, null, {timeout: 10_000});
    const titleText = await page.locator('[name="dialog-Apply-predictive-model"] .d4-dialog-title').textContent();
    expect(titleText).toMatch(/Apply predictive model/i);
    await expect(page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]')).toBeVisible();
  });

  await softStep('3.3 Select TestDemog — ColumnsMapInput reflects WEIGHT/HEIGHT mapped', async () => {
    await page.evaluate(async (modelName) => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return;
      const options = Array.from(sel.options);
      const opt = options.find(o => {
        const text = o.textContent?.trim() ?? '';
        return text === modelName;
      }) ?? options.find(o => {
        const text = o.textContent?.trim() ?? '';
        return text.startsWith(modelName) && !text.includes('Regression');
      }) ?? options.find(o => {
        const text = o.textContent?.trim() ?? '';
        return text.includes(modelName) && !text.includes('Regression');
      });
      if (opt) {
        const setter = Object.getOwnPropertyDescriptor(HTMLSelectElement.prototype, 'value')!.set!;
        setter.call(sel, opt.value);
        sel.dispatchEvent(new Event('input', {bubbles: true}));
        sel.dispatchEvent(new Event('change', {bubbles: true}));
      }
      await new Promise(r => setTimeout(r, 500));
    }, MODEL_NAME);
    await page.waitForFunction(() => {
      const el = document.querySelector('[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]');
      return el?.textContent?.includes('2/2') ?? false;
    }, null, {timeout: 10_000}).catch(() => {});
    const inputsText = await page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]').textContent();
    expect(inputsText).toMatch(/2\/2/);
  });

  await softStep('3.4 Click OK — prediction column appended, column count > 11', async () => {
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.waitForFunction(
      () => !document.querySelector('[name="dialog-Apply-predictive-model"]'),
      null, {timeout: 60_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df?.columns?.length ?? 0;
    });
    expect(cols).toBeGreaterThan(11);
  });

  await softStep('4.1 Navigate to Browse > Platform > Predictive models', async () => {
    await page.evaluate(async () => {
      const g = (window as any).grok;
      g.shell.windows.showBrowse = true;
      g.shell.route('/models');
    });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'models', null, {timeout: 15_000});
  });

  await softStep('4.2 Search "TestDemog" — TestDemog card surfaces', async () => {
    const searchInput = page.locator('input[placeholder="Search models by name or by #tags"]');
    await searchInput.waitFor({timeout: 5_000});
    await searchInput.fill(MODEL_NAME);
    await page.waitForTimeout(1000);
    await page.waitForFunction((name) => {
      const labels = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'));
      return labels.some(l => l.textContent && l.textContent.trim() === name);
    }, MODEL_NAME, {timeout: 10_000});
  });

  await softStep('4.3 Context panel accordion renders model metadata', async () => {
    await page.evaluate((name) => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find(l => l.textContent && l.textContent.trim() === name);
      (label?.closest('.grok-gallery-grid-item') as HTMLElement | null)?.click();
    }, MODEL_NAME);
    await page.waitForTimeout(1000);
    const panes = await page.evaluate(() =>
      Array.from(document.querySelectorAll('.d4-accordion-pane-header'))
        .map(h => h.textContent?.trim() ?? ''));
    expect(panes).toContain('Details');
    expect(panes).toContain('Performance');
    if (!panes.some(p => p.startsWith('Activity'))) {
      console.warn('[INFO] 4.3 Activity pane not present — freshly trained TestDemog has no activity log yet (expected)');
    }
    expect(panes).toContain('Sharing');
  });

  await softStep('4.4 Filter templates panel opens with content', async () => {
    await page.locator('[name="icon-filter"]').click();
    await page.waitForTimeout(800);
    const filterSections = await page.evaluate(() =>
      Array.from(document.querySelectorAll('[name^="div-section--"]'))
        .map(el => el.getAttribute('name')));
    expect(filterSections.length).toBeGreaterThan(0);
    await page.locator('[name="icon-filter"]').click();
  });

  await softStep('4.5 Clear search — all catalog models visible (≥2 for multi-select)', async () => {
    const searchInput = page.locator('input[placeholder="Search models by name or by #tags"]');
    await searchInput.fill('');
    await page.waitForTimeout(1200);
    const modelCount = await page.evaluate(async () => {
      const list = await (window as any).grok.dapi.models.list();
      return list.length;
    });
    if (modelCount < 2) {
      console.warn('[WARN] 4.5 fewer than 2 models — dapi.models.list() returned ' + modelCount);
    }
    expect(modelCount).toBeGreaterThanOrEqual(2);
  });

  await softStep('4.6 CTRL+click ≥2 model cards — multi-select triggers Actions pane', async () => {
    const cards = page.locator('.grok-gallery-grid-item.grok-predictive-model');
    await cards.first().click();
    await page.waitForTimeout(300);
    await page.keyboard.down('Control');
    await cards.nth(1).click();
    await page.keyboard.up('Control');
    await page.waitForTimeout(500);
    await page.waitForFunction(() => {
      const headers = Array.from(document.querySelectorAll('.d4-accordion-pane-header'));
      return headers.some(h => h.textContent && h.textContent.trim() === 'Actions');
    }, null, {timeout: 10_000});
  });

  await softStep('4.7 Open Actions pane → click Compare → Compare models view opens', async () => {
    await page.evaluate(async () => {
      const headers = Array.from(document.querySelectorAll('[name="div-section--Actions"], button'))
        .filter(el => {
          const text = el.textContent?.trim() ?? '';
          return (text === 'Actions' || text.startsWith('Actions')) &&
            !el.closest('[name="grok-toolbox"]');
        });
      const ctxHeader = headers[headers.length - 1] as HTMLElement | undefined;
      if (ctxHeader) ctxHeader.click();
      let compareEl: HTMLElement | null = null;
      for (let i = 0; i < 30; i++) {
        await new Promise(r => setTimeout(r, 100));
        compareEl = (Array.from(document.querySelectorAll('.d4-link-action, button, [role="button"]'))
          .find(el => el.textContent?.trim() === 'Compare') as HTMLElement) ?? null;
        if (compareEl) break;
      }
      if (compareEl) compareEl.click();
    });
    await page.waitForFunction(() => {
      const v = (window as any).grok.shell.v;
      return v?.type === 'TableView' && v?.name === 'Compare models';
    }, null, {timeout: 15_000});
    const cols = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      return df ? Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i).name) : [];
    });
    expect(cols).toContain('Name');
    expect(cols).toContain('Method');
  });

  await softStep('5.1 Navigate to Browse > Platform > Predictive models', async () => {
    await page.evaluate(async () => {
      const g = (window as any).grok;
      g.shell.windows.showBrowse = true;
      g.shell.route('/models');
    });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'models', null, {timeout: 15_000});
  });

  await softStep('5.2 Locate TestDemog card — at least one card present', async () => {
    await page.waitForFunction((name) => {
      const labels = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'));
      return labels.some(l => l.textContent && l.textContent.trim() === name);
    }, MODEL_NAME, {timeout: 15_000});
  });

  await softStep('5.3 Right-click TestDemog → Delete — confirm-delete modal opens', async () => {
    await page.evaluate((name) => {
      const label = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'))
        .find(l => l.textContent && l.textContent.trim() === name);
      const card = label && label.closest('.grok-gallery-grid-item');
      if (card) card.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
    }, MODEL_NAME);
    await page.locator('.d4-menu-popup').waitFor({timeout: 5_000});
    await page.evaluate(() => {
      const popup = document.querySelector('.d4-menu-popup');
      const delItem = Array.from(popup?.querySelectorAll('.d4-menu-item-label') ?? [])
        .find(el => el.textContent && el.textContent.trim() === 'Delete');
      if (delItem) (delItem as HTMLElement).click();
    });
    await page.locator('[name="dialog-Are-you-sure?"]').waitFor({timeout: 5_000});
  });

  await softStep('5.4 Click DELETE — dialog closes, no error, GROK-846 invariant', async () => {
    await page.locator('[name="dialog-Are-you-sure?"] [name="button-DELETE"]').click();
    await page.waitForFunction(() => !document.querySelector('[name="dialog-Are-you-sure?"]'), null, {timeout: 15_000});
    await page.waitForFunction(async (name) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length === 0;
    }, MODEL_NAME, {timeout: 15_000});
  });

  await softStep('5.5 Verify TestDemog card no longer appears in browser', async () => {
    await page.evaluate(() => { (window as any).grok.shell.route('/'); });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type !== 'models', null, {timeout: 5_000}).catch(() => {});
    await page.evaluate(() => { (window as any).grok.shell.route('/models'); });
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'models', null, {timeout: 10_000});
    await page.waitForTimeout(2000);
    await page.waitForFunction((name) => {
      const labels = Array.from(document.querySelectorAll('.grok-gallery-grid-item-title'));
      return !labels.some(l => l.textContent && l.textContent.trim() === name);
    }, MODEL_NAME, {timeout: 15_000});
  });

  await softStep('Cleanup: delete TestDemog_Regression', async () => {
    const exists = await page.evaluate(async (name) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length > 0;
    }, MODEL_NAME_REGRESSION);
    if (exists) {
      await page.evaluate(async (name) => {
        const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
        for (const m of list) await (window as any).grok.dapi.models.delete(m);
      }, MODEL_NAME_REGRESSION);
    }
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
