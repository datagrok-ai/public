/* ---
sub_features_covered: [models.validators.class-imbalance, models.validators.highly-correlated, models.validators.string-features, models.validators.too-many-unique-categories]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

async function openTrainView(page: Page) {
  await page.locator('[name="div-ML"]').click();
  await page.evaluate(() => {
    const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
    if (!models) throw new Error('ML > Models submenu not found');
    const r = models.getBoundingClientRect();
    const ev = (type: string) => new MouseEvent(type, {
      bubbles: true, cancelable: true, view: window,
      clientX: r.left + 5, clientY: r.top + 5,
    });
    models.dispatchEvent(ev('mouseover'));
    models.dispatchEvent(ev('mouseenter'));
    models.dispatchEvent(ev('mousemove'));
  });
  await page.locator('[name="div-ML---Models---Train-Model..."]').click();
  await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel',
    null, {timeout: 15_000});
  await page.waitForFunction(() => {
    const root = (window as any).grok.shell.v?.root;
    if (!root) return false;
    const predictHost = root.querySelector('[name="input-host-Predict"] .d4-column-selector');
    const tipsWidget = Array.from(root.querySelectorAll('.d4-pm-model-widget'))
      .find((w: Element) => w.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Insights & Tips');
    return !!predictHost && !!tipsWidget;
  }, null, {timeout: 15_000});
  await page.waitForTimeout(500);
}

async function setPredict(page: Page, columnName: string) {
  const current = await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim());
  if (current === columnName) return;
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
    sel.dispatchEvent(new MouseEvent('mousedown', {
      bubbles: true, cancelable: true, view: window, button: 0,
    }));
  });
  await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
    null, {timeout: 5_000});
  await page.evaluate(() =>
    (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
  await page.keyboard.type(columnName);
  await page.waitForTimeout(400);
  await page.keyboard.press('Enter').catch(() => {});
  await page.waitForTimeout(700);
  const stillOpen = await page.evaluate(() => !!document.querySelector('.d4-column-selector-backdrop'));
  if (stillOpen) {
    await page.keyboard.press('Escape').catch(() => {});
    await page.waitForTimeout(300);
  }
  await expect.poll(async () => await page.evaluate(() =>
    (window as any).grok.shell.v.root
      .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
    {timeout: 10_000}).toBe(columnName);
}

async function selectFeatures(page: Page, columnNames: string[]) {
  await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const editor = root.querySelector('[name="div-Features"]') as HTMLElement;
    const r = editor.getBoundingClientRect();
    const opts = (t: string) => new MouseEvent(t, {
      bubbles: true, cancelable: true, view: window, button: 0,
      clientX: r.left + 10, clientY: r.top + 5,
    });
    editor.dispatchEvent(opts('mousedown'));
    editor.dispatchEvent(opts('mouseup'));
    editor.dispatchEvent(opts('click'));
  });
  const dlg = page.locator('[name="dialog-Select-columns..."]');
  await dlg.waitFor({timeout: 10_000});
  await dlg.locator('[name="label-None"]').click();
  const readChecked = async () => {
    return await dlg.evaluate((d) => {
      const lbl = Array.from(d.querySelectorAll('label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      const m = lbl?.textContent?.match(/^(\d+)/);
      return m ? parseInt(m[1], 10) : 0;
    });
  };
  await expect.poll(readChecked, {timeout: 5_000, intervals: [100, 200, 300],
    message: 'Select-columns dialog never reached 0 checked after label-None'}).toBe(0);
  const search = dlg.locator('input[placeholder="Search"]');
  for (const name of columnNames) {
    const before = await readChecked();
    await search.fill('');
    await page.waitForTimeout(200);
    await search.fill(name);
    await page.waitForTimeout(500);
    const overlay = dlg.locator('[name="viewer-Grid"] [name="overlay"]');
    const box = await overlay.boundingBox();
    if (!box) throw new Error('column-picker overlay canvas not measurable');
    let toggled = false;
    for (const [dx, dy] of [[40, 28], [40, 34], [42, 28], [35, 28]]) {
      const clientX = box.x + box.width - dx;
      const clientY = box.y + dy;
      const synthetic = await dlg.evaluate((d, {cx, cy}) => {
        const ov = d.querySelector('[name="viewer-Grid"] [name="overlay"]') as HTMLElement;
        if (!ov) return 'no-overlay';
        const opts = (t: string) => new MouseEvent(t, {
          bubbles: true, cancelable: true, view: window,
          button: 0, clientX: cx, clientY: cy,
        });
        ov.dispatchEvent(opts('mousedown'));
        ov.dispatchEvent(opts('mouseup'));
        ov.dispatchEvent(opts('click'));
        return 'dispatched';
      }, {cx: clientX, cy: clientY});
      if (synthetic !== 'dispatched') break;
      await page.waitForTimeout(300);
      if ((await readChecked()) === before + 1) {
        toggled = true;
        break;
      }
    }
    if (!toggled) {
      outer: for (const dx of [40, 35, 45, 30, 50]) {
        for (const dy of [28, 32, 24, 36, 20, 40]) {
          await page.mouse.click(box.x + box.width - dx, box.y + dy);
          await page.waitForTimeout(300);
          if ((await readChecked()) === before + 1) {
            toggled = true;
            break outer;
          }
        }
      }
    }
    if (!toggled)
      throw new Error(`Failed to toggle "${name}" via synthetic dispatchEvent + page.mouse.click fallback on overlay (search-filtered to "${name}")`);
  }
  await search.fill('');
  await page.waitForTimeout(200);
  const finalChecked = await readChecked();
  if (finalChecked !== columnNames.length) {
    throw new Error(
      `Expected ${columnNames.length} checked after toggles, got ${finalChecked}`);
  }
  await dlg.locator('[name="button-OK"]').click();
  await page.waitForFunction(() => !document.querySelector('[name="dialog-Select-columns..."]'),
    null, {timeout: 10_000});
  await expect.poll(async () => await page.evaluate(() => {
    const inp = document.querySelector('[name="input-host-Features"] .ui-input-column-names');
    return inp?.className.includes('d4-invalid') ? 'invalid' : 'valid';
  }), {timeout: 5_000, message: 'Features input should not be d4-invalid after selection'})
    .toBe('valid');
  await page.waitForFunction(() => {
    const host = document.querySelector('[name="input-host-Impute-missing"]') as HTMLElement | null;
    return !!host && host.offsetParent !== null;
  }, null, {timeout: 10_000}).catch(() => {});
  await page.waitForTimeout(1_000);
}

async function getTipsText(page: Page): Promise<string> {
  return await page.evaluate(() => {
    const root = (window as any).grok.shell.v.root;
    const widgets = Array.from(root.querySelectorAll('.d4-pm-model-widget')) as HTMLElement[];
    const tips = widgets.find((w) =>
      w.querySelector('.d4-dialog-title')?.textContent?.trim() === 'Insights & Tips');
    return tips ? (tips as HTMLElement).innerText : '';
  });
}

async function waitForTipsContaining(page: Page, phrase: string, timeout = 30_000) {
  await expect.poll(async () => (await getTipsText(page)).toLowerCase(),
    {timeout, intervals: [500, 1000, 1500], message:
      `waiting for "Insights & Tips" widget to contain "${phrase}"`})
    .toContain(phrase.toLowerCase());
}

async function assertNonBlocking(page: Page) {
  await expect(page.locator('[name="button-Save"]')).toBeAttached({timeout: 5_000});
  const errorBalloons = await page.locator('.d4-balloon-error').count();
  expect(errorBalloons).toBe(0);
  const errCount = await page.evaluate(() => {
    const w: any[] = (window as any).grok.shell.warnings ?? [];
    return w.filter((x) =>
      x && typeof x === 'object' && (x.isError === true || x.level === 'error')).length;
  });
  expect(errCount).toBe(0);
}

async function openInMemoryDataFrame(page: Page, builder: () => unknown) {
  await page.evaluate(async (builderSrc: string) => {
    const g: any = (window as any).grok;
    g.shell.closeAll();
    await new Promise<void>((resolve) => {
      const start = Date.now();
      const poll = () => {
        if (!g.shell.v || g.shell.v.type === undefined || Date.now() - start > 3000) resolve();
        else setTimeout(poll, 100);
      };
      poll();
    });
    await new Promise((r) => setTimeout(r, 400));
    const df = (new Function('DG', `return (${builderSrc})(DG)`))((window as any).DG);
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  }, builder.toString());
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
  await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'TableView',
    null, {timeout: 10_000});
}

test('Models validators — pre-train edge surfaces (class-imbalance, string-features, too-many-unique, highly-correlated)', async ({page}) => {
  // Pure pre-train validator spec: NO model is trained (each block stops at the Insights & Tips
  // widget). Four blocks, each opening the train view + polling a tip up to 30s. 240s is ample.
  test.setTimeout(240_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
  });

  await softStep('1.1 Build 200-row class-imbalance dataframe (target = 180×A / 20×B)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 200;
      const f1 = new Float32Array(n);
      const f2 = new Float32Array(n);
      const tgt: string[] = [];
      for (let i = 0; i < n; i++) {
        f1[i] = ((i * 9301 + 49297) % 233280) / 233280;
        f2[i] = ((i * 1597 + 12345) % 65536) / 65536;
        tgt.push(i < 180 ? 'A' : 'B');
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('feature1', f1),
        DG.Column.fromFloat32Array('feature2', f2),
        DG.Column.fromStrings('target', tgt),
      ]);
    });
  });

  await softStep('1.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('1.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('1.4 Set Features = (feature1, feature2) — exclude target to avoid predict-conflict d4-invalid', async () => {
    await selectFeatures(page, ['feature1', 'feature2']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature1', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature2');
  });

  await softStep('1.5 Insights & Tips surfaces class-imbalance warning citing target', async () => {
    await waitForTipsContaining(page, 'class imbalance');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('target');
    expect(tipsText).toMatch(/A \(\d+\.\d+\)|B \(\d+\.\d+\)/);
  });

  await softStep('1.6 Validator is non-blocking (no error-balloon, no error warnings, SAVE button mounted)', async () => {
    await assertNonBlocking(page);
  });

  await softStep('2.1 Build 60-row string-features dataframe (cat_feature = red/green/blue)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 60;
      const num = new Float32Array(n);
      const tgt = new Float32Array(n);
      const cats: string[] = [];
      for (let i = 0; i < n; i++) {
        num[i] = ((i * 9301 + 49297) % 233280) / 233280;
        tgt[i] = (((i * 1597 + 12345) % 65536) / 65536) * 10;
        cats.push(['red', 'green', 'blue'][i % 3]);
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromStrings('cat_feature', cats),
        DG.Column.fromFloat32Array('num_feature', num),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('2.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('2.3 Set Predict = target (numerical regression target)', async () => {
    await setPredict(page, 'target');
  });

  await softStep('2.4 Set Features = (cat_feature, num_feature) — exclude target', async () => {
    await selectFeatures(page, ['cat_feature', 'num_feature']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('cat_feature', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('num_feature');
  });

  await softStep('2.5 Insights & Tips surfaces string-features hint citing cat_feature', async () => {
    await waitForTipsContaining(page, 'categorical');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('cat_feature');
    expect(tipsText.toLowerCase()).toContain('convert');
    expect(tipsText.toLowerCase()).toContain('numerical');
  });

  await softStep('2.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  await softStep('3.1 Build 50-row dataframe with id_like StringColumn (categories/length = 1.0)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 50;
      const ids: string[] = [];
      const f1 = new Float32Array(n);
      const tgt = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        ids.push(`row_${i}`);
        f1[i] = ((i * 9301 + 49297) % 233280) / 233280;
        tgt[i] = (((i * 1597 + 12345) % 65536) / 65536) * 10;
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromStrings('id_like', ids),
        DG.Column.fromFloat32Array('feature1', f1),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('3.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('3.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('3.4 Set Features = (id_like, feature1) — exclude target', async () => {
    await selectFeatures(page, ['id_like', 'feature1']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('id_like', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feature1');
  });

  await softStep('3.5 Insights & Tips surfaces too-many-unique-categories warning citing id_like', async () => {
    await waitForTipsContaining(page, 'too many unique categories');
    const tipsText = await getTipsText(page);
    expect(tipsText).toContain('id_like');
  });

  await softStep('3.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  await softStep('4.1 Build 30-row dataframe with feat_a + feat_b (Pearson > 0.9)', async () => {
    await openInMemoryDataFrame(page, (DG: any) => {
      const n = 30;
      const a = new Float32Array(n);
      const b = new Float32Array(n);
      const tgt = new Float32Array(n);
      for (let i = 0; i < n; i++) {
        const base = i;
        a[i] = base;
        b[i] = base + (((i * 7919) % 100) / 100) * 0.5;
        tgt[i] = base * 0.3 + (((i * 1597 + 12345) % 17) / 17) * 0.1;
      }
      return DG.DataFrame.fromColumns([
        DG.Column.fromFloat32Array('feat_a', a),
        DG.Column.fromFloat32Array('feat_b', b),
        DG.Column.fromFloat32Array('target', tgt),
      ]);
    });
  });

  await softStep('4.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView(page);
  });

  await softStep('4.3 Set Predict = target', async () => {
    await setPredict(page, 'target');
  });

  await softStep('4.4 Set Features = (feat_a, feat_b) — exclude target', async () => {
    await selectFeatures(page, ['feat_a', 'feat_b']);
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feat_a', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('feat_b');
  });

  await softStep('4.5 Insights & Tips surfaces highly-correlated warning citing feat_a <> feat_b pair', async () => {
    await waitForTipsContaining(page, 'highly correlated');
    const tipsText = await getTipsText(page);
    expect(tipsText).toMatch(/feat_a\s*<>\s*feat_b|feat_b\s*<>\s*feat_a/);
  });

  await softStep('4.6 Validator is non-blocking', async () => {
    await assertNonBlocking(page);
  });

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
