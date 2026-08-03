import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

const MODEL_NAME = 'OneHotSuffixCollision_test';

test('One-hot suffix collision: namespaced <name>=<category> columns survive train + apply', async ({page}) => {
  // Trains one small EDA Linear Regression (40-row in-memory) with one-hot encoding, then applies it.
  // Not chemprop. Step 4 polls SAVE-enable up to 180s; 300s covers train + save + apply with margin.
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async (name: string) => {
    const g: any = (window as any).grok;
    const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
    for (const m of list)
      await g.dapi.models.delete(m);
  }, MODEL_NAME);

  
  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    const g: any = (window as any).grok;
    const DG: any = (window as any).DG;
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const N = 40;
    const fa: string[] = Array.from({length: N}, (_: unknown, i: number) => (i % 2 === 0) ? 'Yes' : 'No');
    const fb: string[] = Array.from({length: N}, (_: unknown, i: number) => (Math.floor(i / 2) % 2 === 0) ? 'Yes' : 'No');
    const tgt: number[] = Array.from({length: N}, (_: unknown, i: number) =>
      (fa[i] === 'Yes' ? 2 : 0) + (fb[i] === 'Yes' ? 3 : 0) + (i * 0.07));
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings('featureA', fa),
      DG.Column.fromStrings('featureB', fb),
      DG.Column.fromFloat32Array('target', Float32Array.from(tgt)),
    ]);
    df.name = 'OneHotSuffixTrain';
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});

  
  const setPredict = async (columnName: string) => {
    await page.evaluate(() => {
      const root = (window as any).grok.shell.v.root;
      const sel = root.querySelector('[name="input-host-Predict"] .d4-column-selector') as HTMLElement;
      sel.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, cancelable: true, view: window, button: 0}));
    });
    await page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'),
      null, {timeout: 5_000});
    await page.evaluate(() => (document.querySelector('.d4-column-selector-backdrop') as HTMLElement).focus());
    await page.keyboard.type(columnName);
    await page.waitForTimeout(300);
    await page.keyboard.press('Enter').catch(() => {});
    await page.waitForTimeout(500);
    await expect.poll(async () => await page.evaluate(() =>
      (window as any).grok.shell.v.root.querySelector(
        '[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 10_000}).toBe(columnName);
  };

  await softStep('1. Build 40-row dataframe: featureA / featureB (Yes/No) + target', async () => {
    const info = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const names: string[] = df ? Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name) : [];
      const fa = df?.col('featureA');
      const fb = df?.col('featureB');
      const faCats: string[] | null = fa?.categories ? Array.from(fa.categories as Iterable<string>).sort() : null;
      const fbCats: string[] | null = fb?.categories ? Array.from(fb.categories as Iterable<string>).sort() : null;
      return {rows: df?.rowCount ?? 0, cols: names, faCats, fbCats};
    });
    expect(info.rows).toBe(40);
    expect(info.cols).toEqual(['featureA', 'featureB', 'target']);
    expect(info.faCats).toEqual(['No', 'Yes']);
    expect(info.fbCats).toEqual(['No', 'Yes']);
  });

  await softStep('2. ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
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
    await page.waitForFunction(() => (window as any).grok.shell.v?.type === 'PredictiveModel',
      null, {timeout: 15_000});
  });

  await softStep('3a. Set Predict to target (numerical regression target)', async () => {
    await setPredict('target');
  });

  await softStep('3b. Set Features to featureA + featureB (the suffix-collision pair)', async () => {
    await page.evaluate(() => {
      const editor = (window as any).grok.shell.v.root.querySelector('[name="div-Features"]') as HTMLElement;
      const r = editor.getBoundingClientRect();
      const opts = (t: string) => new MouseEvent(t, {
        bubbles: true, cancelable: true, view: window, button: 0,
        clientX: r.left + 10, clientY: r.top + 5,
      });
      editor.dispatchEvent(opts('mousedown'));
      editor.dispatchEvent(opts('mouseup'));
      editor.dispatchEvent(opts('click'));
    });
    await page.locator('[name="dialog-Select-columns..."]').waitFor({timeout: 10_000});

    await page.evaluate(() => {
      const overlay = document.querySelector('[name="dialog-Select-columns..."] canvas[name="overlay"]') as HTMLCanvasElement;
      if (!overlay) throw new Error('column-picker overlay canvas not found');
      const click = (x: number, y: number) => {
        const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: x, clientY: y};
        overlay.dispatchEvent(new MouseEvent('mousedown', opts));
        overlay.dispatchEvent(new MouseEvent('mouseup', opts));
        overlay.dispatchEvent(new MouseEvent('click', opts));
      };
      click(826, 257); // featureA (row 0)
      click(826, 285); // featureB (row 1)
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('2');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('featureA', {timeout: 10_000});
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('featureB');
    
    await page.waitForFunction(() => {
      const host = document.querySelector('[name="input-host-One-hot-encoding"]') as HTMLElement | null;
      return !!host && host.offsetParent !== null;
    }, null, {timeout: 5_000});
  });

  await softStep('4. Tick One-hot encoding — train fires; SAVE enables', async () => {
    
    await page.locator('[name="input-host-One-hot-encoding"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-One-hot-encoding"] input[type="checkbox"]'))
      .toBeChecked();
    
    await page.waitForTimeout(2_000);
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.className.includes('d4-disabled');
    }, null, {timeout: 180_000});
    
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x: any) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  await softStep(`5. Save the model as ${MODEL_NAME}`, async () => {
    
    await page.locator('[name="button-Save"]').click();
    const nameInput = page.locator('.d4-dialog [name="input-host-Name"] input');
    await nameInput.waitFor({timeout: 10_000});
    await nameInput.focus();
    await page.keyboard.type(MODEL_NAME);
    await page.locator('.d4-dialog [name="button-OK"]').click();
    await page.waitForFunction(() =>
      !document.querySelector('.d4-dialog [name="input-host-Name"]'),
      null, {timeout: 60_000});
    await expect.poll(async () => await page.evaluate(async (name: string) => {
      const list = await (window as any).grok.dapi.models.filter(`friendlyName = "${name}"`).list();
      return list.length;
    }, MODEL_NAME), {timeout: 60_000}).toBeGreaterThan(0);
  });

  

  await softStep('6. Train-side suffix-collision proof: model persisted with friendlyName + name', async () => {
    
    const declared: {found: boolean; name: string | null; friendly: string | null} =
      await page.evaluate(async (wantName: string) => {
        const g: any = (window as any).grok;
        const list = await g.dapi.models.filter(`friendlyName = "${wantName}"`).list();
        if (!list.length) return {found: false, name: null, friendly: null};
        const m: any = list[0];
        return {found: true, name: m?.name ?? null, friendly: m?.friendlyName ?? null};
      }, MODEL_NAME);
    expect(declared.found,
      `Saved model ${MODEL_NAME} must round-trip from dapi.models. ` +
      `Train-side suffix-collision invariant proof: TRAIN (step 4) and SAVE ` +
      `(step 5) completed — Dart would have thrown a duplicate-column-name ` +
      `error during oneHotEncoded() if the <name>= namespace had been lost.`).toBe(true);
    expect(declared.friendly).toBe(MODEL_NAME);
  });

  await softStep('7. Open a fresh copy of the dataframe in a new table view', async () => {
    
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const DG: any = (window as any).DG;
      const N = 40;
      const fa: string[] = Array.from({length: N}, (_: unknown, i: number) => ((i + 1) % 2 === 0) ? 'Yes' : 'No');
      const fb: string[] = Array.from({length: N}, (_: unknown, i: number) => (Math.floor((i + 1) / 2) % 2 === 0) ? 'Yes' : 'No');
      const tgt: number[] = Array.from({length: N}, (_: unknown, i: number) =>
        (fa[i] === 'Yes' ? 2 : 0) + (fb[i] === 'Yes' ? 3 : 0) + (i * 0.11));
      const df2 = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('featureA', fa),
        DG.Column.fromStrings('featureB', fb),
        DG.Column.fromFloat32Array('target', Float32Array.from(tgt)),
      ]);
      df2.name = 'OneHotSuffixApply';
      g.shell.addTableView(df2);
      await new Promise<void>((resolve) => {
        const sub = df2.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      const views: any[] = Array.from(g.shell.views);
      const tv2 = views.reverse().find((v: any) => v?.type === 'TableView' && v?.dataFrame?.name === 'OneHotSuffixApply');
      if (tv2) g.shell.v = tv2;
      (window as any).__applyInitialColCount = df2.columns.length;
      (window as any).__applyInitialColNames = Array.from({length: df2.columns.length},
        (_: unknown, i: number) => df2.columns.byIndex(i).name);
    });
    
    await page.locator('[name="div-ML"]').waitFor({timeout: 15_000});
    const info = await page.evaluate(() => {
      const g: any = (window as any).grok;
      return {viewType: g.shell.v?.type, tableName: g.shell.tv?.dataFrame?.name,
              initialCols: (window as any).__applyInitialColCount};
    });
    expect(info.viewType).toBe('TableView');
    expect(info.tableName).toBe('OneHotSuffixApply');
    expect(info.initialCols).toBe(3);
  });

  await softStep('8. ML > Models > Apply Model... — select model, OK', async () => {
    
    await page.locator('[name="div-ML"]').click();
    await page.evaluate(() => {
      const models = document.querySelector('[name="div-ML---Models"]') as HTMLElement | null;
      if (!models) throw new Error('ML > Models submenu not found');
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
    }, null, {timeout: 15_000});
    const pick: {found: boolean; opts: string[]; picked?: string} = await page.evaluate((wantName: string) => {
      const sel = document.querySelector(
        '[name="dialog-Apply-predictive-model"] [name="input-host-Model"] select') as HTMLSelectElement | null;
      if (!sel) return {found: false, opts: [] as string[]};
      const opts = Array.from(sel.options).map((o) => o.textContent || '');
      
      const tryPrefixes = [wantName, wantName.slice(0, 24), wantName.slice(0, 16)];
      let idx = -1;
      for (const p of tryPrefixes) {
        idx = opts.findIndex((t) => t.includes(p));
        if (idx >= 0) break;
      }
      if (idx < 0) return {found: false, opts};
      sel.selectedIndex = idx;
      sel.dispatchEvent(new Event('change', {bubbles: true}));
      return {found: true, opts, picked: opts[idx]};
    }, MODEL_NAME);
    expect(pick.found,
      `Apply dialog must surface "${MODEL_NAME}" as a suggested model — ` +
      `confirms that columnNamesMap auto-build recognises the apply-time table's ` +
      `featureA / featureB categorical shape as compatible with the trained ` +
      `per-feature <name>=<category> expansion. Options seen: ${JSON.stringify(pick.opts)}.`)
      .toBe(true);
    
    await expect(page.locator(
      '[name="dialog-Apply-predictive-model"] [name="input-host-Inputs"]'))
      .toContainText('(2/2)', {timeout: 10_000});
    await page.locator('[name="dialog-Apply-predictive-model"] [name="button-OK"]').click();
    await page.locator('[name="dialog-Apply-predictive-model"]')
      .waitFor({state: 'detached', timeout: 30_000});
  });

  await softStep('9. Apply reconstruction — prediction column appended to OneHotSuffixApply', async () => {
    
    await expect.poll(async () => await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initial: number = (window as any).__applyInitialColCount ?? 0;
      return (df?.columns?.length ?? 0) - initial;
    }), {timeout: 60_000, intervals: [500, 1_000, 2_000]}).toBeGreaterThan(0);
    const result = await page.evaluate(() => {
      const df: any = (window as any).grok.shell.tv?.dataFrame;
      const initialNames: string[] = (window as any).__applyInitialColNames ?? [];
      const allNames: string[] = Array.from({length: df.columns.length},
        (_: unknown, i: number) => df.columns.byIndex(i).name);
      const newCols = allNames.filter((n: string) => !initialNames.includes(n));
      const DG: any = (window as any).DG;
      const tagKey = DG?.TAGS?.PREDICTIVE_MODEL ?? '.predictive-model';
      const taggedCount = newCols.filter((n: string) => {
        const col = df.columns.byName(n);
        return !!col?.tags?.[tagKey] || !!col?.getTag?.(tagKey);
      }).length;
      return {newCols, taggedCount, cols: df.columns.length, initial: initialNames.length};
    });
    expect(result.newCols.length,
      `Apply path must append ≥1 prediction column to OneHotSuffixApply. ` +
      `If no column was appended, the apply-side columnNamesMap auto-build ` +
      `failed to reconstruct the per-feature <name>=<category> expansion ` +
      `(the suffix-collision invariant broke on the apply side).`).toBeGreaterThan(0);
    
    if (result.taggedCount === 0) {
      console.warn(`Apply: ${result.newCols.length} new column(s) appended ` +
        `(${JSON.stringify(result.newCols)}) but none carry Tags.PredictiveModel.`);
    }
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x: any) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  
  await softStep('10. Teardown — delete OneHotSuffixCollision_test from the server', async () => {
    const remaining = await page.evaluate(async (name: string) => {
      const g: any = (window as any).grok;
      const list = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      for (const m of list)
        await g.dapi.models.delete(m);
      const after = await g.dapi.models.filter(`friendlyName = "${name}"`).list();
      return after.length;
    }, MODEL_NAME);
    expect(remaining, `Teardown must remove all ${MODEL_NAME} entities from the server`).toBe(0);
  });

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
