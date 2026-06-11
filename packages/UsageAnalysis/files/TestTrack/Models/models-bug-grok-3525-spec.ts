import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('GROK-3525 regression: target nulls surface in validator tip + Ignore missing recovers', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  const openTrainView = async () => {
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
  };

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
      (window as any).grok.shell.v.root
        .querySelector('[name="input-host-Predict"] .d4-column-selector-column')?.textContent?.trim()),
      {timeout: 10_000}).toBe(columnName);
  };

  const readTipText = async () => await page.evaluate(() => {
    const w = document.querySelector('.d4-pm-model-widget') as HTMLElement | null;
    if (!w || w.offsetParent === null) return '';
    return Array.from(w.querySelectorAll('li')).map((li) => li.textContent || '').join(' || ');
  });

  
  await softStep('1.1 Open demog.csv and inject a single null into RACE row 0', async () => {
    const setupInfo = await page.evaluate(async () => {
      const g: any = (window as any).grok;
      document.body.classList.add('selenium');
      g.shell.settings.showFiltersIconsConstantly = true;
      g.shell.windows.simpleMode = true;
      g.shell.closeAll();
      const df = await g.dapi.files.readCsv('System:DemoFiles/demog.csv');
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      
      const race = df.columns.byName('RACE');
      race.set(0, null);
      return {rows: df.rowCount, raceNullsAfter: race.stats.missingValueCount};
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    expect(setupInfo.rows).toBeGreaterThan(0);
    expect(setupInfo.raceNullsAfter).toBeGreaterThan(0);
  });

  await softStep('1.2 ML > Models > Train Model... opens PredictiveModelingView', async () => {
    await openTrainView();
  });

  await softStep('1.3 Set Predict to RACE (categorical target carrying nulls)', async () => {
    await setPredict('RACE');
  });

  await softStep('1.4 Set Features to AGE, WEIGHT via canvas picker', async () => {
    
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
      click(826, 285); // AGE
      click(826, 397); // WEIGHT
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('2');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('AGE', {timeout: 10_000});
  });

  await softStep('1.5 Validator tip names RACE as a missing-values column (GROK-3525 invariant)', async () => {
    
    await expect.poll(readTipText, {timeout: 15_000})
      .toMatch(/contain.*missing values/i);
    const tip = await readTipText();
    
    expect(tip).toContain('RACE');
    expect(tip).toMatch(/RACE:\s*\d+/);
  });

  await softStep('1.6 SAVE/TRAIN button stays disabled — no model save side-effect fires', async () => {
    
    await page.waitForTimeout(5_000);
    const btnCls = await page.evaluate(() => {
      const b = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return b?.className || '';
    });
    expect(btnCls).toContain('d4-disabled');
    const savedCount = await page.evaluate(async () => {
      const list = await (window as any).grok.dapi.models.filter('friendlyName like "BugGrok3525%"').list();
      return list.length;
    });
    expect(savedCount).toBe(0);
  });
  
  await softStep('2.1 Build small in-memory df (X num feature, Y num target with 5 nulls)', async () => {
    await page.evaluate(async () => {
      const g: any = (window as any).grok;
      const DG: any = (window as any).DG;
      g.shell.closeAll();
      // 30 rows: X = 1..30; Y = X*2 except indices 5,10,15,20,25 are null.
      const xs: number[] = [];
      const ys: (number | null)[] = [];
      for (let i = 1; i <= 30; i++) {
        xs.push(i);
        ys.push([5, 10, 15, 20, 25].includes(i) ? null : i * 2);
      }
      const colX = DG.Column.fromList('double', 'X', xs);
      const colY = DG.Column.fromList('double', 'Y', ys);
      const df = DG.DataFrame.fromColumns([colX, colY]);
      df.name = 'BugGrok3525';
      g.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 2000);
      });
    });
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30_000});
    const info = await page.evaluate(() => {
      const df = (window as any).grok.shell.tv?.dataFrame;
      const y = df?.columns?.byName('Y');
      return {rows: df?.rowCount, yNulls: y?.stats?.missingValueCount, yType: y?.type};
    });
    expect(info.rows).toBe(30);
    expect(info.yNulls).toBe(5);
    expect(info.yType).toBe('double');
  });

  await softStep('2.2 ML > Models > Train Model... re-opens train view on BugGrok3525', async () => {
    await openTrainView();
  });

  await softStep('2.3 Set Predict to Y (numerical target carrying nulls)', async () => {
    await setPredict('Y');
  });

  await softStep('2.4 Set Features to X (only column left after Y is excluded)', async () => {
    // Picker has 1 visible row (X) once Y is the Predict column.
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
      const opts = {bubbles: true, cancelable: true, view: window, button: 0, clientX: 800, clientY: 263};
      overlay.dispatchEvent(new MouseEvent('mousedown', opts));
      overlay.dispatchEvent(new MouseEvent('mouseup', opts));
      overlay.dispatchEvent(new MouseEvent('click', opts));
    });
    await page.waitForFunction(() => {
      const lbl = Array.from(document.querySelectorAll('[name="dialog-Select-columns..."] label'))
        .find((l) => /\d+ checked/.test(l.textContent || ''));
      return lbl?.textContent?.startsWith('1');
    }, null, {timeout: 5_000});
    await page.locator('[name="dialog-Select-columns..."] [name="button-OK"]').click();
    await expect(page.locator('[name="input-host-Features"] .ui-input-column-names'))
      .toContainText('X', {timeout: 10_000});
  });

  await softStep('2.5 Sanity check (mirrors Block 1): tip names Y as a missing-values column', async () => {
    
    await expect.poll(readTipText, {timeout: 15_000})
      .toMatch(/contain.*missing values/i);
    const tip = await readTipText();
    expect(tip).toContain('Y');
    
    expect(tip).toMatch(/Column[s]?\s+['"]?Y['"]?.*missing values/i);
  });

  await softStep('2.6 Tick Ignore missing — re-train fires, SAVE enables', async () => {
    
    await page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]').click();
    await expect(page.locator('[name="input-host-Ignore-missing"] input[type="checkbox"]'))
      .toBeChecked();
    await page.waitForFunction(() => {
      const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
      return !!btn && !btn.className.includes('d4-disabled');
    }, null, {timeout: 180_000});
    const errCount = await page.evaluate(() => {
      const w: any[] = (window as any).grok.shell.warnings ?? [];
      return w.filter((x) => /error|fail/i.test(JSON.stringify(x))).length;
    });
    expect(errCount).toBe(0);
  });

  await softStep('2.7 Post-ignore tip no longer cites Y as a missing-values column', async () => {
    
    const tip = await readTipText();
    
    const missingClauseMatch = tip.match(/Column[s]? '[^']+' contain[s]? missing values/);
    if (missingClauseMatch) {
      expect(missingClauseMatch[0]).not.toContain('Y');
    }
  });

  if (stepErrors.length > 0) {
    throw new Error(`${stepErrors.length} step(s) failed:\n` +
      stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n'));
  }
});
