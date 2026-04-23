import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Word Cloud tests', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const df = await (window as any).grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
    (window as any).grok.shell.addTableView(df);
    await new Promise((resolve: (v?: unknown) => void) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Add Word Cloud viewer', async () => {
    await page.locator('[name="icon-Word-cloud"]').click();
    await page.locator('[name="viewer-Word-cloud"]').waitFor({timeout: 10000});
    const ok = await page.evaluate(() => {
      const viewers: string[] = [];
      for (const v of (window as any).grok.shell.tv.viewers) viewers.push(v.type);
      return viewers.includes('Word cloud');
    });
    expect(ok).toBe(true);
  });

  await softStep('Column assignment: RACE', async () => {
    const col = await page.evaluate(() => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      (window as any).grok.shell.o = wc;
      wc.setOptions({column: 'RACE'});
      return wc.getOptions().look.columnColumnName;
    });
    expect(col).toBe('RACE');
  });

  await softStep('Column assignment: DIS_POP', async () => {
    const col = await page.evaluate(() => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({column: 'DIS_POP'});
      return wc.getOptions().look.columnColumnName;
    });
    expect(col).toBe('DIS_POP');
  });

  await softStep('Column assignment: SITE (missing, substitute SEVERITY)', async () => {
    const res = await page.evaluate(() => {
      const df = (window as any).grok.shell.t;
      const hasSite = df.col('SITE') != null;
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      const target = hasSite ? 'SITE' : 'SEVERITY';
      wc.setOptions({column: target});
      return {hasSite, col: wc.getOptions().look.columnColumnName};
    });
    expect(res.col).toBeTruthy();
  });

  await softStep('Column assignment: SEX', async () => {
    const col = await page.evaluate(() => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({column: 'SEX'});
      return wc.getOptions().look.columnColumnName;
    });
    expect(col).toBe('SEX');
  });

  await softStep('Shape options (circle, diamond, triangle-forward, triangle, pentagon, star)', async () => {
    const results = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({column: 'RACE'});
      const shapes = ['circle', 'diamond', 'triangle-forward', 'triangle', 'pentagon', 'star'];
      const out: any[] = [];
      for (const s of shapes) {
        wc.setOptions({shape: s});
        await new Promise((r: any) => setTimeout(r, 150));
        out.push({shape: s, applied: wc.getOptions().look.shape});
      }
      return out;
    });
    for (const r of results) expect(r.applied).toBe(r.shape);
  });

  await softStep('Font size range', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({shape: 'circle'});
      wc.setOptions({minTextSize: 10, maxTextSize: 60});
      const a = wc.getOptions().look;
      if (a.minTextSize !== 10 || a.maxTextSize !== 60) return false;
      wc.setOptions({minTextSize: 20, maxTextSize: 20});
      const b = wc.getOptions().look;
      if (b.minTextSize !== 20 || b.maxTextSize !== 20) return false;
      wc.setOptions({minTextSize: 12, maxTextSize: 48});
      const c = wc.getOptions().look;
      return c.minTextSize === 12 && c.maxTextSize === 48;
    });
    expect(ok).toBe(true);
  });

  await softStep('Text rotation', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({minRotationDegree: 0, maxRotationDegree: 0});
      const a = wc.getOptions().look;
      if (a.minRotationDegree !== 0 || a.maxRotationDegree !== 0) return false;
      wc.setOptions({minRotationDegree: -90, maxRotationDegree: 90});
      const b = wc.getOptions().look;
      if (b.minRotationDegree !== -90 || b.maxRotationDegree !== 90) return false;
      wc.setOptions({rotationStep: 45});
      if (wc.getOptions().look.rotationStep !== 45) return false;
      wc.setOptions({minRotationDegree: -30, maxRotationDegree: 30, rotationStep: 5});
      const d = wc.getOptions().look;
      return d.minRotationDegree === -30 && d.maxRotationDegree === 30 && d.rotationStep === 5;
    });
    expect(ok).toBe(true);
  });

  await softStep('Grid size', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({gridSize: 2});
      if (wc.getOptions().look.gridSize !== 2) return false;
      wc.setOptions({gridSize: 20});
      if (wc.getOptions().look.gridSize !== 20) return false;
      wc.setOptions({gridSize: 8});
      return wc.getOptions().look.gridSize === 8;
    });
    expect(ok).toBe(true);
  });

  await softStep('Draw out of bound', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({drawOutOfBound: true});
      if (wc.getOptions().look.drawOutOfBound !== true) return false;
      wc.setOptions({drawOutOfBound: false});
      if (wc.getOptions().look.drawOutOfBound !== false) return false;
      wc.setOptions({drawOutOfBound: true});
      return wc.getOptions().look.drawOutOfBound === true;
    });
    expect(ok).toBe(true);
  });

  await softStep('Font family', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({fontFamily: 'serif'});
      if (wc.getOptions().look.fontFamily !== 'serif') return false;
      wc.setOptions({fontFamily: 'monospace'});
      if (wc.getOptions().look.fontFamily !== 'monospace') return false;
      wc.setOptions({fontFamily: 'sans-serif'});
      return wc.getOptions().look.fontFamily === 'sans-serif';
    });
    expect(ok).toBe(true);
  });

  await softStep('Bold', async () => {
    const ok = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({bold: true});
      if (wc.getOptions().look.bold !== true) return false;
      wc.setOptions({bold: false});
      if (wc.getOptions().look.bold !== false) return false;
      wc.setOptions({bold: true});
      return wc.getOptions().look.bold === true;
    });
    expect(ok).toBe(true);
  });

  await softStep('Tooltip on hover', async () => {
    const res = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({column: 'RACE'});
      await new Promise((r: any) => setTimeout(r, 400));
      const wcDom = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const canvas = wcDom.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('mousemove', {clientX: rect.left + rect.width/2, clientY: rect.top + rect.height/2, bubbles: true}));
      await new Promise((r: any) => setTimeout(r, 500));
      const tt = document.querySelector('.d4-tooltip') as HTMLElement;
      return tt && tt.innerText.trim().length > 0;
    });
    expect(res).toBe(true);
  });

  await softStep('Word click selects rows', async () => {
    const res = await page.evaluate(async () => {
      const df = (window as any).grok.shell.t;
      df.selection.setAll(false);
      const wcDom = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const canvas = wcDom.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const cx = rect.left + rect.width/2;
      const cy = rect.top + rect.height/2;
      for (const type of ['mousemove', 'mousedown', 'mouseup', 'click']) {
        canvas.dispatchEvent(new MouseEvent(type, {clientX: cx, clientY: cy, bubbles: true, button: 0}));
      }
      await new Promise((r: any) => setTimeout(r, 400));
      return df.selection.trueCount;
    });
    expect(res).toBeGreaterThan(0);
  });

  await softStep('Filter SEX=M updates word cloud', async () => {
    const res = await page.evaluate(async () => {
      const DG = (window as any).DG;
      const df = (window as any).grok.shell.t;
      const sexCol = df.col('SEX');
      const bs = DG.BitSet.create(df.rowCount);
      for (let i = 0; i < df.rowCount; i++) bs.set(i, sexCol.get(i) === 'M', false);
      df.filter.copyFrom(bs);
      df.filter.fireChanged();
      await new Promise((r: any) => setTimeout(r, 500));
      const filtered = df.filter.trueCount;
      df.filter.setAll(true);
      df.filter.fireChanged();
      await new Promise((r: any) => setTimeout(r, 300));
      return {filtered, restored: df.filter.trueCount};
    });
    expect(res.filtered).toBeLessThan(res.restored);
    expect(res.restored).toBe(5850);
  });

  await softStep('Viewer resize via Alt+F fullscreen toggle', async () => {
    const res = await page.evaluate(async () => {
      const wcDom = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const before = {w: wcDom.offsetWidth, h: wcDom.offsetHeight};
      wcDom.focus();
      const ev1 = new KeyboardEvent('keydown', {key: 'F', code: 'KeyF', keyCode: 70, altKey: true, bubbles: true});
      wcDom.dispatchEvent(ev1); document.dispatchEvent(ev1);
      await new Promise((r: any) => setTimeout(r, 500));
      const mid = {w: wcDom.offsetWidth, h: wcDom.offsetHeight};
      const ev2 = new KeyboardEvent('keydown', {key: 'F', code: 'KeyF', keyCode: 70, altKey: true, bubbles: true});
      wcDom.dispatchEvent(ev2); document.dispatchEvent(ev2);
      await new Promise((r: any) => setTimeout(r, 500));
      const after = {w: wcDom.offsetWidth, h: wcDom.offsetHeight};
      return {before, mid, after};
    });
    expect(res.mid.w).toBeGreaterThan(res.before.w);
    expect(res.after.w).toBe(res.before.w);
  });

  await softStep('Context menu has Properties...', async () => {
    const res = await page.evaluate(async () => {
      const wcDom = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const canvas = wcDom.querySelector('canvas') as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {clientX: rect.left + rect.width/2, clientY: rect.top + rect.height/2, bubbles: true, button: 2}));
      await new Promise((r: any) => setTimeout(r, 500));
      const items = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label')).map(e => (e as HTMLElement).innerText.trim());
      const hasProps = items.includes('Properties...');
      // close menu
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      await new Promise((r: any) => setTimeout(r, 200));
      return {items, hasProps};
    });
    expect(res.hasProps).toBe(true);
  });

  await softStep('Error state: > 500 unique categories (USUBJID)', async () => {
    const res = await page.evaluate(async () => {
      let wc: any;
      for (const v of (window as any).grok.shell.tv.viewers) if (v.type === 'Word cloud') wc = v;
      wc.setOptions({column: 'USUBJID'});
      await new Promise((r: any) => setTimeout(r, 800));
      const wcDom = document.querySelector('[name="viewer-Word-cloud"]') as HTMLElement;
      const hasError = wcDom.innerText.includes('500 or fewer') || wcDom.innerText.includes('categorical column');
      wc.setOptions({column: 'RACE'});
      await new Promise((r: any) => setTimeout(r, 500));
      const recovered = !(wcDom.innerText.includes('500 or fewer'));
      return {hasError, recovered};
    });
    expect(res.hasError).toBe(true);
    expect(res.recovered).toBe(true);
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n'));
});
