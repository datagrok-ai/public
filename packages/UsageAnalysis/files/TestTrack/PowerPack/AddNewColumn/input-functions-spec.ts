import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

test.use(specTestOptions);

test('AddNewColumn: input_functions (+ icon, drag, column-typed params) (SPGI)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Setup + open SPGI.csv + wait for semType and cell rendering
  await page.evaluate(async () => {
    const g: any = (window as any).grok;
    document.body.classList.add('selenium');
    g.shell.settings.showFiltersIconsConstantly = true;
    g.shell.windows.simpleMode = true;
    g.shell.closeAll();
    const df = await g.dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    g.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 3000);
    });
    const hasBioChem = Array.from({length: df.columns.length}, (_, i) => df.columns.byIndex(i))
      .some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').first().waitFor({timeout: 30000});

  await softStep('Step 2: Open Add New Column dialog', async () => {
    await page.evaluate(() => (document.querySelector('[name="icon-add-new-column"]') as HTMLElement).click());
    await page.locator('[name="dialog-Add-New-Column"]').waitFor({timeout: 10000});
    expect(await page.locator('[name="dialog-Add-New-Column"]').count()).toBe(1);
  });

  await softStep('Step 3: Hover Abs, click + icon → function inserted with param type', async () => {
    const result = await page.evaluate(async () => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]')!;
      const rows = Array.from(dialog.querySelectorAll('tr'));
      const absRow = rows.find((tr) => tr.querySelector('span[name="span-Abs"]'))!;
      const plus = absRow.querySelector('[name="icon-plus"]') as HTMLElement;
      const r = plus.getBoundingClientRect();
      for (const t of ['mouseenter', 'mouseover', 'mousemove'])
        plus.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, view: window}));
      await new Promise((res) => setTimeout(res, 200));
      plus.click();
      await new Promise((res) => setTimeout(res, 400));
      return (dialog.querySelector('.cm-content') as HTMLElement).textContent ?? '';
    });
    expect(result).toBe('Abs(num)');
  });

  await softStep('Step 4: Clear, add Abs via drag-and-drop → Abs(num)', async () => {
    // Clear formula
    await page.evaluate(() => {
      const cm = document.querySelector('[name="dialog-Add-New-Column"] .cm-content') as HTMLElement;
      cm.focus();
      document.execCommand('selectAll', false);
      document.execCommand('delete', false);
    });
    // Drag Abs function row to cm-content (uses Playwright real mouse events)
    const absSpan = page.locator('[name="dialog-Add-New-Column"] span[name="span-Abs"]');
    const cm = page.locator('[name="dialog-Add-New-Column"] .cm-content');
    await absSpan.scrollIntoViewIfNeeded();
    await absSpan.hover();
    await page.mouse.down();
    const cmBox = (await cm.boundingBox())!;
    // Move in small increments to trigger Dart's distance>5 threshold dnd
    const absBox = (await absSpan.boundingBox())!;
    const sx = absBox.x + absBox.width / 2, sy = absBox.y + absBox.height / 2;
    const tx = cmBox.x + cmBox.width / 2, ty = cmBox.y + cmBox.height / 2;
    const steps = 20;
    for (let i = 1; i <= steps; i++)
      await page.mouse.move(sx + (tx - sx) * i / steps, sy + (ty - sy) * i / steps, {steps: 1});
    await page.mouse.up();
    await page.waitForTimeout(500);
    const content = await page.evaluate(() => document.querySelector('[name="dialog-Add-New-Column"] .cm-content')?.textContent ?? '');
    expect(content).toBe('Abs(num)');
  });

  await softStep('Step 5: Click Structure column, + on getCLogP → Chem:getCLogP(${Structure})', async () => {
    const result = await page.evaluate(async () => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]')!;
      const cm = dialog.querySelector('.cm-content') as HTMLElement;
      cm.focus();
      document.execCommand('selectAll', false);
      document.execCommand('delete', false);
      await new Promise((r) => setTimeout(r, 150));
      // Click Structure row (row 1) in column grid canvas
      const gridEl = dialog.querySelector('.add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const yClick = r.y + rowHeight * 1 + rowHeight / 2;
      const opts = (x: number, y: number): MouseEventInit => ({bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0, buttons: 1});
      mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, yClick)));
      await new Promise((r) => setTimeout(r, 600));
      // + on getCLogP
      const rows = Array.from(dialog.querySelectorAll('tr'));
      const row = rows.find((tr) => tr.querySelector('span[name="span-getCLogP"]'))!;
      const plus = row.querySelector('[name="icon-plus"]') as HTMLElement;
      const pr = plus.getBoundingClientRect();
      for (const t of ['mouseenter', 'mouseover', 'mousemove'])
        plus.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: pr.x + pr.width / 2, clientY: pr.y + pr.height / 2, view: window}));
      await new Promise((res) => setTimeout(res, 200));
      plus.click();
      await new Promise((res) => setTimeout(res, 400));
      return cm.textContent ?? '';
    });
    expect(result).toBe('Chem:getCLogP(${Structure})');
  });

  await softStep('Step 6: Clear, drag getCLogP onto formula → Chem:getCLogP(${Structure})', async () => {
    await page.evaluate(() => {
      const cm = document.querySelector('[name="dialog-Add-New-Column"] .cm-content') as HTMLElement;
      cm.focus();
      document.execCommand('selectAll', false);
      document.execCommand('delete', false);
    });
    const span = page.locator('[name="dialog-Add-New-Column"] span[name="span-getCLogP"]');
    const cm = page.locator('[name="dialog-Add-New-Column"] .cm-content');
    await span.scrollIntoViewIfNeeded();
    await span.hover();
    await page.mouse.down();
    const sBox = (await span.boundingBox())!;
    const cmBox = (await cm.boundingBox())!;
    const sx = sBox.x + sBox.width / 2, sy = sBox.y + sBox.height / 2;
    const tx = cmBox.x + cmBox.width / 2, ty = cmBox.y + cmBox.height / 2;
    const steps = 20;
    for (let i = 1; i <= steps; i++)
      await page.mouse.move(sx + (tx - sx) * i / steps, sy + (ty - sy) * i / steps, {steps: 1});
    await page.mouse.up();
    await page.waitForTimeout(500);
    const content = await page.evaluate(() => document.querySelector('[name="dialog-Add-New-Column"] .cm-content')?.textContent ?? '');
    expect(content).toBe('Chem:getCLogP(${Structure})');
  });

  await softStep('Step 7: Click numeric column, + on Abs → Abs(${NumericCol})', async () => {
    const result = await page.evaluate(async () => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]')!;
      const cm = dialog.querySelector('.cm-content') as HTMLElement;
      cm.focus();
      document.execCommand('selectAll', false);
      document.execCommand('delete', false);
      await new Promise((r) => setTimeout(r, 150));
      // Click a numeric column (row 18 in SPGI — Chemical Space Y/X)
      const gridEl = dialog.querySelector('.add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const yClick = r.y + rowHeight * 18 + rowHeight / 2;
      const opts = (x: number, y: number): MouseEventInit => ({bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0, buttons: 1});
      mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, yClick)));
      await new Promise((r) => setTimeout(r, 600));
      // + on Abs
      const rows = Array.from(dialog.querySelectorAll('tr'));
      const row = rows.find((tr) => tr.querySelector('span[name="span-Abs"]'))!;
      const plus = row.querySelector('[name="icon-plus"]') as HTMLElement;
      const pr = plus.getBoundingClientRect();
      for (const t of ['mouseenter', 'mouseover', 'mousemove'])
        plus.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: pr.x + pr.width / 2, clientY: pr.y + pr.height / 2, view: window}));
      await new Promise((res) => setTimeout(res, 200));
      plus.click();
      await new Promise((res) => setTimeout(res, 400));
      return cm.textContent ?? '';
    });
    // Formula should reference a numeric column by name (Chemical Space X/Y)
    expect(result).toMatch(/^Abs\(\$\{Chemical Space [XY]\}\)$/);
  });

  await softStep('Step 8: Clear the text field', async () => {
    const cleared = await page.evaluate(async () => {
      const cm = document.querySelector('[name="dialog-Add-New-Column"] .cm-content') as HTMLElement;
      cm.focus();
      document.execCommand('selectAll', false);
      document.execCommand('delete', false);
      await new Promise((r) => setTimeout(r, 200));
      return cm.textContent ?? '';
    });
    expect(cleared).toBe('');
  });

  await softStep('Step 9: Click Id column, sort by name → functions alphabetical', async () => {
    const first10 = await page.evaluate(async () => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]')!;
      const gridEl = dialog.querySelector('.add-new-column-columns-grid')!;
      const canvases = gridEl.querySelectorAll('canvas');
      const mainCanvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const r = mainCanvas.getBoundingClientRect();
      const rowHeight = r.height / 21;
      const xClick = r.x + r.width / 2;
      const yClick = r.y + rowHeight * 0 + rowHeight / 2; // Id is row 0
      const opts = (x: number, y: number): MouseEventInit => ({bubbles: true, cancelable: true, view: window, clientX: x, clientY: y, button: 0, buttons: 1});
      mainCanvas.dispatchEvent(new MouseEvent('mousedown', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('mouseup', opts(xClick, yClick)));
      mainCanvas.dispatchEvent(new MouseEvent('click', opts(xClick, yClick)));
      await new Promise((r) => setTimeout(r, 600));
      // Sort icon → By name
      const sortIcon = dialog.querySelector('[name="icon-sort-alt"]') as HTMLElement;
      sortIcon.dispatchEvent(new MouseEvent('mousedown', {bubbles: true, button: 0}));
      sortIcon.dispatchEvent(new MouseEvent('mouseup', {bubbles: true, button: 0}));
      sortIcon.dispatchEvent(new MouseEvent('click', {bubbles: true, button: 0}));
      await new Promise((r) => setTimeout(r, 500));
      const byName = document.querySelector('[name="div-By-name"]') as HTMLElement | null;
      byName?.dispatchEvent(new MouseEvent('click', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 800));
      return Array.from(dialog.querySelectorAll('#actions span[name^="span-"]'))
        .slice(0, 10)
        .map((s) => s.getAttribute('name')!.replace(/^span-/, ''));
    });
    const sorted = [...first10].sort((a, b) => a.localeCompare(b, undefined, {sensitivity: 'base'}));
    expect(first10).toEqual(sorted);
  });

  await softStep('Step 10: + on Abs with Id selected → Abs(num), column not passed (type mismatch)', async () => {
    const result = await page.evaluate(async () => {
      const dialog = document.querySelector('[name="dialog-Add-New-Column"]')!;
      const rows = Array.from(dialog.querySelectorAll('tr'));
      const row = rows.find((tr) => tr.querySelector('span[name="span-Abs"]'))!;
      const plus = row.querySelector('[name="icon-plus"]') as HTMLElement;
      const pr = plus.getBoundingClientRect();
      for (const t of ['mouseenter', 'mouseover', 'mousemove'])
        plus.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, clientX: pr.x + pr.width / 2, clientY: pr.y + pr.height / 2, view: window}));
      await new Promise((res) => setTimeout(res, 200));
      plus.click();
      await new Promise((res) => setTimeout(res, 400));
      return (dialog.querySelector('.cm-content') as HTMLElement).textContent ?? '';
    });
    // Id is string — Abs wants num — column should NOT be substituted
    expect(result).toBe('Abs(num)');
  });

  // Cleanup: cancel dialog
  await page.evaluate(() => {
    const cancel = document.querySelector('[name="button-Add-New-Column---CANCEL"]') as HTMLElement | null;
    cancel?.click();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
