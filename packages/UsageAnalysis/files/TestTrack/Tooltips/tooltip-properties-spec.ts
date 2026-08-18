import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Viewers: tooltip properties and heuristics for column selection', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (grok as any).shell.settings.showFiltersIconsConstantly = true;
    (grok as any).shell.windows.simpleMode = true;
    (grok as any).shell.closeAll();
    const df = await (grok as any).dapi.files.readCsv('System:DemoFiles/chem/SPGI.csv');
    (grok as any).shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(() => resolve(), 3000);
    });
    const cols = Array.from({length: df.columns.length}, (_, i: number) => df.columns.byIndex(i));
    const hasBioChem = cols.some((c: any) => c.semType === 'Molecule' || c.semType === 'Macromolecule');
    if (hasBioChem) {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise<void>((r) => setTimeout(r, 200));
      }
      await new Promise<void>((r) => setTimeout(r, 5000));
    }
  });
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});

  await softStep('Add a scatter plot and a box plot', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="icon-scatter-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 10000});
    await page.waitForTimeout(500);
    await page.evaluate(() => {
      (document.querySelector('[name="icon-box-plot"]') as HTMLElement).click();
    });
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 10000});
    await page.waitForTimeout(500);
    const types = await page.evaluate(() =>
      Array.from((grok as any).shell.tv.viewers).map((v: any) => v.type));
    expect(types).toContain('Scatter plot');
    expect(types).toContain('Box plot');
  });

  await softStep('Tooltip section: Row Tooltip is disabled, Show Tooltip has 3 options', async () => {

    await page.evaluate(() => {
      const v = document.querySelector('[name="viewer-Scatter-plot"]') as HTMLElement;
      let node: HTMLElement | null = v;
      while (node && !(typeof node.className === 'string' && node.className.includes('panel-base')))
        node = node.parentElement;
      const settings = node!.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      settings.click();
    });
    await page.waitForTimeout(600);

    await page.evaluate(() => {
      const cat = document.querySelector('[name="prop-category-tooltip"]') as HTMLElement | null;
      cat?.click();
    });
    await page.waitForTimeout(400);
    const info = await page.evaluate(() => {
      const showRow = document.querySelector('[name="prop-show-tooltip"]') as HTMLElement | null;
      const rowRow = document.querySelector('[name="prop-row-tooltip"]') as HTMLElement | null;
      const rowOpacity = rowRow ? (rowRow.getAttribute('style') ?? '').includes('opacity') : false;
      const rowInput = rowRow?.querySelector('input.property-grid-ellipsis-editor-input') as HTMLInputElement | null;
      const rowValue = rowInput ? rowInput.value : null;
      const showLabel = showRow?.querySelector('[name="prop-view-show-tooltip"]') as HTMLElement | null;

      showLabel?.click();
      const select = document.querySelector('[name="prop-show-tooltip"] select') as HTMLSelectElement | null;
      const options = select ? Array.from(select.options).map((o) => o.value) : [];
      const defaultValue = showLabel?.textContent ?? select?.value ?? null;
      return {rowOpacity, rowValue, options, defaultValue};
    });
    expect(info.rowOpacity).toBe(true);
    expect(info.rowValue).toBe('');
    expect(info.options).toEqual(['do not show', 'inherit from table', 'show custom tooltip']);
    expect(info.defaultValue).toBe('inherit from table');
  });

  await softStep('Set Show Tooltip = show custom tooltip on scatter plot and box plot; enable Show Visible Columns In Tooltip on grid', async () => {

    const result = await page.evaluate(async () => {
      const tv: any = (grok as any).shell.tv;
      let sp: any = null;
      let bp: any = null;
      for (const v of tv.viewers as any) {
        if (v.type === 'Scatter plot' && !sp) sp = v;
        if (v.type === 'Box plot' && !bp) bp = v;
      }
      sp.setOptions({showTooltip: 'show custom tooltip'});
      bp.setOptions({showTooltip: 'show custom tooltip'});
      tv.grid.setOptions({showVisibleColumnsInTooltip: true});
      await new Promise<void>((r) => setTimeout(r, 300));
      return {
        sp: sp.props.showTooltip,
        bp: bp.props.showTooltip,
        gridShowTooltip: tv.grid.props.showTooltip,
        gridShowVisibleCols: tv.grid.props.showVisibleColumnsInTooltip,
      };
    });
    expect(result.sp).toBe('show custom tooltip');
    expect(result.bp).toBe('show custom tooltip');
    expect(result.gridShowTooltip).toBe('show custom tooltip');
    expect(result.gridShowVisibleCols).toBe(true);
  });

  await softStep('Add reference scatter plot (default tooltip) and compare custom vs default tooltip content', async () => {
    await page.evaluate(() => {
      (document.querySelector('[name="icon-scatter-plot"]') as HTMLElement).click();
    });
    await page.waitForTimeout(800);
    const modes = await page.evaluate(() => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      return {sp1Mode: sps[0].props.showTooltip, sp2Mode: sps[1].props.showTooltip};
    });
    expect(modes.sp1Mode).toBe('show custom tooltip');
    expect(modes.sp2Mode).toBe('inherit from table');

    async function readTooltipByDomIndex(viewerName: string, domIndex: number): Promise<string> {
      const locator = page.locator(`[name="${viewerName}"]`).nth(domIndex).locator('canvas').first();
      const box = await locator.boundingBox();
      if (!box) return '';
      await page.evaluate(() => (window as any).ui?.tooltip?.hide?.());
      await page.mouse.move(0, 0);
      await page.waitForTimeout(200);
      await page.mouse.move(box.x - 30, box.y + box.height / 2);
      await page.waitForTimeout(150);
      await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
      await page.waitForTimeout(1200);
      return page.evaluate(() => {
        const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
        return tt && getComputedStyle(tt).display !== 'none' ? (tt.textContent ?? '') : '';
      });
    }
    const firstTooltip = await readTooltipByDomIndex('viewer-Scatter-plot', 0);
    const secondTooltip = await readTooltipByDomIndex('viewer-Scatter-plot', 1);
    const axesInBoth = (t: string) => t.includes('CAST Idea ID') && t.includes('Chemical Space Y');
    expect(axesInBoth(firstTooltip)).toBe(true);
    expect(axesInBoth(secondTooltip)).toBe(true);
  });

  await softStep('Tooltip > Hide on the custom scatter plot toggles the submenu to Show Custom', async () => {
    const after = await page.evaluate(async () => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      const sp1 = sps[0]; 
      const canvas = sp1.root.querySelector('canvas') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const x = r.left + r.width / 2;
      const y = r.top + r.height / 2;
      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const popup = document.querySelector('.d4-menu-popup')!;
      const tg = popup.querySelector('[name="div-Tooltip"]') as HTMLElement;
      tg.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Tooltip---Hide"]') as HTMLElement).click();
      await new Promise<void>((r) => setTimeout(r, 500));

      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const tg2 = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
      tg2.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg2.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      const sp1Items = Array.from(document.querySelectorAll('[name^="div-Tooltip---"]'))
        .map((e) => e.getAttribute('name'));
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
      document.body.click();
      await new Promise<void>((r) => setTimeout(r, 300));

      const sp2 = sps[1];
      const c2 = sp2.root.querySelector('canvas') as HTMLCanvasElement;
      const rr = c2.getBoundingClientRect();
      const xx = rr.left + rr.width / 2;
      const yy = rr.top + rr.height / 2;
      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        c2.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: xx, clientY: yy}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const tg3 = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
      tg3.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg3.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      const sp2Items = Array.from(document.querySelectorAll('[name^="div-Tooltip---"]'))
        .map((e) => e.getAttribute('name'));
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
      document.body.click();
      return {sp1Items, sp2Items};
    });
    expect(after.sp1Items).toContain('div-Tooltip---Show-Custom');
    expect(after.sp1Items).not.toContain('div-Tooltip---Hide');
    expect(after.sp2Items).toContain('div-Tooltip---Hide');
  });

  await softStep('Tooltip > Show Custom restores the custom tooltip on scatter plot', async () => {
    const after = await page.evaluate(async () => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      const sp1 = sps[0];
      const canvas = sp1.root.querySelector('canvas') as HTMLCanvasElement;
      const r = canvas.getBoundingClientRect();
      const x = r.left + r.width / 2;
      const y = r.top + r.height / 2;
      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const tg = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
      tg.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Tooltip---Show-Custom"]') as HTMLElement).click();
      await new Promise<void>((r) => setTimeout(r, 500));

      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        canvas.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const tg2 = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
      tg2.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg2.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      const items = Array.from(document.querySelectorAll('[name^="div-Tooltip---"]'))
        .map((e) => e.getAttribute('name'));
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
      document.body.click();
      return {items};
    });
    expect(after.items).toContain('div-Tooltip---Hide');
    expect(after.items).not.toContain('div-Tooltip---Show-Custom');
  });

  await softStep('Properties panel toggle: switch show custom tooltip → do not show on scatter plot', async () => {

    await page.evaluate(async () => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      const sp1 = sps[0];
      let node: HTMLElement | null = sp1.root;
      while (node && !(typeof node.className === 'string' && node.className.includes('panel-base')))
        node = node.parentElement;
      const settings = node!.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
      settings.click();
      await new Promise<void>((r) => setTimeout(r, 600));
      (document.querySelector('[name="prop-category-tooltip"]') as HTMLElement | null)?.click();
      await new Promise<void>((r) => setTimeout(r, 300));
      sp1.setOptions({showTooltip: 'do not show'});
      await new Promise<void>((r) => setTimeout(r, 300));
    });
    const sp1Mode = await page.evaluate(() => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      return sps[0].props.showTooltip;
    });
    expect(sp1Mode).toBe('do not show');
  });

  await softStep('Hide default tooltip via reference viewer: hides for all "inherit from table" viewers; custom viewers unaffected', async () => {

    await page.evaluate(() => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      const bp = Array.from((grok as any).shell.tv.viewers).find((v: any) => v.type === 'Box plot') as any;
      sps[0].setOptions({showTooltip: 'show custom tooltip'});
      bp.setOptions({showTooltip: 'inherit from table'});
    });
    await page.waitForTimeout(400);
    const after = await page.evaluate(async () => {
      const sps = Array.from((grok as any).shell.tv.viewers).filter((v: any) => v.type === 'Scatter plot') as any[];
      const bp = Array.from((grok as any).shell.tv.viewers).find((v: any) => v.type === 'Box plot') as any;

      const sp2 = sps[1];
      const c = sp2.root.querySelector('canvas') as HTMLCanvasElement;
      const r = c.getBoundingClientRect();
      const x = r.left + r.width / 2;
      const y = r.top + r.height / 2;
      for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
        c.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: x, clientY: y}));
      }
      await new Promise<void>((r) => setTimeout(r, 500));
      const tg = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
      tg.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
      tg.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
      await new Promise<void>((r) => setTimeout(r, 400));
      (document.querySelector('[name="div-Tooltip---Hide"]') as HTMLElement).click();
      await new Promise<void>((r) => setTimeout(r, 600));
      document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
      document.body.click();
      await new Promise<void>((r) => setTimeout(r, 200));

      async function inspect(viewerRoot: HTMLElement) {
        const c2 = viewerRoot.querySelector('canvas') as HTMLCanvasElement;
        const rr = c2.getBoundingClientRect();
        const xx = rr.left + rr.width / 2;
        const yy = rr.top + rr.height / 2;
        for (const t of ['mousedown', 'mouseup', 'contextmenu']) {
          c2.dispatchEvent(new MouseEvent(t, {bubbles: true, cancelable: true, button: 2, buttons: 2, clientX: xx, clientY: yy}));
        }
        await new Promise<void>((r) => setTimeout(r, 500));
        const tg2 = document.querySelector('.d4-menu-popup [name="div-Tooltip"]') as HTMLElement;
        tg2.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, cancelable: true}));
        tg2.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, cancelable: true}));
        await new Promise<void>((r) => setTimeout(r, 400));
        const items = Array.from(document.querySelectorAll('[name^="div-Tooltip---"]'))
          .map((e) => e.getAttribute('name'));
        document.body.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.body.click();
        await new Promise<void>((r) => setTimeout(r, 300));
        return items;
      }
      const sp1Items = await inspect(sps[0].root);
      const sp2Items = await inspect(sps[1].root);
      const bpItems = await inspect(bp.root);
      return {sp1Items, sp2Items, bpItems};
    });

    expect(after.sp1Items).toContain('div-Tooltip---Hide');
    expect(after.sp1Items).not.toContain('div-Tooltip---Show-Custom');

    expect(after.sp2Items).toContain('div-Tooltip---Show-Custom');
    expect(after.bpItems).toContain('div-Tooltip---Show-Custom');
  });

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((s) => `  - ${s.step}: ${s.error}`).join('\n'));
});
