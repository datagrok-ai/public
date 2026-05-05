import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — query post-process balloon and saved layout', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  await page.evaluate(async () => {
    const w: any = window;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    w.grok.shell.windows.showBrowse = true;
    // Hook grok.shell.info to capture post-process balloons.
    if (!w.__capturedInfo) {
      w.__capturedInfo = [];
      const orig = w.grok.shell.info.bind(w.grok.shell);
      w.grok.shell.info = (msg: any) => { w.__capturedInfo.push(String(msg)); return orig(msg); };
    }
    // Pre-clean any leftover queries from previous runs (server normalizes
    // `Test_Postprocessing` → `TestPostprocessing`; suffixed copies pile up).
    try {
      const list = await w.grok.dapi.queries
        .filter('name in ("Test_Postprocessing", "TestPostprocessing")').list();
      for (const q of list)
        try { await w.grok.dapi.queries.delete(q); } catch (e) {}
    } catch (e) {}
  });

  await page.goto(`${process.env.DATAGROK_URL ?? 'http://localhost:8888'}/browse`);
  await page.waitForFunction(() => {
    return document.querySelectorAll('.d4-tree-view-group-label').length > 5;
  }, undefined, {timeout: 30_000});

  let savedQueryId = '';

  await softStep('Browse → Databases → Postgres → NorthwindTest → right-click → New Query...', async () => {
    const result = await page.evaluate(async () => {
      const findChildGroupByLabel = (parent: ParentNode, label: string): Element | null => {
        const candidates = Array.from(parent.querySelectorAll('.d4-tree-view-group-label'))
          .filter((el) => el.textContent?.trim() === label);
        for (const c of candidates) {
          const g = c.closest('.d4-tree-view-group');
          if (!g) continue;
          if (parent === document || (parent as Element).contains(g))
            return g;
        }
        return null;
      };
      const expandAndWaitFor = async (group: Element, nextLabel: string) => {
        const tri = group.querySelector(':scope > .d4-tree-view-node > .d4-tree-view-tri') as HTMLElement | null;
        if (tri && !tri.classList.contains('d4-tree-view-tri-expanded'))
          tri.click();
        const deadline = Date.now() + 30_000;
        while (Date.now() < deadline) {
          await new Promise((r) => setTimeout(r, 200));
          const found = findChildGroupByLabel(group, nextLabel);
          if (found) return true;
        }
        return false;
      };
      const path = ['Databases', 'Postgres', 'NorthwindTest'];
      let scope: ParentNode = document;
      for (let i = 0; i < path.length - 1; i++) {
        const label = path[i];
        const group = findChildGroupByLabel(scope, label);
        if (!group) return {ok: false, missing: label};
        const next = path[i + 1];
        const ok = await expandAndWaitFor(group, next);
        if (!ok) return {ok: false, expandTimedOutAt: label};
        scope = group;
      }
      const nwGroup = findChildGroupByLabel(scope, 'NorthwindTest');
      if (!nwGroup) return {ok: false, missing: 'NorthwindTest'};
      const node = nwGroup.querySelector(':scope > .d4-tree-view-node') as HTMLElement | null;
      if (!node) return {ok: false, missing: 'NorthwindTest tree-node'};
      (node.querySelector('.d4-tree-view-group-label') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 250));
      node.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      let target: HTMLElement | undefined;
      const menuDeadline = Date.now() + 8_000;
      while (Date.now() < menuDeadline) {
        await new Promise((r) => setTimeout(r, 150));
        target = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'New Query...') as HTMLElement | undefined;
        if (target) break;
      }
      if (!target) return {ok: false, missing: 'New Query... menu item'};
      target.click();
      for (let i = 0; i < 80; i++) {
        await new Promise((r) => setTimeout(r, 200));
        const cm: any = document.querySelector('.CodeMirror');
        if (cm && cm.CodeMirror) {
          const v: any = (window as any).grok.shell.v;
          return {ok: true, viewType: v?.type, body: cm.CodeMirror.getValue()};
        }
      }
      return {ok: false, missing: 'CodeMirror editor'};
    });
    expect(result.ok).toBe(true);
    expect(result.viewType).toBe('DataQueryView');
  });

  await softStep('Set name = Test_Postprocessing and body = select * from products', async () => {
    const result = await page.evaluate(() => {
      const nameInput = document.querySelector('input[name="input-Name"]') as HTMLInputElement | null;
      if (!nameInput) return {ok: false, missing: 'name input'};
      nameInput.focus();
      document.execCommand('selectAll');
      document.execCommand('insertText', false, 'Test_Postprocessing');
      const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
      if (!cm) return {ok: false, missing: 'CodeMirror'};
      cm.setValue('select * from products');
      return {ok: true, name: nameInput.value, body: cm.getValue()};
    });
    expect(result.ok).toBe(true);
    expect(result.name).toBe('Test_Postprocessing');
    expect(result.body).toBe('select * from products');
  });

  await softStep('Run query — inline preview grid appears', async () => {
    const result = await page.evaluate(async () => {
      const playIcon = document.querySelector('[name="icon-play"]') as HTMLElement | null;
      if (!playIcon) return {ok: false, missing: 'icon-play'};
      const ribbonItem = (playIcon.closest('.d4-ribbon-item') ?? playIcon) as HTMLElement;
      const v: any = (window as any).grok.shell.v;
      const root: HTMLElement | undefined = v?.root;
      if (!root) return {ok: false, missing: 'view root'};
      const deadline = Date.now() + 60_000;
      let attempts = 0;
      while (Date.now() < deadline) {
        ribbonItem.click();
        playIcon.click();
        attempts++;
        for (let i = 0; i < 12; i++) {
          await new Promise((r) => setTimeout(r, 250));
          const grid = root.querySelector('[name="viewer-Grid"]');
          const canvas = grid?.querySelector('canvas');
          if (canvas) {
            const r2 = grid!.getBoundingClientRect();
            return {ok: true, width: r2.width, attempts};
          }
        }
      }
      return {ok: false, missing: 'inline grid', attempts};
    });
    expect(result.ok).toBe(true);
    expect(result.width).toBeGreaterThan(100);
  });

  await softStep('Post-Process: replace line 7 with grok.shell.info(result.rowCount);', async () => {
    const result = await page.evaluate(async () => {
      const tab = Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
        .find((el) => el.textContent?.trim() === 'Post-Process') as HTMLElement | undefined;
      if (!tab) return {ok: false, missing: 'Post-Process tab'};
      tab.click();
      let cm: any = null;
      for (let i = 0; i < 30; i++) {
        await new Promise((r) => setTimeout(r, 200));
        cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
        if (cm) break;
      }
      if (!cm) return {ok: false, missing: 'CodeMirror'};
      // Template seeds line 7 (index 6) with `console.log(result.rowCount)`.
      // Replace it with the scenario's grok.shell.info(result.rowCount); call.
      cm.replaceRange(
        'grok.shell.info(result.rowCount);',
        {line: 6, ch: 0},
        {line: 6, ch: cm.getLine(6).length}
      );
      return {ok: true, line7: cm.getLine(6)};
    });
    expect(result.ok).toBe(true);
    expect(result.line7).toBe('grok.shell.info(result.rowCount);');
  });

  await softStep('Layout: add scatter plot and correlation plot', async () => {
    const result = await page.evaluate(async () => {
      const tab = Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
        .find((el) => el.textContent?.trim() === 'Layout') as HTMLElement | undefined;
      if (!tab) return {ok: false, missing: 'Layout tab'};
      tab.click();
      // Layout tab needs the result DF (Run was clicked already). Wait for Grid.
      for (let i = 0; i < 40; i++) {
        await new Promise((r) => setTimeout(r, 200));
        const root = (window as any).grok.shell.v?.root;
        if (root?.querySelector('[name="viewer-Grid"]')) break;
      }
      const sIcon = document.querySelector('[name="icon-scatter-plot"]') as HTMLElement | null;
      if (!sIcon) return {ok: false, missing: 'icon-scatter-plot'};
      sIcon.click();
      for (let i = 0; i < 30; i++) {
        await new Promise((r) => setTimeout(r, 200));
        if (document.querySelector('[name="viewer-Scatter-plot"]')) break;
      }
      const cIcon = document.querySelector('[name="icon-correlation-plot"]') as HTMLElement | null;
      if (!cIcon) return {ok: false, missing: 'icon-correlation-plot'};
      cIcon.click();
      for (let i = 0; i < 30; i++) {
        await new Promise((r) => setTimeout(r, 200));
        if (document.querySelector('[name="viewer-Correlation-plot"]')) break;
      }
      return {
        ok: true,
        sp: !!document.querySelector('[name="viewer-Scatter-plot"]'),
        cp: !!document.querySelector('[name="viewer-Correlation-plot"]'),
      };
    });
    expect(result.ok).toBe(true);
    expect(result.sp).toBe(true);
    expect(result.cp).toBe(true);
  });

  await softStep('Save the query', async () => {
    await page.locator('[name="button-Save"]').first().click();
    const id = await page.evaluate(async () => {
      for (let i = 0; i < 50; i++) {
        await new Promise((r) => setTimeout(r, 300));
        const q = await (window as any).grok.dapi.queries
          .filter('name in ("Test_Postprocessing", "TestPostprocessing")').first().catch(() => null);
        if (q) return q.id;
      }
      return null;
    });
    expect(id).toBeTruthy();
    savedQueryId = id ?? '';
  });

  await softStep('Close all and re-open Browse', async () => {
    await page.evaluate(() => { (window as any).grok.shell.closeAll(); });
    await page.goto(`${process.env.DATAGROK_URL ?? 'http://localhost:8888'}/browse`);
    await page.waitForFunction(() => {
      return document.querySelectorAll('.d4-tree-view-group-label').length > 5;
    }, undefined, {timeout: 30_000});
    // Re-install the info hook on the new document, since closeAll/navigation
    // doesn't necessarily reset window state but the hook needs to survive any
    // grok.shell rebind on view re-init.
    await page.evaluate(() => {
      const w: any = window;
      if (!w.__capturedInfo) w.__capturedInfo = [];
      const orig = w.grok.shell.info.bind(w.grok.shell);
      w.grok.shell.info = (msg: any) => { w.__capturedInfo.push(String(msg)); return orig(msg); };
    });
  });

  await softStep('Preview Test_Postprocessing — both viewers + green balloon 77', async () => {
    const result = await page.evaluate(async () => {
      const findChildGroupByLabel = (parent: ParentNode, label: string): Element | null => {
        const candidates = Array.from(parent.querySelectorAll('.d4-tree-view-group-label'))
          .filter((el) => el.textContent?.trim() === label);
        for (const c of candidates) {
          const g = c.closest('.d4-tree-view-group');
          if (!g) continue;
          if (parent === document || (parent as Element).contains(g))
            return g;
        }
        return null;
      };
      const expandAndWaitFor = async (group: Element, nextLabel: string) => {
        const tri = group.querySelector(':scope > .d4-tree-view-node > .d4-tree-view-tri') as HTMLElement | null;
        if (tri && !tri.classList.contains('d4-tree-view-tri-expanded'))
          tri.click();
        const deadline = Date.now() + 30_000;
        while (Date.now() < deadline) {
          await new Promise((r) => setTimeout(r, 200));
          const found = findChildGroupByLabel(group, nextLabel);
          if (found) return true;
        }
        return false;
      };
      let scope: ParentNode = document;
      for (const label of ['Databases', 'Postgres']) {
        const group = findChildGroupByLabel(scope, label);
        if (!group) return {ok: false, lost: label};
        const next = label === 'Databases' ? 'Postgres' : 'NorthwindTest';
        await expandAndWaitFor(group, next);
        scope = group;
      }
      const nw = findChildGroupByLabel(scope, 'NorthwindTest');
      if (!nw) return {ok: false, missing: 'NorthwindTest'};
      const tri = nw.querySelector(':scope > .d4-tree-view-node > .d4-tree-view-tri') as HTMLElement | null;
      if (tri && !tri.classList.contains('d4-tree-view-tri-expanded')) tri.click();
      let tpGroup: Element | null = null;
      const deadline = Date.now() + 30_000;
      while (Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 250));
        tpGroup = findChildGroupByLabel(nw, 'Test_Postprocessing');
        if (tpGroup) break;
      }
      if (!tpGroup) return {ok: false, missing: 'Test_Postprocessing in tree'};
      (window as any).__capturedInfo = [];
      const lbl = tpGroup.querySelector('.d4-tree-view-group-label') as HTMLElement | null;
      lbl?.click();
      for (let i = 0; i < 80; i++) {
        await new Promise((r) => setTimeout(r, 250));
        const v: any = (window as any).grok.shell.v;
        if (v?.type === 'TableView' && v?.dataFrame?.rowCount > 0) {
          const captured = [...((window as any).__capturedInfo ?? [])];
          const viewers = Array.from((window as any).grok.shell.tv.viewers).map((x: any) => x.type);
          if (captured.length > 0 && viewers.includes('Scatter plot') && viewers.includes('Correlation plot'))
            return {ok: true, viewers, rowCount: v.dataFrame.rowCount, captured};
        }
      }
      const v: any = (window as any).grok.shell.v;
      return {
        ok: false,
        viewers: v?.type === 'TableView' ? Array.from((window as any).grok.shell.tv.viewers).map((x: any) => x.type) : null,
        rowCount: v?.dataFrame?.rowCount,
        captured: [...((window as any).__capturedInfo ?? [])],
      };
    });
    expect(result.ok).toBe(true);
    expect(result.rowCount).toBe(77);
    expect(result.captured).toContain('77');
    expect(result.viewers).toContain('Scatter plot');
    expect(result.viewers).toContain('Correlation plot');
  });

  await softStep('Right-click Test_Postprocessing → Edit...', async () => {
    const result = await page.evaluate(async () => {
      const findChildGroupByLabel = (parent: ParentNode, label: string): Element | null => {
        const candidates = Array.from(parent.querySelectorAll('.d4-tree-view-group-label'))
          .filter((el) => el.textContent?.trim() === label);
        for (const c of candidates) {
          const g = c.closest('.d4-tree-view-group');
          if (!g) continue;
          if (parent === document || (parent as Element).contains(g))
            return g;
        }
        return null;
      };
      let scope: ParentNode = document;
      for (const label of ['Databases', 'Postgres', 'NorthwindTest']) {
        scope = findChildGroupByLabel(scope, label) as Element;
        if (!scope) return {ok: false, lost: label};
      }
      const tpGroup = findChildGroupByLabel(scope, 'Test_Postprocessing');
      if (!tpGroup) return {ok: false, missing: 'Test_Postprocessing'};
      const node = tpGroup.querySelector(':scope > .d4-tree-view-node') as HTMLElement | null;
      if (!node) return {ok: false, missing: 'tree node'};
      node.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      let target: HTMLElement | undefined;
      const deadline = Date.now() + 8_000;
      while (Date.now() < deadline) {
        await new Promise((r) => setTimeout(r, 150));
        target = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
        if (target) break;
      }
      if (!target) return {ok: false, missing: 'Edit... menu item'};
      target.click();
      for (let i = 0; i < 80; i++) {
        await new Promise((r) => setTimeout(r, 200));
        const v: any = (window as any).grok.shell.v;
        const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
        if (v?.type === 'DataQueryView' && cm)
          return {ok: true, viewType: v.type, body: cm.getValue()};
      }
      return {ok: false, missing: 'DataQueryView'};
    });
    expect(result.ok).toBe(true);
    expect(result.viewType).toBe('DataQueryView');
  });

  await softStep('Run from Post-Process tab — green balloon 77', async () => {
    const result = await page.evaluate(async () => {
      const tab = Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
        .find((el) => el.textContent?.trim() === 'Post-Process') as HTMLElement | undefined;
      if (!tab) return {ok: false, missing: 'Post-Process tab'};
      tab.click();
      await new Promise((r) => setTimeout(r, 400));
      (window as any).__capturedInfo = [];
      const playIcon = document.querySelector('[name="icon-play"]') as HTMLElement | null;
      if (!playIcon) return {ok: false, missing: 'icon-play'};
      const ribbonItem = (playIcon.closest('.d4-ribbon-item') ?? playIcon) as HTMLElement;
      const deadline = Date.now() + 30_000;
      let attempts = 0;
      while (Date.now() < deadline) {
        ribbonItem.click();
        playIcon.click();
        attempts++;
        for (let i = 0; i < 12; i++) {
          await new Promise((r) => setTimeout(r, 250));
          const captured = [...((window as any).__capturedInfo ?? [])];
          if (captured.length > 0)
            return {ok: true, attempts, captured};
        }
      }
      return {ok: false, missing: 'green balloon', attempts};
    });
    expect(result.ok).toBe(true);
    expect(result.captured).toContain('77');
  });

  await softStep('Run from Layout tab — green balloon 77 + both viewers', async () => {
    const result = await page.evaluate(async () => {
      const tab = Array.from(document.querySelectorAll('.d4-tab-host .d4-tab-header'))
        .find((el) => el.textContent?.trim() === 'Layout') as HTMLElement | undefined;
      if (!tab) return {ok: false, missing: 'Layout tab'};
      tab.click();
      await new Promise((r) => setTimeout(r, 600));
      (window as any).__capturedInfo = [];
      const playIcon = document.querySelector('[name="icon-play"]') as HTMLElement | null;
      if (!playIcon) return {ok: false, missing: 'icon-play'};
      const ribbonItem = (playIcon.closest('.d4-ribbon-item') ?? playIcon) as HTMLElement;
      const v: any = (window as any).grok.shell.v;
      const root = v?.root;
      const deadline = Date.now() + 30_000;
      let attempts = 0;
      while (Date.now() < deadline) {
        ribbonItem.click();
        playIcon.click();
        attempts++;
        for (let i = 0; i < 12; i++) {
          await new Promise((r) => setTimeout(r, 250));
          const captured = [...((window as any).__capturedInfo ?? [])];
          const sp = !!root?.querySelector('[name="viewer-Scatter-plot"]');
          const cp = !!root?.querySelector('[name="viewer-Correlation-plot"]');
          const grid = !!root?.querySelector('[name="viewer-Grid"]');
          if (captured.length > 0 && sp && cp && grid)
            return {ok: true, attempts, captured, sp, cp, grid};
        }
      }
      return {ok: false, missing: 'balloon + viewers', attempts};
    });
    expect(result.ok).toBe(true);
    expect(result.captured).toContain('77');
    expect(result.sp).toBe(true);
    expect(result.cp).toBe(true);
  });

  // Cleanup — delete the saved query.
  await page.evaluate(async (id) => {
    if (!id) {
      const list = await (window as any).grok.dapi.queries
        .filter('name in ("Test_Postprocessing", "TestPostprocessing")').list();
      for (const q of list)
        try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
      return;
    }
    const q = await (window as any).grok.dapi.queries.find(id).catch(() => null);
    if (q)
      try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
  }, savedQueryId);

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
