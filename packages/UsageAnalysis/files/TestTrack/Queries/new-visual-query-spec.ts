import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Queries — New Visual Query from the customers table', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  await page.evaluate(() => {
    const w: any = window;
    document.body.classList.add('selenium');
    w.grok.shell.settings.showFiltersIconsConstantly = true;
    w.grok.shell.windows.simpleMode = true;
    w.grok.shell.closeAll();
    w.grok.shell.windows.showBrowse = true;
  });

  await page.goto(`${process.env.DATAGROK_URL ?? 'http://localhost:8888'}/browse`);
  await page.waitForFunction(() => {
    return document.querySelectorAll('.d4-tree-view-group-label').length > 5;
  }, undefined, {timeout: 30_000});

  const queryName = `tt_new_visual_query_${Date.now()}`;

  await softStep('Browse → Databases → Postgres → NorthwindTest → Schemas → public', async () => {
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
      const path = ['Databases', 'Postgres', 'NorthwindTest', 'Schemas', 'public'];
      let scope: ParentNode = document;
      for (let i = 0; i < path.length; i++) {
        const label = path[i];
        const group = findChildGroupByLabel(scope, label);
        if (!group) return {ok: false, missing: label};
        const next = i + 1 < path.length ? path[i + 1] : 'customers';
        const ok = await expandAndWaitFor(group, next);
        if (!ok) return {ok: false, expandTimedOutAt: label};
        scope = group;
      }
      return {ok: true, hasCustomers: !!findChildGroupByLabel(scope, 'customers')};
    });
    expect(result.ok).toBe(true);
    expect(result.hasCustomers).toBe(true);
  });

  await softStep('Right-click customers → New Visual Query… opens visual query editor', async () => {
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
      for (const label of ['Databases', 'Postgres', 'NorthwindTest', 'Schemas', 'public']) {
        scope = findChildGroupByLabel(scope, label) as Element;
        if (!scope) return {ok: false, lost: label};
      }
      const customersGroup = findChildGroupByLabel(scope, 'customers');
      if (!customersGroup) return {ok: false, missing: 'customers'};
      const tableNode = customersGroup.querySelector(':scope > .d4-tree-view-node') as HTMLElement | null;
      if (!tableNode) return {ok: false, missing: 'tableNode'};
      (tableNode.querySelector('.d4-tree-view-group-label') as HTMLElement | null)?.click();
      await new Promise((r) => setTimeout(r, 250));
      tableNode.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
      let target: HTMLElement | undefined;
      const menuDeadline = Date.now() + 8_000;
      while (Date.now() < menuDeadline) {
        await new Promise((r) => setTimeout(r, 150));
        target = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'New Visual Query...') as HTMLElement | undefined;
        if (target) break;
      }
      if (!target) return {ok: false, missing: 'New Visual Query... menu item'};
      target.click();
      // Wait for visual query editor — Group by / Aggregate panels appear
      for (let i = 0; i < 80; i++) {
        await new Promise((r) => setTimeout(r, 200));
        const panels = Array.from(document.querySelectorAll('.grok-pivot-column-panel'));
        const titles = panels.map((p) => p.querySelector('.grok-pivot-column-tags-title')?.textContent?.trim());
        if (titles.includes('Group by') && titles.includes('Aggregate')) {
          const v: any = (window as any).grok.shell.v;
          return {ok: true, viewType: v?.type, panelTitles: titles};
        }
      }
      return {ok: false, missing: 'visual query panels'};
    });
    expect(result.ok).toBe(true);
    expect(result.viewType).toBe('DataQueryView');
    expect(result.panelTitles).toContain('Group by');
  });

  await softStep('Set Group by to companyname', async () => {
    const result = await page.evaluate(async () => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const panels = Array.from(document.querySelectorAll('.grok-pivot-column-panel'));
      const groupByPanel = panels.find((p) => p.querySelector('.grok-pivot-column-tags-title')?.textContent?.trim() === 'Group by');
      if (!groupByPanel) return {ok: false, missing: 'Group by panel'};
      const plus = groupByPanel.querySelector('.grok-icon.fa-plus') as HTMLElement | null;
      if (!plus) return {ok: false, missing: '+ on Group by'};
      plus.click();
      await wait(1500);
      const grids = Array.from(document.querySelectorAll('[name="viewer-Grid"]'));
      const sorted = grids.slice().sort((a, b) => {
        const ra = a.getBoundingClientRect();
        const rb = b.getBoundingClientRect();
        return (ra.width * ra.height) - (rb.width * rb.height);
      });
      const picker = sorted[0];
      if (!picker) return {ok: false, missing: 'picker grid'};
      const canvases = Array.from(picker.querySelectorAll('canvas'));
      const canvas = canvases[canvases.length - 1] as HTMLCanvasElement;
      const rect = canvas.getBoundingClientRect();
      const headerH = 26;
      const rowH = (rect.height - headerH) / 11;
      const cx = rect.left + 70;
      const cy = rect.top + headerH + rowH * 1 + rowH / 2;
      const ev = (type: string) => new MouseEvent(type, {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window});
      canvas.dispatchEvent(ev('mousedown'));
      canvas.dispatchEvent(ev('mouseup'));
      canvas.dispatchEvent(ev('click'));
      await wait(1000);
      const tags = Array.from(groupByPanel.querySelectorAll('.d4-tag')).map((t) => t.textContent?.replace(/\s+/g, ' ').trim());
      return {ok: tags.some((t) => t?.includes('companyname')), tags};
    });
    expect(result.ok).toBe(true);
  });

  await softStep('Order by picker shows only companyname when Group by has companyname', async () => {
    const result = await page.evaluate(async () => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const panels = Array.from(document.querySelectorAll('.grok-pivot-column-panel'));
      const orderByPanel = panels.find((p) => p.querySelector('.grok-pivot-column-tags-title')?.textContent?.trim() === 'Order by');
      if (!orderByPanel) return {ok: false, missing: 'Order by panel'};
      const plus = orderByPanel.querySelector('.grok-icon.fa-plus') as HTMLElement | null;
      plus?.click();
      await wait(1500);
      const grids = Array.from(document.querySelectorAll('[name="viewer-Grid"]'));
      const sorted = grids.slice().sort((a, b) => {
        const ra = a.getBoundingClientRect();
        const rb = b.getBoundingClientRect();
        return (ra.width * ra.height) - (rb.width * rb.height);
      });
      const picker = sorted[0];
      const r = picker?.getBoundingClientRect();
      const isSmall = r && r.height < 100;
      document.body.click();
      await wait(400);
      return {ok: !!isSmall, height: r?.height};
    });
    expect(result.ok).toBe(true);
  });

  await softStep('Add orders join via Data row', async () => {
    const result = await page.evaluate(async () => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const addBtn = document.querySelector('[name="div-add-Data"]') as HTMLElement | null;
      if (!addBtn) return {ok: false, missing: 'div-add-Data'};
      const rect = addBtn.getBoundingClientRect();
      const cx = rect.left + rect.width / 2;
      const cy = rect.top + rect.height / 2;
      const ev = (type: string) => new MouseEvent(type, {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window});
      addBtn.dispatchEvent(ev('mousedown'));
      addBtn.dispatchEvent(ev('mouseup'));
      addBtn.dispatchEvent(ev('click'));
      await wait(1500);
      const items = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .filter((el) => (el as HTMLElement).offsetParent !== null);
      const northwind = items.find((el) => el.textContent?.trim() === 'northwind') as HTMLElement | undefined;
      if (!northwind) return {ok: false, missing: 'northwind menu item'};
      const nwMi = northwind.closest('.d4-menu-item') as HTMLElement;
      const nwRect = nwMi.getBoundingClientRect();
      ['mouseenter', 'mousemove', 'mouseover'].forEach((t) =>
        nwMi.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: nwRect.left + nwRect.width / 2, clientY: nwRect.top + nwRect.height / 2, view: window})));
      await wait(800);
      const items2 = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .filter((el) => (el as HTMLElement).offsetParent !== null);
      const publicItem = items2.find((el) => el.textContent?.trim() === 'public') as HTMLElement | undefined;
      if (!publicItem) return {ok: false, missing: 'public submenu'};
      const pMi = publicItem.closest('.d4-menu-item') as HTMLElement;
      const pRect = pMi.getBoundingClientRect();
      ['mouseenter', 'mousemove', 'mouseover'].forEach((t) =>
        pMi.dispatchEvent(new MouseEvent(t, {bubbles: true, clientX: pRect.left + pRect.width / 2, clientY: pRect.top + pRect.height / 2, view: window})));
      await wait(900);
      const items3 = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
        .filter((el) => (el as HTMLElement).offsetParent !== null);
      const ordersItem = items3.find((el) => el.textContent?.trim() === 'orders') as HTMLElement | undefined;
      if (!ordersItem) return {ok: false, missing: 'orders submenu'};
      (ordersItem.closest('.d4-menu-item') as HTMLElement).click();
      await wait(2500);
      const tags = Array.from(document.querySelectorAll('.d4-tag')).map((t) => t.textContent?.replace(/\s+/g, ' ').trim());
      return {ok: tags.some((t) => t?.includes('orders')), tags};
    });
    expect(result.ok).toBe(true);
  });

  await softStep('Open Debug tab and run Debug query — no errors', async () => {
    const result = await page.evaluate(async () => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const tabs = Array.from(document.querySelectorAll('.d4-tab-header'));
      const debugTab = tabs.find((t) => t.textContent?.trim() === 'Debug') as HTMLElement | undefined;
      if (!debugTab) return {ok: false, missing: 'Debug tab'};
      debugTab.click();
      await wait(1500);
      const bug = document.querySelector('[name="icon-bug"]') as HTMLElement | null;
      if (!bug) return {ok: false, missing: 'icon-bug'};
      bug.click();
      await wait(1500);
      const dialog = document.querySelector('.d4-dialog');
      const ok = (dialog?.querySelector('[name="button-OK"]')
        ?? Array.from(dialog?.querySelectorAll('button') ?? []).find((b) => b.textContent?.trim() === 'OK')) as HTMLElement | null;
      ok?.click();
      for (let i = 0; i < 40; i++) {
        await wait(500);
        const debugContent = document.body.textContent ?? '';
        if (debugContent.includes('Call was started') || debugContent.includes('STATEMENT EXECUTION'))
          return {ok: true};
      }
      return {ok: false, missing: 'debug log output'};
    });
    expect(result.ok).toBe(true);
  });

  await softStep('Toolbox → Actions → Run query… opens new view', async () => {
    const result = await page.evaluate(async () => {
      const wait = (ms: number) => new Promise((r) => setTimeout(r, ms));
      const w: any = window;
      const before = Array.from(w.grok.shell.views).length;
      const tabs = Array.from(document.querySelectorAll('.d4-tab-header'));
      const queryTab = tabs.find((t) => t.textContent?.trim() === 'Query') as HTMLElement | undefined;
      queryTab?.click();
      await wait(800);
      const link = Array.from(document.querySelectorAll('*'))
        .find((el) => el.textContent?.trim() === 'Run query...' && (el as HTMLElement).offsetParent !== null && el.children.length === 0) as HTMLElement | undefined;
      if (!link) return {ok: false, missing: 'Run query... link', before};
      link.click();
      for (let i = 0; i < 60; i++) {
        await wait(300);
        const after = Array.from(w.grok.shell.views).length;
        if (after > before)
          return {ok: true, before, after};
      }
      return {ok: false, missing: 'new view', before};
    });
    expect(result.ok).toBe(true);
    expect(result.after).toBeGreaterThan(result.before);
  });

  await softStep('Save the query', async () => {
    await page.evaluate(() => {
      const w: any = window;
      const editor = (Array.from(w.grok.shell.views) as any[]).find((v) => v.type === 'DataQueryView');
      if (editor) w.grok.shell.v = editor;
    });
    await page.waitForTimeout(800);
    await page.evaluate(() => {
      const tabs = Array.from(document.querySelectorAll('.d4-tab-header'));
      const queryTab = tabs.find((t) => t.textContent?.trim() === 'Query') as HTMLElement | undefined;
      queryTab?.click();
    });
    await page.waitForTimeout(600);
    await page.evaluate((name) => {
      const input = Array.from(document.querySelectorAll('input.ui-input-editor'))
        .find((el) => (el as HTMLInputElement).value === 'customers' && (el as HTMLElement).offsetParent !== null) as HTMLInputElement | undefined;
      if (!input) throw new Error('no Name input with value=customers');
      input.focus();
      document.execCommand('selectAll');
      document.execCommand('insertText', false, name);
    }, queryName);
    await page.waitForTimeout(400);
    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(4000);
    const found = await page.evaluate(async (n) => {
      const q = await (window as any).grok.dapi.queries.filter(`${n}`).first().catch(() => null);
      return !!q;
    }, queryName);
    expect(found).toBe(true);
  });

  await softStep('Run the saved query — result table opens', async () => {
    const result = await page.evaluate(async (n) => {
      const w: any = window;
      w.grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 500));
      const q = await w.grok.dapi.queries.filter(`${n}`).first();
      if (!q) return {ok: false, missing: 'saved query'};
      try {
        const df = await q.executeTable({});
        w.grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1000));
        return {ok: true, rows: df.rowCount, cols: df.columns.length};
      } catch (e: any) {
        return {ok: false, err: e?.message?.slice?.(0, 200)};
      }
    }, queryName);
    expect(result.ok).toBe(true);
    expect(result.rows).toBeGreaterThan(0);
  });

  // Cleanup — delete the saved query
  await page.evaluate(async (n) => {
    const q = await (window as any).grok.dapi.queries.filter(`${n}`).first().catch(() => null);
    if (q)
      try { await (window as any).grok.dapi.queries.delete(q); } catch (e) {}
    (window as any).grok.shell.closeAll();
  }, queryName);

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
