import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors, baseUrl} from '../spec-login';

test.use(specTestOptions);

test('Queries — query layout (Layout tab, preview, project save/restore)', async ({page}) => {
  test.setTimeout(420_000);

  await loginToDatagrok(page);

  // Setup: Tabs mode + selenium class + closeAll. Resolve query/connection IDs
  // up front. The PostgresAll query lives on connection PostgresTest
  // (friendlyName: NorthwindTest, dataSource: Postgres, namespace: Dbtests).
  const ids = await page.evaluate(async () => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.settings.showFiltersIconsConstantly = true;
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
    const q = await (window as any).grok.dapi.queries
      .filter('name = "PostgresAll" and connection.friendlyName = "NorthwindTest"').first()
      .catch(() => null);
    return {
      queryId: q?.id ?? 'd625c7b0-015b-5b8d-b06d-59ffe8c0e3c4',
      connId: q?.connection?.id ?? 'a2d74603-7594-56ea-a2bd-844b2fd16ee7',
    };
  });
  expect(ids.queryId).toBeTruthy();

  // Pre-clean any leftover project from prior runs.
  await page.evaluate(async () => {
    try {
      const list = await (window as any).grok.dapi.projects
        .filter('name in ("query-layout-test-claude", "QueryLayoutTestClaude")').list();
      for (const p of list)
        try { await (window as any).grok.dapi.projects.delete(p); } catch (e) {}
    } catch (e) {}
  });

  await softStep('Open PostgresAll editor (Browse → Databases → Postgres → NorthwindTest → Edit...)', async () => {
    // Navigate to the connection's queries gallery (`browse=db`), then dispatch
    // contextmenu on the PostgresAll card and click Edit... — same DataQueryView
    // and Toolbox state as the Browse-tree right-click flow, but reachable via URL.
    await page.goto(`${baseUrl}/queries/Dbtests.PostgresTest?browse=db`);
    await page.locator('[name="div-PostgresAll"]').waitFor({timeout: 30_000});
    const opened = await page.evaluate(async () => {
      const card = document.querySelector('[name="div-PostgresAll"]') as HTMLElement | null;
      if (!card) return {ok: false, stage: 'no-card'};
      const rect = card.getBoundingClientRect();
      card.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2, buttons: 2,
        clientX: rect.left + rect.width / 2, clientY: rect.top + rect.height / 2,
      }));
      let editLabel: HTMLElement | undefined;
      for (let i = 0; i < 30; i++) {
        editLabel = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
          .find((el) => el.textContent?.trim() === 'Edit...') as HTMLElement | undefined;
        if (editLabel) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      if (!editLabel) return {ok: false, stage: 'no-edit-item'};
      (editLabel.closest('.d4-menu-item') as HTMLElement).click();
      // Wait for editor view AND CodeMirror to mount on the Query tab — this
      // signal indicates Dart has finished wiring the editor state. Without it,
      // a subsequent Layout tab click renders blank in fresh contexts.
      for (let i = 0; i < 80; i++) {
        const v = (window as any).grok.shell.v;
        const cm = (document.querySelector('.CodeMirror') as any)?.CodeMirror;
        if (v?.type === 'DataQueryView' && v?.name === 'PostgresAll'
            && document.querySelector('.d4-tab-header[name="Layout"]')
            && cm && cm.getValue?.().length > 0)
          return {ok: true, viewName: v?.name, viewType: v?.type};
        await new Promise((r) => setTimeout(r, 300));
      }
      const v = (window as any).grok.shell.v;
      return {ok: false, stage: 'no-editor-view', viewType: v?.type, viewName: v?.name};
    });
    expect(opened.ok, JSON.stringify(opened)).toBe(true);
  });

  await softStep('Open the Layout tab', async () => {
    await page.locator('.d4-tab-header[name="Layout"]').first().click({force: true});
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('.d4-tab-header[name="Layout"].selected')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      for (let i = 0; i < 80; i++) {
        const link = Array.from(document.querySelectorAll('label.d4-link-label'))
          .find((el) => el.textContent?.trim() === 'Run query'
            && (el as HTMLElement).offsetParent !== null);
        if (link) return {ok: true};
        if (document.querySelector('[name="viewer-Grid"] canvas'))
          return {ok: true};
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-run-link-after-tab'};
    });
    expect(result.ok, JSON.stringify(result)).toBe(true);
  });

  await softStep('Click Run query — populate the Layout preview grid', async () => {
    const result = await page.evaluate(async () => {
      // Layout tab shows a "Run query" link as label.d4-link-label
      let link: HTMLElement | undefined;
      for (let i = 0; i < 30; i++) {
        link = Array.from(document.querySelectorAll('label.d4-link-label'))
          .find((el) => el.textContent?.trim() === 'Run query'
            && (el as HTMLElement).offsetParent !== null) as HTMLElement | undefined;
        if (link) break;
        await new Promise((r) => setTimeout(r, 300));
      }
      if (!link) return {ok: false, stage: 'no-run-link'};
      // Synthesize a mouse click — labels need bubbling click event
      const rect = link.getBoundingClientRect();
      const eventInit: MouseEventInit = {
        bubbles: true, cancelable: true, view: window, button: 0,
        clientX: rect.left + 5, clientY: rect.top + 5,
      };
      link.dispatchEvent(new MouseEvent('mousedown', eventInit));
      link.dispatchEvent(new MouseEvent('mouseup', eventInit));
      link.dispatchEvent(new MouseEvent('click', eventInit));
      for (let i = 0; i < 80; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) {
          await new Promise((r) => setTimeout(r, 1500));
          return {ok: true, rows: (window as any).grok.shell.v?.dataFrame?.rowCount};
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-grid'};
    });
    expect(result.ok, JSON.stringify(result)).toBe(true);
    expect(result.rows).toBeGreaterThan(0);
  });

  await softStep('Add four viewers (Histogram, Bar chart, Scatter plot, Pie chart)', async () => {
    const viewers = await page.evaluate(async () => {
      const click = (sel: string) => (document.querySelector(sel) as HTMLElement | null)?.click();
      click('[name="icon-histogram"]');
      await new Promise((r) => setTimeout(r, 1200));
      click('[name="icon-bar-chart"]');
      await new Promise((r) => setTimeout(r, 1200));
      click('[name="icon-scatter-plot"]');
      await new Promise((r) => setTimeout(r, 1200));
      click('[name="icon-pie-chart"]');
      await new Promise((r) => setTimeout(r, 1200));
      return Array.from(document.querySelectorAll('[name^="viewer-"]')).map((e) => e.getAttribute('name'));
    });
    expect(viewers).toEqual(expect.arrayContaining([
      'viewer-Grid', 'viewer-Histogram', 'viewer-Bar-chart', 'viewer-Scatter-plot', 'viewer-Pie-chart',
    ]));
  });

  await softStep('Save the query (with layout)', async () => {
    await page.locator('[name="button-Save"]').first().click();
    await page.waitForTimeout(2500);
  });

  await softStep('Close All', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(1500);
    const home = await page.evaluate(() => (window as any).grok.shell.v?.name);
    expect(home).toBe('Home');
  });

  await softStep('Open PostgresAll preview — saved layout applied', async () => {
    // Substitute for "click in Browse → preview" — `/func/<nqName>` is the
    // direct URL that opens a query as a TableView preview with its saved
    // layout applied (same code path as Browse-tree click).
    await page.goto(`${baseUrl}/func/Dbtests.PostgresAll`);
    // Page reload — wait for grok API to remount and Browse sidebar tab to appear
    await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
    await page.waitForFunction(
      () => typeof (window as any).grok?.shell?.v !== 'undefined',
      null, {timeout: 30_000});
    const result = await page.evaluate(async () => {
      for (let i = 0; i < 120; i++) {
        const v = (window as any).grok?.shell?.v;
        if (v?.name === 'PostgresAll' && v?.type === 'TableView') {
          if (document.querySelector('[name="viewer-Histogram"]')
              && document.querySelector('[name="viewer-Bar-chart"]')
              && document.querySelector('[name="viewer-Scatter-plot"]')
              && document.querySelector('[name="viewer-Pie-chart"]')) {
            await new Promise((r) => setTimeout(r, 2000));
            const viewers = Array.from(document.querySelectorAll('[name^="viewer-"]'))
              .map((e) => e.getAttribute('name'));
            return {ok: true, viewers};
          }
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      const v = (window as any).grok?.shell?.v;
      const viewers = Array.from(document.querySelectorAll('[name^="viewer-"]'))
        .map((e) => e.getAttribute('name'));
      return {ok: false, stage: 'no-layout-viewers', viewType: v?.type, viewName: v?.name, viewers};
    });
    expect(result.ok, JSON.stringify(result)).toBe(true);
    expect(result.viewers).toEqual(expect.arrayContaining([
      'viewer-Grid', 'viewer-Histogram', 'viewer-Bar-chart', 'viewer-Scatter-plot', 'viewer-Pie-chart',
    ]));
  });

  await softStep('Run the query (open in new view) — saved layout applied', async () => {
    // Substitute for "right-click PostgresAll → Run" — execute the query and
    // add a fresh TableView, then apply the saved layout.
    const result = await page.evaluate(async (queryId) => {
      const q = await (window as any).grok.dapi.queries.find(queryId);
      if (!q) return {ok: false, stage: 'no-query'};
      const df = await q.executeTable();
      const tv = (window as any).grok.shell.addTableView(df);
      // Apply saved query layout
      const layout = q.layout ?? (await (window as any).grok.dapi.layouts.getApplicable?.(df))?.[0];
      if (layout) tv.loadLayout(layout);
      // Wait for layout viewers to mount
      for (let i = 0; i < 80; i++) {
        if (document.querySelector('[name="viewer-Histogram"]')
            && document.querySelector('[name="viewer-Bar-chart"]')
            && document.querySelector('[name="viewer-Scatter-plot"]')
            && document.querySelector('[name="viewer-Pie-chart"]')) {
          await new Promise((r) => setTimeout(r, 2000));
          const viewers = Array.from(document.querySelectorAll('[name^="viewer-"]'))
            .map((e) => e.getAttribute('name'));
          return {ok: true, viewers};
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      const viewers = Array.from(document.querySelectorAll('[name^="viewer-"]'))
        .map((e) => e.getAttribute('name'));
      return {ok: false, stage: 'no-layout-viewers', viewers};
    }, ids.queryId);
    expect(result.ok, JSON.stringify(result)).toBe(true);
    expect(result.viewers).toEqual(expect.arrayContaining([
      'viewer-Grid', 'viewer-Histogram', 'viewer-Bar-chart', 'viewer-Scatter-plot', 'viewer-Pie-chart',
    ]));
  });

  await softStep('Add two more viewers (Box plot, Tree map)', async () => {
    const viewers = await page.evaluate(async () => {
      (document.querySelector('[name="icon-box-plot"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1500));
      (document.querySelector('[name="icon-tree-map"]') as HTMLElement)?.click();
      await new Promise((r) => setTimeout(r, 1500));
      return Array.from(document.querySelectorAll('[name^="viewer-"]')).map((e) => e.getAttribute('name'));
    });
    expect(viewers).toEqual(expect.arrayContaining(['viewer-Box-plot', 'viewer-Tree-map']));
  });

  let projectId: string | null = null;

  await softStep('Save the project — SAVE → name → OK (Cancel Share dialog)', async () => {
    await page.locator('[name="button-Save"]').first().click();
    await page.waitForSelector('.d4-dialog', {timeout: 10_000});

    // Diagnose what dialog appeared and find the Name input
    const dialogInfo = await page.evaluate(() => {
      const dlg = document.querySelector('.d4-dialog');
      const title = dlg?.querySelector('.d4-dialog-title, .d4-dialog-header')?.textContent?.trim();
      const inputs = Array.from(dlg?.querySelectorAll('input, textarea') ?? [])
        .map((i: any) => ({type: i.type, placeholder: i.placeholder, value: i.value, name: i.name}));
      return {title, inputs};
    });

    // The Save Project dialog has the project name input. Find by current value
    // matching "PostgresAll" (the query name auto-fills the project name).
    const nameInputIdx = dialogInfo.inputs.findIndex((i) => i.value === 'PostgresAll');
    if (nameInputIdx < 0) {
      // Diagnostic — surface dialog state on failure
      expect(nameInputIdx, JSON.stringify(dialogInfo)).toBeGreaterThanOrEqual(0);
    }

    const nameInput = page.locator('.d4-dialog input, .d4-dialog textarea').nth(nameInputIdx);
    await nameInput.click();
    await page.keyboard.press('Control+a');
    await page.keyboard.type('query-layout-test-claude');

    await page.locator('.d4-dialog [name="button-OK"]').click();

    // The OK click can leave a Share dialog or other modal. Cancel any extra dialogs.
    for (let i = 0; i < 5; i++) {
      await page.waitForTimeout(800);
      const dlg = page.locator('.d4-dialog');
      if (!(await dlg.isVisible({timeout: 500}).catch(() => false))) break;
      const cancelBtn = page.locator('.d4-dialog [name="button-CANCEL"]');
      if (await cancelBtn.isVisible({timeout: 500}).catch(() => false))
        await cancelBtn.click().catch(() => {});
      else
        break;
    }

    const lookup = await page.evaluate(async () => {
      for (let i = 0; i < 80; i++) {
        // Try multiple name variants the server may normalize to
        for (const candidate of ['QueryLayoutTestClaude', 'query-layout-test-claude', 'Query-layout-test-claude']) {
          const p = await (window as any).grok.dapi.projects
            .filter(`name = "${candidate}"`).first().catch(() => null);
          if (p) return {id: p.id, name: p.name};
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      // Fallback: list recent projects with similar names
      const list = await (window as any).grok.dapi.projects
        .filter('query-layout-test-claude').list().catch(() => []);
      return {id: null, recentNames: list.slice(0, 10).map((p: any) => p.name)};
    });
    projectId = lookup.id;

    expect(projectId, JSON.stringify(lookup)).toBeTruthy();
  });

  await softStep('Close All again', async () => {
    await page.evaluate(() => (window as any).grok.shell.closeAll());
    await page.waitForTimeout(1500);
    const home = await page.evaluate(() => (window as any).grok.shell.v?.name);
    expect(home).toBe('Home');
  });

  await softStep('Open the project — layout restored', async () => {
    const result = await page.evaluate(async (id) => {
      if (!id) return {ok: false, stage: 'no-project-id'};
      const proj = await (window as any).grok.dapi.projects.find(id);
      if (!proj) return {ok: false, stage: 'no-project'};
      await proj.open();
      for (let i = 0; i < 80; i++) {
        const v = (window as any).grok.shell.v;
        if (v?.type === 'TableView' && v.name === 'PostgresAll') {
          await new Promise((r) => setTimeout(r, 2000));
          const viewers = Array.from(document.querySelectorAll('[name^="viewer-"]'))
            .map((e) => e.getAttribute('name'));
          return {ok: true, viewers};
        }
        await new Promise((r) => setTimeout(r, 300));
      }
      return {ok: false, stage: 'no-table-view-after-open'};
    }, projectId);
    expect(result.ok).toBe(true);
    expect(result.viewers).toEqual(expect.arrayContaining([
      'viewer-Grid', 'viewer-Histogram', 'viewer-Bar-chart', 'viewer-Scatter-plot',
      'viewer-Pie-chart', 'viewer-Box-plot', 'viewer-Tree-map',
    ]));
  });

  await softStep('Toolbox → Source → Refresh — layout unchanged', async () => {
    // Try the UI-driven Refresh button first; fall back to tv.reloadData()
    // when the Source toolbox section isn't visible (project-opened-via-JS-API
    // case, where the sidebar may default to Browse and not expose the Source
    // section's REFRESH button).
    await page.locator('.d4-tab-header[name="Toolbox"]').first().click({force: true})
      .catch(() => {});
    await page.waitForTimeout(1500);

    const result = await page.evaluate(async () => {
      const tv = (window as any).grok.shell.tv;
      if (!tv) return {ok: false, stage: 'no-tv'};
      const before = Array.from(tv.viewers).map((v: any) => v.type).sort();

      // Expand Source section if visible
      const sourceSection = document.querySelector('[name="div-section--Source"]') as HTMLElement | null;
      if (sourceSection) {
        const header = sourceSection.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
        const content = sourceSection.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
        if (header && content && getComputedStyle(content).display === 'none')
          header.click();
        await new Promise((r) => setTimeout(r, 500));
      }

      let refreshBtn: HTMLElement | undefined;
      for (let i = 0; i < 10; i++) {
        refreshBtn = Array.from(document.querySelectorAll('button'))
          .find((b) => b.textContent?.trim().toUpperCase() === 'REFRESH'
            && (b as HTMLElement).offsetParent !== null) as HTMLElement | undefined;
        if (refreshBtn) break;
        await new Promise((r) => setTimeout(r, 300));
      }
      let refreshMethod = 'ui-button';
      if (refreshBtn) {
        refreshBtn.click();
      } else {
        // Fallback: call reloadData() which is what the REFRESH button binds to
        refreshMethod = 'tv.reloadData()';
        if (typeof tv.reloadData === 'function') await tv.reloadData();
        else return {ok: false, stage: 'no-refresh-btn-or-method'};
      }
      await new Promise((r) => setTimeout(r, 8000));
      const after = Array.from((window as any).grok.shell.tv.viewers).map((v: any) => v.type).sort();
      return {ok: true, before, after, refreshMethod,
        unchanged: JSON.stringify(before) === JSON.stringify(after)};
    });
    expect(result.ok, JSON.stringify(result)).toBe(true);
    expect(result.unchanged, JSON.stringify(result)).toBe(true);
  });

  // Cleanup: delete the project we created.
  await page.evaluate(async (id) => {
    if (!id) return;
    try {
      const p = await (window as any).grok.dapi.projects.find(id);
      if (p) await (window as any).grok.dapi.projects.delete(p);
    } catch (e) {}
  }, projectId);

  if (stepErrors.length > 0)
    throw new Error('Soft step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
