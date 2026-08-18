import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Notebooks / Create (UI Smoke): open demog in Notebook -> view as HTML', async ({page}) => {

  test.setTimeout(240_000);
  stepErrors.length = 0;

  let notebookId: string | null = null;

  await loginToDatagrok(page);
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
    (window as any).grok.shell.closeAll();
  });

  await page.waitForFunction(() => {
    try { return ((window as any).grok.shell.tableViews?.length ?? 0) === 0; }
    catch (e) { return false; }
  }, null, {timeout: 15_000}).catch(() => {});

  await softStep('Setup: open demog.csv and make it the active table', async () => {
    const info = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      const tv = grok.shell.addTableView(df);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 3000);
      });
      return {rows: df.rowCount, type: tv.type};
    });
    await page.locator('[name="viewer-Grid"]').waitFor({timeout: 60_000});
    expect(info.rows).toBeGreaterThan(0);
  });

  await softStep('Scenario 1: ML | Notebooks | Open in Notebook -> notebook view opens + persists', async () => {

    await page.locator('[name="div-ML"]').first().click({timeout: 10_000});

    const nbGroup = page.locator('[name="div-ML---Notebooks"]').first();
    await nbGroup.waitFor({state: 'attached', timeout: 8_000});
    await nbGroup.hover({timeout: 8_000});
    const leaf = page.locator('[name="div-ML---Notebooks---Open-in-Notebook"]').first();
    await leaf.waitFor({state: 'visible', timeout: 15_000});
    await leaf.click();

    const res = await page.waitForFunction(() => {
      const grok = (window as any).grok;
      let v; let o;
      try { v = grok.shell.v; o = grok.shell.o; } catch (e) { return false; } 
      if (v && v.type === 'Notebook' && o && o.id)
        return {viewName: v.name, oId: o.id};
      return false;
    }, null, {timeout: 90_000, polling: 250});
    const opened = await res.jsonValue() as {viewName: string; oId: string | null};

    expect(opened.viewName).toBeTruthy();
    const viewType = await page.evaluate(() => (window as any).grok.shell.v?.type);
    expect(viewType).toBe('Notebook');

    notebookId = opened.oId;
    const persisted = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      if (!id) return {found: false};
      for (let i = 0; i < 10; i++) {
        try {
          const nb = await grok.dapi.notebooks.find(id);
          if (nb && nb.id) return {found: true};
        } catch (e) {  }
        await new Promise((r) => setTimeout(r, 500));
      }
      return {found: false};
    }, notebookId);
    expect(persisted.found).toBe(true);

    const initErrors = await page.evaluate(() => {

      try { return (window as any).grok?.shell?.lastError ?? null; } catch (e) { return null; }
    });
    if (initErrors)
      console.warn(`[OBSERVE GROK-16296] notebook-init surfaced an error (expected while unfixed): ${initErrors}`);
  });

  await softStep('Scenario 2: click HTML in the notebook ribbon -> HTML mode is reachable', async () => {

    const htmlBtn = page.locator('[name="button-HTML"]').first();
    await htmlBtn.waitFor({timeout: 60_000, state: 'attached'});
    await htmlBtn.waitFor({timeout: 60_000, state: 'visible'});
    await htmlBtn.click();

    await expect(page.locator('[name="button-EDIT"]').first()).toBeVisible();
    const viewType = await page.evaluate(() => (window as any).grok.shell.v?.type);
    expect(viewType).toBe('Notebook');

    const contentChildren = await page.evaluate(() => {
      const el = document.querySelector('[data-source="Notebooks:Notebook"]');
      return el ? el.children.length : -1;
    });
    if (contentChildren <= 0)
      console.warn('[OBSERVE GROK-13999] HTML-mode content did not render (expected while unfixed); ' +
        'Download (As HTML / As PDF) is consequently unavailable — covered manually in create-ui.md.');
  });

  await softStep('Scenario 3: delete the notebook via Context Panel Actions > Delete -> removed from server', async () => {
    await page.evaluate(() => {
      const grok = (window as any).grok;
      try { grok.shell.windows.showContextPanel = true; } catch (e) {  }
    });

    await page.locator('[name="div-section--Actions"]').first().waitFor({state: 'attached', timeout: 30_000});
    await page.locator('[name="div-section--Actions"]').first().click().catch(() => {});

    await page.waitForFunction(() => {
      const sec = document.querySelector('[name="div-section--Actions"]');
      const pane = sec?.closest('.d4-accordion-pane');
      if (!pane) return false;
      return Array.from(pane.querySelectorAll('.d4-link-action'))
        .some((e) => (e.textContent ?? '').trim() === 'Delete');
    }, null, {timeout: 30_000}).catch(() => {});
    const clickedDelete = await page.evaluate(() => {
      const sec = document.querySelector('[name="div-section--Actions"]');
      const pane = sec?.closest('.d4-accordion-pane');
      const del = pane
        ? Array.from(pane.querySelectorAll('.d4-link-action')).find((e) => (e.textContent ?? '').trim() === 'Delete') as HTMLElement | undefined
        : undefined;
      if (del) { del.click(); return true; }
      return false;
    });
    expect(clickedDelete, 'Context Panel Actions > "Delete" link should be present and clickable').toBe(true);

    await page.evaluate(async () => {
      for (let i = 0; i < 60; i++) {
        const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => /are you sure|delete notebook/i.test(d.textContent || ''));
        if (dlg) { const yes = dlg.querySelector('[name="button-YES"]') as HTMLElement | null; if (yes) { yes.click(); return; } }
        await new Promise((r) => setTimeout(r, 100));
      }
    });

    const gone = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      if (!id) return false;
      for (let i = 0; i < 30; i++) {
        const ent: any = await grok.dapi.notebooks.find(id).catch(() => undefined);
        if (!ent || !ent.id) return true;
        await new Promise((r) => setTimeout(r, 500));
      }
      return false;
    }, notebookId);
    expect(gone, 'the deleted notebook should no longer resolve on the server').toBe(true);
  });

  await page.evaluate(async (id) => {
    const grok = (window as any).grok;
    try {
      if (id) {
        const nb = await grok.dapi.notebooks.find(id);
        if (nb) await grok.dapi.notebooks.delete(nb);
      }
      grok.shell.closeAll();
    } catch (e) {  }
  }, notebookId).catch(() => {});

  if (stepErrors.length > 0)
    throw new Error('Step failures:\n' + stepErrors.map((e) => `- ${e.step}: ${e.error}`).join('\n'));
});
