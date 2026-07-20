/* ---
sub_features_covered: [charts.echart-base.table, charts.radar]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {saveAllTablesWithProvenance, deleteProjectWithCleanup} from '@datagrok-libraries/test/src/playwright/projects';

test.use(specTestOptions);

test('Charts / Radar — table-rebind on project save/reopen (GROK-18085)', async ({page}) => {
  test.setTimeout(120_000);

  const consoleErrors: string[] = [];
  const isBenignError = (text: string) =>
    /Failed to load resource/.test(text) || /404 \(\)/.test(text) || /favicon/.test(text);
  page.on('console', (msg) => {
    if (msg.type() === 'error' && !isBenignError(msg.text())) consoleErrors.push(msg.text());
  });
  page.on('pageerror', (err) => consoleErrors.push(`pageerror: ${err.message}`));

  await loginToDatagrok(page);
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    (window as any).grok.shell.closeAll();
    document.body.classList.add('selenium');
    (window as any).grok.shell.windows.simpleMode = true;
  });

  let projectId: string | null = null;
  let savedTableInfoIds: string[] = [];
  let savedLayoutId: string | null = null;
  let demogName: string | null = null;
  const projectName = `radar-rebind-${Date.now()}`;

  await softStep('Steps 1-3: Open demog + SPGI; add Radar to SPGI; rebind table to demog', async () => {
    const result = await page.evaluate(async () => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 30_000, step = 150) => {
        const deadline = Date.now() + timeout;
        while (Date.now() < deadline) {
          try { if (pred()) return true; } catch (_) { /* not ready */ }
          await new Promise((r) => setTimeout(r, step));
        }
        return false;
      };
      const norm = (t: any) => (typeof t === 'string' ? t : (t?.name ?? null));
      const hasTable = (name: string) => { for (const t of grok.shell.tables) if (t.name === name) return true; return false; };
      const dfDemog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(dfDemog);
      await poll(() => hasTable(dfDemog.name));
      const dfSpgi = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tvSpgi = grok.shell.addTableView(dfSpgi);
      await poll(() => hasTable(dfSpgi.name));
      const radar = tvSpgi.addViewer('Radar');
      await poll(() => { for (const v of tvSpgi.viewers) if (v.type === 'Radar') return true; return false; });
      let rebindOk = false;
      let readBackTable = null;
      try {
        radar.setOptions({table: dfDemog.name});
        await poll(() => norm(radar.props.get('table')) === dfDemog.name);
        readBackTable = norm(radar.props.get('table'));
        rebindOk = true;
      } catch (e) {
        return {rebindOk: false, err: String(e).substring(0, 200)};
      }
      return {rebindOk, readBackTable, demogName: dfDemog.name, spgiName: dfSpgi.name};
    });
    expect(result.rebindOk, result.rebindOk ? '' : `Radar table rebind failed: ${(result as any).err}`).toBe(true);
    // The rebind must actually take effect: the Radar's table prop reads back the demog table (not SPGI).
    expect(result.readBackTable, `Radar table prop did not rebind to demog (got ${JSON.stringify(result.readBackTable)})`)
      .toBe(result.demogName);
    demogName = result.demogName;
  });

  await softStep('Steps 4-6: Save project, close all, reopen by id (GROK-18085 invariant)', async () => {
    const errorsBefore = consoleErrors.length;
    // Make the Radar's TableView active so its layout (carrying the SPGI->demog rebind under test)
    // is the one persisted by saveAllTablesWithProvenance.
    await page.evaluate(() => {
      const grok = (window as any).grok;
      for (const tv of grok.shell.tableViews) {
        for (const v of tv.viewers) if (v.type === 'Radar') { grok.shell.v = tv; return; }
      }
    });
    // Persist via the platform JS-API contract (upload each in-memory dataframe + addChild the
    // TableInfo + save the active view's layout, then projects.save). Saving grok.shell.project
    // directly fails server-side ("Unable to add entity") because readCsv tables are not uploaded.
    // Both demog and SPGI are persisted so the reopened project can re-materialize the rebind.
    let saved: {projectId: string; tableInfoIds: string[]} | null = null;
    try {
      saved = await saveAllTablesWithProvenance(page, projectName);
    } catch (e) {
      expect(false, `Project save failed: ${String(e).substring(0, 200)}`).toBe(true);
    }
    projectId = saved!.projectId;
    savedTableInfoIds = saved!.tableInfoIds;
    savedLayoutId = saved!.layoutId;
    expect(typeof projectId).toBe('string');
    expect((projectId as string).length).toBeGreaterThan(0);
    // Step 4 persists BOTH demog and SPGI so the reopen can re-materialize the rebind.
    expect(savedTableInfoIds.length).toBe(2);

    const reopenResult = await page.evaluate(async ({id, layoutId}) => {
      const grok = (window as any).grok;
      const poll = async (pred: () => boolean, timeout = 30_000, step = 200) => {
        const deadline = Date.now() + timeout;
        while (Date.now() < deadline) {
          try { if (pred()) return true; } catch (_) { /* not ready */ }
          await new Promise((r) => setTimeout(r, step));
        }
        return false;
      };
      const norm = (t: any) => (typeof t === 'string' ? t : (t?.name ?? null));
      const viewCount = () => { let n = 0; for (const _ of grok.shell.tableViews) n++; return n; };
      const hasRadar = () => { for (const view of grok.shell.tableViews) for (const v of view.viewers) if (v.type === 'Radar') return true; return false; };
      grok.shell.closeAll();
      await poll(() => viewCount() === 0, 15_000);
      try {
        const proj = await grok.dapi.projects.find(id);
        if (proj.open) await proj.open();
        // proj.open() re-materializes the tables + default Grid but does not auto-apply the
        // attached custom-viewer layout, so the Radar (carrying the rebind under test) is absent
        // until its saved layout is loaded explicitly onto the active TableView.
        if (!hasRadar() && layoutId) {
          try {
            const lay = await grok.dapi.layouts.find(layoutId);
            if (lay && grok.shell.tv?.loadLayout) grok.shell.tv.loadLayout(lay);
          } catch (_) { /* fall through to poll */ }
        }
        await poll(hasRadar, 30_000);
        const tv = grok.shell.tv;
        const types: string[] = [];
        let radarTable: string | null = null;
        let radarFound = false;
        for (const view of grok.shell.tableViews)
          for (const v of view.viewers) {
            types.push(v.type);
            if (v.type === 'Radar' && !radarFound) { radarFound = true; radarTable = norm(v.props.get('table')); }
          }
        return {ok: true, viewerTypes: types, hasTv: !!tv, radarFound, radarTable};
      } catch (e) {
        return {ok: false, err: String(e).substring(0, 300)};
      }
    }, {id: projectId as string, layoutId: savedLayoutId});

    // The reopen is the operation under test — a thrown reopen must fail, not just warn.
    expect(reopenResult.ok, reopenResult.ok ? '' : `Project reopen failed: ${(reopenResult as any).err}`).toBe(true);
    expect(reopenResult.hasTv, 'reopened project produced no active TableView').toBe(true);
    // Primary GROK-18085 invariant: the SPGI->demog Radar rebind survives save/reopen.
    expect(reopenResult.viewerTypes, 'reopened project has no Radar viewer').toContain('Radar');
    expect((reopenResult as any).radarTable, `reopened Radar not rebound to demog (got ${JSON.stringify((reopenResult as any).radarTable)})`)
      .toBe(demogName);
    // Secondary guard: the reported bug was a console error on reopen — assert none.
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Cleanup: delete saved project and uploaded tables', async () => {
    if (!projectId) return;
    await deleteProjectWithCleanup(page, {projectId});
    for (const tableInfoId of savedTableInfoIds)
      await deleteProjectWithCleanup(page, {tableInfoId});
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
