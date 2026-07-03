/* ---
sub_features_covered: [charts.echart-base.table, charts.radar]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {saveAllTablesWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(specTestOptions);

test('Charts / Radar — table-rebind on project save/reopen (GROK-18085)', async ({page}) => {
  test.setTimeout(300_000);

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
  const projectName = `radar-rebind-${Date.now()}`;

  await softStep('Steps 1-3: Open demog + SPGI; add Radar to SPGI; rebind table to demog', async () => {
    const result = await page.evaluate(async (pName) => {
      const grok = (window as any).grok;
      const dfDemog = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(dfDemog);
      await new Promise((r) => setTimeout(r, 1500));
      const dfSpgi = await grok.dapi.files.readCsv('System:DemoFiles/SPGI.csv');
      const tvSpgi = grok.shell.addTableView(dfSpgi);
      await new Promise((r) => setTimeout(r, 1500));
      const radar = tvSpgi.addViewer('Radar');
      await new Promise((r) => setTimeout(r, 3000));
      let rebindOk = false;
      let readBackTable = null;
      try {
        radar.setOptions({table: dfDemog.name});
        await new Promise((r) => setTimeout(r, 1500));
        readBackTable = radar.props.get('table');
        rebindOk = true;
      } catch (e) {
        return {rebindOk: false, err: String(e).substring(0, 200)};
      }
      return {rebindOk, readBackTable, demogName: dfDemog.name, spgiName: dfSpgi.name};
    }, projectName);
    expect(result.rebindOk, result.rebindOk ? '' : `Radar table rebind failed: ${(result as any).err}`).toBe(true);
    // The rebind must actually take effect: the Radar's table prop reads back the demog table
    // (prop may surface the name or the Table object — assert it resolved to the demog binding, not SPGI).
    const readBack = result.readBackTable as any;
    const readBackName = typeof readBack === 'string' ? readBack : (readBack?.name ?? readBack);
    expect(readBackName, `Radar table prop did not rebind to demog (got ${JSON.stringify(readBack)})`)
      .toBe(result.demogName);
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
    expect(projectId).toBeTruthy();

    const reopenResult = await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1500));
      try {
        const proj = await grok.dapi.projects.find(id);
        if (proj.open) await proj.open();
        await new Promise((r) => setTimeout(r, 3000));
        const tv = grok.shell.tv;
        const types: string[] = [];
        if (tv) for (const v of tv.viewers) types.push(v.type);
        return {ok: true, viewerTypes: types, hasTv: !!tv};
      } catch (e) {
        return {ok: false, err: String(e).substring(0, 300)};
      }
    }, projectId);

    // The reopen is the operation under test — a thrown reopen must fail, not just warn.
    expect(reopenResult.ok, reopenResult.ok ? '' : `Project reopen failed: ${(reopenResult as any).err}`).toBe(true);
    expect(reopenResult.hasTv, 'reopened project produced no active TableView').toBe(true);
    // Primary GROK-18085 invariant: no console error during reopen.
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
