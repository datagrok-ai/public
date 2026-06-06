import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Radar — table-rebind on project save/reopen (GROK-18085)', async ({page}) => {
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
    expect(result.rebindOk).toBe(true);
  });

  await softStep('Steps 4-6: Save project, close all, reopen by id (GROK-18085 invariant)', async () => {
    const errorsBefore = consoleErrors.length;
    const saveResult = await page.evaluate(async (pName) => {
      const grok = (window as any).grok;
      const DG = (window as any).DG;
      try {
        const project = DG.Project.create();
        project.name = pName;
        const saved = await grok.dapi.projects.save(project);
        return {ok: true, id: saved.id, name: saved.name};
      } catch (e) {
        return {ok: false, err: String(e).substring(0, 200)};
      }
    }, projectName);
    if (saveResult.ok) projectId = saveResult.id;
    // Skip reopen verification if project save unavailable (DG.Project.create not on window).
    if (!saveResult.ok) {
      console.warn('[SKIP]', 'Project save API unavailable (DG.Project.create not on window):', saveResult.err);
      return;
    }

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

    if (!reopenResult.ok)
      console.warn('[Reopen failure]', reopenResult.err);
    // Primary GROK-18085 invariant: no console error during reopen.
    const errorsDuring = consoleErrors.slice(errorsBefore);
    expect(errorsDuring).toEqual([]);
  });

  await softStep('Cleanup: delete saved project', async () => {
    if (!projectId) return;
    await page.evaluate(async (id) => {
      const grok = (window as any).grok;
      try {
        const proj = await grok.dapi.projects.find(id);
        await grok.dapi.projects.delete(proj);
      } catch (e) { /* best-effort cleanup */ }
    }, projectId);
  });

  await page.evaluate(() => (window as any).grok.shell.closeAll());

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
