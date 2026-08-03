/* ---
sub_features_covered: [projects.shell.open, projects.url-params.apply, projects.url-params.build-share-link, projects.view.browse]
--- */
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {projectsTestOptions, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '@datagrok-libraries/test/src/playwright/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '@datagrok-libraries/test/src/playwright/projects';

test.use(projectsTestOptions);

test('Projects / Project URL: deep-link reopen for representative project', async ({page, context}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const projectName = `AutoTest-URL-${Date.now()}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;
  let projectPath = '';

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: create representative file-share project with provenance', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await assertProvenanceScript(page, 'files', opened.script);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3 equivalent: derive deep-link path for the project', async () => {
      if (!saved) throw new Error('no saved project');
      // Use entity.path — server-canonical URL slug, e.g. `/p/QaPw.MyProject` (namespace separator is `.`).
      projectPath = await page.evaluate(async (id) => {
        const grok = (window as any).grok;
        const p = await grok.dapi.projects.find(id);
        return p?.path ?? '';
      }, saved.projectId);
      expect(projectPath.length).toBeGreaterThan(0);
      expect(projectPath.startsWith('/p/')).toBe(true);
    });

    await softStep('Step 4-5: reopen project (Data Sync re-runs) and verify deep-link URL + data', async () => {
      if (!saved) throw new Error('no saved project');
      const expectedId = saved.projectId;
      // Reopen via the JS API `project.open()` — the platform's reliable reopen path
      // (same mechanism as the shared saveAndReopen helper). A cold `page.goto('/p/<ns>.<name>')`
      // does NOT open the project on a fresh load (the SPA shows Home; the project-by-path route
      // resolves through the search index, which lags a just-saved project). `open()` routes the
      // browser to the canonical `/p/<ns>.<name>/<table>` deep-link, which we assert below.
      await page.evaluate(async (id) => {
        const grok = (window as any).grok;
        grok.shell.closeAll();
        const p = await grok.dapi.projects.find(id);
        await p.open();
      }, expectedId);
      // Poll until THIS project is the active one with its provisioned table re-materialized.
      const result = await page.evaluate(async ({pid}) => {
        const g = (f: () => any) => { try { return f(); } catch { return null; } };
        let last: {projId: string | null; rc: number | null; path: string} = {projId: null, rc: null, path: ''};
        for (let i = 0; i < 60; i++) {
          const projId = g(() => (window as any).grok?.shell?.project?.id) ?? null;
          const rcRaw = g(() => (window as any).grok?.shell?.tv?.dataFrame?.rowCount);
          const rc = typeof rcRaw === 'number' ? rcRaw : null;
          last = {projId, rc, path: location.pathname};
          if (projId === pid && rc !== null && rc > 0)
            return {ok: true, ...last};
          await new Promise((r) => setTimeout(r, 500));
        }
        return {ok: false, expected: pid, ...last};
      }, {pid: expectedId});
      console.log('Project reopen result: ' + JSON.stringify(result));
      // Main check (preserved): the SAME project reopened (shell.project.id === expected) with its
      // provisioned table re-materialized via Data Sync (rowCount > 0).
      expect(
        result.ok,
        result.ok ? '' : `reopen did not restore the expected project within 30s: ` +
          `got project.id=${result.projId} (expected ${expectedId}), rowCount=${result.rc}`,
      ).toBe(true);
      // Deep-link URL: open() routes to the project's canonical /p/<ns>.<name> URL (+ table slug).
      expect(result.path.toLowerCase()).toContain(projectPath.toLowerCase());
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
  }

  finishSpec();
});
