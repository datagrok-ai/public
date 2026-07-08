/* ---
sub_features_covered: [projects.shell.open, projects.url-params.apply, projects.url-params.build-share-link, projects.view.browse]
--- */
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {projectsTestOptions, BASE_URL, evalJs, gotoApp, setupSession} from './_helpers';
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

    await softStep('Step 4-5: navigate to project URL and verify project loads (Data Sync re-runs)', async () => {
      if (!saved) throw new Error('no saved project');
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(500);
      await page.goto(`${BASE_URL}${projectPath}`);
      await page.locator('[name="Browse"]').waitFor({timeout: 60_000});
      // Verify via project.id match (primary) or TableView rowCount > 0 (fallback). Avoid grok.shell.tables —
      // Dart-side Tn.grok_TableNames throws on dev for URL-opened projects.
      const expectedId = saved.projectId;
      const result = await page.evaluate(async ({pid}) => {
        const grok = (window as any).grok;
        let lastProjId: string | null = null;
        let lastRc: number | null = null;
        for (let i = 0; i < 90; i++) {
          const projId = grok?.shell?.project?.id;
          const rc = grok?.shell?.tv?.dataFrame?.rowCount;
          lastProjId = projId ?? null;
          lastRc = typeof rc === 'number' ? rc : null;
          if (projId === pid && typeof rc === 'number' && rc > 0)
            return {ok: true, signal: 'matched-id+rowCount', projId, rc};
          await new Promise((r) => setTimeout(r, 500));
        }
        return {ok: false, signal: 'timeout', projId: lastProjId, rc: lastRc, expected: pid};
      }, {pid: expectedId});
      console.log('Project URL load result: ' + JSON.stringify(result));
      // Deep-link contract: navigating to the project URL must open THIS project
      // (shell.project.id === expected) with its table re-materialized (rowCount > 0).
      expect(
        result.ok,
        result.ok ? '' : `deep-link did not open the expected project within 45s: ` +
          `got project.id=${result.projId} (expected ${expectedId}), rowCount=${result.rc}`,
      ).toBe(true);
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
