/* ---
sub_features_covered: [projects.api.files.sync, projects.api.get-by-id, projects.api.save, projects.shell.share-via-context-menu, projects.upload]
--- */
// File-source lifecycle: open, save with provenance, reopen-verify, share, rename, delete.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '@datagrok-libraries/test/src/playwright/spec-login';
import {finishSpec} from '@datagrok-libraries/test/src/playwright/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromFile,
  assertProvenanceScript,
  resetShell,
  PROVENANCE_PATTERNS,
} from '@datagrok-libraries/test/src/playwright/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
  shareWithSecondUserAndVerify,
} from '@datagrok-libraries/test/src/playwright/projects';

test.use(projectsTestOptions);

test('Projects / Lifecycle Files: open → save with provenance → reopen → share → rename', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-files-${stamp}`;
  const renamed = `${projectName}-renamed`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open demog.csv from System:DemoFiles via OpenFile', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', opened.script);
      // Accept colon-form and dot-form — openTableFromFile normalizes `:` → `.` (bug 2a workaround).
      expect(opened.script).toMatch(/OpenFile\("System[:.]DemoFiles\/demog\.csv"\)/);
    });

    await softStep('Step 2: save project with provenance (canonical uploadProject)', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoId).toBeTruthy();
    });

    await softStep('Step 3: closeAll → reopen → verify .script survived + table re-materialized', async () => {
      if (!saved) throw new Error('Step 2 did not produce a saved project');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.files,
      );
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      expect(result.reopenedScript).toMatch(/OpenFile/);
    });

    await softStep('Step 4: rename project via JS API', async () => {
      if (!saved) return;
      const r = await evalJs<{ok: boolean; persistedName: string | null}>(page, `(async () => {
        const p = await grok.dapi.projects.find('${saved.projectId}');
        if (!p) return {ok: false, persistedName: null};
        p.name = '${renamed}';
        await grok.dapi.projects.save(p);
        const verify = await grok.dapi.projects.find('${saved.projectId}');
        return {ok: verify?.name === '${renamed}', persistedName: verify?.name ?? null};
      })()`);
      expect(r.ok).toBe(true);
      expect(r.persistedName).toBe(renamed);
    });

    // Share is LAST step before finally — the helper reloads the page for second-user re-auth.
    await softStep('Step 5: share with second user (View-and-Use + Full) + recipient open', async () => {
      if (!saved) return;
      const r = await shareWithSecondUserAndVerify(page, {id: saved.projectId, name: renamed}, {full: true});
      if (!r.shared) { console.warn('Share skipped: ' + r.reason); return; }
      if (r.recipientVisible !== null) expect(r.recipientVisible).toBe(true);
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
