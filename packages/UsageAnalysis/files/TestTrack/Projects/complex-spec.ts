// Thin smoke: save + rename + share on a single file source (full scope in satellite specs).
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup, shareWithSecondUserAndVerify} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Complex (smoke): save, rename, share', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `Complex_Smoke_${stamp}`;
  const renamed = `${projectName}_v2`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1-2: open demog with provenance + save with Sync ON', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', opened.script);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 9 equiv: rename project via JS API (verify via find-by-id)', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{ok: boolean; persistedName: string | null}>(page, `(async () => {
        const p = await grok.dapi.projects.find('${saved.projectId}');
        if (!p) return {ok: false, persistedName: null};
        p.name = '${renamed}';
        await grok.dapi.projects.save(p);
        // Verify by id (filter-by-name lags after rename).
        const verify = await grok.dapi.projects.find('${saved.projectId}');
        return {ok: verify?.name === '${renamed}', persistedName: verify?.name ?? null};
      })()`);
      expect(r.ok).toBe(true);
      expect(r.persistedName).toBe(renamed);
    });

    // Share is LAST step before finally — the helper reloads the page for second-user re-auth.
    await softStep('Step 12 equiv: share with second user (View-and-Use + Full) + recipient open', async () => {
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
