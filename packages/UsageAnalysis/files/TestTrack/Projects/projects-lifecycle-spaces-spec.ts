// Spaces-source lifecycle: provision a transient root Space with demog.csv, save with Sync ON, reopen.
// GROK-18345: share to the second user's GROUP at View-and-Use + Full and verify recipient can open it.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {
  projectsTestOptions,
  gotoApp,
  setupSession,
  provisionSpaceFixture,
  releaseSpaceFixture,
  SpaceFixture,
} from './_helpers';
import {
  openTableFromFile,
  assertProvenanceScript,
  resetShell,
  PROVENANCE_PATTERNS,
} from '../helpers/openers';
import {
  saveAllTablesWithProvenance,
  reopenAndAssertProvenance,
  shareWithSecondUserAndVerify,
  deleteProjectWithCleanup,
  SavedAllTables,
} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Lifecycle Spaces: open from Space + sync + reopen', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `Test_LifecycleSpaces_${stamp}`;
  let saved: SavedAllTables | null = null;
  let fixture: SpaceFixture | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  // Provision a Space with demog.csv copied in; it stays alive so Sync-ON reopen can re-run OpenFile against it.
  const probe = await provisionSpaceFixture(page, {
    namePrefix: 'lifecycle-spaces',
    fileName: 'demog.csv',
  });
  if (probe.blocked || !probe.fixture) {
    console.warn('Spaces lifecycle env-skip — ' + probe.reason);
    test.skip(true, probe.reason);
    return;
  }
  fixture = probe.fixture;

  try {
    await softStep('open demog.csv from Space', async () => {
      const t1 = await openTableFromFile(page, fixture!.filePath);
      expect(t1.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t1.script);
    });

    await softStep('save with Data Sync ON', async () => {
      saved = await saveAllTablesWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoIds.length).toBeGreaterThanOrEqual(1);
    });

    await softStep('reopen re-runs OpenFile against the still-alive Space', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId, PROVENANCE_PATTERNS.files);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(1);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });

    // GROK-18345 recipient-open leg. MUST be last before cleanup — the helper restores the owner session before delete.
    await softStep('GROK-18345: share with second user + recipient-open verification', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await shareWithSecondUserAndVerify(page, {id: saved.projectId, name: projectName}, {full: true});
      if (!r.shared) {
        console.warn('GROK-18345 share skipped — ' + r.reason);
        return;
      }
      if (r.recipientVisible !== null)
        expect(r.recipientVisible).toBe(true);
    });
  } finally {
    if (saved) {
      await deleteProjectWithCleanup(page, {projectId: saved.projectId});
      for (const tableInfoId of saved.tableInfoIds)
        await deleteProjectWithCleanup(page, {tableInfoId});
    }
    if (fixture) await releaseSpaceFixture(page, fixture);
  }

  finishSpec();
});
