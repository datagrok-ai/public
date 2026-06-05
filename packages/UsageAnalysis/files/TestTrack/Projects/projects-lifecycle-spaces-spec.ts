/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.namespaces, projects.api.get-by-id]
related_bugs: [GROK-18345]
generated_from: projects-lifecycle-spaces.md
--- */
// Spaces-source lifecycle. Provisions a transient root Space, copies
// demog.csv into it (so `Spaces:<spaceName>/demog.csv` becomes openable),
// then runs the canonical upload + reopen flow with Data Sync ON. On
// reopen the persisted `.script` re-runs `OpenFile("Spaces.<...>/demog.csv")`,
// re-reading from the still-alive Space — the proactive coverage cell
// `source_class=spaces × dep_lifecycle_op=save_with_sync_on`.
//
// GROK-18345 recipient-open invariant (share + datasync under a different
// user identity) is now wired via the `shareWithSecondUserAndVerify` helper
// (token-based, needs DATAGROK_AUTH_TOKEN_2). It shares the project to the
// second user's GROUP at both View-and-Use + Full access and verifies the
// recipient can see it after re-auth; when no second-user token is configured
// it degrades to recipientVisible: null (no hard fail). Runs as the last step
// before cleanup so the owner session is restored before delete.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
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

  // Provision a Space with demog.csv physically copied into its storage.
  // The Space stays alive for the duration of the test — Sync-ON reopen
  // re-runs `OpenFile("Spaces.<spaceName>/demog.csv")`, so the file must
  // still be there at reopen time.
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

    // GROK-18345 recipient-open leg (View-and-Use + Full). MUST be the last
    // action before the finally cleanup — the helper reloads the page during
    // re-auth and restores the owner session before returning, so the delete
    // in finally runs as the owner. Defensive skip when no second-user token.
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

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
