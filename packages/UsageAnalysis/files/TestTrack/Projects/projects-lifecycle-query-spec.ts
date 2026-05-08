/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.relations.list, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
related_bugs: [github-3550]
generated_from: projects-lifecycle-query.md
--- */
// Query-source lifecycle. Provisions a saved query on System:Datagrok
// (via helpers/openers.ts:provisionSystemDatagrokQuery) — no Samples
// package dependency. Step 4 is the full github-3550 reproduction:
// since the test owns the query, we actually rename it (no permission
// fallback) and assert the platform's reference-resolution invariant.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromDbQuery,
  provisionSystemDatagrokQuery,
  resetShell,
  assertProvenanceScript,
  PROVENANCE_PATTERNS,
  ProvisionedQuery,
  SYSTEM_DATAGROK_QUERIES,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Lifecycle Query: provisioned System:Datagrok query source', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const projectName = `lifecycle-query-${Date.now()}`;
  let provisioned: ProvisionedQuery | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 0: provision saved query on System:Datagrok', async () => {
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'lifecycle_query',
        sql: SYSTEM_DATAGROK_QUERIES.GROUPS_SAMPLE,
      });
      expect(provisioned.queryId).toBeTruthy();
      expect(provisioned.queryNqName).toMatch(/^[\w]+:[\w]+$/);
    });

    await softStep('Step 1: run provisioned query via canonical recorder', async () => {
      if (!provisioned) throw new Error('no provisioned query');
      const opened = await openTableFromDbQuery(page, provisioned.queryNqName);
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_query', opened.script);
      const callPattern = new RegExp(`${provisioned.queryNqName}\\(\\)`);
      expect(opened.script).toMatch(callPattern);
    });

    await softStep('Step 2: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: reopen → query re-executes server-side', async () => {
      if (!saved) throw new Error('no saved project');
      if (!provisioned) throw new Error('no provisioned query');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.db_query,
      );
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      expect(result.reopenedScript).toMatch(new RegExp(provisioned.resolvedName));
    });

    await softStep('Step 4: github-3550 reproduction — rename query, reopen, verify resolution', async () => {
      if (!saved) throw new Error('no saved project');
      if (!provisioned) throw new Error('no provisioned query');

      const renamedName = `${provisioned.resolvedName}_renamed`;
      const renameOk = await evalJs<boolean>(page, `(async () => {
        try {
          const q = await grok.dapi.queries.find('${provisioned.queryId}');
          if (!q) return false;
          q.name = '${renamedName}';
          await grok.dapi.queries.save(q);
          return true;
        } catch (_) { return false; }
      })()`);
      expect(renameOk).toBe(true);

      // github-3550 invariant: after the external query is renamed, reopening
      // the project must EITHER auto-resolve (happy path: relations.list
      // updated to new name, table re-materializes) OR fail with an explicit
      // error referencing the missing query (graceful failure). Silent null
      // references — relations.list has entries but table doesn't load and
      // no error surfaces — is the regression this test guards against.
      const result = await evalJs<{
        tablesAfter: number;
        relationsCount: number;
        loadedOk: boolean;
        errorMessage: string | null;
      }>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        let errorMessage = null;
        try { await p.open(); } catch (e) { errorMessage = String(e).slice(0, 300); }
        for (let i = 0; i < 60; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        const fresh = await grok.dapi.projects.find('${saved.projectId}');
        return {
          tablesAfter: grok.shell.tables.length,
          relationsCount: fresh.relations ? fresh.relations.length : 0,
          loadedOk: grok.shell.tables.length > 0,
          errorMessage,
        };
      })()`);

      const isHappyPath = result.loadedOk;
      const isGracefulFailure = !result.loadedOk && result.errorMessage !== null;
      expect(isHappyPath || isGracefulFailure).toBe(true);
    });

    await softStep('Step 5: rename project itself (verify via find-by-id)', async () => {
      if (!saved) return;
      const renamedProj = `${projectName}-renamed`;
      const r = await evalJs<{ok: boolean; persistedName: string | null}>(page, `(async () => {
        const p = await grok.dapi.projects.find('${saved.projectId}');
        if (!p) return {ok: false, persistedName: null};
        p.name = '${renamedProj}';
        await grok.dapi.projects.save(p);
        const verify = await grok.dapi.projects.find('${saved.projectId}');
        return {ok: verify?.name === '${renamedProj}', persistedName: verify?.name ?? null};
      })()`);
      expect(r.ok).toBe(true);
      expect(r.persistedName).toBe(renamedProj);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (provisioned) await provisioned.cleanup();
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
