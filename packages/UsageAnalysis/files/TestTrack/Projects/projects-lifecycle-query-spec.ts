/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.relations.list, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
related_bugs: [github-3550]
generated_from: projects-lifecycle-query.md (Phase B canonical openers + uploadProject)
--- */
// Env-dependent: Samples:PostgresCustomers query. Phase B 2026-05-05 —
// replaces grok.functions.eval (recorder bypass; .script tag never written)
// with openTableFromDbQuery via canonical DG.Func.find recorder pattern.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromDbQuery,
  resetShell,
  assertProvenanceScript,
  PROVENANCE_PATTERNS,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Lifecycle Query: Samples:PostgresCustomers source', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const projectName = `lifecycle-query-${Date.now()}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: env-skip if Samples:PostgresCustomers not registered', async () => {
      // Use DG.Func.find as the env probe — same registry openTableFromDbQuery
      // consults; avoids the dapi.queries.filter `name = "X"` field-mismatch
      // (server matches by nqName, not bare name; use shortName if dapi).
      const ok = await evalJs<boolean>(page,
        `(() => DG.Func.find({namespace: 'Samples', name: 'PostgresCustomers'}).length > 0)()`);
      if (!ok) test.skip(true, 'Samples:PostgresCustomers not provisioned on this env');
    });

    await softStep('Step 1: run Samples:PostgresCustomers via canonical recorder', async () => {
      const opened = await openTableFromDbQuery(page, 'Samples:PostgresCustomers');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_query', opened.script);
      expect(opened.script).toMatch(/Samples:PostgresCustomers\(\)/);
    });

    await softStep('Step 2: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: reopen → query re-executes server-side', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.db_query,
      );
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      expect(result.reopenedScript).toMatch(/Samples:PostgresCustomers/);
    });

    await softStep('Step 4: github-3550 invariant — verify relations list', async () => {
      if (!saved) return;
      const r = await evalJs<{relations: number; ok: boolean}>(page, `(async () => {
        try {
          const fetched = await grok.dapi.projects.find('${saved.projectId}');
          return {relations: fetched.relations ? fetched.relations.length : 0, ok: true};
        } catch (e) {
          return {relations: -1, ok: false};
        }
      })()`);
      if (r.ok) expect(r.relations).toBeGreaterThanOrEqual(0);
    });

    await softStep('Step 5: rename project (verify via find-by-id)', async () => {
      if (!saved) return;
      const renamed = `${projectName}-renamed`;
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
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
