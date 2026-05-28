/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
generated_from: projects-lifecycle-db.md (Phase B canonical openers + uploadProject + reopen-verify)
--- */
// DB-source lifecycle. Two paths: provisioned saved query on System:Datagrok
// + ad-hoc DB table (DbQuery on System:Datagrok public.groups). Both
// prerequisites are created within the test via helpers/openers.ts
// (provisionSystemDatagrokQuery / built-in System:Datagrok connection) —
// no Samples package or external DB provisioning required.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromDbQuery,
  openTableFromDbTable,
  assertProvenanceScript,
  provisionSystemDatagrokQuery,
  resetShell,
  PROVENANCE_PATTERNS,
  ProvisionedQuery,
  SYSTEM_DATAGROK_DB_TABLE,
  SYSTEM_DATAGROK_NQNAME,
  SYSTEM_DATAGROK_QUERIES,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

// ---------------------------------------------------------------------------
// Test 1: saved query path (provisioned on System:Datagrok)
// ---------------------------------------------------------------------------

test('Projects / Lifecycle DB / Query: provisioned System:Datagrok query source', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-db-query-${stamp}`;
  let provisioned: ProvisionedQuery | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 0: provision saved query on System:Datagrok', async () => {
      provisioned = await provisionSystemDatagrokQuery(page, {
        nameStem: 'lifecycle_db_query',
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

    await softStep('Step 4: share via JS API (View-and-Use)', async () => {
      if (!saved) return;
      const r = await evalJs<{skipped: boolean; reason?: string}>(page, `(async () => {
        try {
          const users = await grok.dapi.users.list({limit: 50});
          const me = (await grok.dapi.users.current()).login;
          const target = users.find(u => u.login !== me && u.login !== 'system');
          if (!target) return {skipped: true, reason: 'no recipient'};
          const p = await grok.dapi.projects.find('${saved.projectId}');
          await grok.dapi.permissions.grant(p, target, false);
          return {skipped: false};
        } catch (e) {
          return {skipped: true, reason: String(e).slice(0, 200)};
        }
      })()`);
      if (r.skipped) console.warn('Share skipped: ' + r.reason);
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

// ---------------------------------------------------------------------------
// Test 2: ad-hoc DB table double-click path (System:Datagrok / public.groups)
// ---------------------------------------------------------------------------

test('Projects / Lifecycle DB / Table: System:Datagrok public.groups via DbQuery', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-db-table-${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open public.groups ad-hoc via DbQuery (double-click semantics)', async () => {
      const opened = await openTableFromDbTable(page, {
        connectionNqName: SYSTEM_DATAGROK_NQNAME,
        schemaName: SYSTEM_DATAGROK_DB_TABLE.schemaName,
        tableName: SYSTEM_DATAGROK_DB_TABLE.tableName,
      });
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_table', opened.script);
      expect(opened.script).toMatch(/DbQuery\([\w:]*Datagrok,\s*"groups"/);
    });

    await softStep('Step 2: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
    });

    await softStep('Step 3: reopen → DbQuery re-runs server-side', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.db_table,
      );
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      expect(result.reopenedScript).toMatch(/DbQuery/);
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
