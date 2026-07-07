// DB-source lifecycle: provisioned saved query + ad-hoc DB table on System:Datagrok public.groups.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, gotoApp, setupSession} from './_helpers';
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
  shareWithSecondUserAndVerify,
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

    await softStep('Step 4: share with second user (View-and-Use) + recipient open', async () => {
      if (!saved) return;
      const r = await shareWithSecondUserAndVerify(page, {id: saved.projectId, name: projectName});
      if (!r.shared) { console.warn('Share skipped: ' + r.reason); return; }
      // When a second-user token is configured, the recipient must see it.
      if (r.recipientVisible !== null) expect(r.recipientVisible).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (provisioned) await provisioned.cleanup();
  }

  finishSpec();
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

  finishSpec();
});
