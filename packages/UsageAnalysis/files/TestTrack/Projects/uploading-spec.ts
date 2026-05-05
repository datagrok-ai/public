/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync]
generated_from: uploading.md (Phase B canonical openers + uploadProject + reopen-verify)
--- */
// Source-matrix scenario. Original uploading.md has 8 cases × 2 sync states
// = 16 saved projects post Phase B cleanup (Case 7 split out to
// ../Queries/get-all-get-top-100.md). This spec covers a representative
// subset (Cases 1, 8, 9 in Sync ON) — all three exercising the canonical
// recorder + uploadProject pattern verified live 2026-05-05.
//
// Phase B 2026-05-05: replaced grok.dapi.files.readCsv (no .script) +
// df.groupBy().aggregate() (no .script for derived) + saveProjectViaDialog
// with helpers.playwright.openers.* + helpers.playwright.projects.*.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromFile,
  openTableFromDbTable,
  addAggregateToWorkspace,
  assertProvenanceScript,
  resetShell,
  PROVENANCE_PATTERNS,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

// ---------------------------------------------------------------------------
// Case 1 — two file-share tables (Link Tables UI delegated to projects-ui-smoke).
// ---------------------------------------------------------------------------

test('Projects / Uploading / Case 1: Files + Files (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `Test_Case1_Sync_${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;
  let secondTableInfoId: string | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open demog.csv and cars.csv via OpenFile', async () => {
      const t1 = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(t1.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t1.script);

      const t2 = await openTableFromFile(page, 'System:DemoFiles/cars.csv');
      expect(t2.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', t2.script);
      secondTableInfoId = t2.tableInfoId;

      const tableCount = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      expect(tableCount).toBe(2);
    });

    await softStep('Step 2: save project with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: reopen → verify provenance survived', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.files,
      );
      expect(result.tablesAfter).toBeGreaterThanOrEqual(1);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (secondTableInfoId)
      await deleteProjectWithCleanup(page, {tableInfoId: secondTableInfoId});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------------
// Case 8 — Files + Pivot Table (Add to workspace).
// ---------------------------------------------------------------------------

test('Projects / Uploading / Case 8: Files + Pivot Table > Add (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `Test_Case8_Sync_${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;
  let baseTableInfoId: string | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open demog.csv as base', async () => {
      const base = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(base.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', base.script);
      baseTableInfoId = base.tableInfoId;
    });

    await softStep('Step 2: add Pivot Table viewer + click ADD', async () => {
      const tablesBefore = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      const derived = await addAggregateToWorkspace(page, {via: 'pivot-viewer'});
      expect(derived.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'derived', derived.script);
      const tablesAfter = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      expect(tablesAfter).toBe(tablesBefore + 1);
    });

    await softStep('Step 3: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
    });

    await softStep('Step 4: reopen → derived rematerializes from base', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(1);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (baseTableInfoId)
      await deleteProjectWithCleanup(page, {tableInfoId: baseTableInfoId});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------------
// Case 9 — DB table + Aggregate Rows (Add to workspace).
// ---------------------------------------------------------------------------

test('Projects / Uploading / Case 9: DB + Aggregate Rows > Add (Sync ON)', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `Test_Case9_Sync_${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;
  let baseTableInfoId: string | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: env-skip if Samples:PostgresNorthwind not registered', async () => {
      const ok = await evalJs<boolean>(page, `(async () => {
        const list = await grok.dapi.connections.filter('shortName = "PostgresNorthwind"').list();
        return list.length > 0;
      })()`);
      if (!ok) test.skip(true, 'Samples:PostgresNorthwind not provisioned on this env');
    });

    await softStep('Step 1: open public.products via DbQuery (double-click semantics)', async () => {
      const base = await openTableFromDbTable(page, {
        connectionNqName: 'Samples:PostgresNorthwind',
        schemaName: 'public',
        tableName: 'products',
      });
      expect(base.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_table', base.script);
      baseTableInfoId = base.tableInfoId;
    });

    await softStep('Step 2: Aggregate Rows via Top menu + ADD', async () => {
      const tablesBefore = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      const derived = await addAggregateToWorkspace(page, {via: 'menu'});
      expect(derived.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'derived', derived.script);
      const tablesAfter = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      expect(tablesAfter).toBe(tablesBefore + 1);
    });

    await softStep('Step 3: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
    });

    await softStep('Step 4: reopen → DbQuery + Aggregate re-execute server-side', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(1);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (baseTableInfoId)
      await deleteProjectWithCleanup(page, {tableInfoId: baseTableInfoId});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
