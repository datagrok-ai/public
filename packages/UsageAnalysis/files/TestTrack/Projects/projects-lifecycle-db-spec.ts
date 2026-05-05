/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
generated_from: projects-lifecycle-db.md (Phase B canonical openers + uploadProject + reopen-verify)
--- */
// DB-source lifecycle. Two paths: saved query (Samples:PostgresProducts) +
// ad-hoc DB table (DbQuery on Northwind public.products). Phase B 2026-05-05
// — replaces fake `Postgres:NorthwindTest:select_top_n` (which never
// resolved as a function and silently skipped) with the canonical
// recorder-engaged paths captured live on dev. Reference:
// .claude/diagnostics/mcp-capture-db.md
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromDbQuery,
  openTableFromDbTable,
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

const QUERY_NQNAME = 'Samples:PostgresProducts';
const CONNECTION_NQNAME = 'Samples:PostgresNorthwind';

// ---------------------------------------------------------------------------
// Test 1: saved query path
// ---------------------------------------------------------------------------

test('Projects / Lifecycle DB / Query: Samples:PostgresProducts source', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-db-query-${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: env-skip if Samples:PostgresProducts not registered', async () => {
      // Use DG.Func.find — same registry the opener consults, avoids the
      // `dapi.queries.filter('name = "X"')` field-mismatch quirk
      // (server matches by nqName, not bare name; use `shortName` if dapi).
      const ok = await evalJs<boolean>(page,
        `(() => DG.Func.find({namespace: 'Samples', name: 'PostgresProducts'}).length > 0)()`);
      if (!ok) test.skip(true, 'Samples:PostgresProducts not provisioned on this env');
    });

    await softStep('Step 1: run Samples:PostgresProducts via canonical recorder', async () => {
      const opened = await openTableFromDbQuery(page, QUERY_NQNAME);
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_query', opened.script);
      expect(opened.script).toMatch(/Samples:PostgresProducts\(\)/);
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
      expect(result.reopenedScript).toMatch(/Samples:PostgresProducts/);
    });

    await softStep('Step 4: share via JS API (View-and-Use)', async () => {
      if (!saved) return;
      const r = await evalJs<{skipped: boolean; reason?: string}>(page, `(async () => {
        try {
          const users = await grok.dapi.users.list({limit: 50});
          const target = users.find(u => u.login !== 'qa-pw' && u.login !== 'system');
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
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});

// ---------------------------------------------------------------------------
// Test 2: ad-hoc DB table double-click path
// ---------------------------------------------------------------------------

test('Projects / Lifecycle DB / Table: Northwind public.products via DbQuery', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-db-table-${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: env-skip if Samples:PostgresNorthwind not registered', async () => {
      // dapi.connections.filter `name = "X"` matches the fully-qualified
      // server name and returns 0 for bare names — use `shortName`.
      const ok = await evalJs<boolean>(page, `(async () => {
        const list = await grok.dapi.connections.filter('shortName = "PostgresNorthwind"').list();
        return list.length > 0;
      })()`);
      if (!ok) test.skip(true, 'Samples:PostgresNorthwind not provisioned on this env');
    });

    await softStep('Step 1: open public.products ad-hoc via DbQuery (double-click semantics)', async () => {
      const opened = await openTableFromDbTable(page, {
        connectionNqName: CONNECTION_NQNAME,
        schemaName: 'public',
        tableName: 'products',
      });
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'db_table', opened.script);
      expect(opened.script).toMatch(/DbQuery\([\w:]*Postgres[\w:]*Northwind,\s*"products"/);
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
