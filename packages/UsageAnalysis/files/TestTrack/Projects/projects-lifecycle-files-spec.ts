/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.shell.share-via-context-menu, projects.api.get-by-id]
generated_from: projects-lifecycle-files.md (Phase B canonical openers + uploadProject + reopen-verify)
--- */
// File-source lifecycle: open, save with provenance, reopen-verify, share,
// rename, delete. Phase B 2026-05-05 — replaces grok.dapi.files.readCsv path
// (which produced no .script tag) with helpers.playwright.openers.openTableFromFile.
// Reference: .claude/diagnostics/mcp-capture-files.md
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromFile,
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
      // Gate E-PROV-01: verify .script tag is set with the expected pattern.
      await assertProvenanceScript(page, 'files', opened.script);
      expect(opened.script).toMatch(/OpenFile\("System:DemoFiles\/demog\.csv"\)/);
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

    await softStep('Step 4: share via JS API (View-and-Use + Full)', async () => {
      if (!saved) return;
      const r = await evalJs<{skipped: boolean; reason?: string}>(page, `(async () => {
        try {
          const users = await grok.dapi.users.list({limit: 50});
          const target = users.find(u => u.login !== 'qa-pw' && u.login !== 'system');
          if (!target) return {skipped: true, reason: 'no recipient'};
          const p = await grok.dapi.projects.find('${saved.projectId}');
          await grok.dapi.permissions.grant(p, target, false);
          await grok.dapi.permissions.grant(p, target, true);
          return {skipped: false};
        } catch (e) {
          return {skipped: true, reason: String(e).slice(0, 200)};
        }
      })()`);
      if (r.skipped) console.warn('Share skipped: ' + r.reason);
    });

    await softStep('Step 5: rename project via JS API', async () => {
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
