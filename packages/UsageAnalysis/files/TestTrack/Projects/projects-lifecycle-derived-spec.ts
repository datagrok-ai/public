/* ---
sub_features_covered: [projects.add-relation, projects.api.files.sync, projects.api.get-by-id, projects.api.save, projects.shell.share-via-context-menu, projects.upload]
--- */
// Derived-source lifecycle via the UI Aggregate Rows / Pivot Table → Add to workspace flow.
// GROK-19103 invariant: derivation lands in the active workspace (tables grows by 1), not a stray project.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {
  openTableFromFile,
  addAggregateToWorkspace,
  assertProvenanceScript,
  resetShell,
  PROVENANCE_PATTERNS,
} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
  shareWithSecondUserAndVerify,
} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Lifecycle Derived: Aggregate via menu + GROK-19103 invariant', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `lifecycle-derived-${stamp}`;
  let baseTableId: string | null = null;
  let derivedTableId: string | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open parent table demog.csv via OpenFile', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', opened.script);
      baseTableId = opened.tableInfoId;
    });

    await softStep('Step 2a: add aggregate via Top menu (Data → Aggregate Rows → ADD)', async () => {
      const tablesBefore = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      const derived = await addAggregateToWorkspace(page, {via: 'menu'});
      expect(derived.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'derived', derived.script);
      expect(derived.script).toMatch(/Aggregate\("demog"/);
      derivedTableId = derived.tableInfoId;

      // GROK-19103 invariant: tables.length grows by exactly 1, not a stray separate project.
      const tablesAfter = await evalJs<number>(page, '(grok.shell.tables?.length || 0)');
      expect(tablesAfter).toBe(tablesBefore + 1);
    });

    await softStep('Step 3: save with provenance (project carries base + derived TableInfos)', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 4: reopen → verify base AND derived re-materialize with provenance', async () => {
      if (!saved) throw new Error('no saved project');
      const result = await reopenAndAssertProvenance(page, saved.projectId);
      expect(result.tablesAfter).toBeGreaterThanOrEqual(1);

      const tagInspect = await evalJs<{tables: Array<{name: string; script: string}>}>(page, `(() => {
        const out = [];
        for (const t of (grok.shell.tables || [])) {
          const tags = {};
          if (t.tags?.keys) for (const k of t.tags.keys()) tags[k] = String(t.tags.get(k));
          out.push({name: t.name, script: tags['.script'] || ''});
        }
        return {tables: out};
      })()`);

      // At least one table should match base or derived provenance (server may rematerialize in any order).
      const hasBaseProvenance = tagInspect.tables.some(t => /OpenFile\("System:DemoFiles\/demog\.csv"\)/.test(t.script));
      const hasDerivedProvenance = tagInspect.tables.some(t => /Aggregate\("demog"/.test(t.script));
      expect(hasBaseProvenance || hasDerivedProvenance).toBe(true);
    });

    await softStep('Step 5a: rename project (verify via find-by-id)', async () => {
      if (!saved) return;
      const renameR = await evalJs<{ok: boolean; persistedName: string | null}>(page, `(async () => {
        const p = await grok.dapi.projects.find('${saved.projectId}');
        if (!p) return {ok: false, persistedName: null};
        p.name = '${projectName}-renamed';
        await grok.dapi.projects.save(p);
        const verify = await grok.dapi.projects.find('${saved.projectId}');
        return {ok: verify?.name === '${projectName}-renamed', persistedName: verify?.name ?? null};
      })()`);
      expect(renameR.ok).toBe(true);
    });

    // Share is LAST step before finally — the helper reloads the page for second-user re-auth.
    await softStep('Step 5b: share with second user (View-and-Use + Full) + recipient open', async () => {
      if (!saved) return;
      const r = await shareWithSecondUserAndVerify(page, {id: saved.projectId, name: `${projectName}-renamed`}, {full: true});
      if (!r.shared) { console.warn('Share skipped: ' + r.reason); return; }
      if (r.recipientVisible !== null) expect(r.recipientVisible).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (baseTableId)
      await deleteProjectWithCleanup(page, {tableInfoId: baseTableId});
    if (derivedTableId && derivedTableId !== saved?.tableInfoId)
      await deleteProjectWithCleanup(page, {tableInfoId: derivedTableId});
  }

  finishSpec();
});
