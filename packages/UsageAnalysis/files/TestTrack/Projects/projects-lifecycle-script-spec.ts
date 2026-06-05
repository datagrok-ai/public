/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.relations.list, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
related_bugs: [GROK-19403, GROK-19728]
generated_from: projects-lifecycle-script.md (Phase B canonical openers + uploadProject + reopen-verify)
--- */
// Script-source lifecycle. Phase B 2026-05-05 — replaces the broken
// `grok.functions.eval('Samples:Cars()')` (Cars is scalar-output and
// rejects no-input call → spec always silently skipped) with a
// fixture script provisioned via DG.Script.create that wraps
// grok.data.getDemoTable('demog.csv') and has output: dataframe df.
// This mirrors the test_Layout fixture in scripts-layout.test.ts:24-29.
// Reference: .claude/diagnostics/mcp-capture-scripts.md
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, gotoApp, setupSession} from './_helpers';
import {
  openTableFromScript,
  provisionDataframeScript,
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

test('Projects / Lifecycle Script: provisioned df-output script source', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const scriptName = `lifecycleScript${stamp}`;
  const projectName = `lifecycle-script-${stamp}`;
  let provisioned: {scriptId: string; resolvedName: string; resolvedNqName: string} | null = null;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: provision a JS script with output: dataframe', async () => {
      provisioned = await provisionDataframeScript(page, {
        name: scriptName,
        body: `df = await grok.data.getDemoTable('demog.csv');`,
      });
      expect(provisioned.scriptId).toBeTruthy();
      expect(provisioned.resolvedNqName).toMatch(/^[\w]+:[\w]+$/);
    });

    await softStep('Step 1: run provisioned script via canonical recorder', async () => {
      if (!provisioned) throw new Error('no provisioned script');
      const opened = await openTableFromScript(page, provisioned.resolvedNqName, {idx: 0});
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'script', opened.script);
      const scriptCallPattern = new RegExp(`${provisioned.resolvedName}\\s*\\(`);
      expect(opened.script).toMatch(scriptCallPattern);
    });

    await softStep('Step 2: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: reopen → script re-executes + table re-materializes', async () => {
      if (!saved) throw new Error('no saved project');
      if (!provisioned) throw new Error('no provisioned script');
      const result = await reopenAndAssertProvenance(
        page, saved.projectId, PROVENANCE_PATTERNS.script,
      );
      expect(result.tablesAfter).toBeGreaterThan(0);
      expect(result.reopenedRowCount).toBeGreaterThan(0);
      expect(result.reopenedScript).toMatch(new RegExp(provisioned.resolvedName));
    });

    // Share is LAST step before finally — the helper reloads the page for
    // second-user re-auth, so nothing UI/JS-state-dependent may follow it.
    // GROK-19403 recipient re-auth is now WIRED via shareWithSecondUserAndVerify
    // (token2-based second-user round-trip); recipient visibility is asserted
    // when DATAGROK_AUTH_TOKEN_2 is configured.
    await softStep('Step 4: GROK-19403 — share with second user (View-and-Use + Full) + recipient open', async () => {
      if (!saved) return;
      const r = await shareWithSecondUserAndVerify(page, {id: saved.projectId, name: projectName}, {full: true});
      if (!r.shared) { console.warn('Share skipped: ' + r.reason); return; }
      if (r.recipientVisible !== null) expect(r.recipientVisible).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (provisioned)
      await deleteProjectWithCleanup(page, {scriptId: provisioned.scriptId});
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
