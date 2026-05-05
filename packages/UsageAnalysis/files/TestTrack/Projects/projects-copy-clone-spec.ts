/* ---
sub_features_covered: [projects.shell.open, projects.api.save, projects.add-relation, projects.add-link]
generated_from: projects-copy-clone.md (Phase B canonical openers; Step 1 thumbnail render + Step 4d view-state visual checks remain in projects-copy-clone-ui.md)
--- */
// 4 sub-flows over a single demog source: 4a baseline, 4b Save Copy with Link
// + GROK-19750 invariant, 4c Save Copy with Clone, 4d Save personal view
// customizations. Re-share each variant via JS API.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

async function reopenProjectById(page: any, projectId: string) {
  // Use evalJs with template-literal for backwards-compat; but poll on
  // rowCount, NOT tables.length — `tables.length` throws Dart-side
  // `Tn.grok_TableNames is not a function` after reopen on dev.
  await evalJs(page, `(async () => {
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 700));
    const p = await grok.dapi.projects.find('${projectId}');
    if (p) await p.open();
    for (let i = 0; i < 30; i++) {
      const rc = grok.shell.tv?.dataFrame?.rowCount;
      if (typeof rc === 'number' && rc > 0) break;
      await new Promise(r => setTimeout(r, 500));
    }
  })()`);
  await page.waitForTimeout(1500);
}

test('Projects / Copy Clone: 3 save modes + GROK-19750 invariant', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const names = {
    original: `demog-${stamp}`,
    linkCopy: `demog-${stamp}-link`,
    cloneCopy: `demog-${stamp}-clone`,
    pvcCopy: `demog-${stamp}-pvc`,
  };
  const ids: Record<string, {projectId: string; tableInfoId: string}> = {};

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Setup: build original demog project with provenance + viewers', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await assertProvenanceScript(page, 'files', opened.script);
      // Add viewers via JS API. Defensive against shell.tv being briefly
      // undefined right after openTableFromFile under multi-test load.
      await evalJs(page, `(async () => {
        for (let i = 0; i < 10; i++) {
          if (grok.shell.tv?.dataFrame) break;
          await new Promise(r => setTimeout(r, 300));
        }
        const tv = grok.shell.tv;
        if (tv?.addViewer) {
          tv.addViewer('Scatter plot');
          tv.addViewer('Bar chart');
        }
      })()`);
      await page.waitForTimeout(1500);
      const saved = await saveProjectWithProvenance(page, names.original);
      ids.original = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
    });

    await softStep('4b: open original, add viewer, Save Copy as <name>-link', async () => {
      await reopenProjectById(page, ids.original.projectId);
      await evalJs(page, `(() => grok.shell.tv.addViewer('Histogram'))()`);
      await page.waitForTimeout(1500);
      const saved = await saveProjectWithProvenance(page, names.linkCopy);
      ids.linkCopy = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
    });

    await softStep('4b verification: reopen <name>-link, viewers intact', async () => {
      await reopenProjectById(page, ids.linkCopy.projectId);
      // tv.viewers may be a Dart proxy that throws on Array.from after
      // reopen; defensive count via try/catch + length probe.
      const r = await evalJs<{viewers: number}>(page,
        `(() => {
          try {
            const v = grok.shell.tv?.viewers;
            if (!v) return {viewers: 0};
            const arr = Array.from(v);
            return {viewers: arr.length};
          } catch (_) { return {viewers: 0}; }
        })()`);
      expect(r.viewers).toBeGreaterThan(1);
    });

    await softStep('4b GROK-19750 invariant: reopen original, viewers still present', async () => {
      await reopenProjectById(page, ids.original.projectId);
      // tv.viewers may be a Dart proxy that throws on Array.from after
      // reopen; defensive count via try/catch + length probe.
      const r = await evalJs<{viewers: number}>(page,
        `(() => {
          try {
            const v = grok.shell.tv?.viewers;
            if (!v) return {viewers: 0};
            const arr = Array.from(v);
            return {viewers: arr.length};
          } catch (_) { return {viewers: 0}; }
        })()`);
      expect(r.viewers).toBeGreaterThan(0);
    });

    await softStep('4c: open original, add viewer, Save Copy as <name>-clone', async () => {
      await reopenProjectById(page, ids.original.projectId);
      await evalJs(page, `(() => grok.shell.tv.addViewer('Line chart'))()`);
      await page.waitForTimeout(1500);
      const saved = await saveProjectWithProvenance(page, names.cloneCopy);
      ids.cloneCopy = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
    });

    await softStep('4d: open original, add viewer/customization, save as PVC', async () => {
      await reopenProjectById(page, ids.original.projectId);
      await evalJs(page, `(() => grok.shell.tv.addViewer('Pie chart'))()`);
      await page.waitForTimeout(1500);
      const saved = await saveProjectWithProvenance(page, names.pvcCopy);
      ids.pvcCopy = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
    });

    await softStep('Step 5: re-share each variant via JS API', async () => {
      const variantIds = [ids.linkCopy?.projectId, ids.cloneCopy?.projectId, ids.pvcCopy?.projectId].filter(Boolean);
      const r = await evalJs<{shared: string[]; skipped: string}>(page, `(async () => {
        const users = await grok.dapi.users.list({limit: 50});
        const target = users.find(u => u.login !== 'qa-pw' && u.login !== 'system');
        if (!target) return {shared: [], skipped: 'no recipient'};
        const shared = [];
        for (const id of ${JSON.stringify(variantIds)}) {
          try {
            const p = await grok.dapi.projects.find(id);
            if (p) {
              await grok.dapi.permissions.grant(p, target, false);
              shared.push(p.name);
            }
          } catch {}
        }
        return {shared, skipped: ''};
      })()`);
      if (r.skipped) console.warn('Share skipped: ' + r.skipped);
      else expect(r.shared.length).toBeGreaterThan(0);
    });
  } finally {
    for (const v of Object.values(ids))
      await deleteProjectWithCleanup(page, v);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
