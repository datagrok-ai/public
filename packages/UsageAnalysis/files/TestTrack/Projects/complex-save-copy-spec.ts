/* ---
sub_features_covered: [projects.api.save, projects.api.files.sync, projects.add-relation]
generated_from: complex-save-copy.md (Phase B canonical openers)
--- */
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Complex Save Copy: round-trip Save Copy with sync OFF/ON', async ({page}) => {
  test.setTimeout(420_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const original = `save-copy-${stamp}-Original`;
  const noSync = `save-copy-${stamp}-NoSync`;
  const sync = `save-copy-${stamp}-Sync`;
  const savedIds: Array<{projectId: string; tableInfoId: string}> = [];

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1-2: open demog with provenance + save baseline (Sync ON)', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await assertProvenanceScript(page, 'files', opened.script);
      const saved = await saveProjectWithProvenance(page, original);
      savedIds.push({projectId: saved.projectId, tableInfoId: saved.tableInfoId});
    });

    await softStep('Step 3: Save Copy under NoSync name (re-open + save under new name)', async () => {
      await evalJs(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 600));
        const p = await grok.dapi.projects.find('${savedIds[0].projectId}');
        await p.open();
        for (let i = 0; i < 30; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
      })()`);
      await page.waitForTimeout(1000);
      const saved = await saveProjectWithProvenance(page, noSync);
      savedIds.push({projectId: saved.projectId, tableInfoId: saved.tableInfoId});
    });

    await softStep('Step 4: close and reopen NoSync copy', async () => {
      const r = await evalJs<{tables: number}>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${savedIds[1].projectId}');
        await p.open();
        for (let i = 0; i < 30; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        return {tables: grok.shell.tables.length};
      })()`);
      expect(r.tables).toBeGreaterThan(0);
    });

    await softStep('Step 5-6: re-save the copy under Sync name', async () => {
      const saved = await saveProjectWithProvenance(page, sync);
      savedIds.push({projectId: saved.projectId, tableInfoId: saved.tableInfoId});
    });
  } finally {
    for (const ids of savedIds)
      await deleteProjectWithCleanup(page, ids);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
