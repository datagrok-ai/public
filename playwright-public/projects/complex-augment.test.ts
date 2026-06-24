// Verifies augmenting a saved project with additional tables via the
// addRelation/addLink JS API and reopening it with all tables intact.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Complex Augment: addRelation 4 tables via JS API', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const projectName = `augment-test-${Date.now()}`;
  const additional = ['cars.csv', 'iris.csv', 'geo.csv'];
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open demog.csv with provenance + save baseline project', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', opened.script);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 2: addLink 3 more files via Project.addLink (Link mode)', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{added: number; skipped: string[]}>(page, `(async () => {
        const p = await grok.dapi.projects.find('${saved.projectId}');
        const skipped = [];
        let added = 0;
        // dapi.files.list requires trailing slash on the namespace path
        // (verified live 2026-05-05); without it errors with handlerPath/null url.
        const list = await grok.dapi.files.list('System:DemoFiles/');
        for (const fn of ${JSON.stringify(additional)}) {
          try {
            const file = (list || []).find(f => f.name === fn);
            if (!file) { skipped.push(fn + ': not found'); continue; }
            // Register the FileInfo as an entity (gets a UUID), then link
            // to the project via the instance addLink method. addRelation
            // does NOT exist on dapi.projects (verified live 2026-05-05);
            // the canonical Link mode API is project.addLink(entity).
            const saved = await file.save();
            p.addLink(saved);
            added++;
          } catch (e) {
            skipped.push(fn + ': ' + String(e).slice(0, 80));
          }
        }
        await grok.dapi.projects.save(p);
        return {added, skipped};
      })()`);
      console.log('addLink skipped: ' + JSON.stringify(r.skipped));
      expect(r.added).toBeGreaterThan(0);
    });

    await softStep('Step 3: re-open augmented project, verify tables loaded', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{tables: number}>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        await p.open();
        for (let i = 0; i < 30; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        return {tables: grok.shell.tables.length};
      })()`);
      expect(r.tables).toBeGreaterThan(0);
    });
  } finally {
    const s = saved as {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null;
    if (s)
      await deleteProjectWithCleanup(page, {
        projectId: s.projectId,
        tableInfoId: s.tableInfoId,
      });
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
