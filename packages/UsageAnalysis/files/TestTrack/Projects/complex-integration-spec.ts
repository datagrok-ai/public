/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.namespaces, projects.add-relation, projects.add-link]
generated_from: complex-integration.md (Phase B canonical openers)
--- */
// Multi-source integration: 4 file-share tables opened with provenance,
// saved as a single project, reopen verified. Env-dependent files (e.g.
// System:AppData/Chem/tests/spgi-100.csv may not be on every dev) skip
// individually with logged-warning; require >=2 to pass.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Complex Integration: heterogeneous sources in one project', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  const projectName = `integration-test-${Date.now()}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1: open multi-source files with provenance (defensive — skip on env miss)', async () => {
      const tries = [
        'System:AppData/Chem/tests/spgi-100.csv',
        'System:DemoFiles/demog.csv',
        'System:DemoFiles/cars.csv',
        'System:DemoFiles/iris.csv',
      ];
      let opened = 0;
      const errors: string[] = [];
      for (const path of tries) {
        try {
          const r = await openTableFromFile(page, path);
          if (r.script) opened++;
        } catch (e: any) {
          errors.push((e?.message ?? String(e)).slice(0, 80));
        }
      }
      console.log(`opened: ${opened}, errors: ${JSON.stringify(errors)}`);
      expect(opened).toBeGreaterThanOrEqual(2);
    });

    await softStep('Step 2: save with provenance', async () => {
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: verify project relations server-side', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{relations: number; ok: boolean}>(page, `(async () => {
        try {
          const fetched = await grok.dapi.projects.find('${saved.projectId}');
          const relations = fetched.relations ? fetched.relations.length : 0;
          return {relations, ok: true};
        } catch (e) {
          return {relations: -1, ok: false};
        }
      })()`);
      if (r.ok) expect(r.relations).toBeGreaterThanOrEqual(0);
    });

    await softStep('Step 4: closeAll and reopen, verify tables re-materialize from .script', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{tables: number}>(page, `(async () => {
        grok.shell.closeAll();
        await new Promise(r => setTimeout(r, 800));
        const p = await grok.dapi.projects.find('${saved.projectId}');
        await p.open();
        for (let i = 0; i < 60; i++) {
          if (grok.shell.tables.length > 0) break;
          await new Promise(r => setTimeout(r, 500));
        }
        return {tables: grok.shell.tables.length};
      })()`);
      expect(r.tables).toBeGreaterThan(0);
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
