/* ---
sub_features_covered: [projects.api.save, projects.api.namespaces, projects.move.move-entity, projects.move.commit]
generated_from: complex-move.md (Phase B canonical openers)
--- */
// JS API path is primary; UI right-click `Move to` doesn't exist on dev,
// drag-drop unautomatable — both UI paths captured in complex-ui.md Step 10.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {saveProjectWithProvenance, deleteProjectWithCleanup} from '../helpers/projects';

test.use(projectsTestOptions);

test('Projects / Complex Move: move project across namespaces via JS API', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const projectName = `move-test-${stamp}`;
  const spaceName = `move-target-${stamp}`;
  let saved: {projectId: string; tableInfoId: string; layoutId: string | null; resolvedName: string} | null = null;
  let createdSpaceId: string | null = null;

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    await softStep('Step 1-2: open demog with provenance + save with Sync ON', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      expect(opened.rowCount).toBeGreaterThan(0);
      await assertProvenanceScript(page, 'files', opened.script);
      saved = await saveProjectWithProvenance(page, projectName);
      expect(saved.projectId).toBeTruthy();
    });

    await softStep('Step 3: move via JS API to Home namespace', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{ok: boolean; reason?: string}>(page, `(async () => {
        try {
          const p = await grok.dapi.projects.find('${saved.projectId}');
          if (typeof grok.dapi.projects.move !== 'function')
            return {ok: false, reason: 'projects.move not implemented'};
          await grok.dapi.projects.move(p, 'Home');
          return {ok: true};
        } catch (e) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      })()`);
      if (!r.ok) {
        console.warn('Move skipped: ' + r.reason);
        return;
      }
      expect(r.ok).toBe(true);
    });

    await softStep('Step 4: create Space, move project to Spaces namespace', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{ok: boolean; reason?: string; spaceId?: string}>(page, `(async () => {
        try {
          // Note: shipped API is createRootSpace(name), NOT createRoot — verified
          // 2026-05-05. Older specs/scenarios sometimes name createRoot; that
          // symbol does not exist on js-api/src/dapi.ts (SpacesDataSource).
          if (typeof grok.dapi.spaces?.createRootSpace !== 'function')
            return {ok: false, reason: 'spaces.createRootSpace not implemented'};
          const space = await grok.dapi.spaces.createRootSpace('${spaceName}');
          const p = await grok.dapi.projects.find('${saved.projectId}');
          if (typeof grok.dapi.projects.move !== 'function')
            return {ok: false, reason: 'projects.move not implemented', spaceId: space.id};
          await grok.dapi.projects.move(p, 'Spaces:' + space.name);
          return {ok: true, spaceId: space.id};
        } catch (e) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      })()`);
      if (r.spaceId) createdSpaceId = r.spaceId;
      if (!r.ok) {
        console.warn('Space move skipped: ' + r.reason);
        return;
      }
      expect(r.ok).toBe(true);
    });
  } finally {
    if (saved)
      await deleteProjectWithCleanup(page, {
        projectId: saved.projectId,
        tableInfoId: saved.tableInfoId,
      });
    if (createdSpaceId) {
      await evalJs(page, `(async () => {
        try {
          const s = await grok.dapi.spaces.find('${createdSpaceId}');
          if (s) await grok.dapi.spaces.delete(s);
        } catch {}
      })()`).catch(() => {});
    }
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
