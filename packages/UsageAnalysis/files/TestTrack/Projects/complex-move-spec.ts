import {expect} from '@playwright/test';
import {test} from '../shared-page';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {deleteProjectWithCleanup} from '../helpers/projects';
import {saveProjectWithProvenance} from './projects-shared';

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

  // The capability probe belongs here, not inside a softStep: test.skip() throws, softStep
  // catches it as a step error, and finishSpec then reported the env-skip as a failure.
  const api = await evalJs<{move: boolean; spaces: boolean}>(page, `(() => ({
    move: typeof grok.dapi.projects.move === 'function',
    spaces: typeof grok.dapi.spaces?.createRootSpace === 'function',
  }))()`);
  test.skip(!api.move || !api.spaces,
    'no namespace move on this build: grok.dapi.projects.move ' +
    (api.move ? 'exists' : 'is absent (public/js-api/src/dapi.ts declares no ProjectsDataSource.move)') +
    ', grok.dapi.spaces.createRootSpace ' + (api.spaces ? 'exists' : 'is absent'));

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
          await grok.dapi.projects.move(p, 'Home');
          return {ok: true};
        } catch (e) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      })()`);

      expect(r.ok, r.ok ? '' : `project move to Home failed: ${r.reason}`).toBe(true);
    });

    await softStep('Step 4: create Space, move project to Spaces namespace', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{ok: boolean; reason?: string; spaceId?: string}>(page, `(async () => {
        try {
          // Shipped API is createRootSpace(name), NOT createRoot.
          const space = await grok.dapi.spaces.createRootSpace('${spaceName}');
          const p = await grok.dapi.projects.find('${saved.projectId}');
          await grok.dapi.projects.move(p, 'Spaces:' + space.name);
          return {ok: true, spaceId: space.id};
        } catch (e) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      })()`);
      if (r.spaceId) createdSpaceId = r.spaceId;

      expect(r.ok, r.ok ? '' : `space create + move failed: ${r.reason}`).toBe(true);
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

  finishSpec();
});
