import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
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
      // Legitimate skip ONLY when the API isn't shipped on this build; any other
      // reason (a thrown move error) must fail so a real regression is caught.
      if (!r.ok && r.reason === 'projects.move not implemented') {
        test.skip(true, 'projects.move not implemented on this build');
        return;
      }
      expect(r.ok, r.ok ? '' : `project move to Home failed: ${r.reason}`).toBe(true);
    });

    await softStep('Step 4: create Space, move project to Spaces namespace', async () => {
      if (!saved) throw new Error('no saved project');
      const r = await evalJs<{ok: boolean; reason?: string; spaceId?: string}>(page, `(async () => {
        try {
          // Shipped API is createRootSpace(name), NOT createRoot.
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
      // Legitimate skip ONLY when an API isn't shipped on this build; any other
      // reason (a thrown create/move error) must fail so a real regression is caught.
      if (!r.ok && (r.reason === 'spaces.createRootSpace not implemented' || r.reason === 'projects.move not implemented')) {
        test.skip(true, r.reason);
        return;
      }
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
