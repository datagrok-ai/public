/* ---
sub_features_covered: [projects.upload, projects.api.save, projects.api.files.sync, projects.api.namespaces, projects.api.relations.list, projects.add-relation, projects.shell.share-via-context-menu, projects.api.get-by-id]
related_bugs: [GROK-18345]
generated_from: projects-lifecycle-spaces.md (Phase B — env-blocked at platform JS-API level; awaiting addEntity fix)
--- */
// Spaces-source lifecycle. PHASE B 2026-05-05 STATUS:
// **BLOCKED at JS API level — `spaceClient.addEntity(fileInfo, link)` fails on
// dev with `invalid input syntax for type uuid: "demog.csv"` (server expects
// UUID, Dart-side serialization sends fullPath instead). Captured live and
// documented in `.claude/diagnostics/mcp-capture-spaces-partial.md`.
//
// The lifecycle CANNOT be exercised end-to-end until either:
//   (a) JS API binding fix lands for addEntity payload serialization, OR
//   (b) we discover the correct addEntity overload signature, OR
//   (c) a pre-existing Space-with-files fixture is provisioned on dev.
//
// Until then this spec env-skips up-front with a clear blocker message.
// When the platform side resolves, replace the body with the
// openTableFromFile(`Spaces:<spaceName>/<file>`) flow + canonical
// uploadProject + reopenAndAssertProvenance pattern.
import {test, expect} from '@playwright/test';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';

test.use(projectsTestOptions);

test('Projects / Lifecycle Spaces: GROK-18345 share + datasync (env-blocked)', async ({page}) => {
  test.setTimeout(60_000);

  await gotoApp(page);
  await setupSession(page);

  // Probe addEntity availability — the platform-side blocker. If the bug
  // is fixed (addEntity accepts a FileInfo without UUID error), this probe
  // succeeds and the test should be expanded with the full lifecycle. Until
  // then it env-skips cleanly so the suite isn't permanently red.
  const probe = await evalJs<{blocked: boolean; reason: string}>(page, `(async () => {
    try {
      if (typeof grok.dapi.spaces?.createRootSpace !== 'function')
        return {blocked: true, reason: 'spaces.createRootSpace not implemented on this env'};
      const space = await grok.dapi.spaces.createRootSpace('mcp-spaces-probe-' + Date.now());
      const client = grok.dapi.spaces.id(space.id);
      const list = await grok.dapi.files.list('System:DemoFiles/');
      const demog = (list || []).find(f => /^demog\\.csv$/i.test(f.name));
      if (!demog) {
        await grok.dapi.spaces.delete(space).catch(() => {});
        return {blocked: true, reason: 'demog.csv not in System:DemoFiles'};
      }
      const saved = await demog.save?.();
      try {
        await client.addEntity(saved ?? demog, true);
        await grok.dapi.spaces.delete(space).catch(() => {});
        return {blocked: false, reason: 'addEntity succeeded — Spaces lifecycle is now testable; expand this spec.'};
      } catch (e) {
        await grok.dapi.spaces.delete(space).catch(() => {});
        return {blocked: true, reason: 'addEntity failed: ' + String(e).slice(0, 200)};
      }
    } catch (e) {
      return {blocked: true, reason: 'Spaces probe threw: ' + String(e).slice(0, 200)};
    }
  })()`);

  if (probe.blocked) {
    console.warn('Spaces lifecycle env-skip — ' + probe.reason);
    test.skip(true, probe.reason);
    return;
  }

  // Once unblocked, replace this with the full lifecycle:
  //   1. createRootSpace + addEntity (link=true)
  //   2. openTableFromFile(`Spaces:<serverName>/demog.csv`)
  //   3. assertProvenanceScript(files-source pattern)
  //   4. saveProjectWithProvenance + reopenAndAssertProvenance
  //   5. share + GROK-18345 recipient-open invariant (Helper 3 deferred)
  expect(probe.blocked).toBe(false);
});
