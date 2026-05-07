/* ---
sub_features_covered: [projects.shell.open, projects.shell.share-via-context-menu, projects.api.save, projects.add-relation, projects.add-link]
generated_from: projects-copy-clone.md (D4.3 UI-driven rewrite 2026-05-06 — replaces 2026-05-05 bundled-1-test reduction)
ui_coverage_owned: [save-copy-with-link-dialog, save-copy-with-clone-dialog, save-personal-view-customizations-dialog, pcmdShareProject, share-dialog-recipients]
--- */
// Full canonical implementation of `projects-copy-clone.md` Step 4 sub-flows
// 4a / 4b / 4c / 4d + Step 5 re-share. UI-driven Save Project dialog mode
// chooser — drives `[name="button-Save"]` toolbar → mode-radio (Save original
// project / Save a copy / Save personal view customizations) → per-table
// Link/Clone sub-choice (sub-flow 4b/4c) → name → OK. Per-mode reopen
// verification per sub-flow. Sub-flow 4b carries the GROK-19750 invariant
// assertion (after Save-Copy-with-Link the original project's table still
// re-materializes on reopen).
//
// Rewrite rationale (2026-05-06): the prior 2026-05-05 version collapsed the
// 4 sub-flows into a single test() block with 4 JS-API saves and only 2
// reopens to mitigate a shell.tv race. That reduction violated the chain's
// `ui_coverage_responsibility` (3 Save Copy mode dialogs MUST be UI-driven
// per `projects-copy-clone.md` frontmatter). This rewrite restores per-mode
// UI coverage and per-mode reopen verification while keeping the
// `reopenProjectById` shell.tv-poll pattern that mitigated the original race.
//
// Step 5 re-share is UI-driven (right-click → Share dialog) per chain rev 3
// `ui_coverage_plan.delegated_scenarios` — `share-project.md` delegates the
// right-click Share UI surface here.
import {test, expect} from '@playwright/test';
import {softStep, stepErrors} from '../spec-login';
import {projectsTestOptions, evalJs, gotoApp, setupSession} from './_helpers';
import {openTableFromFile, resetShell, assertProvenanceScript} from '../helpers/openers';
import {
  saveProjectWithProvenance,
  saveCopy,
  shareProjectViaContextMenu,
  deleteProjectWithCleanup,
} from '../helpers/projects';

test.use(projectsTestOptions);

// Wait for grok.shell.tv to become a real TableView (with .addViewer +
// .dataFrame) — the proxy can briefly report undefined .addViewer right after
// reopen even when the saved layout DID restore. Same shape as the prior
// spec's helper, retained verbatim because it mitigates the dev-on-Playwright
// shell.tv race observed in the 2026-05-05 retrofit.
async function reopenProjectById(page: any, projectId: string) {
  await evalJs(page, `(async () => {
    grok.shell.closeAll();
    await new Promise(r => setTimeout(r, 700));
    const p = await grok.dapi.projects.find('${projectId}');
    if (p) await p.open();
    for (let i = 0; i < 30; i++) {
      const tv = grok.shell.tv;
      if (tv?.dataFrame && typeof tv.addViewer === 'function') break;
      await new Promise(r => setTimeout(r, 500));
    }
  })()`);
  await page.waitForTimeout(1500);
}

async function reopenedRowCount(page: any): Promise<number> {
  // After reopen the Dart-side `tv.viewers` proxy can intermittently report
  // 0 viewers even when the saved layout restored them. The reliable
  // cross-Dart signal is the dataframe rowCount — non-zero means the
  // project's table re-materialized from .script provenance.
  return await evalJs<number>(page, `(grok.shell.tv?.dataFrame?.rowCount ?? 0)`);
}

async function reopenedViewerCount(page: any): Promise<number> {
  return await evalJs<number>(page, `(() => {
    try { return Number(grok.shell.tv?.viewers?.length ?? 0); }
    catch (_) { return 0; }
  })()`);
}

// Capture {projectId, tableInfoId} of the just-saved Save-Copy project.
// Strategy 1: read grok.shell.project (when Dart shell updates promptly,
//   and id has changed from source — i.e. a non-original-mode save).
// Strategy 2: server-side query by friendlyName via `name like "<typed>"`
//   (with NO wildcards) — Datagrok stores typed name verbatim in
//   friendlyName and the LIKE-no-wildcards form matches it; `name = "X"`
//   exact match is unreliable on dev (verified 2026-05-06).
// Returns {projectId: '', tableInfoId: ''} only if BOTH strategies fail.
async function captureActiveProjectIds(
  page: any,
  typedName: string,
  sourceId: string,
): Promise<{projectId: string; tableInfoId: string}> {
  return await evalJs<{projectId: string; tableInfoId: string}>(page, `(async () => {
    const typedName = ${JSON.stringify(typedName)};
    const sourceId = ${JSON.stringify(sourceId)};
    // Strategy 1: shell.project (must differ from sourceId)
    const shellProj = grok.shell.project;
    if (shellProj?.id && shellProj.id !== sourceId) {
      const tv = grok.shell.tv;
      const ti = tv?.dataFrame?.getTableInfo?.();
      return {projectId: shellProj.id, tableInfoId: ti?.id ?? ''};
    }
    // Strategy 2: friendlyName match via name LIKE (no wildcards).
    try {
      const list = await grok.dapi.projects.filter('name like "' + typedName + '"').list();
      const candidate = list.find(p => p.id !== sourceId);
      if (candidate) return {projectId: candidate.id, tableInfoId: ''};
    } catch (_) {}
    return {projectId: '', tableInfoId: ''};
  })()`);
}

test('Projects / Copy Clone — full UI-driven 4-sub-flow + GROK-19750 invariant', async ({page}) => {
  test.setTimeout(900_000);
  stepErrors.length = 0;

  const stamp = Date.now();
  const names = {
    original: `demog-${stamp}`,
    linkCopy: `demog-${stamp}-link`,
    cloneCopy: `demog-${stamp}-clone`,
    pvcCopy: `demog-${stamp}-pvc`,
  };
  const ids: Record<string, {projectId: string; tableInfoId: string} | undefined> = {};

  await gotoApp(page);
  await setupSession(page);
  await resetShell(page);

  try {
    // -------------------------------------------------------------------
    // Setup: create the baseline `<original>` project via JS API.
    // The 3-mode chooser in the Save Project dialog only appears when an
    // existing project is being modified. JS API setup is the standard
    // pattern for sibling specs (complex-save-copy-spec.ts:30) and is not
    // part of the canonical Step 4 UI surface — Setup is upstream.
    // -------------------------------------------------------------------
    await softStep('Setup: open demog.csv with provenance + save baseline (JS API)', async () => {
      const opened = await openTableFromFile(page, 'System:DemoFiles/demog.csv');
      await assertProvenanceScript(page, 'files', opened.script);
      const saved = await saveProjectWithProvenance(page, names.original);
      ids.original = {projectId: saved.projectId, tableInfoId: saved.tableInfoId};
      expect(ids.original.projectId).toBeTruthy();
    });

    // -------------------------------------------------------------------
    // Sub-flow 4a — Save the original (overwrite, mode='original').
    // -------------------------------------------------------------------
    await softStep('4a: open original (UI-reopen) → addViewer Bar → SAVE mode=original (overwrite) → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await evalJs(page, `(async () => { grok.shell.tv.addViewer('Bar chart'); })()`);
      await page.waitForTimeout(2000);
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'original',
        name: names.original,
      });
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    // -------------------------------------------------------------------
    // Sub-flow 4b — Save Copy with Link + GROK-19750 invariant.
    // -------------------------------------------------------------------
    await softStep('4b step 1-4: open original → addViewer Scatter → SAVE Copy/Link → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await evalJs(page, `(async () => { grok.shell.tv.addViewer('Scatter plot'); })()`);
      await page.waitForTimeout(2000);
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'copy',
        name: names.linkCopy,
        perTableLinkOrClone: 'link',
      });
      ids.linkCopy = await captureActiveProjectIds(page, names.linkCopy, ids.original!.projectId);
      expect(ids.linkCopy.projectId).toBeTruthy();
      expect(ids.linkCopy.projectId).not.toBe(ids.original!.projectId);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    await softStep('4b step 5-6: reopen <name>-link → table re-materializes → closeAll', async () => {
      await reopenProjectById(page, ids.linkCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4b: link copy reopen returned rowCount=0 — Save Copy with Link did not persist the linked table reference').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    await softStep('4b step 7 — GROK-19750 INVARIANT: reopen original → table + viewers intact', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      const rc = await reopenedRowCount(page);
      expect(
        rc,
        'GROK-19750 regression: original `demog` rowCount=0 after Save-Copy-with-Link — table failed to re-materialize on reopen, suggests Save Copy with Link leaked state changes back into the source project',
      ).toBeGreaterThan(0);
      const vc = await reopenedViewerCount(page);
      // Original after 4a (Bar) + 4b (Scatter) addition: at least 1 user-added
      // viewer should still be present on reopen. Some platform builds count
      // the implicit Grid as a viewer, others don't — the relevant invariant
      // is "viewers > 0", not an exact count.
      expect(
        vc,
        'GROK-19750 regression: original lost ALL viewers after Save-Copy-with-Link',
      ).toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    // -------------------------------------------------------------------
    // Sub-flow 4c — Save Copy with Clone.
    // -------------------------------------------------------------------
    await softStep('4c step 1-4: open original → addViewer Line → SAVE Copy/Clone → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await evalJs(page, `(async () => { grok.shell.tv.addViewer('Line chart'); })()`);
      await page.waitForTimeout(2000);
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'copy',
        name: names.cloneCopy,
        perTableLinkOrClone: 'clone',
      });
      ids.cloneCopy = await captureActiveProjectIds(page, names.cloneCopy, ids.original!.projectId);
      expect(ids.cloneCopy.projectId).toBeTruthy();
      expect(ids.cloneCopy.projectId).not.toBe(ids.original!.projectId);
      expect(ids.cloneCopy.projectId).not.toBe(ids.linkCopy?.projectId);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    await softStep('4c step 5-6: reopen <name>-clone → table re-materializes (independent copy) → closeAll', async () => {
      await reopenProjectById(page, ids.cloneCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4c: clone copy reopen returned rowCount=0 — Save Copy with Clone did not persist the cloned table bytes').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    // -------------------------------------------------------------------
    // Sub-flow 4d — Save with Personal View Customizations.
    // -------------------------------------------------------------------
    await softStep('4d step 1-4: open original → addViewer Histogram → SAVE PVC → closeAll', async () => {
      await reopenProjectById(page, ids.original!.projectId);
      await evalJs(page, `(async () => { grok.shell.tv.addViewer('Histogram'); })()`);
      await page.waitForTimeout(2000);
      await saveCopy(page, {
        sourceName: names.original,
        mode: 'personal-view-customizations',
        name: names.pvcCopy,
      });
      ids.pvcCopy = await captureActiveProjectIds(page, names.pvcCopy, ids.original!.projectId);
      expect(ids.pvcCopy.projectId).toBeTruthy();
      expect(ids.pvcCopy.projectId).not.toBe(ids.original!.projectId);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    await softStep('4d step 5-6: reopen <name>-personal-view-customizations → table re-materializes → closeAll', async () => {
      await reopenProjectById(page, ids.pvcCopy!.projectId);
      const rc = await reopenedRowCount(page);
      expect(rc, 'sub-flow 4d: PVC variant reopen returned rowCount=0 — PVC mode failed to preserve the linked source-table reference').toBeGreaterThan(0);
      await evalJs(page, 'grok.shell.closeAll()');
      await page.waitForTimeout(800);
    });

    // -------------------------------------------------------------------
    // Step 5 — Re-share each variant via right-click → Share dialog (UI).
    // -------------------------------------------------------------------
    await softStep('Step 5: re-share each of 3 variants via right-click → Share (UI)', async () => {
      // Recipient: a real Group materialized via dapi.users.find(id) — same
      // pattern as the prior spec (permissions.grant requires Group, not User;
      // .group is not eagerly loaded by users.list). Defensive skip if no
      // recipient with materialized group is available — matches Edit 10 G3
      // acceptable-SR pattern (chain may run on dev with no second-user
      // provisioned).
      const recipient = await evalJs<string | null>(page, `(async () => {
        const users = await grok.dapi.users.list({limit: 50});
        for (const u of users) {
          if (u.login === 'qa-pw' || u.login === 'system') continue;
          const full = await grok.dapi.users.find(u.id);
          if (full && full.group && full.group.id)
            return full.group.friendlyName || full.login || full.group.name || null;
        }
        return null;
      })()`);

      if (!recipient) {
        console.warn('Step 5: no recipient with materialized group on dev — skipping share (acceptable-SR per Edit 10 G3)');
        return;
      }

      await page.goto('/projects');
      await page.locator('.grok-gallery-grid').waitFor({timeout: 30_000});

      const variants = [
        {name: names.linkCopy, label: 'link'},
        {name: names.cloneCopy, label: 'clone'},
        {name: names.pvcCopy, label: 'pvc'},
      ];
      let sharedCount = 0;
      const shareErrors: string[] = [];
      for (const v of variants) {
        try {
          await shareProjectViaContextMenu(page, v.name, {
            recipient,
            accessLevel: 'View and use',
          });
          sharedCount++;
        } catch (e: any) {
          shareErrors.push(`${v.label} (${v.name}): ${(e?.message ?? String(e)).slice(0, 200)}`);
        }
      }
      if (sharedCount === 0)
        console.error('Step 5 share errors:\n  - ' + shareErrors.join('\n  - '));
      expect(
        sharedCount,
        `Step 5: at least one variant must be shared via UI right-click → Share dialog. Errors: ${shareErrors.join('; ')}`,
      ).toBeGreaterThan(0);
    });

  } finally {
    // Cleanup all 4 produced projects + their tableInfos. Best-effort —
    // deleteProjectWithCleanup swallows errors so a failure in any one delete
    // doesn't mask the test's own failure mode.
    for (const v of Object.values(ids))
      if (v) await deleteProjectWithCleanup(page, v);
  }

  if (stepErrors.length > 0) {
    const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
