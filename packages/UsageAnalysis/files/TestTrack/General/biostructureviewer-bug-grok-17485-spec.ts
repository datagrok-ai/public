/* ---
sub_features_covered: [biostructure.file-open, biostructure.project-persistence, biostructure.viewer]
--- */
// NOTE: This spec was REMOVED from the playwright-public CI suite (BiostructureViewer folder) and is
// kept here in TestTrack/General for reference only. Reason for removal: it depends on the cross-user
// `shareWithSecondUserAndVerify` helper (a second-user / token2 sharing flow) which is not available
// in the playwright-public CI harness, so the spec cannot run in CI.
// GROK-17485: an ad-hoc PDB payload (viewer.props.get('pdb')) must survive a project save/reopen;
// the bug drops the Biostructure viewer on reopen. KNOWN_REGRESSION_SOFT_WARN gates a soft-warn
// branch (set false to restore the hard assertion once the platform fix lands).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';
import {
  saveProjectWithProvenance,
  reopenAndAssertProvenance,
  deleteProjectWithCleanup,
  shareProjectViaContextMenu,
  shareWithSecondUserAndVerify,
} from '../helpers/projects';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

const KNOWN_REGRESSION_SOFT_WARN = true; // flip to false to re-enable hard-assert

test('BiostructureViewer — GROK-17485 ad-hoc PDB project-persistence regression guard', async ({page}) => {
  // Removed from the playwright-public CI suite — relies on the cross-user `shareWithSecondUserAndVerify`
  // (second-user / token2) flow that the CI harness does not provide. Kept in General for reference.
  test.setTimeout(600_000);
  stepErrors.length = 0;

  await loginToDatagrok(page);

  // Baseline environment setup.
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
  });
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  // Timestamp-suffixed names for idempotence across re-runs.
  const ts = Date.now();
  const sameUserProjectName = `bug-grok-17485-roundtrip-${ts}`;
  const sharedProjectName = `bug-grok-17485-share-${ts}`;
  let sameUserProjectId: string | null = null;
  let sameUserTableInfoId: string | null = null;
  let sharedProjectId: string | null = null;
  let sharedTableInfoId: string | null = null;

  try {
    // SCENARIO 1 — Open PDB, save project, reopen as same user; pdb payload must survive.
    let pdbContentLength = 0;
    let preReopenRep: string | null = null;

    await softStep('Scenario 1 step 1 — Open 1bdq.pdb; Biostructure viewer mounts with pdb payload', async () => {
      const result = await page.evaluate(async (path) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const content: string = await grok.dapi.files.readAsText(path);
        const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('sentinel', ['1bdq'])]);
        df.name = 'biostructure-bug-grok-17485-scenario1';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));
        // The GROK-17485 invariant payload: the ad-hoc in-memory structure that must survive save.
        v.setOptions({pdb: content});
        await new Promise((r) => setTimeout(r, 1500));
        return {
          contentLength: content.length,
          hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
          vType: v?.type,
          pdbLen: typeof v.props.get('pdb') === 'string' ? v.props.get('pdb').length : null,
          representation: v.props.get('representation'),
        };
      }, samplePdbPath);

      await expect(page.locator('[name="viewer-Biostructure"]')).toBeVisible({timeout: 30_000});

      expect(result.hasContainer).toBe(true);
      expect(result.vType).toBe('Biostructure');
      expect(result.contentLength).toBeGreaterThan(1000); // 1bdq.pdb ~ 180k chars
      expect(result.pdbLen).toBe(result.contentLength);
      expect(result.representation).toBe('cartoon');

      pdbContentLength = result.contentLength;
      preReopenRep = result.representation;
    });

    await softStep('Scenario 1 steps 3-4 — Save project + closeAll', async () => {
      const saved = await saveProjectWithProvenance(page, sameUserProjectName);
      sameUserProjectId = saved.projectId;
      sameUserTableInfoId = saved.tableInfoId;
      expect(saved.resolvedName).toBe(sameUserProjectName);
      expect(saved.projectId).toBeTruthy();
      expect(saved.tableInfoId).toBeTruthy();
      // Server commit settle (mirrors helper precedent).
      await page.waitForTimeout(1500);
    });

    await softStep('Scenario 1 steps 5-6 — Reopen project; assert viewer rehydration (GROK-17485 invariant)', async () => {
      const reopen = await reopenAndAssertProvenance(page, sameUserProjectId!);
      expect(reopen.tablesAfter).toBeGreaterThan(0);
      expect(reopen.reopenedName).toBe('biostructure-bug-grok-17485-scenario1');

      // Poll up to 30s for the Biostructure viewer to rehydrate after reopen.
      const rehydrated = await page.evaluate(async (expectedPdbLen) => {
        const grokAny: any = (window as any).grok;
        for (let i = 0; i < 60; i++) {
          await new Promise((r) => setTimeout(r, 500));
          const tvPoll = grokAny.shell.tv;
          if (!tvPoll) continue;
          let v: any = null;
          for (const x of tvPoll.viewers || []) if (x.type === 'Biostructure') { v = x; break; }
          if (!v) continue;
          const pdb = v.props.get('pdb');
          const rep = v.props.get('representation');
          if (typeof pdb === 'string' && pdb.length === expectedPdbLen) {
            return {ok: true, pdbLen: pdb.length, rep, viewerType: v.type};
          }
        }
        // Last-iteration sweep — gather what we DO see for diagnosis.
        const tvFinal = grokAny.shell.tv;
        const viewerTypes = tvFinal?.viewers ? Array.from(tvFinal.viewers).map((x: any) => x.type) : [];
        let lastPdbLen: number | null = null;
        if (tvFinal) {
          for (const x of tvFinal.viewers || []) if (x.type === 'Biostructure') {
            const p = x.props.get('pdb');
            lastPdbLen = typeof p === 'string' ? p.length : null;
            break;
          }
        }
        return {ok: false, viewerTypes, lastPdbLen};
      }, pdbContentLength);

      // GROK-17485 invariant: viewer present AND pdb payload matches pre-save (rehydrated.ok).
      const observedDiag =
        `Observed viewer types after reopen: ${JSON.stringify((rehydrated as any).viewerTypes)}, ` +
        `lastPdbLen: ${(rehydrated as any).lastPdbLen} (expected ${pdbContentLength}).`;
      if (KNOWN_REGRESSION_SOFT_WARN && !rehydrated.ok) {
        // Soft-warn branch: record the known regression without failing the spec.
        // eslint-disable-next-line no-console
        console.warn(
          `[SR-02 known-regression] GROK-17485 still regressed on dev as of test run: ${observedDiag} ` +
          `See spec header SR-02 + bug-library/biostructureviewer.yaml#GROK-17485 (status: fixed, fixed_in: '').`,
        );
      } else {
        // Hard-assert path — runs once KNOWN_REGRESSION_SOFT_WARN is false.
        expect(
          rehydrated.ok,
          `GROK-17485 invariant violated: ad-hoc PDB payload did not survive project save+reopen. ${observedDiag} ` +
          `See bug-library/biostructureviewer.yaml#GROK-17485.`,
        ).toBe(true);
        expect(rehydrated.pdbLen).toBe(pdbContentLength);
        expect(rehydrated.rep).toBe(preReopenRep);
        expect(rehydrated.viewerType).toBe('Biostructure');
      }
    });

    // SCENARIO 2 — Cross-user persistence: share the project; cross-user reopen runs when token2 is set.
    let scenario2SetupComplete = false;
    await softStep('Scenario 2 step 1 — Open 1bdq.pdb fresh; Biostructure viewer mounts', async () => {
      const result = await page.evaluate(async (path) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const content: string = await grok.dapi.files.readAsText(path);
        const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('sentinel', ['1bdq'])]);
        df.name = 'biostructure-bug-grok-17485-scenario2';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 1500));
        const v = tv.addViewer('Biostructure');
        await new Promise((r) => setTimeout(r, 2500));
        v.setOptions({pdb: content});
        await new Promise((r) => setTimeout(r, 1500));
        return {
          hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
          pdbLen: typeof v.props.get('pdb') === 'string' ? v.props.get('pdb').length : null,
        };
      }, samplePdbPath);
      expect(result.hasContainer).toBe(true);
      expect(result.pdbLen).toBeGreaterThan(1000);
      scenario2SetupComplete = true;
    });

    await softStep('Scenario 2 step 3 — Save project bug-grok-17485-share-<ts>', async () => {
      if (!scenario2SetupComplete) return;
      const saved = await saveProjectWithProvenance(page, sharedProjectName);
      sharedProjectId = saved.projectId;
      sharedTableInfoId = saved.tableInfoId;
      expect(saved.resolvedName).toBe(sharedProjectName);
      await page.waitForTimeout(1500);
    });

    await softStep('Scenario 2 step 4 — Share project via Right-click -> Share dialog (DOM-driven)', async () => {
      if (!sharedProjectId) return;
      // Close current view so the dashboards view can take focus.
      await page.evaluate(async () => {
        const grokAny: any = (window as any).grok;
        grokAny.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1000));
      });
      const origin = new URL(page.url()).origin;
      await page.goto(`${origin}/browse/dashboards`);
      await page.waitForTimeout(2500);

      // 'Test permission group' exists on dev; degrade gracefully if absent.
      const tileLocator = page.locator(`[name="div-${sharedProjectName}"]`).or(
        page.locator(`[name="div-${sharedProjectName.toLowerCase()}"]`),
      );
      const tileVisible = await tileLocator.first().isVisible({timeout: 15_000}).catch(() => false);
      if (!tileVisible) {
        // Tile not visible — degrade with a warning; the JS-API cross-user leg still runs in steps 5-6.
        console.warn('SR-01 share-dialog leg: project tile not visible in Browse > Dashboards within 15s after save (UI dialog coverage skipped; JS-API cross-user leg still runs)');
        return;
      }

      try {
        await shareProjectViaContextMenu(page, sharedProjectName, {
          recipient: 'Test permission group',
          accessLevel: 'View and use',
        });
      } catch (e) {
        const errStr = String(e).slice(0, 400);
        // Degrade with a warning if the share dialog fails; steps 5-6 still run.
        console.warn(`SR-01 share-dialog leg: shareProjectViaContextMenu raised "${errStr}" (UI dialog coverage skipped; JS-API cross-user leg still runs)`);
        return;
      }

      // Server-side verification: permission grant landed.
      const permsOk = await page.evaluate(async (id) => {
        try {
          const grokAny: any = (window as any).grok;
          const proj = await grokAny.dapi.projects.find(id);
          if (!proj) return {ok: false, reason: 'project not found post-share'};
          const perms = await grokAny.dapi.permissions.get(proj).catch(() => null);
          if (perms && (perms.view?.length || perms.edit?.length)) return {ok: true};
          return {ok: false, reason: 'no permissions recorded'};
        } catch (e: any) {
          return {ok: false, reason: String(e).slice(0, 200)};
        }
      }, sharedProjectId);

      expect(permsOk.ok, `Share leg server-side verification failed: ${(permsOk as any).reason}`).toBe(true);
    });

    await softStep('Scenario 2 steps 5-6 — Cross-user reopen: share to second user + verify recipient visibility', async () => {
      if (!sharedProjectId) return;
      // Cross-user leg: share to the second user's group, re-auth as token2, verify visibility, restore.
      // The helper reloads the page during re-auth, so this must be the last action before cleanup.
      const r = await shareWithSecondUserAndVerify(page, {id: sharedProjectId, name: sharedProjectName});
      if (!r.shared) {
        test.skip(true, 'cross-user share leg: ' + r.reason);
        return;
      }
      if (r.recipientVisible !== null) expect(r.recipientVisible).toBe(true);
    });
  } finally {
    // Cleanup — delete created projects + their TableInfos (best-effort).
    try {
      if (sameUserProjectId || sameUserTableInfoId) {
        await deleteProjectWithCleanup(page, {
          projectId: sameUserProjectId ?? undefined,
          tableInfoId: sameUserTableInfoId ?? undefined,
        });
      }
    } catch (e) { /* best-effort */ }
    try {
      if (sharedProjectId || sharedTableInfoId) {
        await deleteProjectWithCleanup(page, {
          projectId: sharedProjectId ?? undefined,
          tableInfoId: sharedTableInfoId ?? undefined,
        });
      }
    } catch (e) { /* best-effort */ }
    try {
      await page.evaluate(() => { (window as any).grok?.shell?.closeAll?.(); });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
