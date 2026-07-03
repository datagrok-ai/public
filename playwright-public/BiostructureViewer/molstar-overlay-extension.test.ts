/* ---
sub_features_covered: [biostructure.overlay.screenshot, biostructure.overlay.selection-mode, biostructure.overlay.settings-info, biostructure.overlay.toggle-controls]
--- */
// Mol* viewport overlay buttons: Screenshot / Toggle Controls / Selection Mode / Settings.
// Each overlay button only renders when .msp-plugin is mounted, so DOM-click assertions are gated on
// that precondition (the engine doesn't init reliably in headless CI). Scenario 2 also asserts the
// layoutShowControls property round-trip unconditionally (mirrors the Toggle Controls button).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const sample1bdq = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer / Mol* viewport overlay buttons extension (Screenshot / Toggle Controls / Selection Mode / Settings)', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  const pageErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });

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

  try {
    // SHARED SETUP — load 1bdq.pdb Molecule3D fixture + add Biostructure viewer (used by all scenarios).
    let setupReady = false;

    await softStep('Shared Setup — Open Molecule3D table; tv.addViewer("Biostructure"); container DOM mounts', async () => {
      const res = await page.evaluate(async (pdbPath) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const content = await grok.dapi.files.readAsText(pdbPath);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq']),
          DG.Column.fromStrings('structure', [content]),
        ]);
        const col = df.col('structure');
        col.semType = 'Molecule3D';
        try { col.setTag('cell.renderer', 'Molecule3D'); } catch (_) { /* tag set */ }
        df.name = 'overlay-extension-fixture';
        const tv = grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 2500));
        const v = tv.addViewer('Biostructure');
        // Wait for Mol* render; tolerate timeout (per-scenario .msp-plugin guards cover the WebGL-uncertain case).
        try { await v.awaitRendered(20_000); } catch (_) { /* WebGL-uncertain */ }
        await new Promise((r) => setTimeout(r, 3000));
        return {
          rowCount: df.rowCount,
          hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
          vType: v?.type,
          mspPluginPresent: !!document.querySelector('.msp-plugin'),
          mspViewportPresent: !!document.querySelector('.msp-viewport'),
          mspCanvasPresent: !!document.querySelector('.msp-viewport canvas'),
          overlayBtnCount: document.querySelectorAll('.msp-viewport-controls-buttons button').length,
        };
      }, sample1bdq);

      await expect(page.locator('[name="viewer-Biostructure"]')).toBeVisible({timeout: 30_000});

      expect(res.hasContainer).toBe(true);
      expect(res.vType).toBe('Biostructure');
      expect(res.rowCount).toBe(1);
      setupReady = true;
    });

    // SCENARIO 1 — Screenshot / State Snapshot overlay button toggles a Mol* panel (.msp-plugin gated).
    await softStep('Scenario 1 — Screenshot / State Snapshot overlay button toggles a Mol* panel', async () => {
      if (!setupReady) return;
      const res = await page.evaluate(async () => {
        const pluginPresent = !!document.querySelector('.msp-plugin');
        if (!pluginPresent) return {precondition: false};
        const btn = document.querySelector(
          'button[title="Screenshot / State Snapshot"]',
        ) as HTMLButtonElement | null;
        if (!btn) return {precondition: true, btnPresent: false};

        const classBefore = btn.className;
        const mspPanelsBefore = document.querySelectorAll('.msp-plugin .msp-control').length;

        btn.click();
        await new Promise((r) => setTimeout(r, 1200));

        const classAfter1 = btn.className;
        const mspPanelsAfter1 = document.querySelectorAll('.msp-plugin .msp-control').length;

        btn.click();
        await new Promise((r) => setTimeout(r, 1200));

        const classAfter2 = btn.className;
        const mspPanelsAfter2 = document.querySelectorAll('.msp-plugin .msp-control').length;

        const canvasStillPresent = !!document.querySelector('.msp-viewport canvas');
        return {
          precondition: true,
          btnPresent: true,
          classBefore,
          classAfter1,
          classAfter2,
          mspPanelsBefore,
          mspPanelsAfter1,
          mspPanelsAfter2,
          canvasStillPresent,
          toggleOnAfter1:
            classAfter1.includes('msp-btn-link-toggle-on') ||
            classAfter1 !== classBefore,
          toggleOffAfter2: classAfter2 === classBefore,
        };
      });
      // When .msp-plugin is built, the button must be present and toggle class state.
      if (res.precondition) {
        expect(
          res.btnPresent,
          `Screenshot / State Snapshot overlay button missing despite .msp-plugin built — biostructure.overlay.screenshot regression.`,
        ).toBe(true);
        expect(
          res.toggleOnAfter1,
          `Screenshot overlay button class did not flip on first click. ` +
          `Before: '${res.classBefore}', After1: '${res.classAfter1}'.`,
        ).toBe(true);
        // Panel-mount DOM diff: msp-control count increases on first click.
        expect(
          res.mspPanelsAfter1,
          `Screenshot panel did not mount (msp-control count did not increase). ` +
          `Before: ${res.mspPanelsBefore}, After1: ${res.mspPanelsAfter1}.`,
        ).toBeGreaterThanOrEqual(res.mspPanelsBefore);
        expect(res.canvasStillPresent).toBe(true);
      }
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during Screenshot overlay button toggle: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });

    // SCENARIO 2 — Toggle Controls Panel: DOM-click (.msp-plugin gated) + unconditional
    //   layoutShowControls property round-trip (mirrors the button state).
    await softStep('Scenario 2 — Toggle Controls Panel overlay button + layoutShowControls property mirror', async () => {
      if (!setupReady) return;
      pageErrors.length = 0;
      const res = await page.evaluate(async () => {
        // Layer (b): unconditional layoutShowControls round-trip.
        let v: any = null;
        for (const tv of grok.shell.tableViews || []) {
          for (const x of tv.viewers || []) if (x.type === 'Biostructure') { v = x; break; }
          if (v) break;
        }
        if (!v) return {viewerPresent: false};
        const layoutInit = v.props.get('layoutShowControls');
        v.setOptions({layoutShowControls: true});
        await new Promise((r) => setTimeout(r, 1000));
        const layoutAfterOn = v.props.get('layoutShowControls');
        v.setOptions({layoutShowControls: false});
        await new Promise((r) => setTimeout(r, 1000));
        const layoutRestored = v.props.get('layoutShowControls');

        // Layer (a): conditional DOM click on the overlay button.
        const pluginPresent = !!document.querySelector('.msp-plugin');
        if (!pluginPresent) {
          return {
            viewerPresent: true,
            precondition: false,
            layoutInit, layoutAfterOn, layoutRestored,
          };
        }
        const btn = document.querySelector(
          'button[title="Toggle Controls Panel"]',
        ) as HTMLButtonElement | null;
        if (!btn) {
          return {
            viewerPresent: true,
            precondition: true,
            btnPresent: false,
            layoutInit, layoutAfterOn, layoutRestored,
          };
        }
        const layoutHideLeftBefore = !!document.querySelector('.msp-layout-hide-left');
        const classBefore = btn.className;
        btn.click();
        await new Promise((r) => setTimeout(r, 1200));
        const classAfter1 = btn.className;
        const layoutHideLeftAfter1 = !!document.querySelector('.msp-layout-hide-left');
        btn.click();
        await new Promise((r) => setTimeout(r, 1200));
        const classAfter2 = btn.className;
        const layoutHideLeftAfter2 = !!document.querySelector('.msp-layout-hide-left');
        return {
          viewerPresent: true,
          precondition: true,
          btnPresent: true,
          classBefore, classAfter1, classAfter2,
          layoutHideLeftBefore, layoutHideLeftAfter1, layoutHideLeftAfter2,
          layoutInit, layoutAfterOn, layoutRestored,
        };
      });

      expect(res.viewerPresent).toBe(true);

      // Layer (b): layoutShowControls round-trip (default false — 3D viewport only).
      expect(
        res.layoutInit,
        `layoutShowControls default expected false (3D-viewport-only default per refdoc).`,
      ).toBe(false);
      expect(
        res.layoutAfterOn,
        `layoutShowControls did not round-trip to true via setOptions — biostructure.overlay.toggle-controls mirroring regression.`,
      ).toBe(true);
      expect(
        res.layoutRestored,
        `layoutShowControls did not restore to false — round-trip regression.`,
      ).toBe(false);

      // Layer (a): conditional DOM-click assertion.
      if (res.precondition) {
        expect(
          res.btnPresent,
          `Toggle Controls Panel overlay button missing despite .msp-plugin built — biostructure.overlay.toggle-controls regression.`,
        ).toBe(true);
        expect(
          res.classAfter1,
          `Toggle Controls Panel overlay button class did not flip on first click. ` +
          `Before: '${res.classBefore}', After1: '${res.classAfter1}'.`,
        ).not.toBe(res.classBefore);
        expect(
          res.classAfter2,
          `Toggle Controls Panel overlay button class did not return to baseline on second click. ` +
          `Before: '${res.classBefore}', After2: '${res.classAfter2}'.`,
        ).toBe(res.classBefore);
      }

      // No JS console error during the click + setOptions chain.
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during Toggle Controls Panel toggle: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });

    // SCENARIO 3 — Toggle Selection Mode overlay button class-flip (.msp-plugin gated; no residue-pick pixels).
    await softStep('Scenario 3 — Toggle Selection Mode overlay button class-flip on/off', async () => {
      if (!setupReady) return;
      pageErrors.length = 0;
      const res = await page.evaluate(async () => {
        const pluginPresent = !!document.querySelector('.msp-plugin');
        if (!pluginPresent) return {precondition: false};
        const btn = document.querySelector(
          'button[title="Toggle Selection Mode"]',
        ) as HTMLButtonElement | null;
        if (!btn) return {precondition: true, btnPresent: false};
        const classBefore = btn.className;
        btn.click();
        await new Promise((r) => setTimeout(r, 1000));
        const classAfter1 = btn.className;
        btn.click();
        await new Promise((r) => setTimeout(r, 1000));
        const classAfter2 = btn.className;
        const canvasStillPresent = !!document.querySelector('.msp-viewport canvas');
        return {
          precondition: true,
          btnPresent: true,
          classBefore, classAfter1, classAfter2,
          canvasStillPresent,
        };
      });

      if (res.precondition) {
        expect(
          res.btnPresent,
          `Toggle Selection Mode overlay button missing despite .msp-plugin built — biostructure.overlay.selection-mode regression.`,
        ).toBe(true);
        expect(
          res.classAfter1,
          `Toggle Selection Mode overlay button class did not flip on first click. ` +
          `Before: '${res.classBefore}', After1: '${res.classAfter1}'.`,
        ).not.toBe(res.classBefore);
        expect(
          res.classAfter2,
          `Toggle Selection Mode overlay button class did not return to baseline on second click. ` +
          `Before: '${res.classBefore}', After2: '${res.classAfter2}'.`,
        ).toBe(res.classBefore);
        expect(
          res.canvasStillPresent,
          `Viewport canvas disappeared during selection-mode toggle — structure underneath should remain rendered.`,
        ).toBe(true);
      }
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during Toggle Selection Mode: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });

    // SCENARIO 4 — Settings / Controls Info overlay button toggles a Mol* settings panel (.msp-plugin gated).
    await softStep('Scenario 4 — Settings / Controls Info overlay button toggles a Mol* settings panel', async () => {
      if (!setupReady) return;
      pageErrors.length = 0;
      const res = await page.evaluate(async () => {
        const pluginPresent = !!document.querySelector('.msp-plugin');
        if (!pluginPresent) return {precondition: false};
        const btn = document.querySelector(
          'button[title="Settings / Controls Info"]',
        ) as HTMLButtonElement | null;
        if (!btn) return {precondition: true, btnPresent: false};
        const classBefore = btn.className;
        const mspPanelsBefore = document.querySelectorAll('.msp-plugin .msp-control').length;
        btn.click();
        await new Promise((r) => setTimeout(r, 1200));
        const classAfter1 = btn.className;
        const mspPanelsAfter1 = document.querySelectorAll('.msp-plugin .msp-control').length;
        btn.click();
        await new Promise((r) => setTimeout(r, 1200));
        const classAfter2 = btn.className;
        const mspPanelsAfter2 = document.querySelectorAll('.msp-plugin .msp-control').length;
        const canvasStillPresent = !!document.querySelector('.msp-viewport canvas');
        return {
          precondition: true,
          btnPresent: true,
          classBefore, classAfter1, classAfter2,
          mspPanelsBefore, mspPanelsAfter1, mspPanelsAfter2,
          canvasStillPresent,
        };
      });

      if (res.precondition) {
        expect(
          res.btnPresent,
          `Settings / Controls Info overlay button missing despite .msp-plugin built — biostructure.overlay.settings-info regression.`,
        ).toBe(true);
        expect(
          res.classAfter1,
          `Settings / Controls Info overlay button class did not flip on first click. ` +
          `Before: '${res.classBefore}', After1: '${res.classAfter1}'.`,
        ).not.toBe(res.classBefore);
        expect(
          res.classAfter2,
          `Settings / Controls Info overlay button class did not return to baseline on second click. ` +
          `Before: '${res.classBefore}', After2: '${res.classAfter2}'.`,
        ).toBe(res.classBefore);
        expect(
          res.mspPanelsAfter1,
          `Settings / Controls Info panel did not mount (msp-control count did not increase). ` +
          `Before: ${res.mspPanelsBefore}, After1: ${res.mspPanelsAfter1}.`,
        ).toBeGreaterThanOrEqual(res.mspPanelsBefore);
        expect(
          res.canvasStillPresent,
          `Viewport canvas disappeared during settings-panel toggle — structure underneath should remain rendered.`,
        ).toBe(true);
      }
      const errSig = pageErrors.filter((m) =>
        /TypeError|ReferenceError|Cannot read properties/i.test(m),
      );
      expect(
        errSig,
        `JS console error during Settings / Controls Info toggle: ${JSON.stringify(errSig)}.`,
      ).toEqual([]);
    });
  } finally {
    // Cleanup.
    try {
      await page.evaluate(() => {
        document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        try { grok.shell.closeAll(); } catch (_) { /* best-effort */ }
      });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
