/* ---
sub_features_covered: [biostructure.cell-renderer.molecule3d, biostructure.grid-context-menu.copy-raw, biostructure.grid-context-menu.download-raw, biostructure.grid-context-menu.show-biostructure-viewer]
--- */
// GROK-14552: detectors.js context-menu hook must null-guard the cell (no "reading 'semType'" on
// whitespace right-click) while still injecting Copy/Download/Show items on populated Molecule3D cells.
// In-package sentinel: window.$biostructureViewer.contextMenuError. Grid overlay = canvas[2] of tv.grid.root.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

test('BiostructureViewer / GROK-14552 grid null-cell right-click safety regression guard', async ({page}) => {
  test.setTimeout(120_000);
  stepErrors.length = 0;

  // pageerror capture supplementing the in-package contextMenuError sentinel.
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
    // SETUP — Build a Molecule3D-column DataFrame from 1bdq.pdb + force-call autostart to wire the hook.
    let setupDiag: any = null;

    await softStep('Setup — Stage in-memory DataFrame with Molecule3D column from 1bdq.pdb + force BiostructureViewer autostart', async () => {
      setupDiag = await page.evaluate(async (path) => {
        const w: any = window;
        const content = await grok.dapi.files.readAsText(path);
        const df = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['1bdq', '1bdq-clone']),
          DG.Column.fromStrings('structure', [content, content]),
        ]);
        // Stage semType + cell.renderer + units before addTableView so the hook's gate is satisfied.
        const col = df.col('structure');
        col.semType = DG.SEMTYPE.MOLECULE3D || 'Molecule3D';
        col.setTag('cell.renderer', 'Molecule3D');
        try { col.meta.units = 'pdb'; } catch (_) { /* meta API variants */ }
        df.name = 'biostructure-bug-grok-14552-fixture';
        const tv = grok.shell.addTableView(df);
        // Wait for semantic-type detection (best-effort; semType already set explicitly).
        await new Promise((resolve) => {
          try {
            const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
            setTimeout(resolve, 5000);
          } catch (_) { resolve(undefined); }
        });
        // Poll for the grid + overlay canvases to render (canvas[2] is the overlay).
        for (let i = 0; i < 60; i++) {
          const gr = tv && tv.grid ? tv.grid.root : null;
          if (gr && gr.querySelectorAll('canvas').length >= 3) break;
          await new Promise((r) => setTimeout(r, 300));
        }
        // Force-call autostart to wire the onContextMenu hook (idempotent on a cold session).
        let autostartCalled = false;
        let autostartErr: string | null = null;
        try {
          const autoFns = DG.Func.find({package: 'BiostructureViewer', name: 'autostart'});
          if (autoFns && autoFns.length > 0) {
            await autoFns[0].apply({}, {processed: true});
            autostartCalled = true;
          }
        } catch (e: any) { autostartErr = String(e && e.message ? e.message : e); }
        // Hook readiness is polled authoritatively by the next step (probe loop); no blind wait here.
        // Initialize the in-package error sentinel that detectors.js writes to in its outer catch.
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        const gridRoot = tv && tv.grid ? tv.grid.root : null;
        const canvases = gridRoot ? Array.from(gridRoot.querySelectorAll('canvas')) : [];
        return {
          rowCount: df.rowCount,
          colCount: df.columns.length,
          structureSemType: col.semType,
          structureCellRenderer: col.getTag('cell.renderer'),
          hasGridDom: !!document.querySelector('[name="viewer-Grid"]'),
          canvasCount: canvases.length,
          sentinelInitialized: w.$biostructureViewer && w.$biostructureViewer.contextMenuError === null,
          autostartCalled,
          autostartErr,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Grid"]').waitFor({timeout: 30_000});

      expect(setupDiag.hasGridDom).toBe(true);
      expect(setupDiag.rowCount).toBe(2);
      expect(setupDiag.colCount).toBe(2);
      expect(setupDiag.structureSemType).toBe('Molecule3D');
      expect(setupDiag.structureCellRenderer).toBe('Molecule3D');
      expect(setupDiag.canvasCount).toBeGreaterThanOrEqual(3);
      expect(setupDiag.sentinelInitialized).toBe(true);
      expect(
        setupDiag.autostartCalled,
        `BiostructureViewer:autostart was not invocable via DG.Func.find. ` +
        `autostartErr=${JSON.stringify(setupDiag.autostartErr)}. ` +
        `The cold-start race makes the contextmenu hook readiness ` +
        `non-deterministic; we require force-autostart to close the race.`,
      ).toBe(true);
    });

    // SETUP — Readiness poll: probe contextmenu up to 8 times until the BSV "Copy" label appears.
    let probeDiag: any = null;

    await softStep('Setup — Readiness poll: probe contextmenu until BSV hook injects "Copy"', async () => {
      probeDiag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid' };
        const gridRoot = tv.grid.root;
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const gridRect = gridRoot.getBoundingClientRect();
        const probeX = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const probeY = gridRect.top + sb.y + sb.height / 2;

        const observations: any[] = [];
        let ready = false;
        for (let i = 0; i < 8; i++) {
          // Reset sentinel + dismiss any prior menu before each probe.
          if (!w.$biostructureViewer) w.$biostructureViewer = {};
          w.$biostructureViewer.contextMenuError = null;
          document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 300));

          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: probeX, clientY: probeY,
          });
          overlay.dispatchEvent(evt);
          await new Promise((r) => setTimeout(r, 1200));

          const labels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          const hasCopy = labels.includes('Copy');
          observations.push({attempt: i + 1, labelCount: labels.length, hasCopy});
          if (hasCopy) { ready = true; break; }
        }

        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return { ready, observations };
      });

      expect(
        probeDiag.ready,
        `BiostructureViewer hook readiness probe failed: after 8 probe ` +
        `contextmenu dispatches on a populated Molecule3D cell (1s waits), ` +
        `the "Copy" label never appeared in the rendered menu. This ` +
        `indicates the package's autostart subscribe to grok.events.onContextMenu ` +
        `is not wired despite the explicit force-call in the prior step. ` +
        `Observations: ${JSON.stringify(probeDiag.observations)}. ` +
        `Regression-direction: equivalent to the over-application class of ` +
        `the GROK-14552 fix (the null-cell guard accidentally suppresses ` +
        `populated-cell menu injection too). See ` +
        `bug-library/biostructureviewer.yaml#GROK-14552.`,
      ).toBe(true);
    });

    // SCENARIO 1 — Right-click grid whitespace (null cell) must not throw "reading 'semType'".
    let scenario1Diag: any = null;

    await softStep('Scenario 1 step 4 — Right-click row whitespace (past last column); hook MUST NOT throw', async () => {
      // Reset pageErrors and the in-package sentinel just before the
      // load-bearing dispatch.
      pageErrors.length = 0;
      await page.evaluate(() => {
        const w: any = window;
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
      });

      scenario1Diag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid in scenario1' };
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };

        // Whitespace coordinates: aim ~200px past the last column.
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const whitespaceX = gridRect.left + sb.x + sb.width + 200;
        const whitespaceY = gridRect.top + sb.y + Math.min(15, sb.height / 4);

        const evt = new PointerEvent('contextmenu', {
          cancelable: true, bubbles: true, view: window, button: 2,
          clientX: whitespaceX, clientY: whitespaceY,
        });
        overlay.dispatchEvent(evt);

        // Poll for a positive completion signal (menu rendered) or a late async error,
        // so a deferred contextMenuError write is still caught within the budget.
        for (let i = 0; i < 15; i++) {
          await new Promise((r) => setTimeout(r, 200));
          if (w.$biostructureViewer && w.$biostructureViewer.contextMenuError) break;
          if (document.querySelector('.d4-menu-popup')) break;
        }

        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;
        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          whitespaceCoords: [whitespaceX, whitespaceY],
          gridRectShape: [gridRect.left, gridRect.top, gridRect.width, gridRect.height],
        };
      });

      // Assertion #1: the in-package error sentinel must be null (pre-fix it held the semType TypeError).
      expect(
        scenario1Diag.contextMenuErrorIsNull,
        'GROK-14552 grid null-cell right-click crash regressed: ' +
        'window.$biostructureViewer.contextMenuError populated after a ' +
        'whitespace right-click. Captured: ' +
        `${JSON.stringify(scenario1Diag.contextMenuErrorMessage)}. ` +
        'addContextMenuForCell did not null-guard the cell argument. ' +
        'See bug-library/biostructureviewer.yaml#GROK-14552.',
      ).toBe(true);

      // Assertion #2 (independent): pageerror capture must not contain the bug signature.
      const semTypeRegressionSignatures = pageErrors.filter((m) =>
        /Cannot read properties of null.*semType/i.test(m) ||
        /(null|undefined).*semType/i.test(m) ||
        /semType.*(null|undefined)/i.test(m),
      );
      expect(
        semTypeRegressionSignatures,
        'GROK-14552 regression signature surfaced via page.on(pageerror): ' +
        `${JSON.stringify(semTypeRegressionSignatures)}. ` +
        'The BiostructureViewer grid-cell context-menu hook raised a ' +
        '"Cannot read properties of null (reading semType)" TypeError on ' +
        'a row-whitespace right-click — the exact pre-fix bug shape.',
      ).toEqual([]);
    });

    // SCENARIO 2 — Right-click a populated Molecule3D cell must still inject Copy/Download/Show/Biostructure/NGL.
    let scenario2Diag: any = null;

    await softStep('Scenario 2 step 3 — Right-click populated Molecule3D cell; hook MUST inject menu items', async () => {
      // Reset state before the populated-cell dispatch.
      pageErrors.length = 0;
      await page.evaluate(() => {
        const w: any = window;
        if (!w.$biostructureViewer) w.$biostructureViewer = {};
        w.$biostructureViewer.contextMenuError = null;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
      });

      scenario2Diag = await page.evaluate(async () => {
        const w: any = window;
        const tv = w.grok.shell.tv;
        if (!tv || !tv.grid) return { err: 'no grid in scenario2' };
        const gridRoot = tv.grid.root;
        const gridRect = gridRoot.getBoundingClientRect();
        const canvases = Array.from(gridRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
        const overlay = canvases[2];
        if (!overlay) return { err: 'no overlay canvas' };

        // Aim for the centre of the populated structure cell on row 0.
        const cell = tv.grid.cell('structure', 0);
        const sb = cell.bounds;
        const populatedX = gridRect.left + sb.x + Math.min(20, sb.width / 4);
        const populatedY = gridRect.top + sb.y + sb.height / 2;

        // Retry up to 3 times — a single dispatch can miss the render commit.
        let menuLabels: string[] = [];
        let hasCopy = false, hasDownload = false, hasShow = false, hasBiostructure = false, hasNgl = false;
        let attemptCount = 0;
        for (let attempt = 0; attempt < 3; attempt++) {
          attemptCount = attempt + 1;
          document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
          document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
          await new Promise((r) => setTimeout(r, 400));

          const evt = new PointerEvent('contextmenu', {
            cancelable: true, bubbles: true, view: window, button: 2,
            clientX: populatedX, clientY: populatedY,
          });
          overlay.dispatchEvent(evt);
          // Advance as soon as the menu renders instead of a blind 2s wait.
          for (let k = 0; k < 12; k++) {
            await new Promise((r) => setTimeout(r, 200));
            if (document.querySelectorAll('.d4-menu-popup .d4-menu-item-label').length > 0) break;
          }

          menuLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .map((el) => (el.textContent || '').trim());
          hasCopy = menuLabels.includes('Copy');
          hasDownload = menuLabels.includes('Download');
          hasShow = menuLabels.includes('Show');
          hasBiostructure = menuLabels.includes('Biostructure');
          hasNgl = menuLabels.includes('NGL');
          if (hasCopy && hasDownload && hasShow && hasBiostructure) break;
        }

        // Hover the "Show" group to expand a deferred submenu if needed (usually inline).
        let showSubmenuExpanded = false;
        if (!hasBiostructure || !hasNgl) {
          const showLabelEl = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
            .find((el) => (el.textContent || '').trim() === 'Show');
          const showGroupEl = showLabelEl ? showLabelEl.closest('.d4-menu-item') : null;
          if (showGroupEl) {
            const r = (showGroupEl as HTMLElement).getBoundingClientRect();
            showGroupEl.dispatchEvent(new MouseEvent('mouseover', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            showGroupEl.dispatchEvent(new MouseEvent('mouseenter', {
              bubbles: true, clientX: r.left + 10, clientY: r.top + 10,
            }));
            // Poll for the submenu to expand rather than a fixed 800ms wait.
            for (let k = 0; k < 8; k++) {
              await new Promise((rr) => setTimeout(rr, 100));
              if (Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
                .some((el) => ['Biostructure', 'NGL'].includes((el.textContent || '').trim()))) break;
            }
            const afterHoverLabels = Array.from(document.querySelectorAll('.d4-menu-popup .d4-menu-item-label'))
              .map((el) => (el.textContent || '').trim());
            showSubmenuExpanded = afterHoverLabels.includes('Biostructure') ||
                                  afterHoverLabels.includes('NGL');
            if (showSubmenuExpanded) {
              hasBiostructure = hasBiostructure || afterHoverLabels.includes('Biostructure');
              hasNgl = hasNgl || afterHoverLabels.includes('NGL');
            }
          }
        }

        const captured = w.$biostructureViewer && w.$biostructureViewer.contextMenuError;

        // Dismiss menus before returning.
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });

        return {
          contextMenuErrorIsNull: captured === null || captured === undefined,
          contextMenuErrorMessage: captured ? String(captured.message || captured) : null,
          hasCopy, hasDownload, hasShow, hasBiostructure, hasNgl,
          menuItemCount: menuLabels.length,
          menuItemsSample: menuLabels.filter((t) => t.length > 0 && t.length < 50).slice(0, 30),
          populatedCoords: [populatedX, populatedY],
          attemptCount,
        };
      });

      // Assertion #1: the hook sentinel must remain clean on a populated cell.
      expect(
        scenario2Diag.contextMenuErrorIsNull,
        'GROK-14552 hook threw on a POPULATED Molecule3D cell — ' +
        'window.$biostructureViewer.contextMenuError populated. Captured: ' +
        `${JSON.stringify(scenario2Diag.contextMenuErrorMessage)}.`,
      ).toBe(true);

      // Assertion #2 (inverse guard): the four BiostructureViewer menu entries must be injected.
      expect(
        scenario2Diag.hasCopy,
        'GROK-14552 fix over-applied: the BiostructureViewer "Copy" ' +
        'context-menu item is missing on a populated Molecule3D cell. ' +
        `Menu item count: ${scenario2Diag.menuItemCount}. ` +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}. ` +
        `Attempts: ${scenario2Diag.attemptCount}. ` +
        'A correct null-cell guard MUST NOT regress the populated-cell ' +
        'menu injection (detectors.js L76).',
      ).toBe(true);
      expect(
        scenario2Diag.hasDownload,
        'GROK-14552 fix over-applied: the BiostructureViewer "Download" ' +
        'context-menu item is missing on a populated Molecule3D cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario2Diag.hasShow,
        'GROK-14552 fix over-applied: the BiostructureViewer "Show" ' +
        'context-menu group is missing on a populated Molecule3D cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);
      expect(
        scenario2Diag.hasBiostructure,
        'GROK-14552 fix over-applied: the BiostructureViewer "Show > ' +
        'Biostructure" submenu item is missing on a populated Molecule3D ' +
        'cell. ' +
        `Menu items sample: ${JSON.stringify(scenario2Diag.menuItemsSample)}.`,
      ).toBe(true);

      const populatedRegressionSignatures = pageErrors.filter((m) =>
        /Cannot read properties of null.*semType/i.test(m) ||
        /semType.*(null|undefined)/i.test(m),
      );
      expect(
        populatedRegressionSignatures,
        'GROK-14552 regression signature surfaced during populated-cell ' +
        `right-click: ${JSON.stringify(populatedRegressionSignatures)}.`,
      ).toEqual([]);
    });

    // SCENARIO 2 step 4 — Joint invariant cross-check (log-only summary).
    await softStep('Scenario 2 step 4 — Joint invariant cross-check (GROK-14552)', async () => {
      const summary = {
        scenario1: {
          contextMenuErrorIsNull: scenario1Diag ? scenario1Diag.contextMenuErrorIsNull : null,
          whitespaceCoords: scenario1Diag ? scenario1Diag.whitespaceCoords : null,
        },
        scenario2: {
          contextMenuErrorIsNull: scenario2Diag ? scenario2Diag.contextMenuErrorIsNull : null,
          hasCopy: scenario2Diag ? scenario2Diag.hasCopy : null,
          hasDownload: scenario2Diag ? scenario2Diag.hasDownload : null,
          hasShow: scenario2Diag ? scenario2Diag.hasShow : null,
          hasBiostructure: scenario2Diag ? scenario2Diag.hasBiostructure : null,
          populatedCoords: scenario2Diag ? scenario2Diag.populatedCoords : null,
        },
        jointInvariantHolds: !!(
          scenario1Diag && scenario1Diag.contextMenuErrorIsNull &&
          scenario2Diag && scenario2Diag.contextMenuErrorIsNull &&
          scenario2Diag.hasCopy && scenario2Diag.hasDownload &&
          scenario2Diag.hasShow && scenario2Diag.hasBiostructure
        ),
      };
      // eslint-disable-next-line no-console
      console.log(`[GROK-14552 joint-invariant summary] ${JSON.stringify(summary)}`);
      expect(summary.jointInvariantHolds).toBe(true);
    });
  } finally {
    // Cleanup — close menus/dialogs, reset the sentinel, close views.
    try {
      await page.evaluate(() => {
        const w: any = window;
        document.dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape' }));
        document.querySelectorAll('.d4-menu-popup').forEach((m) => { try { m.remove(); } catch (_) {} });
        document.querySelectorAll('.d4-dialog').forEach((d) => {
          const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
          if (cancel) cancel.click();
        });
        if (w.$biostructureViewer) w.$biostructureViewer.contextMenuError = null;
        try { w.grok.shell.closeAll(); } catch (_) { /* best-effort */ }
      });
    } catch (e) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
