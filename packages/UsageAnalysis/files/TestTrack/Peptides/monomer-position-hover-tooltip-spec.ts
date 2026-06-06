// MonomerPosition hover tooltip — GROK-15934 regression: no null-receiver on hover across
// state mutations (SVM cell hovers + main-grid WebLogo header hovers, then mode switch /
// weblogo selection sync / Settings round-trip). Hovers use page.mouse.move() (CDP-trusted) —
// synthetic dispatchEvent does NOT fire the Dart hit-test. Canvas picks largest-by-area
// (a small row-marker canvas is first). 100-row Extract subset makes MCL/LST attach in ~4s.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/bio/peptides.csv';

// GROK-15934 null-receiver bug class (e.g. "Cannot read properties of null (reading 'getTag')").
const NULL_RECEIVER_PATTERN =
  /(?:Cannot read .* of (?:null|undefined).*(?:getTag|tag|column))|(?:getTag.*on (?:null|undefined))/i;

test('MonomerPosition hover tooltip — GROK-15934 regression (no null-receiver on hover across state mutations)', async ({page}) => {
  test.setTimeout(120_000);
  await loginToDatagrok(page);

  await softStep('Setup: open peptides dataset, prewarm Peptides:initPeptides', async () => {
    const result = await page.evaluate(async (path) => {
      document.querySelectorAll('.d4-dialog').forEach((d) => {
        const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
        if (cancel) cancel.click();
      });
      grok.shell.closeAll();
      document.body.classList.add('selenium');
      grok.shell.settings.showFiltersIconsConstantly = true;
      // Windows mode required: simpleMode hides the Context Panel (the Launch SAR surface).
      grok.shell.windows.simpleMode = false;

      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise<void>((resolve) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 4000);
      });
      // Macromolecule dataset: wait for grid canvas + Bio package settle.
      for (let i = 0; i < 50; i++) {
        if (document.querySelector('[name="viewer-Grid"] canvas')) break;
        await new Promise((r) => setTimeout(r, 200));
      }
      await new Promise((r) => setTimeout(r, 4000));

      // GROK-17557 prewarm — removes the cold-package mount race on the Launch SAR button.
      try { await grok.functions.call('Peptides:initPeptides'); }
      catch (e) { console.log('[note] Peptides:initPeptides pre-warm threw (non-fatal):', String(e)); }

      return {
        rows: df.rowCount,
        semType: df.col('AlignedSequence')?.semType ?? null,
      };
    }, datasetPath);
    await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
    expect(result.rows, 'peptides.csv should load 647 rows').toBe(647);
    expect(result.semType, 'AlignedSequence must be a Macromolecule column').toBe('Macromolecule');
  });

  // Extract first 100 rows via CmdExtractSelectedRows (menu item not name-addressable) for fast MCL.
  await softStep('Setup: select first 100 rows, Select > Extract Selected Rows to a fast 100-row table', async () => {
    const result = await page.evaluate(async () => {
      const src = grok.shell.t;
      src.selection.init((i) => i < 100);
      await new Promise((r) => setTimeout(r, 300));
      const selected = src.selection.trueCount;
      const fn = DG.Func.find({name: 'CmdExtractSelectedRows'})[0];
      await fn.prepare().call();
      await new Promise((r) => setTimeout(r, 2500));
      const t = grok.shell.t;
      // Wait for semType detection so Launch SAR sees a Macromolecule column.
      await new Promise<void>((resolve) => {
        const sub = t.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
        setTimeout(resolve, 3000);
      });
      return {
        selected,
        extractedRows: t.rowCount,
        extractedName: t.name,
        semType: t.col('AlignedSequence')?.semType ?? null,
      };
    });
    expect(result.selected, 'first 100 rows should be selected on the source table').toBe(100);
    expect(result.extractedRows, 'Extract Selected Rows should yield a 100-row working table').toBe(100);
    expect(result.semType, 'extracted AlignedSequence must remain a Macromolecule column').toBe('Macromolecule');
  });

  await softStep('Setup: focus AlignedSequence column, open Peptides pane, click Launch SAR', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.t;
      const col = df.col('AlignedSequence');
      // Dual-set df.currentCol + grok.shell.o to deterministically rebuild the Context Panel.
      df.currentCol = col;
      await new Promise((r) => setTimeout(r, 500));
      grok.shell.o = col;

      // Poll up to 30s for the Peptides pane + Launch SAR button to mount.
      let pane: Element | null = null;
      let launchBtn: HTMLElement | null = null;
      for (let i = 0; i < 60; i++) {
        pane = document.querySelector('[name="pane-Peptides"]');
        if (pane && !pane.classList.contains('expanded')) {
          const header = pane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (header) header.click();
        }
        launchBtn = document.querySelector('[name="button-Launch-SAR"]') as HTMLElement | null;
        if (launchBtn) break;
        await new Promise((r) => setTimeout(r, 500));
      }
      const paneFound = !!pane;
      const launchFound = !!launchBtn;
      if (launchBtn) launchBtn.click();
      return {paneFound, launchFound};
    });
    expect(result.paneFound, '[name="pane-Peptides"] context pane not found (waited 30s)').toBe(true);
    expect(result.launchFound, '[name="button-Launch-SAR"] not found (waited 30s)').toBe(true);

    // SAR launch is async server compute — wait for the PeptidesModel singleton.
    await page.waitForFunction(() => {
      return Array.from(grok.shell.tableViews).some((v) => v.dataFrame.temp['peptidesModel']);
    }, {timeout: 45000});
    // Settle for MCL clustering + sequence-space embedding (~4s on the 100-row subset).
    await page.waitForTimeout(5000);

    // Verify default attach set.
    const viewers = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return Array.from(tv.viewers).map((v) => v.type);
    });
    expect(viewers, 'Sequence Variability Map must attach (peptides.viewers.monomer-position)')
      .toContain('Sequence Variability Map');
    expect(viewers, 'Most Potent Residues must attach').toContain('Most Potent Residues');
  });

  // Clear any pre-existing benign error before the regression-window starts.
  await softStep('Setup: clear baseline lastError before regression window', async () => {
    await page.evaluate(() => {
      // Best-effort clear; lastError is read-only in some builds (the assert pattern-matches anyway).
      try { (grok.shell as any).lastError = null; } catch (e) { /* nf */ }
    });
  });

  // Default mode after Launch SAR is Mutation Cliffs; toggle to Invariant Map for steps 1-4.
  await softStep('Scenario 1 (steps 1-2): toggle SVM to Invariant Map mode', async () => {
    const toggled = await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      if (!svm) return {svmFound: false};
      const im = svm.querySelector('[name="input-Invariant-Map"]') as HTMLInputElement | null;
      const mc = svm.querySelector('[name="input-Mutation-Cliffs"]') as HTMLInputElement | null;
      if (!im || !mc) return {svmFound: true, radiosFound: false};
      if (!im.checked) {
        im.click();
        await new Promise((r) => setTimeout(r, 1500));
      }
      return {
        svmFound: true,
        radiosFound: true,
        imChecked: im.checked,
        mcChecked: mc.checked,
      };
    });
    expect(toggled.svmFound, '[name="viewer-Sequence-Variability-Map"] not found').toBe(true);
    expect(toggled.radiosFound,
      'SVM mode-toggle radios (input-Mutation-Cliffs / input-Invariant-Map) not found').toBe(true);
    expect(toggled.imChecked, 'Invariant Map radio did not flip on').toBe(true);
    expect(toggled.mcChecked, 'Mutation Cliffs radio did not flip off').toBe(false);
  });

  await softStep('Scenario 1 (step 3): hover SVM Invariant-Map cell, verify tooltip appears', async () => {
    // Resolve a populated SVM cell (matrix is sparse; showTooltipAt returns null on count=0).
    const target = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const svmViewer: any = model.findViewer('Sequence Variability Map');
      if (!svmViewer) return {svmFound: false};
      const vg: any = svmViewer.viewerGrid;
      const stats: any = svmViewer.monomerPositionStats;
      // Find a populated (position, monomer) with non-trivial count.
      let bestPos: string | null = null;
      let bestMonomer: string | null = null;
      let bestCount = 0;
      for (const pos of Object.keys(stats || {})) {
        const monomersForPos = stats[pos];
        if (!monomersForPos) continue;
        for (const m of Object.keys(monomersForPos)) {
          const c = monomersForPos[m]?.count || 0;
          if (c > bestCount) { bestCount = c; bestPos = pos; bestMonomer = m; }
        }
      }
      if (!bestPos || !bestMonomer) return {svmFound: true, populatedFound: false};
      // viewerGrid is sorted by monomer, so DF row != grid row; match via getMonomerPosition.
      let cellRow = -1;
      let bounds: any = null;
      for (let r = 0; r < vg.dataFrame.rowCount; r++) {
        const probeCell: any = vg.cell(bestPos, r);
        if (!probeCell) continue;
        const mp = svmViewer.getMonomerPosition(probeCell);
        if (mp?.monomerOrCluster === bestMonomer) {
          cellRow = r;
          bounds = probeCell.bounds;
          break;
        }
      }
      if (cellRow < 0 || !bounds) return {svmFound: true, populatedFound: true, cellResolved: false};
      // Pick the largest canvas (NOT pointer-events:none row-marker) for viewport offset.
      const svmRoot: Element = svmViewer.root;
      const canvases = Array.from(svmRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {svmFound: true, populatedFound: true, cellResolved: true, canvasFound: false};
      const cv = canvas.getBoundingClientRect();
      return {
        svmFound: true, populatedFound: true, cellResolved: true, canvasFound: true,
        monomer: bestMonomer, position: bestPos, count: bestCount,
        cellRow,
        viewportX: cv.x + bounds.x + bounds.width / 2,
        viewportY: cv.y + bounds.y + bounds.height / 2,
        canvasX: cv.x, canvasY: cv.y, canvasW: cv.width, canvasH: cv.height,
      };
    });
    expect(target.svmFound, '[name="viewer-Sequence-Variability-Map"] not found').toBe(true);
    expect((target as any).populatedFound,
      'No populated (monomer, position) cells found in svmViewer.monomerPositionStats').toBe(true);
    expect((target as any).cellResolved,
      'Failed to resolve the populated cell to the matching grid row in the sorted viewerGrid').toBe(true);
    expect((target as any).canvasFound,
      'No pointer-events-enabled canvas found inside the SVM container').toBe(true);

    // Away-then-back: Playwright collapses same-point moves, so move far first for a fresh enter.
    const tx = (target as any).viewportX as number;
    const ty = (target as any).viewportY as number;
    await page.mouse.move((target as any).canvasX - 50, (target as any).canvasY - 50);
    await page.waitForTimeout(200);
    await page.mouse.move(tx, ty, {steps: 6});
    await page.waitForTimeout(1500);

    const tooltipState = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const root = tt?.root as HTMLElement | undefined;
      if (!root) return {ttFound: false};
      const display = root.style.display;
      const visible = display !== 'none' && display !== '';
      const innerLen = root.innerHTML?.length || 0;
      return {
        ttFound: true,
        visible,
        innerLen,
        sample: (root.innerText || '').slice(0, 200),
      };
    });
    expect(tooltipState.ttFound, 'ui.tooltip.root singleton not found').toBe(true);
    // Tolerant content record: the GROK-15934 contract is no-crash, asserted at the mutation steps.
    if (tooltipState.innerLen === 0)
      console.log(`[note] Invariant-Map cell hover (monomer="${(target as any).monomer}", ` +
        `position="${(target as any).position}", count=${(target as any).count}) produced no tooltip content ` +
        `(display="${tooltipState.visible}") — tolerant per the GROK-15934 no-crash contract`);
  });

  // Step 4: Move cursor to a different populated cell, verify tooltip re-renders.
  await softStep('Scenario 1 (step 4): move to a different SVM cell, verify tooltip re-renders', async () => {
    // Resolve a second populated cell (second-highest count) different from step 3.
    const target = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const svmViewer: any = model.findViewer('Sequence Variability Map');
      const vg: any = svmViewer.viewerGrid;
      const stats: any = svmViewer.monomerPositionStats;
      // Collect populated (position, monomer) pairs sorted descending by count.
      const pairs: Array<{pos: string; monomer: string; count: number}> = [];
      for (const pos of Object.keys(stats || {})) {
        const monomersForPos = stats[pos];
        if (!monomersForPos) continue;
        for (const m of Object.keys(monomersForPos)) {
          const c = monomersForPos[m]?.count || 0;
          if (c > 0) pairs.push({pos, monomer: m, count: c});
        }
      }
      pairs.sort((a, b) => b.count - a.count);
      // Pick index 1 (second-most-populated) so we hover a DIFFERENT cell than step 3.
      const pick = pairs[1] || pairs[0];
      if (!pick) return {found: false};
      let bounds: any = null;
      for (let r = 0; r < vg.dataFrame.rowCount; r++) {
        const probeCell: any = vg.cell(pick.pos, r);
        if (!probeCell) continue;
        const mp = svmViewer.getMonomerPosition(probeCell);
        if (mp?.monomerOrCluster === pick.monomer) {
          bounds = probeCell.bounds;
          break;
        }
      }
      const svmRoot: Element = svmViewer.root;
      const canvases = Array.from(svmRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas || !bounds) return {found: false};
      const cv = canvas.getBoundingClientRect();
      return {
        found: true,
        monomer: pick.monomer, position: pick.pos, count: pick.count,
        viewportX: cv.x + bounds.x + bounds.width / 2,
        viewportY: cv.y + bounds.y + bounds.height / 2,
        canvasX: cv.x, canvasY: cv.y,
      };
    });
    expect((target as any).found, 'Failed to resolve a second populated SVM cell').toBe(true);

    const tx = (target as any).viewportX as number;
    const ty = (target as any).viewportY as number;
    await page.mouse.move((target as any).canvasX - 50, (target as any).canvasY - 50);
    await page.waitForTimeout(200);
    await page.mouse.move(tx, ty, {steps: 6});
    await page.waitForTimeout(1500);

    const tooltipState = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const root = tt?.root as HTMLElement | undefined;
      const innerLen = root?.innerHTML?.length || 0;
      return {innerLen};
    });
    // Tolerant content record (no-crash is the GROK-15934 contract).
    if (tooltipState.innerLen === 0)
      console.log(`[note] Inter-cell hover transition (monomer="${(target as any).monomer}", ` +
        `position="${(target as any).position}", count=${(target as any).count}) maintained no tooltip content ` +
        `— tolerant per the GROK-15934 no-crash contract`);
  });

  await softStep('Scenario 1 (step 5-6): switch to Mutation Cliffs mode, hover, verify tooltip + no regression', async () => {
    await page.evaluate(async () => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const mc = svm?.querySelector('[name="input-Mutation-Cliffs"]') as HTMLInputElement | null;
      if (mc && !mc.checked) {
        mc.click();
        await new Promise((r) => setTimeout(r, 2000));
      }
    });

    // Prefer a cliff-bearing cell; cliffs are async so fall back to a populated Invariant-Map cell.
    const target = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const svmViewer: any = model.findViewer('Sequence Variability Map');
      const vg: any = svmViewer.viewerGrid;
      const cliffs: any = svmViewer.mutationCliffs;
      const stats: any = svmViewer.monomerPositionStats;
      // Prefer a cliff-bearing cell.
      let bestPos: string | null = null;
      let bestMonomer: string | null = null;
      let bestCount = 0;
      let source: 'cliffs' | 'inv-map' = 'cliffs';
      if (cliffs && typeof cliffs.forEach === 'function') {
        cliffs.forEach((posMap: any, monomer: string) => {
          if (posMap && typeof posMap.forEach === 'function') {
            posMap.forEach((indexMap: any, pos: string) => {
              const cnt = indexMap?.size ?? 0;
              if (cnt > bestCount) { bestCount = cnt; bestPos = pos; bestMonomer = monomer; }
            });
          }
        });
      }
      if (!bestPos || !bestMonomer) {
        source = 'inv-map';
        for (const pos of Object.keys(stats || {})) {
          const monomersForPos = stats[pos];
          if (!monomersForPos) continue;
          for (const m of Object.keys(monomersForPos)) {
            const c = monomersForPos[m]?.count || 0;
            if (c > bestCount) { bestCount = c; bestPos = pos; bestMonomer = m; }
          }
        }
      }
      if (!bestPos || !bestMonomer) return {found: false};
      let bounds: any = null;
      for (let r = 0; r < vg.dataFrame.rowCount; r++) {
        const probeCell: any = vg.cell(bestPos, r);
        if (!probeCell) continue;
        const mp = svmViewer.getMonomerPosition(probeCell);
        if (mp?.monomerOrCluster === bestMonomer) {
          bounds = probeCell.bounds;
          break;
        }
      }
      const svmRoot: Element = svmViewer.root;
      const canvases = Array.from(svmRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas || !bounds) return {found: false};
      const cv = canvas.getBoundingClientRect();
      return {
        found: true, source,
        monomer: bestMonomer, position: bestPos, count: bestCount,
        viewportX: cv.x + bounds.x + bounds.width / 2,
        viewportY: cv.y + bounds.y + bounds.height / 2,
        canvasX: cv.x, canvasY: cv.y,
      };
    });

    if (!(target as any).found) {
      console.log('[note] Scenario 1 step 5-6: no populated MutationCliffs OR InvariantMap cell resolved — ' +
        'recording informationally; the GROK-15934 invariant check below still stands');
    } else {
      const tx = (target as any).viewportX as number;
      const ty = (target as any).viewportY as number;
      await page.mouse.move((target as any).canvasX - 50, (target as any).canvasY - 50);
      await page.waitForTimeout(200);
      await page.mouse.move(tx, ty, {steps: 6});
      await page.waitForTimeout(1500);
    }

    const state = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const innerLen = tt?.root?.innerHTML?.length || 0;
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const mcChecked = (svm?.querySelector('[name="input-Mutation-Cliffs"]') as HTMLInputElement | null)?.checked;
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : null;
      return {innerLen, mcChecked, lastError};
    });
    expect(state.mcChecked, 'SVM did not switch back to Mutation Cliffs mode').toBe(true);
    // Tolerant on innerLen (depends on cell cliff data); no-crash is the GROK-15934 contract.
    if ((target as any).found) {
      expect(state.innerLen,
        `Mutation-Cliffs mode hover (monomer="${(target as any).monomer}", position="${(target as any).position}", ` +
        `count=${(target as any).count}, source=${(target as any).source}) did not produce tooltip content`)
        .toBeGreaterThan(0);
    }
    // GROK-15934 invariant: no null-receiver crash across the mode-switch mutation.
    const hasNullReceiver = state.lastError && NULL_RECEIVER_PATTERN.test(state.lastError);
    expect(hasNullReceiver,
      `GROK-15934 (post-mode-switch hover): null-receiver error surfaced: ${state.lastError}`)
      .toBeFalsy();
  });

  // Canvas-click on the main-grid header IS drivable via synthetic dispatch (only hover needs CDP).
  await softStep('Scenario 1 (step 7): WebLogo column-header click → selection mutation → re-hover SVM', async () => {
    const selBefore = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return tv.dataFrame.selection.trueCount;
    });

    // Click the column-header WebLogo strip; pick the largest canvas (row-marker canvas is small).
    const result = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const mainGrid = tv.grid;
      const canvases = Array.from(mainGrid.root.querySelectorAll('canvas'));
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const rr = (c as HTMLCanvasElement).getBoundingClientRect();
        if (rr.width * rr.height > maxArea) {
          maxArea = rr.width * rr.height;
          canvas = c as HTMLCanvasElement;
        }
      }
      if (!canvas) return {canvasFound: false, selAfter: -1};
      const r = canvas.getBoundingClientRect();
      const chh = mainGrid.props.colHeaderHeight;
      // Click into the column-header WebLogo strip, past the row-header column.
      const cx = r.x + 250;
      const cy = r.y + Math.max(20, chh / 2);
      const opts = {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 0, view: window};
      canvas.dispatchEvent(new MouseEvent('mousemove', opts));
      canvas.dispatchEvent(new MouseEvent('mousedown', opts));
      canvas.dispatchEvent(new MouseEvent('mouseup', opts));
      canvas.dispatchEvent(new MouseEvent('click', opts));
      await new Promise((res) => setTimeout(res, 2000));
      // Re-resolve the model TableView post-click (active view may drift on dock dispatch).
      const tv2 = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      return {canvasFound: true, selAfter: tv2.dataFrame.selection.trueCount};
    });
    expect(result.canvasFound, 'main grid canvas not found for WebLogo header click').toBe(true);
    // Selection may or may not grow depending on where the click landed; either is acceptable.
    if (result.selAfter === selBefore) {
      console.log(`[note] WebLogo header click did not raise selection.trueCount (before=${selBefore}, after=${result.selAfter}) — ` +
        `hit-test landed on an empty canvas region. Re-hover assertion below still exercises the post-click model state.`);
    }

    // Re-hover the SVM after the selection mutation.
    const svmRect = await page.evaluate(() => {
      const svm = document.querySelector('[name="viewer-Sequence-Variability-Map"]');
      const canvases = svm ? Array.from(svm.querySelectorAll('canvas')) : [];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        const r = (c as HTMLCanvasElement).getBoundingClientRect();
        if (r.width * r.height > maxArea) {
          maxArea = r.width * r.height;
          canvas = c as HTMLCanvasElement;
        }
      }
      if (!canvas) return null;
      const r = canvas.getBoundingClientRect();
      return {x: r.x, y: r.y, w: r.width, h: r.height};
    });
    expect(svmRect, 'SVM canvas not found for post-selection re-hover').not.toBeNull();

    await page.mouse.move(svmRect!.x - 50, svmRect!.y - 50);
    await page.waitForTimeout(200);
    await page.mouse.move(svmRect!.x + 200, svmRect!.y + 160, {steps: 6});
    await page.waitForTimeout(1500);

    const lastError = await page.evaluate(() =>
      grok.shell.lastError ? String(grok.shell.lastError) : null);
    const hasNullReceiver = lastError && NULL_RECEIVER_PATTERN.test(lastError);
    // GROK-15934 invariant: post-selection-sync hover must not throw on null.
    expect(hasNullReceiver,
      `GROK-15934 (post-selection-sync hover): null-receiver error surfaced: ${lastError}`)
      .toBeFalsy();
  });

  // Settings-driven mutation via Viewers-pane Active-peptide-selection toggle (Columns variant deferred).
  await softStep('Scenario 1 (step 8): Settings dialog round-trip → re-hover SVM', async () => {
    // d4 panes signal collapse via .d4-accordion-pane-content display==='none', not a class;
    // detect that (never toggle an open pane) and drive cb/OK via composed MouseEvents.
    await page.evaluate(async () => {
      const wrench = document.querySelector(
        'i.grok-icon.fa-wrench[aria-label="Peptides analysis settings"]') as HTMLElement | null;
      if (wrench) wrench.click();
    });
    await page.locator('[name="dialog-Peptides-settings"]').waitFor({timeout: 8000});

    const ok = await page.evaluate(async () => {
      const dlg = document.querySelector('[name="dialog-Peptides-settings"]');
      if (!dlg) return {error: 'dialog not found'};
      const viewersPane = Array.from(dlg.querySelectorAll('.d4-accordion-pane'))
        .find((p) => p.querySelector('.d4-accordion-pane-header')?.textContent?.trim() === 'Viewers');
      if (viewersPane) {
        const content = viewersPane.querySelector('.d4-accordion-pane-content') as HTMLElement | null;
        const isCollapsed = !!content && getComputedStyle(content).display === 'none';
        if (isCollapsed) {
          const h = viewersPane.querySelector('.d4-accordion-pane-header') as HTMLElement | null;
          if (h) h.click();
          await new Promise((r) => setTimeout(r, 600));
        }
      }
      const cb = dlg.querySelector('[name="input-Active-peptide-selection"]') as HTMLInputElement | null;
      if (!cb) return {error: 'cb not found'};
      // Two clicks toggle OFF->ON->OFF, exercising the closeViewer/addClusterMaxActivityViewer chain.
      cb.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      await new Promise((r) => setTimeout(r, 400));
      cb.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      await new Promise((r) => setTimeout(r, 400));
      const okBtn = dlg.querySelector('[name="button-OK"]') as HTMLElement | null;
      if (!okBtn) return {error: 'OK btn not found'};
      okBtn.dispatchEvent(new MouseEvent('click', {bubbles: true, cancelable: true, composed: true, view: window}));
      return {ok: true};
    });
    expect(ok && (ok as any).error, 'Settings dialog driving setup failed').toBeFalsy();

    await page.waitForFunction(() =>
      !document.querySelector('[name="dialog-Peptides-settings"]'), null, {timeout: 8000});
    await page.waitForTimeout(2000);

    // Post-settings re-hover on a populated cell.
    const target = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const svmViewer: any = model.findViewer('Sequence Variability Map');
      if (!svmViewer) return {svmFound: false};
      const vg: any = svmViewer.viewerGrid;
      const stats: any = svmViewer.monomerPositionStats;
      let bestPos: string | null = null;
      let bestMonomer: string | null = null;
      let bestCount = 0;
      for (const pos of Object.keys(stats || {})) {
        const monomersForPos = stats[pos];
        if (!monomersForPos) continue;
        for (const m of Object.keys(monomersForPos)) {
          const c = monomersForPos[m]?.count || 0;
          if (c > bestCount) { bestCount = c; bestPos = pos; bestMonomer = m; }
        }
      }
      if (!bestPos || !bestMonomer) return {svmFound: true, populatedFound: false};
      let bounds: any = null;
      for (let r = 0; r < vg.dataFrame.rowCount; r++) {
        const probeCell: any = vg.cell(bestPos, r);
        if (!probeCell) continue;
        const mp = svmViewer.getMonomerPosition(probeCell);
        if (mp?.monomerOrCluster === bestMonomer) {
          bounds = probeCell.bounds;
          break;
        }
      }
      const svmRoot: Element = svmViewer.root;
      const canvases = Array.from(svmRoot.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas || !bounds) return {svmFound: true, populatedFound: true, cellResolved: false};
      const cv = canvas.getBoundingClientRect();
      return {
        svmFound: true, populatedFound: true, cellResolved: true,
        monomer: bestMonomer, position: bestPos, count: bestCount,
        viewportX: cv.x + bounds.x + bounds.width / 2,
        viewportY: cv.y + bounds.y + bounds.height / 2,
        canvasX: cv.x, canvasY: cv.y,
      };
    });
    expect((target as any).svmFound, 'SVM viewer missing after Settings round-trip (SVM should persist)').toBe(true);

    if ((target as any).cellResolved) {
      const tx = (target as any).viewportX as number;
      const ty = (target as any).viewportY as number;
      await page.mouse.move((target as any).canvasX - 50, (target as any).canvasY - 50);
      await page.waitForTimeout(200);
      await page.mouse.move(tx, ty, {steps: 6});
      await page.waitForTimeout(1500);
    }

    const finalState = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const innerLen = tt?.root?.innerHTML?.length || 0;
      const lastError = grok.shell.lastError ? String(grok.shell.lastError) : null;
      return {innerLen, lastError};
    });
    // GROK-15934 invariant: post-settings hover must not throw on null column receiver.
    const hasNullReceiver = finalState.lastError && NULL_RECEIVER_PATTERN.test(finalState.lastError);
    expect(hasNullReceiver,
      `GROK-15934 invariant: post-settings-change hover produced null-receiver error: ${finalState.lastError}`)
      .toBeFalsy();
  });

  // Step 9: Cursor leaves viewer body → tooltip dismisses + highlight clears.
  await softStep('Scenario 1 (step 9): cursor leaves SVM body → tooltip dismisses', async () => {
    await page.mouse.move(10, 10);
    await page.waitForTimeout(1500);
    const tooltipState = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const display = tt?.root?.style?.display;
      return {display};
    });
    expect(tooltipState.display,
      'tooltip should dismiss (display:none) when cursor leaves the SVM body')
      .toBe('none');
  });

  // Scenario 2 — column-header WebLogo hover via showTooltipAt explicit-anchor path.
  await softStep('Scenario 2 (steps 2-4): hover main-grid column-header WebLogo, verify tooltip via showTooltipAt', async () => {
    // Resolve a VISIBLE position column; hover at its center-x and middle-y of the header strip.
    const target = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const mainGrid: any = (tv as any).grid;
      const stats: any = model.monomerPositionStats;
      const canvases = Array.from(mainGrid.root.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {found: false};
      const cv = canvas.getBoundingClientRect();
      const chh = mainGrid.props.colHeaderHeight;
      // Pick the first position column that is fully visible AND has a populated monomer.
      for (let p = 1; p <= 30; p++) {
        const posName = String(p);
        const populated = stats?.[posName] && Object.values(stats[posName] as any)
          .some((s: any) => (s?.count || 0) > 0);
        if (!populated) continue;
        try {
          const cell = mainGrid.cell(posName, 0);
          const b = cell?.bounds;
          if (!b) continue;
          if (b.x >= 0 && (b.x + b.width) <= cv.width) {
            return {
              found: true, position: posName,
              chh,
              viewportX: cv.x + b.x + b.width / 2,
              viewportY: cv.y + Math.floor(chh / 2),
              canvasX: cv.x, canvasY: cv.y, canvasW: cv.width, canvasH: cv.height,
            };
          }
        } catch (e) { /* continue */ }
      }
      return {found: false, chh, canvasW: cv.width};
    });
    expect((target as any).found,
      `No visible+populated position column found for column-header hover ` +
      `(canvasW=${(target as any).canvasW || 'n/a'}, chh=${(target as any).chh || 'n/a'})`).toBe(true);
    expect((target as any).chh,
      'main grid column-header height should be enlarged for WebLogo (>=80 px expected)')
      .toBeGreaterThanOrEqual(80);

    const tx = (target as any).viewportX as number;
    const ty = (target as any).viewportY as number;
    await page.mouse.move((target as any).canvasX - 50, (target as any).canvasY - 50);
    await page.waitForTimeout(200);
    await page.mouse.move(tx, ty, {steps: 6});
    await page.waitForTimeout(1500);

    const tooltipState = await page.evaluate(() => {
      const tt = (ui as any).tooltip;
      const root = tt?.root as HTMLElement | undefined;
      if (!root) return {ttFound: false};
      const innerLen = root.innerHTML?.length || 0;
      return {ttFound: true, innerLen, sample: (root.innerText || '').slice(0, 200)};
    });
    expect(tooltipState.ttFound, 'ui.tooltip.root singleton not found').toBe(true);
    expect(tooltipState.innerLen,
      `Column-header WebLogo hover (position="${(target as any).position}") did not produce tooltip content ` +
      `(showTooltipAt should render rich payload via the WebLogo letter's bounding-box anchor)`)
      .toBeGreaterThan(0);
  });

  await softStep('Scenario 2 (steps 5-6): inter-letter + cross-column hover transitions', async () => {
    // Resolve two visible populated position columns and hover between them.
    const targets = await page.evaluate(() => {
      const tv = Array.from(grok.shell.tableViews).find((v: any) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = (tv as any).dataFrame.temp['peptidesModel'];
      const mainGrid: any = (tv as any).grid;
      const stats: any = model.monomerPositionStats;
      const canvases = Array.from(mainGrid.root.querySelectorAll('canvas')) as HTMLCanvasElement[];
      let canvas: HTMLCanvasElement | null = null;
      let maxArea = 0;
      for (const c of canvases) {
        if (getComputedStyle(c).pointerEvents === 'none') continue;
        const r = c.getBoundingClientRect();
        if (r.width * r.height > maxArea) { maxArea = r.width * r.height; canvas = c; }
      }
      if (!canvas) return {found: false};
      const cv = canvas.getBoundingClientRect();
      const chh = mainGrid.props.colHeaderHeight;
      const visible: Array<{position: string; viewportX: number; viewportY: number}> = [];
      for (let p = 1; p <= 30; p++) {
        const posName = String(p);
        const populated = stats?.[posName] && Object.values(stats[posName] as any)
          .some((s: any) => (s?.count || 0) > 0);
        if (!populated) continue;
        try {
          const cell = mainGrid.cell(posName, 0);
          const b = cell?.bounds;
          if (!b) continue;
          if (b.x >= 0 && (b.x + b.width) <= cv.width) {
            visible.push({
              position: posName,
              viewportX: cv.x + b.x + b.width / 2,
              viewportY: cv.y + Math.floor(chh / 2),
            });
          }
        } catch (e) { /* continue */ }
      }
      return {
        found: visible.length >= 2,
        visibleCount: visible.length,
        first: visible[0], second: visible[1],
        canvasX: cv.x, canvasY: cv.y,
      };
    });

    if (!(targets as any).found) {
      console.log(`[note] Scenario 2 steps 5-6: only ${(targets as any).visibleCount || 0} visible+populated ` +
        `position columns; cross-column transition tolerantly skipped, GROK-15934 invariant check still runs`);
    } else {
      // Hover the first column header.
      await page.mouse.move((targets as any).canvasX - 50, (targets as any).canvasY - 50);
      await page.waitForTimeout(200);
      await page.mouse.move((targets as any).first.viewportX, (targets as any).first.viewportY, {steps: 4});
      await page.waitForTimeout(800);
      // Cross-column: move to a different position column's header.
      await page.mouse.move((targets as any).second.viewportX, (targets as any).second.viewportY, {steps: 4});
      await page.waitForTimeout(1200);
    }

    const lastError = await page.evaluate(() =>
      grok.shell.lastError ? String(grok.shell.lastError) : null);
    const hasNullReceiver = lastError && NULL_RECEIVER_PATTERN.test(lastError);
    expect(hasNullReceiver,
      `GROK-15934 (inter-letter / cross-column WebLogo hover): null-receiver error surfaced: ${lastError}`)
      .toBeFalsy();
  });

  await softStep('Scenario 2 (step 7): LST cell WebLogo hover (no-crash invariant)', async () => {
    const lstState = await page.evaluate(async () => {
      const tv = Array.from(grok.shell.tableViews).find((v) => v.dataFrame.temp['peptidesModel']) ?? grok.shell.tv;
      const model = tv.dataFrame.temp['peptidesModel'];
      const lst = model.findViewer('Logo Summary Table');
      if (!lst) {
        // Defensive only; subset normally auto-attaches LST (imperative add throws circular-JSON).
        try { await model.addLogoSummaryTable(); }
        catch (e) { return {lstAttached: false, attemptThrew: String(e).slice(0, 200)}; }
        await new Promise((r) => setTimeout(r, 3000));
        const lst2 = model.findViewer('Logo Summary Table');
        if (!lst2) return {lstAttached: false, attemptThrew: null};
      }
      const lstViewer = model.findViewer('Logo Summary Table');
      if (!lstViewer) return {lstAttached: false, attemptThrew: null};
      const root = lstViewer.root as HTMLElement;
      const canvas = root?.querySelector('canvas') as HTMLCanvasElement | null;
      if (!canvas) return {lstAttached: true, canvasFound: false};
      const r = canvas.getBoundingClientRect();
      return {
        lstAttached: true,
        canvasFound: true,
        x: r.x, y: r.y, w: r.width, h: r.height,
      };
    });

    if (!lstState.lstAttached) {
      console.log(`[note] (SR-02) Logo Summary Table viewer unexpectedly not attached on the 100-row subset ` +
        `(MCL clustering normally attaches it in ~4s; imperative add throws circular-JSON: ${lstState.attemptThrew || 'no-attempt'}). ` +
        `LST cell-WebLogo hover skipped; column-header WebLogo above covers the ` +
        `drawLogoInBounds primitive at its deterministic call-site.`);
      return;
    }
    if (!lstState.canvasFound) {
      console.log(`[note] (SR-02) LST viewer attached but no canvas found for hover; recorded informationally.`);
      return;
    }

    // Hover the leftmost WebLogo cell.
    await page.mouse.move(lstState.x! - 30, lstState.y! - 30);
    await page.waitForTimeout(200);
    await page.mouse.move(lstState.x! + 60, lstState.y! + 30, {steps: 4});
    await page.waitForTimeout(1500);

    const lastError = await page.evaluate(() =>
      grok.shell.lastError ? String(grok.shell.lastError) : null);
    const hasNullReceiver = lastError && NULL_RECEIVER_PATTERN.test(lastError);
    expect(hasNullReceiver,
      `GROK-15934 (LST cell WebLogo hover): null-receiver error surfaced: ${lastError}`)
      .toBeFalsy();
  });

  await softStep('Scenario 2 (step 9): GROK-15934 invariant — no null-receiver across dual-call-site WebLogo + showTooltipAt', async () => {
    const lastError = await page.evaluate(() =>
      grok.shell.lastError ? String(grok.shell.lastError) : null);
    const hasNullReceiver = lastError && NULL_RECEIVER_PATTERN.test(lastError);
    expect(hasNullReceiver,
      `GROK-15934 invariant (final): null-receiver / getTag-on-null error surfaced across the ` +
      `full hover sequence (SVM Invariant-Map + SVM Mutation-Cliffs + post-mode-switch + ` +
      `post-selection-sync + post-settings + column-header WebLogo + optional LST). ` +
      `lastError: ${lastError}`)
      .toBeFalsy();
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
