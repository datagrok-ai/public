// CLAUDE-33: Molstar onViewRemoved must not throw "reading 'children'" when an unrelated view closes.
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbPath = 'System:AppData/BiostructureViewer/samples/1bdq.pdb';

// Strictly matches the CLAUDE-33 signature, avoiding rcsb-molstar 'props' noise.
const CLAUDE_33_SIGNATURE = /Cannot read properties of undefined \(reading '?children'?\)/i;

function matchesClaude33(text: string): boolean {
  if (!text) return false;
  return CLAUDE_33_SIGNATURE.test(text);
}

test('BiostructureViewer — CLAUDE-33 Molstar onViewRemoved unrelated-view-close safety regression guard', async ({page}) => {
  test.setTimeout(600_000);
  stepErrors.length = 0;

  // Capture both pageerror and console errors; the CLAUDE-33 signature can surface on either.
  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (err) => { pageErrors.push(err.message); });
  page.on('console', (msg) => {
    if (msg.type() === 'error') consoleErrors.push(msg.text());
  });

  await loginToDatagrok(page);

  // Windows mode (simpleMode=false) so view-tabs render with clickable geometry.
  await page.evaluate(() => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = false;
  });
  await page.locator('[name="Browse"]').waitFor({timeout: 30_000});

  try {
    // SETUP — Stage a host table view + mount Biostructure viewer to wire onViewRemoved.
    let setupDiag: any = null;

    await softStep('Setup — Stage host table view + tv.addViewer(Biostructure) + setOptions({pdb}) to wire onViewRemoved subscription', async () => {
      setupDiag = await page.evaluate(async (path) => {
        const content = await grok.dapi.files.readAsText(path);

        const hostDf = DG.DataFrame.fromColumns([
          DG.Column.fromStrings('id', ['host-1', 'host-2']),
        ]);
        hostDf.name = 'host-claude-33';
        const hostTv = grok.shell.addTableView(hostDf);
        await new Promise((r) => setTimeout(r, 1500));
        const hostViewName = hostTv && hostTv.name ? hostTv.name : 'host-claude-33';

        let bioViewerMounted = false;
        let bioSetOptsErr: string | null = null;
        try {
          const bioViewer = hostTv.addViewer('Biostructure');
          await new Promise((r) => setTimeout(r, 1500));
          try { bioViewer.setOptions({pdb: content}); } catch (e: any) {
            bioSetOptsErr = String(e && e.message ? e.message : e);
          }
          bioViewerMounted = true;
        } catch (e: any) {
          bioSetOptsErr = String(e && e.message ? e.message : e);
        }

        // Wait up to 15s for .msp-plugin or the viewer container to mount.
        let mspPluginMounted = false;
        let viewerContainerMounted = false;
        for (let i = 0; i < 75; i++) {
          if (document.querySelector('.msp-plugin')) mspPluginMounted = true;
          if (document.querySelector('[name="viewer-Biostructure"]')) viewerContainerMounted = true;
          if (mspPluginMounted && viewerContainerMounted) break;
          await new Promise((r) => setTimeout(r, 200));
        }
        await new Promise((r) => setTimeout(r, 2000));

        const viewerTypes = (hostTv && hostTv.viewers)
          ? Array.from(hostTv.viewers).map((v: any) => v.type)
          : [];

        return {
          hostViewName,
          hostViewListLength: (grok.shell.tableViews && grok.shell.tableViews.length) || 0,
          bioViewerMounted,
          bioSetOptsErr,
          mspPluginMounted,
          viewerContainerMounted,
          mspViewportPresent: !!document.querySelector('.msp-viewport'),
          contentLen: content.length,
          viewerTypes,
        };
      }, samplePdbPath);

      await page.locator('[name="viewer-Biostructure"]').first().waitFor({timeout: 30_000});
      await page.locator(`[name="view-handle: ${setupDiag.hostViewName}"]`).first().waitFor({timeout: 30_000});

      expect(setupDiag.bioViewerMounted).toBe(true);
      expect(setupDiag.contentLen).toBeGreaterThan(1000);
      expect(setupDiag.hostViewName).toBe('host-claude-33');
      expect(setupDiag.viewerTypes).toEqual(expect.arrayContaining(['Grid', 'Biostructure']));
    });

    // SCENARIO 1 — Closing an unrelated view (3 cycles) must not throw the 'children' signature.
    const claude33Hits: Array<{cycle: number, channel: string, msg: string}> = [];
    const cycleSummaries: Array<any> = [];

    for (let cycle = 1; cycle <= 3; cycle++) {
      await softStep(`Scenario 1 — Open unrelated view, click view-tab Close (cycle ${cycle} of 3); onViewRemoved handler MUST NOT throw 'children' signature`, async () => {
        // Reset capture buffers to isolate this close action from setup noise.
        pageErrors.length = 0;
        consoleErrors.length = 0;

        const unrelatedName = await page.evaluate((cy) => {
          const name = `unrelated-claude-33-${cy}-${Date.now()}`;
          const v = DG.View.create();
          v.name = name;
          grok.shell.addView(v);
          return name;
        }, cycle);

        const tabLocator = page.locator(`[name="view-handle: ${unrelatedName}"]`).first();
        await tabLocator.waitFor({timeout: 15_000});

        // Settle so the platform commits the addView before dispatching close.
        await page.waitForTimeout(800);

        // Click the Close icon on the unrelated view's tab.
        const closeIcon = tabLocator.locator('[name="Close"]').first();
        let clickErr: string | null = null;
        try {
          await closeIcon.waitFor({timeout: 5_000});
          await closeIcon.click({timeout: 10_000});
        } catch (e: any) {
          clickErr = String(e && e.message ? e.message : e);
        }

        await page.waitForTimeout(1500);

        // Fallback: if the UI click did not remove the tab, close it programmatically
        // so onViewRemoved fires and the assertion below is not vacuous.
        let usedProgrammaticFallback = false;
        const tabStillPresent = (await page.locator(`[name="view-handle: ${unrelatedName}"]`).count()) > 0;
        if (tabStillPresent) {
          await page.evaluate((name) => {
            const v = grok.shell.views && Array.from(grok.shell.views).find((vw: any) => vw && vw.name === name);
            if (v && typeof v.close === 'function') v.close();
          }, unrelatedName);
          usedProgrammaticFallback = true;
          await page.waitForTimeout(1500);
        }

        const tabGone = (await page.locator(`[name="view-handle: ${unrelatedName}"]`).count()) === 0;

        const pageErrSig = pageErrors.filter(matchesClaude33);
        const consoleErrSig = consoleErrors.filter(matchesClaude33);

        for (const m of pageErrSig) claude33Hits.push({cycle, channel: 'pageerror', msg: m});
        for (const m of consoleErrSig) claude33Hits.push({cycle, channel: 'console', msg: m});

        cycleSummaries.push({
          cycle,
          unrelatedName,
          clickErr,
          usedProgrammaticFallback,
          tabGone,
          pageErrCount: pageErrors.length,
          consoleErrCount: consoleErrors.length,
          pageErrSigCount: pageErrSig.length,
          consoleErrSigCount: consoleErrSig.length,
          consoleErrSample: consoleErrors.slice(0, 5).map((s) => s.slice(0, 200)),
        });

        // The CLAUDE-33 signature must not surface on either channel after the close.
        expect(
          pageErrSig,
          'CLAUDE-33 Molstar onViewRemoved unrelated-view-close crash ' +
          `regressed (cycle ${cycle}, pageerror channel): ` +
          `${JSON.stringify(pageErrSig)}. The onViewRemoved handler ` +
          'dereferenced `evtView.root.children[0].children[0]` without ' +
          'null-guard before the view-ID equality check. See ' +
          'bug-library/biostructureviewer.yaml#CLAUDE-33; defect site ' +
          'public/packages/BiostructureViewer/src/viewers/molstar-viewer/utils.ts#L155-L157.',
        ).toEqual([]);
        expect(
          consoleErrSig,
          'CLAUDE-33 Molstar onViewRemoved unrelated-view-close crash ' +
          `regressed (cycle ${cycle}, console channel): ` +
          `${JSON.stringify(consoleErrSig)}. Same defect as above (see ` +
          'bug-library/biostructureviewer.yaml#CLAUDE-33); the handler ' +
          'throw surfaced as a console.error instead of a re-thrown ' +
          'page-level exception.',
        ).toEqual([]);
        expect(
          tabGone,
          `Scenario 1 cycle ${cycle}: unrelated view ${unrelatedName} was ` +
          'NOT removed from the DOM after the close trigger (UI click + ' +
          'programmatic fallback). The bug-invariant assertion above is ' +
          'meaningless if the close trigger did not actually fire ' +
          'onViewRemoved with the unrelated view as evtView. cycleSummary: ' +
          `${JSON.stringify(cycleSummaries[cycleSummaries.length - 1])}.`,
        ).toBe(true);
      });
    }

    // SCENARIO 1 step 6 — Cross-cycle invariant: zero CLAUDE-33 hits across the three cycles.
    await softStep('Scenario 1 step 6 — Cross-cycle invariant: zero CLAUDE-33 signature hits across three open/close cycles', async () => {
      // eslint-disable-next-line no-console
      console.log(`[CLAUDE-33 cycle summaries] ${JSON.stringify(cycleSummaries)}`);
      expect(
        claude33Hits,
        'CLAUDE-33 signature surfaced across the three-cycle invariant ' +
        `check: ${JSON.stringify(claude33Hits)}. The handler must no-op ` +
        'for non-matching views across a range of `evtView.root` shapes.',
      ).toEqual([]);
      // At least two of three cycles must have closed their tabs (robust to one timing hiccup).
      const tabsClosedCount = cycleSummaries.filter((s) => s.tabGone === true).length;
      expect(
        tabsClosedCount,
        'Fewer than two of three unrelated view-tabs were observed as ' +
        'removed after the close trigger — the bug-invariant assertion ' +
        'would be vacuously satisfied. cycleSummaries: ' +
        `${JSON.stringify(cycleSummaries)}.`,
      ).toBeGreaterThanOrEqual(2);
    });

    // SCENARIO 2 — Closing the Molstar host view itself must still tear down cleanly (inverse guard).
    let scenario2Diag: any = null;

    await softStep('Scenario 2 step 3 — Close the Molstar host view; teardown must be clean (no children-signature error)', async () => {
      // Reset capture buffers before the host-view close.
      pageErrors.length = 0;
      consoleErrors.length = 0;

      const hostViewName = setupDiag.hostViewName;

      // Ensure a live host with the Biostructure viewer exists before closing.
      const hostTabPresent = (await page.locator(`[name="view-handle: ${hostViewName}"]`).count()) > 0;
      if (!hostTabPresent) {
        await page.evaluate(async (path) => {
          const content = await grok.dapi.files.readAsText(path);
          const hostDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('id', ['rh-1'])]);
          hostDf.name = 'host-claude-33-rehydrate';
          const hostTv = grok.shell.addTableView(hostDf);
          await new Promise((r) => setTimeout(r, 1500));
          try {
            const bioViewer = hostTv.addViewer('Biostructure');
            await new Promise((r) => setTimeout(r, 1000));
            try { bioViewer.setOptions({pdb: content}); } catch (_) { /* best-effort */ }
          } catch (_) { /* best-effort */ }
          for (let i = 0; i < 75; i++) {
            if (document.querySelector('.msp-plugin')) break;
            await new Promise((r) => setTimeout(r, 200));
          }
          await new Promise((r) => setTimeout(r, 1500));
        }, samplePdbPath);
        await page.locator(`[name="view-handle: host-claude-33-rehydrate"]`).first().waitFor({timeout: 30_000});
      }

      // Pick whichever host-view-tab is currently present.
      const candidateNames = ['host-claude-33', 'host-claude-33-rehydrate'];
      let liveHostName: string | null = null;
      for (const n of candidateNames) {
        if ((await page.locator(`[name="view-handle: ${n}"]`).count()) > 0) { liveHostName = n; break; }
      }

      if (!liveHostName) {
        throw new Error(
          'Scenario 2 setup: no Molstar host view-tab is present in the ' +
          'DOM (neither host-claude-33 nor host-claude-33-rehydrate). The ' +
          'inverse-regression invariant cannot be exercised without a ' +
          'live host to close.',
        );
      }

      const liveHostLocator = page.locator(`[name="view-handle: ${liveHostName}"]`).first();
      const liveHostClose = liveHostLocator.locator('[name="Close"]').first();
      let clickErr: string | null = null;
      try {
        await liveHostClose.waitFor({timeout: 5_000});
        await liveHostClose.click({timeout: 10_000});
      } catch (e: any) {
        clickErr = String(e && e.message ? e.message : e);
      }

      await page.waitForTimeout(2000);

      // Fallback: if the UI click did not remove the host tab, close it programmatically.
      let usedProgrammaticFallback = false;
      const hostStillPresent = (await page.locator(`[name="view-handle: ${liveHostName}"]`).count()) > 0;
      if (hostStillPresent) {
        await page.evaluate((name) => {
          const v = grok.shell.tableViews
            && Array.from(grok.shell.tableViews).find((vw: any) => vw && vw.name === name);
          if (v && typeof v.close === 'function') v.close();
        }, liveHostName);
        usedProgrammaticFallback = true;
        await page.waitForTimeout(1500);
      }

      const hostGone = (await page.locator(`[name="view-handle: ${liveHostName}"]`).count()) === 0;
      const mspViewportGone = (await page.locator('.msp-viewport').count()) === 0;

      const pageErrSig = pageErrors.filter(matchesClaude33);
      const consoleErrSig = consoleErrors.filter(matchesClaude33);

      scenario2Diag = {
        liveHostName,
        clickErr,
        usedProgrammaticFallback,
        hostGone,
        mspViewportGone,
        pageErrSig,
        consoleErrSig,
        pageErrCount: pageErrors.length,
        consoleErrCount: consoleErrors.length,
        consoleErrSample: consoleErrors.slice(0, 5).map((s) => s.slice(0, 200)),
      };

      // Inverse #1: the CLAUDE-33 signature must not surface on the host close either.
      expect(
        pageErrSig,
        'CLAUDE-33 inverse-regression: closing the Molstar host view ' +
        'itself raised the children-signature TypeError ' +
        `(pageerror channel): ${JSON.stringify(pageErrSig)}.`,
      ).toEqual([]);
      expect(
        consoleErrSig,
        'CLAUDE-33 inverse-regression: closing the Molstar host view ' +
        'itself raised the children-signature TypeError ' +
        `(console channel): ${JSON.stringify(consoleErrSig)}.`,
      ).toEqual([]);

      // Inverse #2: the host close must actually tear down (hostGone true).
      expect(
        scenario2Diag.hostGone,
        'CLAUDE-33 inverse-regression: clicking the Close icon on the ' +
        `Molstar host view-tab (${liveHostName}) AND the programmatic ` +
        'fallback both failed to remove the tab from the DOM. This is ' +
        'the regressed-fix-over-applied shape: the handler no-ops for ' +
        'the host view\'s own close event too, breaking the legitimate ' +
        'teardown branch. See scenario .md Scenario 2 inverse-regression ' +
        'signature. scenario2Diag: ' + JSON.stringify(scenario2Diag),
      ).toBe(true);
    });

    // SCENARIO 2 step 4 — Joint invariant cross-check (CLAUDE-33 onViewRemoved safety).
    await softStep('Scenario 2 step 4 — Joint invariant cross-check (CLAUDE-33 onViewRemoved safety)', async () => {
      const summary = {
        scenario1: {
          cycleSummaries: cycleSummaries.map((s) => ({
            cycle: s.cycle,
            tabGone: s.tabGone,
            usedProgrammaticFallback: s.usedProgrammaticFallback,
            pageErrSigCount: s.pageErrSigCount,
            consoleErrSigCount: s.consoleErrSigCount,
          })),
          totalClaude33Hits: claude33Hits.length,
        },
        scenario2: {
          liveHostName: scenario2Diag ? scenario2Diag.liveHostName : null,
          hostGone: scenario2Diag ? scenario2Diag.hostGone : null,
          mspViewportGone: scenario2Diag ? scenario2Diag.mspViewportGone : null,
          usedProgrammaticFallback: scenario2Diag ? scenario2Diag.usedProgrammaticFallback : null,
          pageErrSigCount: scenario2Diag ? scenario2Diag.pageErrSig.length : null,
          consoleErrSigCount: scenario2Diag ? scenario2Diag.consoleErrSig.length : null,
        },
        jointInvariantHolds: !!(
          claude33Hits.length === 0 &&
          scenario2Diag &&
          scenario2Diag.hostGone === true &&
          scenario2Diag.pageErrSig.length === 0 &&
          scenario2Diag.consoleErrSig.length === 0
        ),
      };
      // eslint-disable-next-line no-console
      console.log(`[CLAUDE-33 joint-invariant summary] ${JSON.stringify(summary)}`);
      expect(summary.jointInvariantHolds).toBe(true);
    });
  } finally {
    // Cleanup — close menus/dialogs and reset shell state.
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
    } catch (_) { /* best-effort */ }
  }

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
