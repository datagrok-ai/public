/* ---
realizes: []
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {openChemMenuItem} from '../helpers/chem';

declare const grok: any;
declare const DG: any;

// Selector recon-notes (class-2: live-MCP-observed, not yet in grok-browser chem.md):
//   [name="button-Save"] — table-view toolbar Save button; opens the Save project dialog.
//   [name="dialog-Save-project"] — Save project dialog (name text input + OK/CANCEL).
//   [name="input-host-Data-sync"] — Data sync toggle host inside the save dialog; a
//     ui-input-bool-switch that is style="display:none" in the default table-view save
//     flow (no input[type=checkbox]; state on the .ui-input-switch decorator). The host is
//     hidden for any table without Tags.CreationScript — see project_entity_move.dart:513.
//   All reached via [name="button-Save"] on an open smiles.csv table view; observed live
//   2026-08-19 via chrome-devtools MCP (evaluate_script). Not in chem.md.
test.use(specTestOptions);

// spgi-100, not smiles.csv: the Scaffold Tree does not generate for a 1000-row table — the
// Generate action is disabled above a small-dataset bound, so the tree stays empty and the
// round-trip has no scaffold state to restore. 100 rows is what the green neighbour
// chem-filters-spec.ts drives a scaffold tree on.
const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const DATASET_ROWS = 100;
// The column the substructure search must run on. spgi-100 exposes seven Molecule columns;
// the R-group ones hold fragment SMILES and searching them is not what this scenario is about.
const MOL_COL = 'Structure';
// Structure, Core, R1, R2, R3, R100, R101 — the seven the picker branch depends on. Pinned
// rather than bounded below: a fixture that degraded to one Molecule column would take the
// no-picker branch and silently drop the preselection assertions while a `> 0` guard stayed green.
const MOL_COL_COUNT = 7;
// Scaffold Tree Generate is disabled at or above this many DISTINCT molecules — the icon gains
// `inactive` and a capture-phase blocker swallows the click (widgets/scaffold-tree.ts#L2340, #L86).
const MAX_DISTINCT_MOLECULES = 500;

const roundtripProject = 'ChemStateRoundtripAuto' + Date.now();
const scaffoldOnlyProject = 'ChemScaffoldOnlyAuto' + Date.now();

// A bare [name="viewer-Scaffold-Tree"] also matches a Scaffold Tree living in another, possibly
// hidden, view, and reading that one reports a foreign dataframe (chem.md:378-395). Both the
// in-page reads and the Playwright locators go through this resolver, which scopes by the viewer's
// own dataFrame as that recipe prescribes and stamps the element so a locator aims at exactly one.
const SCAFFOLD_TARGET = '[name="viewer-Scaffold-Tree"][data-spec-scaffold="target"]';

function installScaffoldResolver() {
  (window as any).__scaffoldEl = () => {
    const view = grok.shell.tv;
    if (!view || !view.dataFrame) return null;
    const v = Array.from(view.viewers)
      .find((x: any) => /Scaffold Tree/i.test(x.type || '') && x.dataFrame === view.dataFrame) as any;
    if (!v) return null;
    const root = v.root as HTMLElement;
    const el = (root.matches('[name="viewer-Scaffold-Tree"]') ? root :
      root.closest('[name="viewer-Scaffold-Tree"]') ??
      root.querySelector('[name="viewer-Scaffold-Tree"]')) as HTMLElement | null;
    return el && el.isConnected ? el : null;
  };
  (window as any).__stampScaffold = () => {
    document.querySelectorAll('[data-spec-scaffold="target"]')
      .forEach((e) => e.removeAttribute('data-spec-scaffold'));
    const el = (window as any).__scaffoldEl() as HTMLElement | null;
    if (el) el.setAttribute('data-spec-scaffold', 'target');
    return {found: !!el, domMatches: document.querySelectorAll('[name="viewer-Scaffold-Tree"]').length};
  };
}

async function aimAtScaffoldViewer(page: any, where: string) {
  const aim = await page.evaluate(() => (window as any).__stampScaffold());
  console.log(`[probe] ${where} scaffold viewer aiming: ${JSON.stringify(aim)}`);
  expect(aim.found, `${where}: no [name="viewer-Scaffold-Tree"] element belongs to the Scaffold Tree viewer ` +
    `over the CURRENT table view (${aim.domMatches} such element(s) in the document — a match that is not ` +
    'this view\'s viewer would read a foreign dataframe, chem.md:378-395)').toBe(true);
  return page.locator(SCAFFOLD_TARGET);
}

async function openDatasetWithChem(page: any) {
  await page.evaluate(async (path: string) => {
    document.querySelectorAll('.d4-dialog').forEach((d) => {
      const cancel = d.querySelector('[name="button-CANCEL"]') as HTMLElement | null;
      if (cancel) cancel.click();
    });
    grok.shell.closeAll();
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    const df = await grok.dapi.files.readCsv(path);
    (window as any).__df = df;
    grok.shell.addTableView(df);
    await new Promise<void>((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(); });
      setTimeout(resolve, 4000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000)); // missing-event: grid-first-paint-settled
    (window as any).__errors = [];
    if (!(window as any).__errorTrapInstalled) {
      (window as any).__errorTrapInstalled = true;
      const orig = console.error;
      console.error = function (...args: any[]) {
        (window as any).__errors.push(args.map((a: any) => String(a)).join(' '));
        orig.apply(console, args as any);
      };
    }
  }, datasetPath);
  await page.evaluate(installScaffoldResolver);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 30000});
  await waitForChemMenu(page);
  await waitForMolecule(page);
}

// Generation is driven through the viewer's own magic-wand control, not through the
// viewer.generateTree() JS API: this is a ui-smoke scenario and the control is the documented
// actuation (grok-browser chem.md:329). The icon is scoped to the viewer whose dataframe is the
// one under test, per the aiming trap in chem.md:378-395.
async function generateTreeViaMagicWand(page: any, expectedRows: number) {
  return page.evaluate(async ({rows, limit, cap}: {rows: number; limit: number; cap: number}) => {
    const viewer: any = Array.from(grok.shell.tv.viewers)
      .find((v: any) => /Scaffold Tree/i.test(v.type || '') && v.dataFrame && v.dataFrame.rowCount === rows);
    if (!viewer) return {viewerFound: false, distinct: -1, iconFound: false, inactive: null, clicked: false,
      threw: false, reason: 'no Scaffold Tree viewer over a ' + rows + '-row dataframe', nodes: 0};
    const molCol = viewer.molColumn ?? viewer.dataFrame.columns.bySemType('Molecule');
    const distinct = molCol ? molCol.categories.length : -1;
    const icon = viewer.root.querySelector('[aria-label="Generate"]') as HTMLElement | null;
    const inactive = icon ? icon.classList.contains('inactive') : null;
    if (!icon || inactive || distinct <= 0 || distinct >= limit)
      return {viewerFound: true, distinct, iconFound: !!icon, inactive, clicked: false, threw: false,
        reason: 'generation precondition not met', nodes: 0};
    // `clicked` records a click event actually reaching the magic-wand element, not the absence of
    // an exception: the product's capture-phase blocker cancels the gesture upstream of the icon
    // without throwing, and a listener on the icon itself is what separates the two.
    let dispatched = false;
    const onClick = () => { dispatched = true; };
    icon.addEventListener('click', onClick);
    let threw = false;
    let reason = '';
    try { icon.click(); }
    catch (e) { threw = true; reason = String(e); }
    icon.removeEventListener('click', onClick);
    let nodes = 0;
    const deadline = Date.now() + cap;
    while (Date.now() < deadline) {
      nodes = Array.from(viewer.root.querySelectorAll('.d4-tree-view-node'))
        .filter((n: any) => n.querySelector('canvas.chem-canvas')).length;
      if (nodes > 0) break;
      await new Promise((r) => setTimeout(r, 500));
    }
    return {viewerFound: true, distinct, iconFound: true, inactive, clicked: dispatched, threw, reason, nodes};
  }, {rows: expectedRows, limit: MAX_DISTINCT_MOLECULES, cap: 60_000});
}

// Reaching a row count is not the claim — HOLDING it is. The substructure leg restores through an
// RDKit worker, so the count can pass through a value on its way somewhere else. Once it hits, dwell:
// require the armed onRowsFiltered counter and the count itself to be unchanged `dwell` ms later,
// otherwise keep waiting. Never reaching it burns the cap and fails on the caller's assertion, which
// is what should happen.
async function waitForSettledTrueCount(page: any, expected: number, capMs: number) {
  return page.evaluate(async ({cap, dwell, want}: {cap: number; dwell: number; want: number}) => {
    const d = grok.shell.tv.dataFrame;
    const started = Date.now();
    let dwellHeld = false;
    let eventsAtHit = -1;
    while (Date.now() - started < cap) {
      if (d.filter.trueCount !== want) {
        await new Promise((r) => setTimeout(r, 500));
        continue;
      }
      eventsAtHit = (window as any).__filterEvents ?? -1;
      await new Promise((r) => setTimeout(r, dwell));
      if (d.filter.trueCount === want && ((window as any).__filterEvents ?? -1) === eventsAtHit) {
        dwellHeld = true;
        break;
      }
    }
    return {trueCount: d.filter.trueCount, rowCount: d.rowCount, dwellHeld, eventsAtHit,
      waitedMs: Date.now() - started, events: (window as any).__filterEvents ?? -1};
  }, {cap: capMs, dwell: 2000, want: expected});
}

function expectTreeGenerated(gen: any, where: string) {
  console.log(`[probe] ${where} scaffold generation: ${JSON.stringify(gen)}`);
  expect(gen.viewerFound, `${where}: no Scaffold Tree viewer over the table under test (${gen.reason})`).toBe(true);
  expect(gen.distinct, `${where}: the fixture must carry fewer than ${MAX_DISTINCT_MOLECULES} DISTINCT molecules ` +
    `for Generate to be enabled at all — measured ${gen.distinct}`).toBeLessThan(MAX_DISTINCT_MOLECULES);
  expect(gen.distinct, `${where}: the molecule column reported ${gen.distinct} distinct values`).toBeGreaterThan(0);
  expect(gen.iconFound, `${where}: the magic-wand [aria-label="Generate"] control is not in the viewer`).toBe(true);
  expect(gen.inactive, `${where}: Generate carries the "inactive" class, so the click is swallowed by the ` +
    `product's own capture-phase blocker and nothing is generated (distinct=${gen.distinct})`).toBe(false);
  expect(gen.clicked, `${where}: no click event reached the magic-wand element — the icon was detached ` +
    `mid-step, or an ANCESTOR cancelled the gesture in capture phase before it arrived (${gen.reason}). ` +
    'This does NOT cover the product\'s own disabled-state blocker: that one is registered on the icon ' +
    'itself and only calls stopPropagation, which leaves listeners on the same node running, so the flag ' +
    'still reads true with it armed (chem.md:363-367). The blocker is already covered on both sides — by ' +
    'the "inactive" check above and the node count below — so do not add a third check for it. ' +
    '"Did not throw" cannot stand in for "did run".').toBe(true);
  expect(gen.threw, `${where}: actuating Generate threw: ${gen.reason}`).toBe(false);
  expect(gen.nodes, `${where}: Generate produced no scaffold nodes carrying canvas.chem-canvas — an empty tree ` +
    'has no criterion to arm and nothing downstream can observe its round-trip').toBeGreaterThan(0);
}

function nullErrorsExpr() {
  return () => {
    const ambient = [
      /ProjectMeta\.publish|Project\.upload|GrokClientBase\._postItem|HttpDataSource\.save/,
      /^\s*Stack trace\b/,
      /Unable to find element in cloned iframe/i,
      /Permissions policy violation/i,
      /Failed to load resource/i,
    ];
    const errs = ((window as any).__errors as string[] | undefined) ?? [];
    return errs
      .filter((e) => !ambient.some((a) => a.test(e)))
      .filter((e) => /(^|\s)error:\s*null(\s|$)/i.test(e) ||
        /^\s*null\s*$/i.test(e) ||
        /Cannot read propert(?:y|ies)[^]*of (null|undefined)/i.test(e));
  };
}

// console.error is a real error channel for the reopen path: shell_project.dart:412-414 calls
// log.error(x, s) on a failed project open, and the client logger writes every LogLevel.ERROR to
// window.console.error (features/logging/logger.dart:35). It is not the ONLY one — Balloon.error
// logs to the console only when called with logError: true (d4 balloon.dart:131-134), so a
// user-visible error balloon can coexist with a clean console. Both channels are therefore read.
// Each channel carries its own sentinel round-trip as a positive control: without one, an empty
// list is equally consistent with "no errors fired" and "the collector was never installed, was
// overwritten, or is watching the wrong thing".
async function readErrorChannels(page: any) {
  const nullErrors = await page.evaluate(nullErrorsExpr());
  const rest = await page.evaluate(async () => {
    const sentinel = 'consoleChannelProbe-' + Date.now();
    console.error(sentinel);
    const errs = ((window as any).__errors as string[] | undefined) ?? [];
    const trapLive = errs.some((e) => e.includes(sentinel));
    (window as any).__errors = errs.filter((e) => !e.includes(sentinel));
    // Same treatment for the balloon arm: an installed observer is presence, not function — it can
    // be observing a detached node or filtering the class away. Put a real .d4-balloon.error into
    // the tree and require the collector to have recorded its text before its silence is trusted.
    const balloonSentinel = 'balloonChannelProbe-' + Date.now();
    const probe = document.createElement('div');
    probe.className = 'd4-balloon error';
    probe.textContent = balloonSentinel;
    document.body.appendChild(probe);
    await new Promise((r) => setTimeout(r, 300));
    const balloonWatcherLive = (((window as any).__errorBalloons as string[] | undefined) ?? [])
      .includes(balloonSentinel);
    probe.remove();
    await new Promise((r) => setTimeout(r, 100));
    const errorBalloons = (((window as any).__errorBalloons as string[] | undefined) ?? [])
      .filter((t) => t !== balloonSentinel);
    (window as any).__errorBalloons = errorBalloons;
    return {trapLive, balloonWatcherLive, errorBalloons};
  });
  return {nullErrors, ...rest};
}

function expectNoReopenErrors(channels: any, where: string) {
  console.log(`[probe] ${where} error channels: ${JSON.stringify(channels)}`);
  expect(channels.trapLive, `${where}: the console.error trap did not capture its own sentinel, so an empty ` +
    'null-error list proves nothing').toBe(true);
  expect(channels.balloonWatcherLive, `${where}: the error-balloon collector did not record a .d4-balloon.error ` +
    'placed in the document under it, so an empty balloon list proves nothing').toBe(true);
  expect(channels.nullErrors, `${where}: null-deserialization console errors fired on reopen: ` +
    `${JSON.stringify(channels.nullErrors)}`).toEqual([]);
  expect(channels.errorBalloons, `${where}: error balloons raised during reopen: ` +
    `${JSON.stringify(channels.errorBalloons)}`).toEqual([]);
}

async function saveOpenProject(page: any, name: string) {
  await page.evaluate(() => {
    const btn = document.querySelector('[name="button-Save"]') as HTMLElement | null;
    if (btn) btn.click();
  });
  await page.locator('[name="dialog-Save-project"]').waitFor({timeout: 12000});
  return page.evaluate(async ({projectName, cap, dlgCap}: {projectName: string; cap: number; dlgCap: number}) => {
    const dlg = document.querySelector('[name="dialog-Save-project"]')!;
    const nameInput = dlg.querySelector('input[type="text"]') as HTMLInputElement | null;
    if (nameInput) {
      nameInput.focus();
      nameInput.value = projectName;
      nameInput.dispatchEvent(new Event('input', {bubbles: true}));
      nameInput.dispatchEvent(new Event('change', {bubbles: true}));
      nameInput.blur();
    }
    // Two independent facts, reported separately: whether the host is in the dialog at all (a
    // renamed selector must not read as "datasync not offered"), and what display it carries.
    const dsHost = dlg.querySelector('[name="input-host-Data-sync"]') as HTMLElement | null;
    const dsHostFound = !!dsHost;
    const dsDisplay = dsHost ? getComputedStyle(dsHost).display : null;
    // The decorator is looked up on the host whether or not the host is displayed. Scoping the
    // lookup to the visible branch would leave the state readable only where it is already known,
    // and a switch sitting "on" underneath display:none — a save carrying datasync nobody actuated
    // — is precisely the state the caller asserts against. Only the CLICK is gated on visibility.
    const dsSwitch = dsHost ? dsHost.querySelector('.ui-input-switch') as HTMLElement | null : null;
    const dsSwitchFound = !!dsSwitch;
    if (dsSwitch && dsDisplay !== 'none') dsSwitch.click();
    await new Promise((r) => setTimeout(r, 500)); // missing-event: dialog-input-committed
    // Read the switch back off its own decorator class rather than trusting the click.
    const datasyncEnabled = !!dsSwitch && dsSwitch.classList.contains('ui-input-switch-on');
    (dlg.querySelector('[name="button-OK"]') as HTMLElement).click();
    // Teardown of the dialog is not synchronous with the OK click, and nothing in the product
    // promises it is prompt. Wait for the element to leave the document instead of sampling once;
    // exhausting the cap leaves dlgClosed false, which the caller's assertion turns red.
    const dlgWaitStarted = Date.now();
    while (Date.now() - dlgWaitStarted < dlgCap) {
      if (!document.querySelector('[name="dialog-Save-project"]')) break;
      await new Promise((r) => setTimeout(r, 250));
    }
    const dlgClosed = !document.querySelector('[name="dialog-Save-project"]');
    const dlgCloseMs = Date.now() - dlgWaitStarted;
    let savedId = null;
    const deadline = Date.now() + cap;
    while (Date.now() < deadline) {
      try {
        const p = await grok.dapi.projects.filter('name = "' + projectName + '"').first();
        if (p) { savedId = p.id; break; }
      } catch (e) { }
      await new Promise((r) => setTimeout(r, 500));
    }
    return {
      nameFound: !!nameInput,
      dlgClosed, dlgCloseMs,
      savedId, dsHostFound, dsDisplay, dsSwitchFound, datasyncEnabled,
    };
  }, {projectName: name, cap: 30000, dlgCap: 20000});
}

async function closeAllAndReopen(page: any, name: string) {
  return page.evaluate(async ({projectName, cap}: {projectName: string; cap: number}) => {
    (window as any).__errors = [];
    // Error balloons are collected as they appear rather than snapshotted at the end: an error
    // raised early in the restore can be gone by the time the row set settles.
    (window as any).__errorBalloons = [];
    if (!(window as any).__balloonObserver) {
      const record = () => {
        const seen = (window as any).__errorBalloons as string[];
        document.querySelectorAll('.d4-balloon.error').forEach((b) => {
          const text = (b.textContent ?? '').trim();
          if (text && !seen.includes(text)) seen.push(text);
        });
      };
      const observer = new MutationObserver(record);
      observer.observe(document.body, {childList: true, subtree: true});
      (window as any).__balloonObserver = observer;
      record();
    }
    grok.shell.closeAll();
    await new Promise((r) => setTimeout(r, 2500)); // missing-event: shell-close-all-settled
    const tvAfterClose = Array.from(grok.shell.tableViews).length;
    const project = await grok.dapi.projects.filter('name = "' + projectName + '"').first();
    const found = !!project;
    if (project) {
      await project.open();
      const deadline = Date.now() + cap;
      while (Date.now() < deadline) {
        const tv = grok.shell.tv;
        if (tv && tv.dataFrame && Array.from(tv.viewers).some((v: any) => /Scaffold Tree/i.test(v.type || '')) &&
            (window as any).__scaffoldEl()) break;
        await new Promise((r) => setTimeout(r, 250));
      }
    }
    const tv = grok.shell.tv;
    const df = tv ? tv.dataFrame : null;
    // Filtering is NOT applied at the moment the layout comes back — the criteria are restored
    // first and the row set is recomputed after. Reading trueCount here without waiting measures
    // the gap, not the round-trip. Race the dataframe's own event against a cap and record which
    // of the two ended the wait, so a silent never-fires is distinguishable from a slow apply.
    // Count the recomputes instead of waiting for the first one. The substructure filter restores
    // ASYNCHRONOUSLY: applyState starts an RDKit worker search and the bitset is applied only when
    // it returns (Chem/src/widgets/chem-substructure-filter.ts:576). A layout restore fires
    // onRowsFiltered well before that, so reading at the first event measures the gap, not the
    // round-trip. The settle wait lives in the step; this only arms the counter.
    if (df) {
      (window as any).__filterEvents = 0;
      df.onRowsFiltered.subscribe(() => { (window as any).__filterEvents++; });
    }
    return {
      found, tvAfterClose,
      tvAfterOpen: Array.from(grok.shell.tableViews).length,
      viewers: tv ? Array.from(tv.viewers).map((v: any) => v.type) : null,
      molSem: (() => {
        if (!df) return null;
        const mc = df.columns.toList().find((c: any) => c.semType === 'Molecule');
        return mc ? mc.semType : null;
      })(),
      filterTrue: df ? df.filter.trueCount : null,
      rowCount: df ? df.rowCount : null,
      scaffoldViewerEl: !!(window as any).__scaffoldEl(),
    };
  }, {projectName: name, cap: 30000});
}

test('Chem: project save and reopen with Chem state (GROK-17595)', async ({page}) => {
  test.setTimeout(360_000);
  await loginToDatagrok(page);

  // Carried across steps: the two filtered row counts the round-trip is judged against.
  let subOnlyTrueCount = 0;
  let preSaveTrueCount = 0;

  try {
    await softStep('Setup: open spgi-100.csv and confirm the Molecule column', async () => {
      await openDatasetWithChem(page);
      const setup = await page.evaluate(() => {
        const df = (window as any).__df;
        const mc = df.columns.toList().find((c: any) => c.semType === 'Molecule');
        return {molSem: mc ? mc.semType : null, molCol: mc ? mc.name : null, rowCount: df.rowCount};
      });
      expect(setup.molSem, `${datasetPath} must expose a Molecule column (found: ${setup.molCol})`).toBe('Molecule');
      expect(setup.rowCount, `${datasetPath} should carry ${DATASET_ROWS} rows`).toBe(DATASET_ROWS);
    });

    await softStep('Scenario 1 (step 1): apply Substructure Search benzene via Chem top-menu', async () => {
      await openChemMenuItem(page, 'Substructure Search...', {delayMs: 700});

      // With a single Molecule column the command goes straight to the filter sketcher; with
      // several it first opens a column picker that has no SMILES input at all
      // (Chem/src/package.ts:744-766, branch on bySemTypeAll length). spgi-100 carries seven —
      // Structure plus the Core/R* R-group columns — so the picker step is not optional here.
      // Whether the picker appears is decided by the product, not by us: the command branches on
      // the Molecule-column count (Chem/src/package.ts:744-766). Read that count from the data and
      // require the matching dialog, rather than probing for a selector and skipping when it does
      // not match — a mistyped selector would then silently no-op instead of failing.
      const molCount = await page.evaluate(() =>
        grok.shell.tv.dataFrame.columns.bySemTypeAll('Molecule').length);
      console.log(`[probe] Molecule columns in ${datasetPath}: ${molCount} (expected ${MOL_COL_COUNT})`);
      expect(molCount, `${datasetPath} must expose exactly ${MOL_COL_COUNT} Molecule columns — that count is what ` +
        'sends the command down the picker branch. A fixture that lost them would take the no-picker branch and ' +
        'delete the preselection assertions below without any assertion turning red.').toBe(MOL_COL_COUNT);
      if (molCount > 1) {
        const picker = page.locator('[name="dialog-Substructure-search"]');
        await picker.waitFor({timeout: 15_000});
        const shown = (await picker.locator('[name="input-host-Molecules"] .d4-column-selector-column')
          .first().textContent() ?? '').trim();
        expect(shown, `with ${molCount} Molecule columns the picker must preselect ${MOL_COL}; another column ` +
          'would run the search over R-group fragments instead of the structures').toBe(MOL_COL);
        await picker.locator('[name="button-OK"]').first().click();
      }
      await page.locator('.d4-dialog input[placeholder*="SMILES" i]').waitFor({timeout: 12000});
      await page.evaluate(async () => {
        const dlg = Array.from(document.querySelectorAll('.d4-dialog'))
          .find((d) => d.querySelector('input[placeholder*="SMILES" i]'))!;
        const smi = dlg.querySelector('input[placeholder*="SMILES" i]') as HTMLInputElement;
        smi.focus();
        smi.value = 'c1ccccc1';
        smi.dispatchEvent(new Event('input', {bubbles: true}));
        smi.dispatchEvent(new Event('change', {bubbles: true}));
        smi.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', code: 'Enter', keyCode: 13, bubbles: true}));
        await new Promise((r) => setTimeout(r, 1000)); // missing-event: sketcher-query-committed
        (dlg.querySelector('[name="button-OK"]') as HTMLElement).click();
      });
      await expect.poll(async () => page.evaluate(() => {
        const df = (window as any).__df;
        return df.filter.trueCount > 0 && df.filter.trueCount < df.rowCount;
      }), {timeout: 15000, intervals: [500, 1000, 1500]}).toBe(true);
      subOnlyTrueCount = await page.evaluate(() => (window as any).__df.filter.trueCount);
      console.log(`[probe] filter.trueCount after benzene substructure search: ${subOnlyTrueCount} of ` +
        `${DATASET_ROWS}`);
      expect(subOnlyTrueCount, `benzene substructure search must filter to a strict subset of the ${DATASET_ROWS} rows`)
        .toBeLessThan(DATASET_ROWS);
      expect(subOnlyTrueCount, 'benzene substructure search must leave a non-empty row set').toBeGreaterThan(0);
    });

    await softStep('Scenario 1 (step 2): add the Scaffold Tree viewer via Chem top-menu', async () => {
      await openChemMenuItem(page, 'Scaffold Tree', {delayMs: 700});
      await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 15000});
      const viewers = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((v: any) => v.type));
      expect(viewers, 'Scaffold Tree viewer must attach to the table view').toContain('Scaffold Tree');
    });

    await softStep('Scenario 1 (step 3): generate the scaffold tree from the magic-wand control', async () => {
      const gen = await generateTreeViaMagicWand(page, DATASET_ROWS);
      expectTreeGenerated(gen, 'Scenario 1');

      // Generating the tree does not filter anything: filtering starts when a scaffold node is
      // CHECKED. The checkbox is read back through the product, not through the DOM — a DOM
      // :checked readback is the test reading back its own click, and grok-browser chem.md:670-674
      // and :1139-1143 record that for these Dart tree checkboxes it can read true while the model never
      // moved. The product-side consequence is the row set: checking a node ANDs a scaffold BitSet
      // into the filter, so trueCount must drop below the substructure-only count
      // (chem.md:441-443 measures 50 -> 47 on one click).
      const viewer = await aimAtScaffoldViewer(page, 'Scenario 1 step 3');
      const node = viewer.locator('.d4-tree-view-node input[type="checkbox"]').first();
      await node.waitFor({timeout: 45_000, state: 'attached'});
      await node.check({force: true});
      await expect.poll(() => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount),
        {timeout: 30_000, intervals: [500, 1000, 1500],
          message: 'checking the first scaffold node must narrow the row set below the ' +
            `${subOnlyTrueCount} rows the substructure search alone passes — an unmoved count means no ` +
            'scaffold criterion was applied, and the round-trip below would then be measuring the ' +
            'substructure filter twice'})
        .toBeLessThan(subOnlyTrueCount);
      preSaveTrueCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
      console.log(`[probe] filter.trueCount with the scaffold node checked: ${preSaveTrueCount} ` +
        `(substructure-only was ${subOnlyTrueCount} of ${DATASET_ROWS})`);
      expect(preSaveTrueCount, 'the scaffold criterion must leave a non-empty row set — an empty one is below ' +
        'the substructure-only count too, and would satisfy the narrowing poll without being a usable ' +
        'round-trip reference').toBeGreaterThan(0);
    });

    await softStep('Scenario 1 (step 4): save the project; record the observed Data-sync state', async () => {
      const saved = await saveOpenProject(page, roundtripProject);
      expect(saved.nameFound, 'Save project dialog name input not found').toBe(true);
      console.log(`[probe] Scenario 1 save dialog teardown: closed=${saved.dlgClosed} after ` +
        `${saved.dlgCloseMs} ms`);
      expect(saved.dlgClosed, `Save project dialog was still in the document ${saved.dlgCloseMs} ms after OK — ` +
        'the dialog was waited out to detachment, not sampled once, so this is a dialog that never closed')
        .toBe(true);
      expect(saved.savedId, `saved project "${roundtripProject}" not found on the server after OK`).toBeTruthy();
      console.log(`[probe] Scenario 1 save datasync: hostFound=${saved.dsHostFound}, display=${saved.dsDisplay}, ` +
        `switchFound=${saved.dsSwitchFound}, enabled=${saved.datasyncEnabled}`);
      expect(saved.dsHostFound,
        '[name="input-host-Data-sync"] is not in the Save-project dialog at all. A renamed or removed host ' +
        'reads identically to "datasync is not offered", which would silently disarm the declared-reduction ' +
        'tripwire below — re-find the selector before trusting either of them.').toBe(true);
      expect(saved.dsDisplay,
        'DECLARED REDUCTION "DATASYNC NOT EXERCISED": the Data-sync switch is display:none for a table ' +
        'opened without a creation script, so the GROK-17595 triple runs here as a PAIR (substructure + ' +
        'scaffold, no datasync). If this assertion fails the switch IS now presented — the reduction no ' +
        'longer holds and both the scenario md and this spec must be re-strengthened to drive datasync.')
        .toBe('none');
      expect(saved.dsSwitchFound,
        'no .ui-input-switch decorator inside [name="input-host-Data-sync"]. The switch state below is read ' +
        'off that class, so a missing decorator reads as "off" rather than as "not read" — re-find it before ' +
        'trusting the assertion under this one.').toBe(true);
      expect(saved.datasyncEnabled,
        'the Data-sync switch carries "ui-input-switch-on" — the project was saved with datasync ON, a state ' +
        'this spec never actuated: the switch is hidden here, so nothing clicked it. The state is read off the ' +
        'decorator regardless of the host\'s display, so this is an observation, not a restatement of the ' +
        'assertion above.').toBe(false);
    });

    await softStep('Scenario 1 (step 5): close all, reopen, and verify Chem state round-tripped (GROK-17595)', async () => {
      const reopened = await closeAllAndReopen(page, roundtripProject);
      expect(reopened.tvAfterClose, 'workspace did not close before reopen').toBe(0);
      expect(reopened.found, `saved project "${roundtripProject}" not found for reopen`).toBe(true);
      expect(reopened.tvAfterOpen, 'reopening the project did not open a TableView').toBeGreaterThan(0);
      expect(reopened.molSem, 'the Molecule semType was lost after reopen (serialization failed)')
        .toBe('Molecule');
      expect(reopened.viewers, 'GROK-17595: Scaffold Tree viewer not restored in the reopened layout')
        .toContain('Scaffold Tree');
      expect(reopened.scaffoldViewerEl, 'GROK-17595: Scaffold Tree viewer DOM element absent after reopen').toBe(true);
      // The wait ends on the exact pre-save count, not on a range: the scenario's claim is that the
      // filters come back UNCHANGED, and that number was measured before the save.
      const settled = await waitForSettledTrueCount(page, preSaveTrueCount, 60_000);
      console.log(`[probe] reopen settle: first read ${reopened.filterTrue} of ${reopened.rowCount}, waited ` +
        `${settled.waitedMs} ms, ${settled.events} onRowsFiltered recompute(s) (${settled.eventsAtHit} when the ` +
        `count first hit the expected value), dwell held: ${settled.dwellHeld}, final ${settled.trueCount} of ` +
        `${settled.rowCount}; pre-save reference ${preSaveTrueCount} (substructure-only ${subOnlyTrueCount})`);
      const settleReport =
        'GROK-17595: both restored criteria must reproduce the exact row set saved with the project — ' +
        `${preSaveTrueCount} rows (substructure search alone passed ${subOnlyTrueCount}, so a count back at ` +
        'that value means the scaffold criterion did not come back, 0 means a BitSet deserialized empty and ' +
        `${settled.rowCount} means nothing deserialized at all). Read ${reopened.filterTrue} of ` +
        `${reopened.rowCount} before the wait (that reading measures the restore gap, not the round-trip), ` +
        `then waited ${settled.waitedMs} ms and saw ${settled.events} onRowsFiltered recompute(s); final ` +
        `${settled.trueCount} of ${settled.rowCount}. The substructure filter restores through an RDKit ` +
        'worker, so a late apply is expected — a wrong steady state is not.';
      expect(settled.rowCount, `${settleReport} (the reopened table itself must still carry ${DATASET_ROWS} rows)`)
        .toBe(DATASET_ROWS);
      expect(settled.trueCount, settleReport).toBe(preSaveTrueCount);
      expect(settled.dwellHeld, `${settleReport} The count reached ${preSaveTrueCount} but never held it: the ` +
        'filter kept recomputing for the rest of the cap, so what was observed is a value the row set passed ' +
        'through, not the state the project restored to.').toBe(true);

      // Both restored criteria, independently observable. The equality above cannot separate them:
      // every row matching the root scaffold also matches benzene, so the scaffold-only row set is a
      // SUBSET of the substructure-only one and preSaveTrueCount is what a reopen that restored the
      // scaffold criterion alone would also land on. Releasing the scaffold criterion here must widen
      // the row set back to exactly the substructure-only count — which it can only do if the
      // Substructure Search came back too. A reopen that dropped it lands on the full table instead.
      const restoredViewer = await aimAtScaffoldViewer(page, 'GROK-17595 reopen');
      const restoredNode = restoredViewer.locator('.d4-tree-view-node input[type="checkbox"]').first();
      await restoredNode.waitFor({timeout: 45_000, state: 'attached'});
      // Actuated with a bare click, not uncheck(): uncheck() verifies the DOM `checked` property it
      // has just clicked, and on these Dart tree checkboxes that property is not the criterion — the
      // row set is (chem.md:441-443). Which way the two diverge on a scaffold node is not recorded
      // anywhere: chem.md:670-674 and :1139-1143 document the descriptor tree, and in the OPPOSITE
      // direction (DOM `checked` true while the model never moved). So the release is verified the
      // way the CHECK leg above is — through the product-side row count.
      const checkedBefore = await restoredNode.isChecked();
      await restoredNode.click({force: true});
      const checkedAfter = await restoredNode.isChecked();
      console.log(`[probe] scaffold node release actuation: DOM checked ${checkedBefore} -> ${checkedAfter} ` +
        '(DOM state is diagnostic only; the row count below is the criterion)');
      const released = await waitForSettledTrueCount(page, subOnlyTrueCount, 30_000);
      console.log(`[probe] scaffold criterion released after reopen: ${released.trueCount} of ` +
        `${released.rowCount} after ${released.waitedMs} ms, dwell held: ${released.dwellHeld} ` +
        `(substructure-only reference ${subOnlyTrueCount}, with scaffold ${preSaveTrueCount})`);
      const releaseReport =
        'GROK-17595: with the restored scaffold node unchecked the row set must widen back to exactly the ' +
        `${subOnlyTrueCount} rows the Substructure Search passes on its own. ${DATASET_ROWS} means the ` +
        'substructure criterion did NOT come back — the reopen restored the scaffold leg only, and because ' +
        `every root-scaffold row also matches benzene that reopen still reproduces ${preSaveTrueCount} above. ` +
        `A value still at ${preSaveTrueCount} means unchecking did not release the scaffold criterion, so the ` +
        `equality above was not measuring a live scaffold filter either. Read ${released.trueCount} of ` +
        `${released.rowCount}.`;
      expect(released.trueCount, releaseReport).toBe(subOnlyTrueCount);
      expect(released.dwellHeld, `${releaseReport} The count reached the value but never held it.`).toBe(true);

      expectNoReopenErrors(await readErrorChannels(page), 'GROK-17595 reopen');
    });

    await softStep('Scenario 2 (steps 1-2): Scaffold Tree alone — add viewer and generate scaffold nodes', async () => {
      await openDatasetWithChem(page);
      await openChemMenuItem(page, 'Scaffold Tree', {delayMs: 700});
      await page.locator('[name="viewer-Scaffold-Tree"]').waitFor({timeout: 15000});
      const viewers = await page.evaluate(() => Array.from(grok.shell.tv.viewers).map((v: any) => v.type));
      expect(viewers, 'Scaffold Tree viewer must attach in Scenario 2').toContain('Scaffold Tree');
      expectTreeGenerated(await generateTreeViaMagicWand(page, DATASET_ROWS), 'Scenario 2');
    });

    await softStep('Scenario 2 (step 3): save, close all, reopen — Scaffold Tree round-trips cleanly', async () => {
      const saved = await saveOpenProject(page, scaffoldOnlyProject);
      expect(saved.nameFound, 'Save project dialog did not open in Scenario 2').toBe(true);
      expect(saved.savedId, `scaffold-only project "${scaffoldOnlyProject}" not saved`).toBeTruthy();
      console.log(`[probe] Scenario 2 save datasync: hostFound=${saved.dsHostFound}, display=${saved.dsDisplay}, ` +
        `switchFound=${saved.dsSwitchFound}, enabled=${saved.datasyncEnabled}`);
      expect(saved.dsHostFound,
        '[name="input-host-Data-sync"] is not in the Save-project dialog at all — see the Scenario 1 note: a ' +
        'renamed host disarms the declared-reduction tripwire instead of firing it').toBe(true);
      expect(saved.dsDisplay,
        'DECLARED REDUCTION "DATASYNC NOT EXERCISED": the Data-sync switch is display:none for a table ' +
        'opened without a creation script, so the scaffold-only baseline also saves without datasync. ' +
        'If this assertion fails the switch IS now presented and the reduction must be revisited.')
        .toBe('none');
      expect(saved.dsSwitchFound,
        'no .ui-input-switch decorator inside [name="input-host-Data-sync"] — see the Scenario 1 note: a ' +
        'missing decorator reads as "off" rather than as "not read"').toBe(true);
      expect(saved.datasyncEnabled,
        'the Data-sync switch carries "ui-input-switch-on" — the scaffold-only baseline was saved with datasync ' +
        'ON, which nothing in this spec actuated (the state is read off the decorator regardless of display)')
        .toBe(false);

      const reopened = await closeAllAndReopen(page, scaffoldOnlyProject);
      expect(reopened.tvAfterOpen, 'scaffold-only project did not open a TableView').toBeGreaterThan(0);
      expect(reopened.viewers, 'Scaffold Tree viewer not restored on scaffold-only reopen').toContain('Scaffold Tree');
      expect(reopened.scaffoldViewerEl, 'Scaffold Tree viewer DOM element absent on scaffold-only reopen').toBe(true);
      expectNoReopenErrors(await readErrorChannels(page), 'scaffold-only reopen');
    });
  } finally {
    await page.evaluate(async (names: string[]) => {
      for (const name of names) {
        try {
          const p = await grok.dapi.projects.filter('name = "' + name + '"').first();
          if (p) await grok.dapi.projects.delete(p);
        } catch (e) { console.log('[cleanup] project delete threw (non-fatal):', String(e)); }
      }
      grok.shell.closeAll();
    }, [roundtripProject, scaffoldOnlyProject]).catch(() => {});
  }

  finishSpec();
});
