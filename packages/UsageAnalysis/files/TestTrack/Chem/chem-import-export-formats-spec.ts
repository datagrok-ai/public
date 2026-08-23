/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

// Recon-notes (run of 2026-08-22, 3 attempts):
//   The control-column search in structureInkReadings never reached its row-number-gutter
//   fallback on any of the three measured tables: it found a real data column each time
//   ("CAS number" on ApprovedDrugs2015.sdf, "Name" on the built V3000 SD file, "StdInChI" on
//   aspirin.mol — SDF and molfile property fields become ordinary columns on import). The
//   gutter branch is therefore untested rather than shown to work; the false-red risk the
//   fallback was flagged for is neither confirmed nor ruled out.

// Clears the console-error buffer and installs the capture hook once. Every step that later
// claims "no console.error calls" must call this at its own start: a buffer carried in from earlier
// steps makes the claim cover something other than the step whose label it wears.
// Armed the way the balloon channel is armed further down: a monkey-patched console.error that
// was overwritten by a later page load, or never installed, is indistinguishable from a clean
// run, so a known error is planted and its capture asserted before the buffer is used as
// evidence of silence. The probe is cleared out of the buffer in the same evaluate, so nothing
// downstream sees it.
const CAPTURE_PROBE = 'capture-probe: console-error channel check (chem-import-export-formats-spec)';

async function beginErrorCapture(page: Page): Promise<void> {
  const probed = await page.evaluate((probe) => {
    const w = window as any;
    w.__errs = [];
    if (!w.__errsHooked) {
      w.__errsHooked = true;
      const orig = console.error;
      console.error = function (...a: any[]) {
        w.__errs.push(a.map((x: any) => String(x)).join(' '));
        orig.apply(console, a as any);
      };
    }
    console.error(probe);
    const caught = (w.__errs as string[]).filter((e) => e.indexOf(probe) >= 0).length;
    w.__errs = [];
    return caught;
  }, CAPTURE_PROBE);
  expect(probed, 'the console-error capture hook must record a planted error, or an empty buffer later in the step is evidence of nothing')
    .toBe(1);
}

async function openChemFile(page: Page, path: string): Promise<string> {
  await page.evaluate(() => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
  });
  await beginErrorCapture(page);
  const table = await page.evaluate(async (p) => {
    const fns = DG.Func.find({name: 'OpenFile'});
    let df: any;
    if (/\.(sdf|mol)$/i.test(p) && fns?.length) {
      const res = await fns[0].apply({fullPath: p});
      df = Array.isArray(res) ? res[0] : res;
      if (!grok.shell.tv || grok.shell.tv.dataFrame !== df) grok.shell.addTableView(df);
    } else {
      df = await grok.dapi.files.readCsv(p);
      grok.shell.addTableView(df);
    }
    (window as any).__df = grok.shell.t;
    await new Promise<void>((res) => {
      const sub = grok.shell.t.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
      setTimeout(res, 8000);
    });
    return String(grok.shell.t.name);
  }, path);
  await waitForChemMenu(page);
  // Scoped to the current view, like gridNonBlankPixels below: '[name="viewer-Grid"]' resolved
  // in document order also matches filter mini-grids and a grid left behind by an earlier view,
  // so the barrier could clear on a grid nothing downstream reads.
  let ready = false;
  for (let i = 0; i < 150 && !ready; i++) {
    ready = await page.evaluate(() => {
      const tv = grok.shell.tv;
      const grid = tv && tv.root ? tv.root.querySelector('.d4-grid[name="viewer-Grid"]') : null;
      return !!grid && grid.querySelectorAll('canvas').length > 0;
    });
    if (!ready) await page.waitForTimeout(200);
  }
  expect(ready, `the current view's grid canvas must appear after opening ${path}`).toBe(true);
  await page.waitForTimeout(3000);
  return table;
}

// Reads the CURRENT VIEW's dataframe, not grok.shell.t / window.__df. waitForMolecule already
// resolves on "shell.t or __df carries a Molecule column", so facts read off those would restate
// the barrier; read here, the semType and the table name are a claim about the table whose grid
// the pixel readings below measure, which the barrier does not establish.
async function moleculeColumnFacts(page: Page):
  Promise<{semType: string; cellRenderer: string; rows: number; molCol: string; table: string}> {
  return await page.evaluate(() => {
    const tv = grok.shell.tv;
    const df = tv && tv.dataFrame ? tv.dataFrame : null;
    const mc = df ? df.columns.toList().find((c: any) => c.semType === 'Molecule') : null;
    return {
      semType: mc ? mc.semType : 'NONE',
      cellRenderer: mc ? mc.meta.cellRenderer : 'NONE',
      rows: df ? df.rowCount : -1,
      molCol: mc ? mc.name : 'NONE',
      table: df ? String(df.name) : 'NONE',
    };
  });
}

async function gridNonBlankPixels(page: Page): Promise<{table: string; pixels: number}> {
  return await page.evaluate(() => {
    // Scoped to the CURRENT table view, not to document order: with a previously opened view
    // still in the document, document.querySelector('.d4-grid[name="viewer-Grid"]') returns
    // that older grid, so a new grid that paints nothing still reads the old table's canvas.
    const tv = grok.shell.tv;
    const table = tv && tv.dataFrame ? String(tv.dataFrame.name) : 'NONE';
    const grid = tv && tv.root ? tv.root.querySelector('.d4-grid[name="viewer-Grid"]') : null;
    if (!grid) return {table, pixels: -1};
    const cs = Array.from(grid.querySelectorAll('canvas'))
      .sort((a: any, b: any) => b.width * b.height - a.width * a.height) as HTMLCanvasElement[];
    if (cs.length === 0) return {table, pixels: -1};
    const c = cs[0];
    const ctx = c.getContext('2d');
    if (!ctx) return {table, pixels: -1};
    const w = Math.min(c.width, 600), h = Math.min(c.height, 400);
    if (w <= 0 || h <= 0) return {table, pixels: -1};
    const d = ctx.getImageData(0, 0, w, h).data;
    let cnt = 0;
    // Opacity is required as well as darkness: an untouched canvas reads (0,0,0,0),
    // so a colour-only rule scores a canvas that was never drawn on at its full area.
    for (let i = 0; i < d.length; i += 4)
      if (d[i + 3] > 0 && (d[i] < 250 || d[i + 1] < 250 || d[i + 2] < 250)) cnt++;
    return {table, pixels: cnt};
  });
}

// Barrier on the current view's grid: polls the count until it repeats, so the read lands
// after the new grid stopped painting rather than after a fixed sleep. The barrier's target
// is stabilization, not the assertion's floor — a grid that settles blank settles here too
// and is then rejected by the caller's threshold.
async function waitForGridPaint(page: Page, timeoutMs = 30_000): Promise<{table: string; pixels: number}> {
  const deadline = Date.now() + timeoutMs;
  let previous = -2;
  let measured = {table: 'NONE', pixels: -1};
  while (Date.now() < deadline) {
    measured = await gridNonBlankPixels(page);
    if (measured.pixels > 0 && measured.pixels === previous) break;
    previous = measured.pixels;
    await page.waitForTimeout(250);
  }
  return measured;
}

// The structure-versus-text discriminator, ported from molecule-rendering-spec.ts
// (moleculeColumnColouredPixels, lines 83-145 there). RDKit's default palette draws
// heteroatoms in saturated colour (O ~ rgb(255,0,0), N ~ rgb(0,0,255), channel spread ~255)
// while canvas text is drawn with one grey/black fill, whose three channels are equal by
// construction (spread 0). The 30/255 margin sits far above antialiasing noise and far below
// real heteroatom ink, and the caller asserts only "> 0", so no tuned pixel count is involved.
// The strip is confined to the column's own screen rectangle and starts below the column-header
// band (grid.colHeaderHeight), so the header, the row-number gutter and the neighbouring
// property columns are outside the measured area altogether.
// The control reading is the load-bearing half: a monochrome column measured on the same grid
// in the same run must score zero, which is what makes "> 0" on the molecule column mean ink
// the molecule renderer drew rather than grid chrome the two columns share. Both readings are
// taken in one pass so they cannot land on different paints of the grid.
// Returns -1 for `coloured` when a region cannot be read, so the caller rejects the read
// instead of comparing against a sentinel.
interface ColumnColourReading {
  column: string;
  coloured: number;
  headerHeight: number;
  x0: number;
  x1: number;
  y0: number;
  h: number;
}

function describeReading(r: ColumnColourReading): string {
  return `column="${r.column}" header=${r.headerHeight} strip x=[${r.x0},${r.x1}) y=[${r.y0},${r.y0 + r.h}) coloured=${r.coloured}`;
}

async function structureInkReadings(page: Page, moleculeColumn: string):
  Promise<{molecule: ColumnColourReading; control: ColumnColourReading}> {
  return await page.evaluate((name) => {
    const fail = (col: string): ColumnColourReading =>
      ({column: col, coloured: -1, headerHeight: -1, x0: -1, x1: -1, y0: -1, h: -1});
    const both = (m: ColumnColourReading, c: ColumnColourReading) => ({molecule: m, control: c});
    const tv = grok.shell.tv;
    const grid = tv && tv.root ? tv.root.querySelector('.d4-grid[name="viewer-Grid"]') : null;
    if (!grid) return both(fail(name), fail('NONE'));
    const cs = Array.from(grid.querySelectorAll('canvas'))
      .sort((a: any, b: any) => b.width * b.height - a.width * a.height) as HTMLCanvasElement[];
    if (cs.length === 0) return both(fail(name), fail('NONE'));
    const c = cs[0];
    const ctx = c.getContext('2d');
    if (!ctx || !c.width || !c.height) return both(fail(name), fail('NONE'));
    const g = tv.grid;
    const headerHeight = g ? g.colHeaderHeight : -1;
    if (!(headerHeight >= 0)) return both(fail(name), fail('NONE'));
    const scale = c.width / (c.clientWidth || c.width);
    const read = (gc: any, label: string): ColumnColourReading => {
      if (!gc) return fail(label);
      const x0 = Math.max(0, Math.round(gc.left * scale));
      const x1 = Math.min(c.width, Math.round(gc.right * scale));
      const y0 = Math.min(c.height, Math.round(headerHeight * scale));
      const h = Math.min(c.height - y0, Math.round(400 * scale));
      if (x1 - x0 <= 0 || h <= 0) return {column: label, coloured: -1, headerHeight, x0, x1, y0, h};
      const d = ctx.getImageData(x0, y0, x1 - x0, h).data;
      let coloured = 0;
      for (let i = 0; i < d.length; i += 4)
        if (Math.max(d[i], d[i + 1], d[i + 2]) - Math.min(d[i], d[i + 1], d[i + 2]) >= 30) coloured++;
      return {column: label, coloured, headerHeight, x0, x1, y0, h};
    };
    // The control is a monochrome column of this same grid. A non-Molecule data column is
    // preferred; the single-structure tables (.mol, the built V3000 SD file) carry no other data
    // column, so the row-number gutter — grid-drawn, monochrome, subject to the same row tint
    // and current-row highlight — is the control there.
    let control: any = null;
    for (let i = 1; i < g.columns.length && !control; i++) {
      const gc = g.columns.byIndex(i);
      if (gc && gc.visible && gc.column && gc.column.semType !== 'Molecule') control = gc;
    }
    const controlLabel = control ? String(control.name) : '(row header)';
    return both(read(g.columns.byName(name), name), read(control ?? g.columns.rowHeader, controlLabel));
  }, moleculeColumn);
}

// No pixel floor: the reading is confined to the molecule column's own cells and the control
// column establishes the zero, so "coloured > 0 here, exactly 0 there" discriminates structures
// from a monochrome raw-text fallback without any threshold to derive.
function expectStructureInk(r: {molecule: ColumnColourReading; control: ColumnColourReading}, tag: string): void {
  console.log(`[${tag}] molecule ${describeReading(r.molecule)}`);
  console.log(`[${tag}] control  ${describeReading(r.control)}`);
  expect(r.molecule.coloured, `${tag}: the molecule column region could not be read from the grid canvas (fault guard)`)
    .toBeGreaterThanOrEqual(0);
  expect(r.molecule.h, `${tag}: the molecule strip below the column header must have height — a zero-height read scores nothing whatever the renderer drew`)
    .toBeGreaterThan(0);
  expect(r.molecule.coloured, `${tag}: the Molecule column must be drawn as structures: RDKit paints heteroatoms in colour (O red, N blue) while a raw-text fallback is monochrome, so zero coloured pixels inside the column's own cells means no structure was drawn`)
    .toBeGreaterThan(0);
  expect(r.control.coloured, `${tag}: the control column could not be read from the grid canvas (fault guard)`)
    .toBeGreaterThanOrEqual(0);
  expect(r.control.h, `${tag}: the control strip must have height — a zero-height control scores zero for free and controls for nothing`)
    .toBeGreaterThan(0);
  expect(r.control.coloured, `${tag}: the monochrome control column must score zero coloured pixels — a non-zero count means the counter is scoring grid chrome (row tint, selection highlight) instead of structure ink, which would let the assertion above pass over raw text`)
    .toBe(0);
}

// Balloons auto-hide after roughly 5 s, so a count taken after an unbounded wait reads empty
// whether or not one appeared. These two record additions as they happen and classify them at
// read time — the class is read off the live node, so a balloon whose class is set after
// insertion is still classified correctly.
async function installBalloonRecorder(page: Page): Promise<void> {
  await page.evaluate(() => {
    const w = window as any;
    if (w.__balloonObs) w.__balloonObs.disconnect();
    w.__balloonNodes = [];
    w.__balloonObs = new MutationObserver((muts: MutationRecord[]) => {
      for (const m of muts)
        for (const n of Array.from(m.addedNodes)) {
          if (n.nodeType !== 1) continue;
          const el = n as Element;
          if (el.parentNode === document.body || String(el.className).indexOf('d4-balloon') >= 0)
            w.__balloonNodes.push(el);
        }
    });
    w.__balloonObs.observe(document.body, {childList: true, subtree: true});
  });
}

async function readRecordedBalloons(page: Page): Promise<Array<{cls: string; text: string}>> {
  return await page.evaluate(() => {
    const nodes = ((window as any).__balloonNodes ?? []) as Element[];
    const out: Array<{cls: string; text: string}> = [];
    for (const n of nodes) {
      const b = n.classList && n.classList.contains('d4-balloon') ? n : n.querySelector('.d4-balloon');
      // textContent, not innerText: a balloon that already auto-hid is detached and has no layout.
      if (b) out.push({cls: String(b.className), text: String(b.textContent || '').slice(0, 200)});
    }
    return out;
  });
}

// The zero-assertions consume the TOTAL captured count minus the named exclusions below; the
// rdkit patterns only classify what survived exclusion, so an error outside them is not invisible
// to the assertion. Widened 2026-08-22 on a measured baseline: Gate B logged captured=0
// rdkit-matched=0 other=0 at all five read points on three consecutive attempts.
// The channel is console.error alone: window.onerror, unhandled promise rejections and
// console.warn are never hooked, so what the callers assert is "no console.error call", not "no
// error". The claim also covers only the windows beginErrorCapture arms — login and platform boot
// precede the first arming, and the .mol2 step arms nothing and claims nothing.
//
// Named exclusion, added 2026-08-22 from an observation and not from expectation. What the rule
// actually covers is the stand's datagrok_reader credential failing server-side, whichever home
// widget happens to surface it — not one named widget. Two have been seen so far:
// PowerPack:RecentlySharedWithMe, in Gate B cycle
// direct-gate-b-2026-08-22-chem-import-export-formats-r7 (attempt 1 only, InChI read point only,
// three console.error entries); and PowerPack:MostRecentEntities alongside it in cycle
// direct-gate-b-2026-08-22-chem-import-export-formats-r8 (attempt 1, InChI read point, captured=4
// counted=0 named-excluded=4). Both failed the same connection with "password authentication
// failed for user datagrok_reader", so the second was excluded by the credential clause rather
// than by name, which is the rule working as written. The other read points read zero on every
// attempt of both cycles. Two occurrences in two cycles make this recurring, not a one-off: the
// exclusion makes this spec silent about the credential failure BY DESIGN, and not because the
// failure is fixed. The run log is the only remaining witness — excluded messages are counted and
// printed in full, and the log line prints the whole captured total before exclusion.
// The two identifying substrings are required narrowly — a known widget's own function name, or
// the datagrok_reader credential named together with the authentication failure — so a chem-side
// database or renderer error is not swallowed with it.
function capturedConsoleErrors(page: Page): Promise<{matched: string[]; other: string[]; excluded: string[]}> {
  return page.evaluate(() => {
    const errs = ((window as any).__errs as string[] | undefined) ?? [];
    const isAmbientNoise = (e: string) =>
      e.indexOf('PowerPack:RecentlySharedWithMe') >= 0 ||
      (/password authentication failed for user/i.test(e) && e.indexOf('datagrok_reader') >= 0);
    // The platform logger emits the stack trace of an error as a separate, content-free
    // continuation ("Translating stack trace... ID = pGCYd", then "Stack trace pGCYd"), and the id
    // is minted per occurrence, so it cannot be written into a rule. A continuation is dropped
    // only when it belongs to an already-excluded event: it either directly follows one, or it
    // carries an id that an already-excluded continuation announced. A stack trace that follows a
    // kept error is kept and counted with it.
    const idOfTrace = /Stack trace ([A-Za-z0-9]+)/;
    const idAnnounced = /Translating stack trace[\s\S]*?ID\s*=\s*([A-Za-z0-9]+)/;
    const isTraceContinuation = (e: string) => /^\s*(Translating stack trace|Stack trace\b)/.test(e);
    const knownIds: {[id: string]: boolean} = {};
    const excluded: string[] = [];
    const kept: string[] = [];
    let previousExcluded = false;
    for (const e of errs) {
      const carried = e.match(idOfTrace);
      const drop: boolean = isAmbientNoise(e) ||
        (isTraceContinuation(e) && (previousExcluded || (!!carried && knownIds[carried[1]] === true)));
      if (drop) {
        const announced = e.match(idAnnounced);
        if (announced) knownIds[announced[1]] = true;
        if (carried) knownIds[carried[1]] = true;
        excluded.push(e);
      }
      else kept.push(e);
      previousExcluded = drop;
    }
    // 'gS' is the minified Dart symbol from the GROK-16870 crash and is only matched inside its
    // quotes — bare, it matches the "gs" in settings/logs/flags and turns any error into a hit.
    const re = /rdkit[-_]cell[-_]renderer|method not found|'gS'|cellRenderer\.render|NullError/i;
    return {matched: kept.filter((e) => re.test(e)), other: kept.filter((e) => !re.test(e)), excluded};
  });
}

function logAndCountConsoleErrors(errs: {matched: string[]; other: string[]; excluded: string[]}, tag: string): number {
  const counted = errs.matched.length + errs.other.length;
  const captured = counted + errs.excluded.length;
  console.log(`[${tag}] console.error calls captured=${captured} counted=${counted} named-excluded=${errs.excluded.length} rdkit-matched=${errs.matched.length} other=${errs.other.length} matched=${JSON.stringify(errs.matched.slice(0, 3))} other=${JSON.stringify(errs.other.slice(0, 3))} excluded=${JSON.stringify(errs.excluded)}`);
  return counted;
}

test.use(specTestOptions);

test('Chem: Import/Export Formats (SDF/mol import, Save as SDF round-trip, InChI-to-Molecule converter)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  // Carried from the export step to the round-trip step: the artefact Save as
  // SDF actually wrote, and the row count it was written from.
  let exportedSdf = '';
  let exportedFromRows = 0;

  await softStep('Scenario 1 fixture note (SR-03): mol1K.sdf, named by the scenario, is absent from the server — the SDF steps below run against ApprovedDrugs2015.sdf instead', async () => {

    // `.catch(() => true)` used to sit on this probe, so an auth loss or an API
    // change produced exactly the value the assertion demanded. The probe now
    // reports its own failure and the assertion rejects it.
    const mol1K = await page.evaluate(() =>
      grok.dapi.files.exists('System:DemoFiles/chem/mol1K.sdf')
        .then((e: boolean) => (e ? 'present' : 'absent'))
        .catch((err: any) => `probe-failed: ${String(err)}`));
    expect(mol1K, `fixture note: mol1K.sdf named in the scenario must be genuinely absent, got "${mol1K}"`)
      .toBe('absent');
  });

  await softStep('Scenario 1 Step 3, Scenario 1 Step 4: import ApprovedDrugs2015.sdf — Molecule semType, renderer bound, structure ink somewhere in the column\'s first screenful against a zero-scoring control (not per cell, SR-07), record count matches "$$$$" terminators, no console.error calls', async () => {
    const sdfTable = await openChemFile(page, 'System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf');
    await waitForMolecule(page);
    const facts = await moleculeColumnFacts(page);
    expect(facts.table, 'the facts below must be read from the table this step opened, not from one left behind by an earlier step')
      .toBe(sdfTable);
    expect(facts.semType, 'SDF import must type the structure column as Molecule so the RDKit cell renderer binds').toBe('Molecule');
    expect(facts.cellRenderer, 'Molecule cell renderer must be bound to the imported SDF column').toBe('Molecule');
    const dollarCount = await page.evaluate(async () => {
      const bytes = await grok.dapi.files.readAsBytes('System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf');
      const txt = new TextDecoder().decode(bytes);
      return (txt.match(/\$\$\$\$/g) || []).length;
    });
    expect(dollarCount, 'fixture sanity: ApprovedDrugs2015.sdf must contain $$$$ record terminators').toBeGreaterThan(0);
    expect(facts.rows, 'imported row count must equal the number of "$$$$"-terminated records in the SDF file').toBe(dollarCount);
    const painted = await waitForGridPaint(page);
    console.log(`[sdf-render] measured table="${painted.table}" grid non-blank pixels = ${painted.pixels}`);
    expect(painted.table, 'the SDF pixel count must be read from the SDF table, not from a table left open by an earlier step')
      .toBe(sdfTable);
    expect(painted.pixels, 'grid canvas read failed (fault guard)').toBeGreaterThanOrEqual(0);
    expectStructureInk(await structureInkReadings(page, facts.molCol), 'sdf-render');
    const errs = await capturedConsoleErrors(page);
    expect(logAndCountConsoleErrors(errs, 'sdf-render'), `console.error calls during SDF Molecule render — rdkit-matched: ${JSON.stringify(errs.matched.slice(0, 3))}, other: ${JSON.stringify(errs.other.slice(0, 3))}`).toBe(0);
  });

  await softStep('Scenario 1 Step 5: V3000 SDF variant — build an SD file from the V3000 molfile, open it, assert one imported row typed as Molecule', async () => {
    // The step's subject is V3000 support, so the fixture must be a V3000 *SD file*.
    // mol/v3000.mol is a bare molfile with no "$$$$" record terminator, and the SD parser
    // ends a record on that terminator regardless of format version: strip "$$$$" from the
    // V2000 aspirin.mol and it yields 0 records too, add one to v3000.mol and it yields 1
    // (openchemlib 7.5.0, checked directly 2026-08-20). Opening the bare molfile therefore
    // measures the fixture, not V3000.
    const SDF_PATH = 'System:AppData/Chem/tmp-v3000-check.sdf';
    await beginErrorCapture(page);
    const v3000 = await page.evaluate(async (path) => {
      const source = await grok.dapi.files.readAsText('System:DemoFiles/chem/mol/v3000.mol');
      const isV3000 = /V3000/.test(source);
      try {
        await grok.dapi.files.writeAsText(path, source.trimEnd() + String.fromCharCode(10) + '$$$$' + String.fromCharCode(10));
        const res = await grok.functions.call('OpenFile', {fullPath: path});
        const df = Array.isArray(res) ? res[0] : res;
        if (!df) return {rows: -1, isV3000, semType: null, molCol: null, smiles: null, smilesIsStructure: false,
          convertError: null, cellIsMolblock: null, dfName: 'NONE'};
        // "Close this table" (scenario step 5): the previous view must go before this one is
        // added, or the grid read below measures the table the previous step already asserted.
        grok.shell.closeAll();
        grok.shell.addTableView(df);
        let col: any = null;
        for (let i = 0; i < 40 && !col; i++) {
          col = df.columns.toList().find((c: any) => c.semType === 'Molecule') ?? null;
          if (!col) await new Promise((r) => setTimeout(r, 250));
        }
        // The count and the semType are both test-controlled: the record ends where this
        // step wrote "$$$$", and _importSdfString stamps semType on the column without ever
        // parsing the molblock (Chem/src/open-chem/sdf-importer.ts:41). Converting the cell
        // is the assertion that actually depends on V3000 being understood.
        let smiles: string | null = null;
        let convertError: string | null = null;
        try { smiles = col ? grok.chem.convert(col.get(0), DG.chem.Notation.Unknown, DG.chem.Notation.Smiles) : null; }
        catch (e) { convertError = String(e); }
        // convert() does not throw on a parse failure: _convertMolNotation pre-seeds its result
        // with 'MALFORMED_INPUT_VALUE' and returns that (Chem/src/utils/convert-notation-utils.ts:60),
        // which is a non-empty string free of "V30" — so the assertions on the value above are all
        // satisfied by a total V3000 parse failure. Same guard the round-trip step uses:
        // checkSmiles runs through Chem:isSmiles (RDKit get_mol, package.ts:2289), a different
        // platform function from Chem:convertMolNotation, so it cannot degenerate along with it.
        const smilesIsStructure = (() => {
          try { return typeof smiles === 'string' && smiles.length > 0 && !!grok.chem.checkSmiles(smiles); }
          catch (e) { return false; }
        })();
        return {rows: df.rowCount, isV3000, semType: col ? col.semType : null, molCol: col ? col.name : null,
          smiles, smilesIsStructure, convertError, cellIsMolblock: col ? /V30/.test(String(col.get(0))) : null,
          dfName: String(df.name)};
      } finally {
        try { await grok.dapi.files.delete(path); } catch (e) {}
      }
    }, SDF_PATH);
    expect(v3000.isV3000, 'fixture sanity: mol/v3000.mol must be a V3000-format molfile').toBe(true);
    expect(v3000.rows, 'a V3000 SD record must import as one row').toBe(1);
    expect(v3000.semType, `the V3000 structure column must be typed as Molecule (column: ${v3000.molCol})`).toBe('Molecule');
    expect(v3000.cellIsMolblock, 'fixture sanity: the imported cell must still hold the V3000 molblock').toBe(true);
    console.log(`[v3000-convert] smiles=${JSON.stringify(v3000.smiles)} isStructure=${v3000.smilesIsStructure} convertError=${v3000.convertError}`);
    expect(v3000.convertError,
      `converting the V3000 cell must not throw — a throw means the V3000 molblock was not parsed; got ${v3000.convertError}`)
      .toBeNull();
    expect(typeof v3000.smiles, `converting the V3000 cell must return a string; got ${JSON.stringify(v3000.smiles)}`).toBe('string');
    expect(v3000.smiles!.length,
      `the V3000 molblock must convert to a non-empty SMILES — an empty result means the structure was never parsed; got ${JSON.stringify(v3000.smiles)}`)
      .toBeGreaterThan(0);
    expect(v3000.smiles,
      'converting the V3000 cell must yield SMILES, not the molblock echoed back unparsed')
      .not.toMatch(/V30/);
    expect(v3000.smilesIsStructure,
      `the conversion result must itself be a parseable structure — convert() returns the sentinel 'MALFORMED_INPUT_VALUE' instead of throwing when the molblock is not parsed, and that sentinel satisfies every assertion above; got ${JSON.stringify(v3000.smiles)}`)
      .toBe(true);
    const painted = await waitForGridPaint(page);
    console.log(`[v3000-render] measured table="${painted.table}" grid non-blank pixels = ${painted.pixels}`);
    expect(painted.table, 'the V3000 pixel count must be read from the V3000 table, not from a table left open by an earlier step')
      .toBe(v3000.dfName);
    // Not >= 0: this is the step's only reading of the grid, and a blank canvas scores 0, so the
    // "no console.error calls during rendering" claim below would cover a grid that never painted. The
    // floor is waitForGridPaint's own exit condition (pixels > 0), not a tuned threshold.
    expect(painted.pixels, 'the V3000 grid must have painted — a blank canvas means nothing was rendered for the "no console.error calls during rendering" claim to be about').toBeGreaterThan(0);
    // Declared gap: this step makes no structure-versus-raw-text claim about the drawn cell.
    // The colour-spread discriminator the SDF and .mol steps use cannot apply here — v3000.mol
    // is a pure hydrocarbon (M V30 COUNTS 19 19, every atom C) and Chem's atom palette paints
    // carbon and hydrogen black (chem-common-rdkit.ts:34), so a correctly drawn structure scores
    // zero coloured pixels and the reading is a constant rather than a discriminator. The
    // fixture is pinned by the scenario, and every threshold-free substitute (ink extent against
    // this grid's text control, a two-state contrast against a forced text fallback) rests on
    // measurements this run could not take. The step's V3000 claim is carried by the checkSmiles
    // guard above — not by the four assertions preceding it, which the 'MALFORMED_INPUT_VALUE'
    // sentinel satisfies with the molblock unparsed.

    const v3000Errs = await capturedConsoleErrors(page);
    expect(logAndCountConsoleErrors(v3000Errs, 'v3000-render'), `console.error calls on V3000 render — rdkit-matched: ${JSON.stringify(v3000Errs.matched.slice(0, 3))}, other: ${JSON.stringify(v3000Errs.other.slice(0, 3))}`).toBe(0);
    const removed = await page.evaluate((path) =>
      grok.dapi.files.exists(path)
        .then((e: boolean) => (e ? 'present' : 'absent'))
        .catch((err: any) => `probe-failed: ${String(err)}`), SDF_PATH);
    expect(removed, `the temporary V3000 fixture must be removed; got "${removed}"`).toBe('absent');
  });

  await softStep('Scenario 1 Step 6: import a single-structure .mol (aspirin.mol, V2000) — one-row table, Molecule semType, structure renders', async () => {
    const molTable = await openChemFile(page, 'System:DemoFiles/chem/mol/aspirin.mol');
    await waitForMolecule(page);
    const facts = await moleculeColumnFacts(page);
    expect(facts.table, 'the facts below must be read from the .mol table this step opened, not from one left behind by an earlier step')
      .toBe(molTable);
    expect(facts.semType, 'a single .mol file must type its structure column as Molecule').toBe('Molecule');
    expect(facts.cellRenderer, 'Molecule cell renderer must bind on a single-mol import').toBe('Molecule');
    expect(facts.rows, 'a single-structure .mol file produces a one-row table').toBe(1);
    const painted = await waitForGridPaint(page);
    console.log(`[mol-render] measured table="${painted.table}" grid non-blank pixels = ${painted.pixels}`);
    expect(painted.table, 'the .mol pixel count must be read from the .mol table, not from a table left open by an earlier step')
      .toBe(molTable);
    expect(painted.pixels, 'mol grid canvas read failed (fault guard)').toBeGreaterThanOrEqual(0);
    expectStructureInk(await structureInkReadings(page, facts.molCol), 'mol-render');
    const errs = await capturedConsoleErrors(page);
    expect(logAndCountConsoleErrors(errs, 'mol-render'), `console.error calls on .mol render — rdkit-matched: ${JSON.stringify(errs.matched.slice(0, 3))}, other: ${JSON.stringify(errs.other.slice(0, 3))}`).toBe(0);
  });

  await softStep('Scenario 1 Step 8: import a single-molecule .mol2 through its own handler — one-row table, Molecule semType', async () => {
    // The Step 4 anchor names four formats; mol2 had no assertion. adrenalin.mol2 carries
    // exactly one @<TRIPOS>MOLECULE block, so the expected record count is derived from the
    // fixture rather than assumed. Routed through OpenFile so importMol2
    // (Chem/src/package.ts:2069, ext: 'mol2') is the handler under test.
    const MOL2_PATH = 'System:DemoFiles/bio/ngl-formats/adrenalin.mol2';
    const mol2 = await page.evaluate(async (path) => {
      const present = await grok.dapi.files.exists(path)
        .then((e: boolean) => (e ? 'present' : 'absent'))
        .catch((err: any) => `probe-failed: ${String(err)}`);
      if (present !== 'present') return {present, blocks: -1, rows: -1, semType: null};
      const source = await grok.dapi.files.readAsText(path);
      const blocks = (source.match(/@<TRIPOS>MOLECULE/g) || []).length;
      const res = await grok.functions.call('OpenFile', {fullPath: path});
      const df = Array.isArray(res) ? res[0] : res;
      if (!df) return {present, blocks, rows: -1, semType: null};
      grok.shell.closeAll();
      grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 3000));
      const col = df.columns.toList().find((c: any) => c.semType === 'Molecule') ?? null;
      return {present, blocks, rows: df.rowCount, semType: col ? col.semType : null};
    }, MOL2_PATH);
    expect(mol2.present, `the mol2 fixture must be readable on the server; got "${mol2.present}"`).toBe('present');
    expect(mol2.blocks, 'fixture sanity: adrenalin.mol2 must hold exactly one TRIPOS molecule block').toBe(1);
    expect(mol2.rows, 'a one-molecule .mol2 must import as one row').toBe(mol2.blocks);
    expect(mol2.semType, 'the .mol2 structure column must be typed as Molecule').toBe('Molecule');
  });

  await softStep('Scenario 1 Step 7: .smi (SMILES) file — write a temporary fixture, open it through the .smi handler, assert three rows + Molecule semType, then delete it', async () => {
    // No .smi ships with the platform, so the fixture is created here and removed
    // in the same step. It must be opened as a FILE: `.smi` has its own importer
    // (Chem/src/package.ts:2055, ext: 'smi'); asserting on smiles.csv instead
    // would run the ordinary CSV parser and miss any regression in importSmi.
    const SMI_PATH = 'System:AppData/Chem/tmp-import-check.smi';
    await beginErrorCapture(page);
    const facts = await page.evaluate(async (path) => {
      await grok.dapi.files.writeAsText(path, ['CCO', 'c1ccccc1', 'CC(=O)O'].join(String.fromCharCode(10)));
      try {
        const opened = await grok.functions.call('OpenFile', {fullPath: path});
        const df = Array.isArray(opened) ? opened[0] : opened;
        if (!df) return {rows: -1, semType: null, nonEmptyValues: -1, dfName: 'NONE'};
        grok.shell.closeAll();
        grok.shell.addTableView(df);
        await new Promise((r) => setTimeout(r, 4000));
        const t = grok.shell.t;
        const molCol = t.columns.toList().find((c: any) => c.semType === 'Molecule');
        return {
          rows: t.rowCount,
          semType: molCol ? 'Molecule' : null,
          nonEmptyValues: molCol ? molCol.toList().filter((v: any) => typeof v === 'string' && v.length > 0).length : -1,
          dfName: String(t.name),
        };
      } finally {
        await grok.dapi.files.delete(path);
      }
    }, SMI_PATH);

    expect(facts.rows, `.smi import must yield one row per SMILES line (3 written), got ${facts.rows}`).toBe(3);
    expect(facts.semType, '.smi import must produce a column detected as Molecule').toBe('Molecule');
    expect(facts.nonEmptyValues, 'every imported .smi row must carry a non-empty molecule value').toBe(3);

    // The three assertions above are all about values this step itself wrote back out of the
    // dataframe; none of them observes the grid. The reachable rendering triple is added here:
    // the renderer binding, a grid that actually painted, and this step's own error buffer.
    // The colour-spread ink helper the SDF and .mol steps use is NOT usable here, on two
    // independently checked facts: _importSmi builds the table by prepending a single "SMILES"
    // header and parsing the result as CSV (Chem/src/file-importers/smi-importer.ts:4-12), so the
    // table has exactly one column and the helper's control search finds no non-Molecule data
    // column, falling through to the row-number-gutter branch the recon note at the top of this
    // file records as never exercised; and the middle fixture line c1ccccc1 is all carbon, which
    // Chem paints [0,0,0] (chem-common-rdkit.ts:34-37) — the same constant-discriminator trap
    // SR-06 documents. What the triple does and does not cover is declared in SR-08.
    const view = await moleculeColumnFacts(page);
    console.log(`[smi-render] view table="${view.table}" molCol="${view.molCol}" cellRenderer="${view.cellRenderer}"`);
    expect(view.table, 'the renderer binding must be read from the .smi table this step opened, not from one left behind by an earlier step')
      .toBe(facts.dfName);
    expect(view.cellRenderer, 'the .smi Molecule column must have the Molecule cell renderer bound — the semType assertion above is satisfied by a column the grid still draws as raw text')
      .toBe('Molecule');
    const painted = await waitForGridPaint(page);
    console.log(`[smi-render] measured table="${painted.table}" grid non-blank pixels = ${painted.pixels}`);
    expect(painted.table, 'the .smi pixel count must be read from the .smi table, not from a table left open by an earlier step')
      .toBe(facts.dfName);
    // Not >= 0: this is the step's only reading of the grid, so a blank canvas would leave the
    // zero-error claim below with nothing rendered to be about. Floor is waitForGridPaint's own
    // exit condition (pixels > 0), not a tuned threshold.
    expect(painted.pixels, 'the .smi grid must have painted — a blank canvas means nothing was drawn for the "no console.error calls while rendering" claim to be about')
      .toBeGreaterThan(0);
    const smiErrs = await capturedConsoleErrors(page);
    expect(logAndCountConsoleErrors(smiErrs, 'smi-render'), `console.error calls on .smi render — rdkit-matched: ${JSON.stringify(smiErrs.matched.slice(0, 3))}, other: ${JSON.stringify(smiErrs.other.slice(0, 3))}`).toBe(0);

    const gone = await page.evaluate((path) =>
      grok.dapi.files.exists(path)
        .then((e: boolean) => (e ? 'present' : 'absent'))
        .catch((err: any) => `probe-failed: ${String(err)}`), SMI_PATH);
    expect(gone, `the temporary .smi fixture must be deleted — this check leaves nothing behind; got "${gone}"`).toBe('absent');
  });

  await softStep('Scenario 2 Step 1-2: open smiles.csv, invoke Save as SDF exporter — the export dialog opens with a Molecules column selector', async () => {
    const csvTable = await openChemFile(page, 'System:DemoFiles/chem/smiles.csv');
    await waitForMolecule(page);
    const facts = await moleculeColumnFacts(page);
    expect(facts.table, 'the 1000-row count below is this file\'s only nailed comparand, so it must be read from the smiles.csv table this step opened, not from one left behind by an earlier step')
      .toBe(csvTable);
    expect(facts.semType, 'smiles.csv must load with a Molecule column so it can be exported as SDF').toBe('Molecule');
    expect(facts.rows, 'smiles.csv fixture must be the full 1000-row table for the round-trip count').toBe(1000);
    exportedFromRows = facts.rows;
    await page.evaluate(() => grok.functions.call('Chem:saveAsSdf', {}));
    const dlg = page.locator('.d4-dialog:has([name="input-host-Molecules"])');
    await dlg.waitFor({timeout: 20_000, state: 'visible'});
    await expect(dlg.locator('[name="input-host-Molecules"]'), 'Save as SDF dialog must expose the Molecules column selector').toBeVisible();
  });

  await softStep('Scenario 2 Step 3: confirm the export (OK) — the dialog closes, the SDF artefact is written, and no error notification appears', async () => {
    const dlg = page.locator('.d4-dialog:has([name="input-host-Molecules"])');
    // There is no grok.shell.warnings collection — shell.warning() only raises a balloon
    // (js-api/src/shell.ts:176), so the earlier read of that property could never be
    // anything but absent. Warnings are observable as .d4-balloon.warning, and the channel
    // is proved live here before it is used as evidence of absence: raise a known warning,
    // watch it appear, wait it out, and only then assert that the export raises none.
    // The probe proves the recorder as well as the channel: the same recorder that reports
    // "no balloon" after the export must be seen reporting one here first.
    await installBalloonRecorder(page);
    const warnBalloon = page.locator('.d4-balloon.warning');
    await page.evaluate(() => grok.shell.warning('export-probe'));
    await expect(warnBalloon.first(),
      'the warning-balloon channel must be observable, or "no warning was raised" proves nothing')
      .toBeVisible({timeout: 10_000});
    const probed = await readRecordedBalloons(page);
    console.log(`[balloon-probe] recorded=${JSON.stringify(probed)}`);
    expect(probed.filter((b) => /warning/.test(b.cls)).length,
      'the balloon recorder must catch the probe warning, or its silence after the export proves nothing')
      .toBeGreaterThan(0);
    await expect(warnBalloon,
      'the probe warning must clear before the export, so any balloon seen afterwards belongs to the export')
      .toHaveCount(0, {timeout: 30_000});

    // Re-armed here, so everything it holds afterwards belongs to the export itself.
    await installBalloonRecorder(page);

    // The exporter's last act is DG.Utils.download (Chem/src/utils/sdf-utils.ts:154),
    // so the artefact leaves through the browser download channel. Reading it here is
    // what makes the round-trip below a check of the file the product wrote rather
    // than of a second in-memory conversion, which a broken writer would survive.
    const downloadPromise = page.waitForEvent('download', {timeout: 120_000});
    await dlg.locator('[name="button-OK"]').click();
    await dlg.waitFor({timeout: 20_000, state: 'detached'});
    const download = await downloadPromise;
    expect(download.suggestedFilename(), 'Save as SDF must write an .sdf artefact').toMatch(/\.sdf$/i);
    const stream = await download.createReadStream();
    const chunks: Buffer[] = [];
    for await (const c of stream) chunks.push(c as Buffer);
    exportedSdf = Buffer.concat(chunks).toString('utf8');
    const dollar = (exportedSdf.match(/\$\$\$\$/g) || []).length;
    console.log(`[export] ${download.suggestedFilename()} bytes=${exportedSdf.length} records=${dollar}`);
    expect(exportedSdf.length, 'the exported SDF file must not be empty').toBeGreaterThan(0);
    expect(dollar, 'the exported SDF file must carry one "$$$$" record per exported row').toBe(exportedFromRows);

    const raised = await readRecordedBalloons(page);
    const errBalloons = raised.filter((b) => /\berror\b/.test(b.cls));
    const warnBalloons = raised.filter((b) => /\bwarning\b/.test(b.cls));
    console.log(`[export-balloons] recorded=${raised.length} error=${errBalloons.length} warning=${warnBalloons.length}`);
    expect(errBalloons.length,
      `confirming Save as SDF must not raise an error balloon; first: ${JSON.stringify(errBalloons[0] ?? null)}`).toBe(0);
    expect(warnBalloons.length,
      `Save as SDF export must not raise a warning balloon; first: ${JSON.stringify(warnBalloons[0] ?? null)}`).toBe(0);
  });

  await softStep('Scenario 2 Step 4: round-trip — re-importing the exported SDF file reproduces the source row count, and three spot-check rows (first, middle, last) survive as matching canonical SMILES; the 997 rows between them are not read (SR-09)', async () => {
    expect(exportedSdf.length,
      'the SDF artefact written by the export must have been captured in Step 3, or nothing below is a round-trip of the exported file')
      .toBeGreaterThan(0);
    const rt = await page.evaluate(async (sdf) => {
      const df = grok.shell.t;
      const mc = df.columns.toList().find((c: any) => c.semType === 'Molecule');
      const dfs = await grok.functions.call('Chem:importSdf', {bytes: new TextEncoder().encode(sdf)});
      const df2 = Array.isArray(dfs) ? dfs[0] : dfs;
      const mc2 = df2.columns.toList().find((c: any) => c.semType === 'Molecule');
      const canon = (s: string) => grok.chem.convert(s, DG.chem.Notation.Unknown, DG.chem.Notation.Smiles);
      // Both sides of the equality below come from `canon`, so a uniformly degenerate
      // conversion channel would make them equal without any structure being compared.
      // checkSmiles runs through a different platform function (Chem:isSmiles, RDKit
      // get_mol) than convert (Chem:convertMolNotation), so it cannot degenerate with it;
      // the length test covers RDKit accepting '' as an empty molecule.
      const isStructure = (s: any) => {
        try { return typeof s === 'string' && s.length > 0 && !!grok.chem.checkSmiles(s); }
        catch (e) { return false; }
      };
      const spot: Array<{orig: string; rt: string; origIsStructure: boolean}> = [];
      for (const i of [0, Math.floor(df.rowCount / 2), df.rowCount - 1]) {
        const orig = canon(mc.get(i));
        spot.push({orig, rt: canon(mc2.get(i)), origIsStructure: isStructure(orig)});
      }
      return {srcRows: df.rowCount, reimportRows: df2.rowCount, mc2: mc2 ? mc2.semType : 'NONE', spot};
    }, exportedSdf);
    expect(rt.reimportRows, 're-importing the exported SDF file must reproduce the source row count').toBe(rt.srcRows);
    expect(rt.mc2, 're-imported SDF must carry a Molecule column').toBe('Molecule');
    for (const p of rt.spot) {
      expect(p.origIsStructure,
        `the source-side canonical SMILES must itself be a parseable structure, or the equality below compares two degenerate values instead of molecules; got "${p.orig}"`)
        .toBe(true);
      expect(p.rt, `round-trip must preserve molecule identity: canonical SMILES "${p.orig}" must survive export-to-file + reimport`).toBe(p.orig);
    }
  });

  await softStep('Scenario 3 Step 2 + Step 3: InChI-to-Molecule converter — its output parameter is typed semType=Molecule, each output is a parseable SMILES differing from the raw InChI input, malformed InChI fails rather than yielding a structure, and the grid paints with no console.error calls', async () => {
    await beginErrorCapture(page);
    const res = await page.evaluate(async () => {
      grok.shell.closeAll();
      const inchis = [
        'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
        'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H',
      ];
      const col = DG.Column.fromStrings('inchi', inchis);
      const df = DG.DataFrame.fromColumns([col]);
      df.name = 'inchi_test';
      grok.shell.addTableView(df);
      const conv = DG.Func.find({name: 'InchiToMol'});
      const smiles: string[] = [];
      // The semantic type is read off the converter's own call, not stamped here: the
      // platform builds the output parameter from the script's header declaration
      // `#output: string smiles { semType: Molecule }` (Chem/scripts/inchiToMol.py:7),
      // so a dropped annotation — or a platform that stops honouring one — changes
      // what this reads. Stamping it on the column, as this step did before
      // 2026-08-20, asserted the test's own write.
      const outSemTypes: string[] = [];
      const outParams: string[] = [];
      for (const v of inchis) {
        const call = await conv[0].prepare({id: v}).call();
        smiles.push(call.getOutputParamValue());
        const names = Object.keys(call.outputParams);
        outParams.push(names.join(',') || 'NONE');
        outSemTypes.push(names.length === 1 ? String(call.outputParams[names[0]].property.semType) : 'NONE');
      }
      const molCol = DG.Column.fromStrings('molecule', smiles);
      molCol.semType = outSemTypes[0];
      df.columns.add(molCol);
      await new Promise<void>((res) => {
        const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
        setTimeout(res, 5000);
      });
      const parseable = smiles.map((s) => {
        try { return !!grok.chem.checkSmiles(s); } catch (e) { return false; }
      });

      let malformedThrew = false;
      let malformedValue: string | null = null;
      let malformedError = '';
      try { malformedValue = await conv[0].apply({id: 'InChI=NONSENSE-BROKEN'}); }
      catch (e) { malformedThrew = true; malformedError = String(e); }
      return {inchis, smiles, parseable, malformedThrew, malformedValue, malformedError, outSemTypes, outParams};
    });
    expect(res.smiles.length, 'the converter must produce one output per valid InChI row').toBe(res.inchis.length);
    for (let i = 0; i < res.smiles.length; i++) {
      expect(res.parseable[i], `InchiToMol output for row ${i} must be a valid, parseable structure (got "${res.smiles[i]}")`).toBe(true);
      expect(res.smiles[i], `the converter must produce a structure DIFFERENT from the raw InChI input for row ${i} — a same-string echo means no conversion happened`).not.toBe(res.inchis[i]);
    }
    for (let i = 0; i < res.outSemTypes.length; i++)
      expect(res.outSemTypes[i],
        `the converter's output parameter must be typed semType=Molecule so its result renders as a structure — read from the call's own output metadata for row ${i} (output params: "${res.outParams[i]}")`)
        .toBe('Molecule');
    expect(res.malformedThrew, `malformed InChI must NOT silently yield a valid structure — the parse failure must surface (got value "${res.malformedValue}" instead of an error)`).toBe(true);
    // The valid conversions above are the live control, but they go through prepare().call()
    // while this goes through apply(), and infrastructure could still fall over in between —
    // either way the throw would be about reaching the converter, not about the input. No
    // exact vendor message is pinned (it would rot); the unreachability vocabulary is excluded
    // instead, so the throw is only accepted as a rejection of the input.
    expect(res.malformedError.length, 'the malformed-InChI failure must carry an error message').toBeGreaterThan(0);
    expect(res.malformedError,
      `malformed InChI must fail because the converter rejected the input, not because it could not be reached; got "${res.malformedError}"`)
      .not.toMatch(/failed to fetch|networkerror|econnrefused|socket hang up|bad gateway|service unavailable|gateway time-?out|\b50[234]\b|unauthorized|forbidden|not found|is not a function|cannot read propert/i);
    // The zero-error claim below is about rendering the converted Molecule column, so the render
    // has to have happened first: read straight after the conversions, the buffer is empty
    // because nothing has been drawn yet, not because the renderer stayed quiet.
    const painted = await waitForGridPaint(page);
    console.log(`[inchi-render] measured table="${painted.table}" grid non-blank pixels = ${painted.pixels}`);
    expect(painted.table, 'the render barrier must settle on the InChI table, not on one left open by an earlier step').toBe('inchi_test');
    // Same as the V3000 step: the only reading of this grid, so >= 0 would let the zero-error
    // claim below stand over a canvas that never painted. Floor is waitForGridPaint's exit condition.
    expect(painted.pixels, 'the InChI grid must have painted — a blank canvas means nothing was rendered for the "no console.error calls while rendering" claim to be about').toBeGreaterThan(0);
    const errs = await capturedConsoleErrors(page);
    expect(logAndCountConsoleErrors(errs, 'inchi-render'), `console.error calls while rendering the converted Molecule column — rdkit-matched: ${JSON.stringify(errs.matched.slice(0, 3))}, other: ${JSON.stringify(errs.other.slice(0, 3))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});
