/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForChemMenu, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';

declare const grok: any;
declare const DG: any;

async function openChemTable(page: Page, path: string): Promise<void> {
  await page.evaluate(async (p) => {
    document.body.classList.add('selenium');
    try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
    try { grok.shell.windows.simpleMode = true; } catch (e) {}
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(p);
    grok.shell.addTableView(df);
    (window as any).__df = df;
    (window as any).__errs = [];
    const orig = console.error;
    console.error = function (...a: any[]) {
      (window as any).__errs.push(a.map((x: any) => String(x)).join(' '));
      orig.apply(console, a as any);
    };
    await new Promise<void>((res) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); res(); });
      setTimeout(res, 8000);
    });
  }, path);
  await waitForChemMenu(page);
  for (let i = 0; i < 50; i++) {
    if (await page.locator('[name="viewer-Grid"] canvas').count() > 0) break;
    await page.waitForTimeout(200);
  }
  await page.locator('[name="viewer-Grid"] canvas').first().waitFor({timeout: 30_000, state: 'attached'});
  await page.waitForTimeout(4000);
}

// Separates a painted canvas from a blank or unreadable one, and nothing finer:
// a raw-text fallback paints plenty of dark pixels too. Callers must phrase the
// message as a paint/fault guard — never as "structures, not text".
// A pixel counts only when it is opaque: an untouched canvas reads (0,0,0,0),
// whose three colour channels satisfy the darkness test, so without the alpha
// term a never-drawn canvas would score its full area.
// With a column name the read is narrowed to that grid column's x-range, so a
// blank column beside painted neighbours no longer satisfies it; without one it
// reads the leading 600x400 of the canvas.
async function gridNonBlankPixels(page: Page, columnName?: string): Promise<number> {
  return await page.evaluate((name) => {
    // Categorical filter cards render mini-grids that carry [name="viewer-Grid"] too,
    // so a first match can be one of those while the x-range comes from
    // grok.shell.tv.grid — an arbitrary strip that still clears the pixel floor,
    // i.e. a silent green rather than a red.
    const grid = Array.from(document.querySelectorAll('[name="viewer-Grid"]'))
      .find((g) => !g.closest('.d4-filter')) ?? null;
    if (!grid) return -1;
    const cs = Array.from(grid.querySelectorAll('canvas'))
      .sort((a: any, b: any) => b.width * b.height - a.width * a.height) as HTMLCanvasElement[];
    if (cs.length === 0) return -1;
    const c = cs[0];
    const ctx = c.getContext('2d');
    if (!ctx || !c.width || !c.height) return -1;
    let x0 = 0, x1 = Math.min(c.width, 600);
    if (name) {
      let gc: any = null;
      try { gc = grok.shell.tv.grid.columns.byName(name); } catch (e) { return -1; }
      if (!gc) return -1;
      const scale = c.width / (c.clientWidth || c.width);
      x0 = Math.max(0, Math.round(gc.left * scale));
      x1 = Math.min(c.width, Math.round(gc.right * scale));
    }
    const h = Math.min(c.height, 400);
    if (x1 - x0 <= 0 || h <= 0) return -1;
    const d = ctx.getImageData(x0, 0, x1 - x0, h).data;
    let cnt = 0;
    for (let i = 0; i < d.length; i += 4)
      if (d[i + 3] > 0 && (d[i] < 250 || d[i + 1] < 250 || d[i + 2] < 250)) cnt++;
    return cnt;
  }, columnName);
}

// The structure-versus-text discriminator. RDKit's default palette draws
// heteroatoms in saturated colour (O ~ rgb(255,0,0), N ~ rgb(0,0,255), channel
// spread ~255) while canvas text is drawn with one grey/black fill, whose three
// channels are equal by construction (spread 0). The 30/255 margin therefore
// sits far above antialiasing noise and far below real heteroatom ink, and the
// caller asserts only "> 0", so no tuned pixel count is involved. Returns -1
// when the column region cannot be read, so the caller can reject the read
// instead of comparing against a sentinel.
// The strip starts below the column-header band (grid.colHeaderHeight) rather
// than at y = 0, so the header and the filter chrome openChemTable pins there
// are outside the measured area altogether. Controlling for that chrome is not
// enough on its own: a negative control on another column only covers chrome
// the two columns share, while ink belonging to the Molecule column's own
// header would satisfy "coloured > 0" with no molecule drawn at all.
// Reports the header height and the strip it settled on so a red run does not
// have to guess which region was read.
interface ColumnColourReading {
  coloured: number;
  headerHeight: number;
  x0: number;
  x1: number;
  y0: number;
  h: number;
}

function describeReading(r: ColumnColourReading): string {
  return `header=${r.headerHeight} strip x=[${r.x0},${r.x1}) y=[${r.y0},${r.y0 + r.h}) coloured=${r.coloured}`;
}

async function moleculeColumnColouredPixels(page: Page, columnName: string): Promise<ColumnColourReading> {
  return await page.evaluate((name) => {
    const fail = {coloured: -1, headerHeight: -1, x0: -1, x1: -1, y0: -1, h: -1};
    // Filter mini-grids carry [name="viewer-Grid"] too — see gridNonBlankPixels.
    const grid = Array.from(document.querySelectorAll('[name="viewer-Grid"]'))
      .find((g) => !g.closest('.d4-filter')) ?? null;
    if (!grid) return fail;
    const cs = Array.from(grid.querySelectorAll('canvas'))
      .sort((a: any, b: any) => b.width * b.height - a.width * a.height) as HTMLCanvasElement[];
    if (cs.length === 0) return fail;
    const c = cs[0];
    const ctx = c.getContext('2d');
    if (!ctx || !c.width || !c.height) return fail;
    let gc: any = null;
    let headerHeight = -1;
    try {
      gc = grok.shell.tv.grid.columns.byName(name);
      headerHeight = grok.shell.tv.grid.colHeaderHeight;
    }
    catch (e) { return fail; }
    if (!gc || !(headerHeight >= 0)) return fail;
    const scale = c.width / (c.clientWidth || c.width);
    const x0 = Math.max(0, Math.round(gc.left * scale));
    const x1 = Math.min(c.width, Math.round(gc.right * scale));
    const y0 = Math.min(c.height, Math.round(headerHeight * scale));
    const h = Math.min(c.height - y0, Math.round(400 * scale));
    if (x1 - x0 <= 0 || h <= 0) return {coloured: -1, headerHeight, x0, x1, y0, h};
    const d = ctx.getImageData(x0, y0, x1 - x0, h).data;
    let coloured = 0;
    for (let i = 0; i < d.length; i += 4)
      if (Math.max(d[i], d[i + 1], d[i + 2]) - Math.min(d[i], d[i + 1], d[i + 2]) >= 30) coloured++;
    return {coloured, headerHeight, x0, x1, y0, h};
  }, columnName);
}

function rdkitErrors(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const errs = ((window as any).__errs as string[] | undefined) ?? [];
    // 'gS' is the minified Dart symbol from the GROK-16870 crash and is only matched inside its
    // quotes — bare, it matches the "gs" in settings/logs/flags and turns any error into a hit.
    return errs.filter((e) => /rdkit[-_]cell[-_]renderer|method not found|'gS'|cellRenderer\.render|NullError/i.test(e));
  });
}

test.use(specTestOptions);

test('Chem: Molecule rendering end-to-end (RDKit cell / cross-viewer tooltip / reaction + mixture renderers)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await softStep('Scenario 1 Step 1-2: open smiles.csv, canonical_smiles detected as Molecule', async () => {
    await openChemTable(page, 'System:DemoFiles/chem/smiles.csv');
    await waitForMolecule(page);
    const detected = await page.evaluate(() => {
      const df = (window as any).__df;
      const c = df.columns.byName('canonical_smiles');
      return {semType: c.semType, cellRenderer: c.meta.cellRenderer, rows: df.rowCount};
    });
    expect(detected.semType, 'canonical_smiles must detect as Molecule so the RDKit cell renderer binds').toBe('Molecule');
    expect(detected.cellRenderer, 'Molecule cell renderer must be bound to the column').toBe('Molecule');
  });

  await softStep('Scenario 1 Step 3: Molecule cells render as structures (coloured heteroatom ink), no renderer errors', async () => {
    const nonBlank = await gridNonBlankPixels(page);
    console.log(`[molecule-render] grid non-blank pixels = ${nonBlank}`);
    expect(nonBlank, 'grid canvas read failed (fault guard)').toBeGreaterThanOrEqual(0);
    // 2000 is a fault floor, not a discriminator. The only measured non-blank
    // count on record for a painted Chem grid is 27274 (single-cell V3000 grid,
    // chem-import-export-formats gate_verdicts.b diagnostic, 2026-08-20), an
    // order of magnitude above it, so the floor separates "canvas painted" from
    // "canvas blank or unreadable" and nothing else.
    expect(nonBlank, 'the grid canvas must be painted at all — a count at the floor means nothing was drawn (this says nothing about structures vs text)').toBeGreaterThan(2000);
    const coloured = await moleculeColumnColouredPixels(page, 'canonical_smiles');
    console.log(`[molecule-render] canonical_smiles ${describeReading(coloured)}`);
    expect(coloured.coloured, 'the canonical_smiles column region could not be read from the grid canvas (fault guard)').toBeGreaterThanOrEqual(0);
    expect(coloured.h, 'the canonical_smiles strip below the column header must have height — a zero-height read scores nothing whatever the renderer drew').toBeGreaterThan(0);
    expect(coloured.coloured, 'the Molecule column must be drawn as structures: RDKit paints heteroatoms in colour (O red, N blue) while a raw-SMILES text fallback is monochrome, so zero coloured pixels in the column region means no structure was drawn').toBeGreaterThan(0);
    // Negative control, measured on the same grid in the same run: molregno is a
    // plain int column, so nothing inside its strip is drawn in colour. Both
    // strips start below the column header, so this covers what remains shared
    // inside the data rows — row tint, selection and current-row highlight. A
    // non-zero reading here means the counter is scoring that rather than
    // heteroatom ink and the assertion above proves nothing.
    const control = await moleculeColumnColouredPixels(page, 'molregno');
    console.log(`[molecule-render] molregno control ${describeReading(control)}`);
    expect(control.coloured, 'the molregno control column could not be read from the grid canvas (fault guard)').toBeGreaterThanOrEqual(0);
    expect(control.h, 'the molregno control strip must have height — a zero-height control scores zero for free and controls for nothing').toBeGreaterThan(0);
    expect(control.coloured, 'the monochrome molregno column must score zero coloured pixels — a non-zero count means the counter is picking up grid chrome (row tint, selection highlight) instead of structure ink, which would let the Molecule assertion pass over raw SMILES text').toBe(0);
    const errs = await rdkitErrors(page);
    expect(errs.length, `rdkit-cell-renderer errors during Molecule render: ${JSON.stringify(errs.slice(0, 3))}`).toBe(0);
  });

  await softStep('Scenario 2 Step 1-2: add Box Plot and Scatter Plot over numeric axes', async () => {
    await page.evaluate(() => {
      const tv = grok.shell.tv;
      tv.addViewer('Box plot');
      tv.addViewer('Scatter plot');
    });
    await page.locator('[name="viewer-Box-plot"]').waitFor({timeout: 30_000, state: 'visible'});
    await page.locator('[name="viewer-Scatter-plot"]').waitFor({timeout: 30_000, state: 'visible'});
  });

  await softStep('Scenario 2 Step 3: hover Box Plot — molecule tooltip renders, no rdkit-cell-renderer NullError (GROK-16870)', async () => {
    const tip = await hoverSweep(page, '[name="viewer-Box-plot"]');
    console.log(`[tooltip] Box Plot seen=${tip.seen} canvasPixels=${tip.canvasPixels}`);
    expectTooltipPainted(tip, 'Box Plot');
    const errs = await rdkitErrors(page);
    expect(errs.length, `GROK-16870 regression on Box Plot tooltip: ${JSON.stringify(errs.slice(0, 3))}`).toBe(0);
    const alive = await viewerAlive(page, '[name="viewer-Box-plot"]');
    expect(alive, 'Box Plot unmounted/zero-size after hover — renderer crash poisoned layout').toBe(true);
  });

  await softStep('Scenario 2 Step 4: hover Scatter Plot — molecule tooltip renders, no rdkit errors', async () => {
    const tip = await hoverSweep(page, '[name="viewer-Scatter-plot"]');
    console.log(`[tooltip] Scatter Plot seen=${tip.seen} canvasPixels=${tip.canvasPixels}`);
    expectTooltipPainted(tip, 'Scatter Plot');
    const errs = await rdkitErrors(page);
    expect(errs.length, `Scatter Plot tooltip rdkit errors: ${JSON.stringify(errs.slice(0, 3))}`).toBe(0);
    const alive = await viewerAlive(page, '[name="viewer-Scatter-plot"]');
    expect(alive, 'Scatter Plot unmounted/zero-size after hover — renderer crash').toBe(true);
  });


  await softStep('Scenario 3 Step 3: reaction column detected as ChemicalReaction, reaction cells render, no errors', async () => {
    await openChemTable(page, 'System:AppData/Chem/test-reactions.csv');
    const detected = await page.evaluate(() => {
      const df = (window as any).__df;
      const c = df.columns.byName('reaction');
      let arrowRows = 0;
      for (let i = 0; i < df.rowCount; i++) if (String(df.get('reaction', i)).includes('>>')) arrowRows++;
      return {semType: c.semType, cellRenderer: c.meta.cellRenderer, rows: df.rowCount, arrowRows};
    });
    expect(detected.semType, 'reaction column must detect as ChemicalReaction so the reaction cell renderer binds').toBe('ChemicalReaction');
    expect(detected.cellRenderer, 'ChemicalReaction cell renderer must be bound').toBe('ChemicalReaction');
    expect(detected.arrowRows, 'the reaction column must be genuine reaction data — most rows carry a reactant>>product arrow (lower bound proves the fixture is reactions, not mislabeled molecules; strict equality would assert fixture uniformity, which the renderer does not depend on)').toBeGreaterThan(detected.rows / 2);
    expect(detected.arrowRows, 'at least 5 genuine reaction rows so the scenario "view at least 5 reaction cells" step is realizable').toBeGreaterThanOrEqual(5);
    const nonBlank = await gridNonBlankPixels(page, 'reaction');
    console.log(`[reaction-render] reaction column non-blank pixels = ${nonBlank}`);
    expect(nonBlank, 'the reaction column region could not be read from the grid canvas (fault guard)').toBeGreaterThanOrEqual(0);
    // 2000 is a fault floor for the column strip, not a discriminator, and it is
    // placed from measurement (standalone node+chromium probe against dev,
    // 2026-08-22): the painted reaction strip reads 6324, the same-width strip of
    // canvas below the last data row reads 400, and the column-header band with
    // its pinned filter chrome contributes at most 937 — so no combination of
    // chrome and empty grid lines reaches the floor.
    expect(nonBlank, 'the reaction column region must be painted — a count at the floor means the reaction renderer drew nothing in its own cells (this says nothing about arrow layout or text fallback)').toBeGreaterThan(2000);
    const errs = await rdkitErrors(page);
    expect(errs.length, `reaction renderer errors: ${JSON.stringify(errs.slice(0, 3))}`).toBe(0);
  });

  await softStep('Scenario 3 Step 6: mixture column detected as ChemicalMixture, multi-component cells render, no errors', async () => {
    await openChemTable(page, 'System:AppData/Chem/test_mixtures.csv');
    const detected = await page.evaluate(() => {
      const df = (window as any).__df;
      const c = df.columns.byName('mixture');
      const compCounts: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        try { compCounts.push((JSON.parse(String(df.get('mixture', i))).contents || []).length); }
        catch (e) { compCounts.push(-1); }
      }
      return {semType: c.semType, cellRenderer: c.meta.cellRenderer, rows: df.rowCount, compCounts};
    });
    expect(detected.semType, 'mixture column must detect as ChemicalMixture so the mixture cell renderer binds').toBe('ChemicalMixture');
    expect(detected.cellRenderer, 'ChemicalMixture cell renderer must be bound').toBe('ChemicalMixture');
    expect(Math.max(...detected.compCounts), 'a genuine mixture cell must carry more than one component').toBeGreaterThan(1);
    expect(Math.min(...detected.compCounts), 'no mixture cell may fail to parse its components').toBeGreaterThanOrEqual(1);
    const nonBlank = await gridNonBlankPixels(page, 'mixture');
    console.log(`[mixture-render] mixture column non-blank pixels = ${nonBlank}`);
    expect(nonBlank, 'the mixture column region could not be read from the grid canvas (fault guard)').toBeGreaterThanOrEqual(0);
    // Same fault floor, same probe: the painted mixture strip reads 3963 against
    // 400 for the identically sized blank strip below the data rows and 522 for
    // its header band.
    expect(nonBlank, 'the mixture column region must be painted — a count at the floor means the mixture renderer drew nothing in its own cells (this says nothing about multi-fragment layout or text fallback)').toBeGreaterThan(2000);
    const errs = await rdkitErrors(page);
    expect(errs.length, `mixture renderer errors: ${JSON.stringify(errs.slice(0, 3))}`).toBe(0);
  });

  await page.evaluate(() => grok.shell.closeAll());

  finishSpec();
});

interface TooltipObservation { seen: boolean; canvasPixels: number; }

// `canvasPixels` is -1 while no readable canvas has been found inside a visible
// tooltip, so "the tooltip painted nothing" and "the tooltip was never read" are
// distinguishable at the assertion instead of collapsing onto 0.
// The pixel test requires opacity as well as darkness. A tooltip canvas the
// renderer never touched reads (0,0,0,0) everywhere, which passes the darkness
// test on all three channels, so a colour-only rule scores a blank transparent
// canvas at its full area and cannot see a renderer that produced nothing.
async function probeTooltip(page: Page): Promise<TooltipObservation> {
  return await page.evaluate(() => {
    const tips = Array.from(document.querySelectorAll('.d4-tooltip'))
      .filter((t) => {
        const r = (t as HTMLElement).getBoundingClientRect();
        return r.width * r.height > 0;
      });
    if (tips.length === 0) return {seen: false, canvasPixels: -1};
    const cv = tips
      .map((t) => t.querySelector('canvas') as HTMLCanvasElement | null)
      .find((c) => !!c && !!c.width && !!c.height) ?? null;
    if (!cv) return {seen: true, canvasPixels: -1};
    const ctx = cv.getContext('2d');
    if (!ctx) return {seen: true, canvasPixels: -1};
    const d = ctx.getImageData(0, 0, cv.width, cv.height).data;
    let painted = 0;
    for (let i = 0; i < d.length; i += 4)
      if (d[i + 3] > 0 && (d[i] < 250 || d[i + 1] < 250 || d[i + 2] < 250)) painted++;
    return {seen: true, canvasPixels: painted};
  });
}

// The painted-pixel bar stays at "> 0" and gets no floor. With opacity required,
// both ways of producing nothing score exactly zero — an untouched canvas is
// transparent, and a canvas the renderer only background-filled is opaque white
// — so zero already separates "drew something" from "drew nothing", and a floor
// above it would only be a threshold tuned to whichever molecules the fixture
// happens to hold.
function expectTooltipPainted(tip: TooltipObservation, viewer: string): void {
  expect(tip.seen, `hovering the ${viewer} must raise a tooltip — with no tooltip at all, the console and viewer-alive checks below pass whether or not the molecule cell renderer ever ran`).toBe(true);
  expect(tip.canvasPixels, `the ${viewer} tooltip carried no readable canvas, so its molecule region could not be inspected`).toBeGreaterThanOrEqual(0);
  expect(tip.canvasPixels, `the ${viewer} tooltip's molecule region must be painted — a blank tooltip canvas means the molecule cell renderer produced nothing`).toBeGreaterThan(0);
}

// Hover positions taken from where the viewer's own canvas carries ink, densest
// cell first. Fractions of the viewer's bounding box are not a substitute: the box
// includes the header and the axis gutters, so a viewer whose ink sits in one
// narrow band would be missed entirely.
async function inkHoverTargets(page: Page, viewerSelector: string): Promise<{x: number; y: number}[]> {
  return await page.evaluate((sel) => {
    const v = document.querySelector(sel) as HTMLElement | null;
    if (!v) return [];
    const cs = Array.from(v.querySelectorAll('canvas'))
      .sort((a: any, b: any) => b.width * b.height - a.width * a.height) as HTMLCanvasElement[];
    const c = cs[0];
    if (!c || !c.width || !c.height) return [];
    const ctx = c.getContext('2d');
    if (!ctx) return [];
    const d = ctx.getImageData(0, 0, c.width, c.height).data;
    const r = c.getBoundingClientRect();
    const n = 16;
    const cells: {x: number; y: number; ink: number}[] = [];
    for (let cy = 0; cy < n; cy++) {
      for (let cx = 0; cx < n; cx++) {
        const x0 = Math.floor(c.width * cx / n), x1 = Math.floor(c.width * (cx + 1) / n);
        const y0 = Math.floor(c.height * cy / n), y1 = Math.floor(c.height * (cy + 1) / n);
        let ink = 0;
        for (let y = y0; y < y1; y += 2)
          for (let x = x0; x < x1; x += 2) {
            const i = (y * c.width + x) * 4;
            if (d[i + 3] > 0 && (d[i] < 250 || d[i + 1] < 250 || d[i + 2] < 250)) ink++;
          }
        if (ink > 0)
          cells.push({x: r.x + r.width * (cx + 0.5) / n, y: r.y + r.height * (cy + 0.5) / n, ink});
      }
    }
    cells.sort((a, b) => b.ink - a.ink);
    return cells.slice(0, 16).map((t) => ({x: t.x, y: t.y}));
  }, viewerSelector);
}

// Sweeps the pointer over the viewer until a tooltip with a painted canvas is
// observed, and reports the best observation. Stops early on success, so the
// full sweep cost is paid only when nothing ever appears.
async function hoverSweep(page: Page, viewerSelector: string): Promise<TooltipObservation> {
  const targets = await inkHoverTargets(page, viewerSelector);
  expect(targets.length,
    `viewer ${viewerSelector} has no readable canvas ink to hover — nothing was drawn, so no hover position could raise a tooltip`)
    .toBeGreaterThan(0);
  const best: TooltipObservation = {seen: false, canvasPixels: -1};
  for (const t of targets) {
    await page.mouse.move(t.x, t.y, {steps: 5});
    for (let i = 0; i < 4; i++) {
      await page.waitForTimeout(250);
      const p = await probeTooltip(page);
      if (p.seen) best.seen = true;
      if (p.canvasPixels > best.canvasPixels) best.canvasPixels = p.canvasPixels;
      if (best.canvasPixels > 0) return best;
    }
  }
  return best;
}

async function viewerAlive(page: Page, viewerSelector: string): Promise<boolean> {
  return await page.evaluate((sel) => {
    const v = document.querySelector(sel) as HTMLElement | null;
    if (!v) return false;
    const r = v.getBoundingClientRect();
    return r.width * r.height > 0 && v.querySelectorAll('canvas').length >= 1;
  }, viewerSelector);
}
