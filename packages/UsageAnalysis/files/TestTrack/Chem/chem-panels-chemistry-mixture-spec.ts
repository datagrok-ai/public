/* ---
realizes: [chem.cp.panels-chemistry-mixture]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, waitForMolecule} from '../spec-login';
import {finishSpec} from '../helpers/viewers';
import {withConsoleErrorCount} from '../helpers/forms';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

// Fixture-path reconciliation (the scenario's Setup names DemoFiles paths for both):
//   test_mixtures.csv ships in Chem's package files (public/packages/Chem/files/test_mixtures.csv),
//   which the server exposes under System:AppData/Chem/ — there is no copy under DemoFiles/chem.
//   mol1K.sdf is the file Scenario 2 Step 4 names, and it likewise lives under System:AppData/Chem/;
//   chem-import-export-formats-spec.ts#L181 records System:DemoFiles/chem/mol1K.sdf as absent, and
//   activity-cliffs.md#L54 / calculate-ui.md#L62 both use the AppData path as their molV2000 fixture.
//   The previous constant here pointed at System:DemoFiles/chem/sdf/aspirin-2d.sdf, a path attested
//   nowhere else in this corpus and not the file the scenario names.
const MIXTURES_PATH = 'System:AppData/Chem/test_mixtures.csv';
const SMILES_PATH = 'System:DemoFiles/chem/smiles.csv';
const MOLBLOCK_SDF_PATH = 'System:AppData/Chem/mol1K.sdf';

// Ambient console noise is excluded BY NAME, not by count — the same list the sibling Chem specs
// carry (chem-calculate-properties-spec.ts#L57-58).
const AMBIENT_CONSOLE = [/favicon/i, /ResizeObserver loop/i, /Permissions policy violation/i,
  /Unable to find element in cloned iframe/i];
const nonAmbient = (texts: string[]) => texts.filter((t) => !AMBIENT_CONSOLE.some((re) => re.test(t)));

// Positive control for the console channel, same shape as chem-filters-spec.ts#L1005-1010: an
// empty capture is what a DETACHED listener produces too, so the channel is proved live by
// making it receive a deliberately emitted error before any step relies on its silence.
const CONSOLE_PROBE = 'chem-panels-chemistry-mixture spec console channel probe';
// Collection window: withConsoleErrorCount only listens while its callback runs plus the settle
// milliseconds passed to it (helpers/forms.ts#L64) — 4000 ms here. Errors arriving later are
// outside every console claim in this spec.
const CONSOLE_SETTLE_MS = 4000;

// Selector recon-notes (class-2):
//   .grok-prop-panel .d4-accordion-pane-header — context (property) panel accordion headers.
//   Chem panels are nested under a "Chemistry" GROUP pane; a collapsed group renders no child
//   headers, so ["Actions","Chemistry"] is the correct state before expansion, not an absent
//   panel. Bind via DG.SemanticValue.fromTableCell(...) — see chem.md.
const CHEM_PANEL_FUNCS = ['mixtureWidget', 'mixtureTreeWidget', 'ChemistryGasteigerPartialCharges'];

// The skip gate asks ONLY whether Chem/RDKit is usable on this build. It must not ask whether
// the mixture panel functions are registered: a Chemistry-group leaf is offered if and only if
// its panel function is registered, so gating on the registrations turned the exact failing
// state of `headers.includes('Mixture')` into a silent skip labelled environmental. The
// registrations are an ASSERTION below, not a precondition.
async function chemRdKitReady(page: Page): Promise<{ok: boolean; detail: string}> {
  return page.evaluate(async () => {
    let detail = '(no probe completed)';
    for (let i = 0; i < 120; i++) {
      let chemFuncs = 0;
      try { chemFuncs = DG.Func.find({package: 'Chem'}).length; }
      catch (e) { chemFuncs = 0; }
      let rdkit = false;
      let err = '';
      try { rdkit = !!(await grok.functions.call('Chem:getRdKitModule')); }
      catch (e: any) { err = e && e.message ? e.message : String(e); }
      detail = `Chem functions registered=${chemFuncs}; Chem:getRdKitModule resolved=${rdkit}` +
        `${err ? ` (threw: ${err})` : ''}`;
      if (chemFuncs > 0 && rdkit) return {ok: true, detail};
      await new Promise((r) => setTimeout(r, 1000));
    }
    return {ok: false, detail};
  });
}

async function missingPanelFuncs(page: Page): Promise<string[]> {
  return page.evaluate((names) => names.filter((n: string) => {
    try { return DG.Func.find({package: 'Chem', name: n}).length === 0; }
    catch (e) { return true; }
  }), CHEM_PANEL_FUNCS);
}

async function paneHeaders(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
      .map((h) => (h as HTMLElement).textContent!.trim()));
}

// Bind to a REAL table cell (Bio/bio-cell-actions-panels-spec.ts#L168, Gate B PASS 3/3).
// A SemanticValue built from a bare string carries no owning cell, and the Chem panels
// resolve from the current object's cell/column context.
async function openContextPanes(page: Page, semType: string,
  requiredPane?: string): Promise<string[]> {
  // Bind, and wait only for the GROUP — a leaf cannot appear while its group is collapsed,
  // so requiring the leaf here would wait out the timeout for something only the expansion
  // below can produce.
  await page.waitForFunction((st) => {
    grok.shell.windows.showContextPanel = true;
    const o = grok.shell.o;
    if (!o || o.semType !== st) {
      grok.functions.call('Chem:getRdKitModule').catch(() => {});
      const df = grok.shell.t;
      const col = df.columns.toList().find((c: any) => c.semType === st);
      if (!col) return false;
      grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, col.name));
      return false;
    }
    return Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
      .some((h) => /^Chemistry$/i.test((h as HTMLElement).textContent!.trim()));
  }, semType, {timeout: 120_000, polling: 1000});

  await expandPane(page, '^Chemistry$');
  if (requiredPane) {
    await page.waitForFunction((req) => Array.from(
      document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
      .some((h) => (h as HTMLElement).textContent!.trim() === req),
    requiredPane, {timeout: 60_000, polling: 500}).catch(() => {});
  }
  const headers = await paneHeaders(page);
  console.log(`[panels] semType=${semType} required=${requiredPane ?? '-'} headers=${JSON.stringify(headers)}`);
  return headers;
}

// A missing header is a failure, not a skip: the previous `if (await header.count())` made
// expandPane a silent no-op whenever the pane it was asked to open did not exist, so every
// downstream assertion measured an unexpanded panel and could still hold.
async function expandPane(page: Page, reSource: string): Promise<void> {
  const re = new RegExp(reSource, 'i');
  const header = page.locator('.grok-prop-panel .d4-accordion-pane-header').filter({hasText: re}).first();
  const count = await header.count();
  const headers = await paneHeaders(page);
  console.log(`[panels] expandPane(${reSource}) headerMatches=${count} headers=${JSON.stringify(headers)}`);
  expect(count,
    `no context-panel header matches /${reSource}/i — expandPane cannot open a pane that is not ` +
    `offered. headers=${JSON.stringify(headers)}`).toBeGreaterThan(0);
  const expanded = await header.evaluate((h) => h.classList.contains('expanded'));
  if (!expanded) await header.click();
  // The settle barrier stays advisory — Gate B green runs rely on it only to reduce flake, and the
  // real claim is made by the caller's assertions. What must not be advisory is the header's
  // existence, asserted above.
  await page.waitForFunction((src) => {
    const rx = new RegExp(src, 'i');
    const h = Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
      .find((x) => rx.test((x as HTMLElement).textContent!.trim())) as HTMLElement | undefined;
    if (!h || !h.classList.contains('expanded')) return false;
    let el = h.nextElementSibling;
    while (el && !el.classList.contains('d4-accordion-pane-content')) el = el.nextElementSibling;
    if (!el) return false;
    const hasText = (el.textContent || '').trim().length > 0;
    const hasVisual = !!el.querySelector('canvas, img, [name="viewer-Grid"], .d4-accordion-pane-header');
    return hasText || hasVisual;
  }, reSource, {timeout: 30_000}).catch(() => {});
}

async function paneContent(page: Page, name: string) {
  return page.evaluate((paneName) => {
    const header = Array.from(document.querySelectorAll('.grok-prop-panel .d4-accordion-pane-header'))
      .find((x) => (x as HTMLElement).textContent!.trim() === paneName) as HTMLElement | undefined;
    let el = header ? header.nextElementSibling : null;
    while (el && !el.classList.contains('d4-accordion-pane-content')) el = el.nextElementSibling;
    if (!el) return {found: false, hasGrid: false, hasMixfileVersion: false, innerPanes: 0, text: ''};
    return {
      found: true,
      hasGrid: !!el.querySelectorAll('[name="viewer-Grid"]').length,
      hasMixfileVersion: /mixfileVersion/.test(el.textContent || ''),
      innerPanes: el.querySelectorAll('.d4-accordion-pane-header').length,
      text: (el.textContent || '').trim(),
    };
  }, name);
}

interface GraphicMetrics {
  found: boolean;
  boxWidth: number;
  boxHeight: number;
  imgWidth: number;
  imgHeight: number;
  payloadLength: number;
  payloadHash: string;
  opaquePixels: number;
  inkPixels: number;
}

// How the selector was established: `Chemistry | Gasteiger Partial Charges` declares
// `#output: graphics charges` (Chem/scripts/gasteiger_charges.py#L13). The client renders a
// GRAPHICS run parameter as `div('grok-scripting-image-container-info-panel')` whose
// backgroundImage is a `data:image/png;base64,…` URL — core/client/d4/lib/src/widgets/ui.dart#L1552-1556.
// So the panel's graphic is neither a <canvas> nor an <img>: it is that div, and the bytes the
// PRODUCT produced are the base64 payload. `.d4-pane-<name>` scoping per chem.md:1877-1879.
async function gasteigerGraphic(page: Page): Promise<GraphicMetrics> {
  return page.evaluate(async () => {
    const empty = {found: false, boxWidth: 0, boxHeight: 0, imgWidth: 0, imgHeight: 0,
      payloadLength: 0, payloadHash: '', opaquePixels: 0, inkPixels: 0};
    const pane = document.querySelector('.d4-pane-gasteiger_partial_charges') as HTMLElement | null;
    const el = pane
      ? pane.querySelector('.grok-scripting-image-container-info-panel') as HTMLElement | null
      : null;
    if (!el) return empty;
    const rect = el.getBoundingClientRect();
    const bg = el.style.backgroundImage || getComputedStyle(el).backgroundImage || '';
    const m = /url\(["']?(data:image\/png;base64,[^"')]+)["']?\)/.exec(bg);
    const url = m ? m[1] : '';
    let hash = 0;
    for (let i = 0; i < url.length; i++) hash = ((hash * 31) + url.charCodeAt(i)) | 0;

    let imgWidth = 0; let imgHeight = 0; let opaquePixels = 0; let inkPixels = 0;
    if (url) {
      const img = new Image();
      await new Promise<void>((res) => { img.onload = () => res(); img.onerror = () => res(); img.src = url; });
      imgWidth = img.naturalWidth;
      imgHeight = img.naturalHeight;
      if (imgWidth > 0 && imgHeight > 0) {
        const c = document.createElement('canvas');
        c.width = imgWidth;
        c.height = imgHeight;
        const ctx = c.getContext('2d')!;
        ctx.clearRect(0, 0, imgWidth, imgHeight);
        ctx.drawImage(img, 0, 0);
        const d = ctx.getImageData(0, 0, imgWidth, imgHeight).data;
        for (let i = 0; i < d.length; i += 4) {
          // Alpha guard first: a fully transparent pixel is not paint, whatever its RGB says.
          if (d[i + 3] === 0) continue;
          opaquePixels++;
          if (d[i] === 255 && d[i + 1] === 255 && d[i + 2] === 255) continue;
          inkPixels++;
        }
      }
    }
    return {found: true, boxWidth: rect.width, boxHeight: rect.height, imgWidth, imgHeight,
      payloadLength: url.length, payloadHash: String(hash), opaquePixels, inkPixels};
  });
}

// Rebind UNCONDITIONALLY. `grok.shell.o` survives closeAll() + addTableView, so the
// `o.semType !== semType` guard inside openContextPanes is a no-op the second time a Molecule cell
// is needed: the panel would keep measuring the FIRST table's cell and the molblock branch would
// never be exercised. Returns the bound object's own value so the caller can assert WHAT is bound,
// rather than re-reading the dataframe string it just supplied.
async function bindFirstCell(page: Page, semType: string): Promise<string | null> {
  return page.evaluate(async (st) => {
    grok.shell.windows.showContextPanel = true;
    const df = grok.shell.t;
    const col = df.columns.toList().find((c: any) => c.semType === st);
    if (!col) return null;
    grok.shell.o = DG.SemanticValue.fromTableCell(df.cell(0, col.name));
    await new Promise((r) => setTimeout(r, 500));
    const bound = grok.shell.o;
    return bound && bound.value != null ? String(bound.value) : null;
  }, semType);
}

test('Chem: Chemistry Mixture / MixtureTree panels differentiate mixture input and Gasteiger panel runs error-free',
  async ({page}) => {
    test.setTimeout(900_000);

    // Carried across Scenario 2's two steps: the molblock render is compared against it, which is
    // what makes the dual-format claim about the PRODUCT'S output rather than about the input string.
    let smilesGraphic: GraphicMetrics | null = null;

    await loginToDatagrok(page);

    const probeSeen: string[] = [];
    await expect.poll(async () => {
      await withConsoleErrorCount(page, async () => {
        await page.evaluate((marker) => console.error(marker), CONSOLE_PROBE);
      }, 1000, probeSeen);
      return probeSeen.filter((t) => t.includes(CONSOLE_PROBE)).length;
    }, {timeout: 15_000, intervals: [100, 200, 300, 500], message:
      'positive control: withConsoleErrorCount must deliver a deliberately emitted console.error ' +
      'before any step relies on it — a detached listener would satisfy every "no new console ' +
      'error" assertion below by receiving nothing'}).toBeGreaterThan(0);
    console.log(`[console] probe capture = ${JSON.stringify(probeSeen)}`);

    await softStep('Setup: open test_mixtures.csv (mixture column, semType ChemicalMixture)', async () => {
      await page.evaluate(async (path) => {
        try { grok.shell.settings.showFiltersIconsConstantly = true; } catch (e) {}
        try { grok.shell.windows.simpleMode = true; } catch (e) {}
        grok.shell.closeAll();
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
      }, MIXTURES_PATH);
      await waitForMolecule(page).catch(() => {});
      const name = await page.evaluate(() => {
        const c = grok.shell.t.columns.toList().find((x: any) => x.semType === 'ChemicalMixture');
        return c ? c.name : null;
      });
      expect(name, 'test_mixtures.csv must expose a ChemicalMixture column').toBe('mixture');
    });

    const chem = await chemRdKitReady(page);
    console.log(`[panels] Chem/RDKit readiness: ${chem.detail}`);
    test.skip(!chem.ok,
      `Chem is not usable on this build — ${chem.detail} — environmental, not a panel defect`);

    const missing = await missingPanelFuncs(page);
    console.log(`[panels] panel-function registrations: expected=${JSON.stringify(CHEM_PANEL_FUNCS)} missing=${JSON.stringify(missing)}`);
    await expect.poll(() => missingPanelFuncs(page),
      {timeout: 60_000, intervals: [500, 1000, 1000, 2000], message:
        `Chem must register the panel functions this test is about — an unregistered function is a ` +
        `panel defect, not an environment fact, and it is the same condition the pane assertions ` +
        `below check. expected=${JSON.stringify(CHEM_PANEL_FUNCS)}`})
      .toEqual([]);

    await softStep('Scenario 1 Step 3: Chemistry context shows Mixture and MixtureTree panes for a mixture cell', async () => {
      const headers = await openContextPanes(page, 'ChemicalMixture', 'Mixture');
      expect(headers.includes('Mixture'),
        `Mixture pane absent from context panel. headers=${JSON.stringify(headers)}`).toBe(true);
      expect(headers.includes('MixtureTree'),
        `MixtureTree pane absent from context panel. headers=${JSON.stringify(headers)}`).toBe(true);
    });

    // Mixture renders a TABLE, not component panes — the per-component panes belong to
    // MixtureTree and are asserted in Step 5. This step previously read MixtureTree's inner
    // headers while being labelled a Mixture check, so it measured a different panel.
    await softStep('Scenario 1 Step 4: Mixture panel renders its component table', async () => {
      await expandPane(page, '^Mixture$');
      const mixture = await paneContent(page, 'Mixture');
      console.log(`[panels] Mixture pane = ${JSON.stringify({found: mixture.found, hasGrid: mixture.hasGrid, innerPanes: mixture.innerPanes})}`);
      expect(mixture.found, 'Mixture pane must be present once the Chemistry group is expanded').toBe(true);
      expect(mixture.hasGrid, 'Mixture panel must render its component table').toBe(true);
    });

    await softStep('Scenario 1 Step 5: MixtureTree renders mixfileVersion header and >1 component pane', async () => {
      await expandPane(page, '^MixtureTree$');
      const tree = await paneContent(page, 'MixtureTree');
      expect(tree.hasMixfileVersion,
        `MixtureTree must show the mixfileVersion header. text=${tree.text.slice(0, 80)}`).toBe(true);
      expect(tree.innerPanes,
        `MixtureTree must render >1 component pane (got ${tree.innerPanes})`).toBeGreaterThan(1);
    });

    // There is no "single-component Mixture" state to observe: the Mixture/MixtureTree panel
    // parameters are typed ChemicalMixture, so for a Molecule cell those panels are not offered
    // at all. The claim is therefore about WHICH panels the Chemistry group offers, and the
    // absence below is only meaningful next to the positive control that molecular panels ARE
    // offered in the same expanded group.
    await softStep('Scenario 1 Step 6: for a molecule cell the Chemistry group offers molecular panels and not Mixture/MixtureTree', async () => {
      await page.evaluate(async (path) => {
        grok.shell.closeAll();
        const df = await grok.dapi.files.readCsv(path);
        grok.shell.addTableView(df);
      }, SMILES_PATH);
      await waitForMolecule(page).catch(() => {});
      const value = await page.evaluate(() => {
        const c = grok.shell.t.columns.toList().find((x: any) => x.semType === 'Molecule');
        return c ? String(grok.shell.t.get(c.name, 0)) : null;
      });
      expect(value, 'smiles.csv must expose a Molecule column').not.toBeNull();
      const headers = await openContextPanes(page, 'Molecule', 'Gasteiger Partial Charges');
      expect(headers.some((h) => /^Chemistry$/.test(h)),
        `Molecule cell must open the Chemistry context group. headers=${JSON.stringify(headers)}`).toBe(true);
      // Positive control: the group is expanded and really offering molecular panels. Without
      // it, "no Mixture pane" is equally satisfied by a group that renders nothing at all.
      expect(headers.includes('Gasteiger Partial Charges'),
        `the expanded Chemistry group must offer molecular panels for a Molecule cell, or the ` +
        `absence asserted below proves nothing. headers=${JSON.stringify(headers)}`).toBe(true);
      expect(headers.includes('Mixture') || headers.includes('MixtureTree'),
        `A molecule cell must NOT be offered the mixture panels — their parameters are typed ` +
        `ChemicalMixture. headers=${JSON.stringify(headers)}`)
        .toBe(false);
    });

    await softStep('Scenario 2 Step 3: Gasteiger Partial Charges panel renders a non-empty charge graphic for a SMILES cell without a script error', async () => {
      const value = await page.evaluate(() => {
        const c = grok.shell.t.columns.toList().find((x: any) => x.semType === 'Molecule');
        return c ? String(grok.shell.t.get(c.name, 0)) : null;
      });
      expect(value, 'smiles.csv must expose a Molecule column').not.toBeNull();
      const smilesConsole: string[] = [];
      await withConsoleErrorCount(page, async () => {
        await openContextPanes(page, 'Molecule', 'Gasteiger Partial Charges');
        await expandPane(page, 'gasteiger|partial charge');
        // The Python script runs server-side; the div only carries a payload once it returns.
        await page.waitForFunction(() => {
          const pane = document.querySelector('.d4-pane-gasteiger_partial_charges');
          const el = pane ? pane.querySelector('.grok-scripting-image-container-info-panel') as HTMLElement | null : null;
          return !!el && /data:image\/png;base64,/.test(el.style.backgroundImage || '');
        }, undefined, {timeout: 180_000, polling: 1000});
      }, CONSOLE_SETTLE_MS, smilesConsole);

      smilesGraphic = await gasteigerGraphic(page);
      console.log(`[gasteiger] SMILES graphic = ${JSON.stringify(smilesGraphic)}`);
      expect(smilesGraphic.found,
        'the Gasteiger pane must hold a rendered graphics host (.grok-scripting-image-container-info-panel)').toBe(true);
      expect(Math.min(smilesGraphic.boxWidth, smilesGraphic.boxHeight),
        `the graphic element must have non-zero dimensions on screen (got ${smilesGraphic.boxWidth}x${smilesGraphic.boxHeight})`)
        .toBeGreaterThan(0);
      expect(Math.min(smilesGraphic.imgWidth, smilesGraphic.imgHeight),
        `the decoded PNG must have non-zero dimensions (got ${smilesGraphic.imgWidth}x${smilesGraphic.imgHeight})`)
        .toBeGreaterThan(0);
      // The alpha guard runs before the white comparison, which is the right order — but on this
      // payload it never fires: the decoded PNG carries no transparent pixels, so opaque equals the
      // full pixel count. It is kept for the day the payload does carry alpha, not because it is
      // doing work here.
      expect(smilesGraphic.inkPixels,
        `the charge graphic must carry non-white paint. This is a NON-EMPTY FLOOR, not proof of a ` +
        `contour — a skeleton drawn without contours would pass it too. The panel draws a matplotlib ` +
        `figure on white, so an unpainted result is all-white and scores ink=0 ` +
        `(got ink=${smilesGraphic.inkPixels}, opaque=${smilesGraphic.opaquePixels}, ` +
        `png=${smilesGraphic.imgWidth}x${smilesGraphic.imgHeight})`)
        .toBeGreaterThan(0);

      const relevant = nonAmbient(smilesConsole);
      console.log(`[console] SMILES step captured ${smilesConsole.length} errors, non-ambient=${JSON.stringify(relevant)}`);
      expect(relevant,
        `Gasteiger Partial Charges panel must run for a SMILES cell without a console script error. ` +
        `Window: the render callback plus ${CONSOLE_SETTLE_MS} ms — later errors are outside this claim. ` +
        `The channel itself was proved live by the CONSOLE_PROBE control at the top of this test. ` +
        `non-ambient=${JSON.stringify(relevant)}`).toEqual([]);
    });

    await softStep('Scenario 2 Step 4: Gasteiger Partial Charges panel renders a DIFFERENT charge graphic for a molV2000 cell without a script error', async () => {
      await page.evaluate(async (path) => {
        grok.shell.closeAll();
        const bytes = await grok.dapi.files.readAsBytes(path);
        const dfs = await grok.functions.call('Chem:importSdf', {bytes});
        const df = Array.isArray(dfs) ? dfs[0] : dfs;
        grok.shell.addTableView(df);
      }, MOLBLOCK_SDF_PATH);
      await waitForMolecule(page).catch(() => {});

      // Bind first and unconditionally, then read back what is ACTUALLY bound. The previous
      // version read the string it had just pulled out of its own dataframe, which says nothing
      // about which cell the panel is rendering.
      const bound = await bindFirstCell(page, 'Molecule');
      console.log(`[gasteiger] bound current object length=${bound ? bound.length : 0} hasMEnd=${!!bound && bound.includes('M  END')}`);
      expect(bound, 'the SDF table must expose a Molecule cell and it must become the current object').not.toBeNull();
      expect(bound!.includes('M  END'),
        `grok.shell.o must be the molV2000 cell the panel will render, not the previous SMILES ` +
        `binding. bound="${bound!.slice(0, 60)}"`).toBe(true);

      const molblockConsole: string[] = [];
      await withConsoleErrorCount(page, async () => {
        await openContextPanes(page, 'Molecule', 'Gasteiger Partial Charges');
        await expandPane(page, 'gasteiger|partial charge');
        await page.waitForFunction((prevHash) => {
          const pane = document.querySelector('.d4-pane-gasteiger_partial_charges');
          const el = pane ? pane.querySelector('.grok-scripting-image-container-info-panel') as HTMLElement | null : null;
          const bg = el ? el.style.backgroundImage || '' : '';
          if (!/data:image\/png;base64,/.test(bg)) return false;
          const url = /url\(["']?(data:image\/png;base64,[^"')]+)["']?\)/.exec(bg)![1];
          let h = 0;
          for (let i = 0; i < url.length; i++) h = ((h * 31) + url.charCodeAt(i)) | 0;
          return String(h) !== prevHash;
        }, smilesGraphic ? smilesGraphic.payloadHash : '', {timeout: 180_000, polling: 1000})
          // The barrier computes the same predicate the assertion below names, and there is no
          // weaker settle signal on this pane (the script's completion is only visible as the new
          // payload). So on expiry, re-measure and rethrow with both hashes — otherwise the
          // diagnosis arrives as a bare 180 s timeout.
          .catch(async () => {
            const stuck = await gasteigerGraphic(page);
            throw new Error(
              `the molV2000 render never produced a payload different from the SMILES one within 180s: ` +
              `smiles hash=${smilesGraphic ? smilesGraphic.payloadHash : '(none)'} ` +
              `len=${smilesGraphic ? smilesGraphic.payloadLength : 0}, ` +
              `current hash=${stuck.payloadHash} len=${stuck.payloadLength} found=${stuck.found}`);
          });
      }, CONSOLE_SETTLE_MS, molblockConsole);

      const molblockGraphic = await gasteigerGraphic(page);
      console.log(`[gasteiger] molblock graphic = ${JSON.stringify(molblockGraphic)}`);
      expect(smilesGraphic,
        'Step 3 must have captured the SMILES graphic — without it the dual-format comparison below is unmade')
        .not.toBeNull();
      expect(molblockGraphic.found,
        'the Gasteiger pane must hold a rendered graphics host for the molV2000 cell').toBe(true);
      expect(Math.min(molblockGraphic.boxWidth, molblockGraphic.boxHeight),
        `the graphic element must have non-zero dimensions on screen (got ${molblockGraphic.boxWidth}x${molblockGraphic.boxHeight})`)
        .toBeGreaterThan(0);
      expect(Math.min(molblockGraphic.imgWidth, molblockGraphic.imgHeight),
        `the decoded PNG must have non-zero dimensions (got ${molblockGraphic.imgWidth}x${molblockGraphic.imgHeight})`)
        .toBeGreaterThan(0);
      expect(molblockGraphic.inkPixels,
        `the molV2000 charge graphic must carry non-white paint. This is a NON-EMPTY FLOOR, not proof ` +
        `of a contour — a skeleton drawn without contours would pass it too; the dual-format claim is ` +
        `carried by the hash difference below, not by this count ` +
        `(got ink=${molblockGraphic.inkPixels}, opaque=${molblockGraphic.opaquePixels}, ` +
        `png=${molblockGraphic.imgWidth}x${molblockGraphic.imgHeight})`)
        .toBeGreaterThan(0);
      // The product-produced discriminator: the PNG bytes the script emitted for the molblock input
      // differ from the ones it emitted for the SMILES input. Identical bytes mean the panel is
      // still rendering the previous binding — the exact false green this step used to hide.
      console.log(`[gasteiger] payloadHash smiles=${smilesGraphic!.payloadHash} molblock=${molblockGraphic.payloadHash} ` +
        `len smiles=${smilesGraphic!.payloadLength} molblock=${molblockGraphic.payloadLength}`);
      expect(molblockGraphic.payloadHash,
        `the molV2000 render must be a NEW graphic produced from the molblock, not the SMILES ` +
        `render still on screen (smiles hash=${smilesGraphic!.payloadHash} len=${smilesGraphic!.payloadLength}, ` +
        `molblock hash=${molblockGraphic.payloadHash} len=${molblockGraphic.payloadLength})`)
        .not.toBe(smilesGraphic!.payloadHash);

      const relevant = nonAmbient(molblockConsole);
      console.log(`[console] molblock step captured ${molblockConsole.length} errors, non-ambient=${JSON.stringify(relevant)}`);
      expect(relevant,
        `Gasteiger Partial Charges panel must run for a molV2000 cell without a console script error. ` +
        `Window: the render callback plus ${CONSOLE_SETTLE_MS} ms — later errors are outside this claim. ` +
        `The channel itself was proved live by the CONSOLE_PROBE control at the top of this test. ` +
        `non-ambient=${JSON.stringify(relevant)}`).toEqual([]);
    });

    await page.evaluate(() => grok.shell.closeAll());
    finishSpec();
  });
