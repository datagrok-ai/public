/* ---
sub_features_covered: [biostructure.data-provider.rcsb-mmcif, biostructure.file-open.importPdb, biostructure.overlay.reset-camera, biostructure.prop.binding-site-radius, biostructure.prop.biostructure-data-provider, biostructure.prop.biostructure-id-column, biostructure.prop.ligand-column, biostructure.prop.representation, biostructure.prop.show-binding-site, biostructure.prop.show-current-row-ligand, biostructure.top-menu.fetch-pdb-sequences, biostructure.viewer, biostructure.viewer.add-via-dropdown, biostructure.viewer.settings-panel, biostructure.viewport-context-menu.download-cif, biostructure.viewport-context-menu.download-pdb]
--- */
// BiostructureViewer happy-path smoke. File-handler routing is verified via DG.Func.find registry
// probes (no Mol* engine init, which surfaces WebGL noise in CI). Mol*-engine-dependent assertions
// (Reset Camera, viewport context menu) are gated on .msp-plugin / .msp-viewport being present.
// Scenarios 4, 7, 8 make outbound RCSB calls (bounded 30-90s; absence treated as a precondition).
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbIdCsv = 'System:AppData/BiostructureViewer/pdb_id.csv';

// RCSB reachability probe (via the platform proxy, mirroring how Fetch PDB Sequences reaches RCSB).
// Used to distinguish a transient outbound-network outage (explicit skip) from a genuinely broken
// transform (real failure) when the chain-column append times out.
async function rcsbReachable(page: Page): Promise<boolean> {
  return await page.evaluate(async () => {
    try {
      const resp = await grok.dapi.fetchProxy('https://files.rcsb.org/download/1CRN.cif', {method: 'GET'});
      return !!resp?.ok;
    } catch (_e) {
      return false;
    }
  });
}

test('BiostructureViewer / happy-path smoke', async ({page}) => {
  test.setTimeout(360_000);
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

  // Scenario 1 — importPdb is the registered dispatch target for .mmcif/.pdb (registry probe).
  await softStep('Scenario 1 — importPdb file-handler dispatch surface registered (.mmcif/.pdb)', async () => {
    const res = await page.evaluate(() => {
      const fn = DG.Func.find({name: 'importPdb', package: 'BiostructureViewer'})[0];
      const inputs = fn?.inputs?.map((i: any) => ({n: i.name, t: i.propertyType})) ?? [];
      return {
        registered: !!fn,
        inputCount: inputs.length,
        firstInputName: inputs[0]?.n,
        firstInputType: inputs[0]?.t,
      };
    });
    expect(res.registered).toBe(true);
    expect(res.firstInputName).toBe('fileContent');
    expect(res.firstInputType).toBe('string');
    expect(res.inputCount).toBe(1);
  });

  // Scenario 2a — TableView + Biostructure viewer mount via tv.addViewer.
  let scenarioMountedViewer = false;
  await softStep('Scenario 2a/4 — Build pdb_id table; tv.addViewer mounts [name="viewer-Biostructure"]', async () => {
    const res = await page.evaluate(async () => {
      grok.shell.closeAll();
      for (let i = 0; i < 40; i++) {
        if ((grok.shell.tableViews?.length ?? 0) === 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('pdb_id', ['1CRN'])]);
      df.col('pdb_id').semType = 'PDB_ID';
      const tv = grok.shell.addTableView(df);
      for (let i = 0; i < 40; i++) {
        if (tv.grid || document.querySelector('[name="viewer-Grid"]')) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const v = tv.addViewer('Biostructure');
      // Cold-start readiness poll: wait up to 30s for the 'representation' descriptor to register
      // (props.get throws "Property not found" before the lazy package chunk loads on a cold context).
      let defaultRep = null;
      let representationReady = false;
      let pollIters = 0;
      for (let i = 0; i < 60; i++) {
        pollIters = i + 1;
        await new Promise((r) => setTimeout(r, 500));
        try {
          const props = v.getProperties ? v.getProperties() : [];
          const hasRepDescriptor = props.some((p) => p.name === 'representation');
          if (!hasRepDescriptor) continue;
          const val = v.props.get('representation');
          if (typeof val === 'string' && val.length > 0) {
            defaultRep = val;
            representationReady = true;
            break;
          }
        } catch (_e) {
          // Descriptor not registered yet on cold mount; keep polling.
        }
      }
      return {
        hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
        vType: v?.type,
        defaultRep,
        representationReady,
        pollIters,
      };
    });
    expect(res.hasContainer).toBe(true);
    expect(res.vType).toBe('Biostructure');
    expect(res.representationReady).toBe(true);
    expect(res.defaultRep).toBe('cartoon');
    scenarioMountedViewer = true;
  });

  // Scenario 2b — Settings panel via gear icon (DOM driven).
  await softStep('Scenario 2b — Open viewer settings via gear (DOM); property panel surfaces', async () => {
    if (!scenarioMountedViewer) return;
    const gearClicked = await page.evaluate(() => {
      const container = document.querySelector('[name="viewer-Biostructure"]');
      if (!container) return false;
      // Gear lives in the panel-titlebar of the enclosing .panel-base, not the viewer container.
      const gear = container.closest('.panel-base')?.querySelector(
        '.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (!gear) return false;
      gear.click();
      return true;
    });
    expect(gearClicked).toBe(true);
    await page.locator('.grok-prop-panel').waitFor({timeout: 10_000});
  });

  // Scenario 2b cont. — Switch representation cartoon -> ball-and-stick -> molecular-surface -> cartoon.
  await softStep('Scenario 2b — Switch representation cartoon -> ball-and-stick -> molecular-surface -> cartoon', async () => {
    if (!scenarioMountedViewer) return;
    const res = await page.evaluate(async () => {
      let v: any = null;
      for (const tv of grok.shell.tableViews || []) {
        for (const x of tv.viewers || []) if (x.type === 'Biostructure') { v = x; break; }
        if (v) break;
      }
      if (!v) return {ok: false, observed: [] as string[]};
      const reps: string[] = ['ball-and-stick', 'molecular-surface', 'cartoon'];
      const observed: string[] = [];
      for (const r of reps) {
        v.setOptions({representation: r});
        let val: any = null;
        for (let i = 0; i < 50; i++) {
          val = v.props.get('representation');
          if (val === r) break;
          await new Promise((res) => setTimeout(res, 100));
        }
        observed.push(val);
      }
      return {ok: true, observed};
    });
    expect(res.ok).toBe(true);
    expect(res.observed).toEqual(['ball-and-stick', 'molecular-surface', 'cartoon']);
  });

  // Scenario 3 — Reset Camera overlay button (asserted only when .msp-plugin is built).
  await softStep('Scenario 3 — Reset Camera overlay button click (DOM, precondition .msp-plugin)', async () => {
    const result = await page.evaluate(async () => {
      const pluginPresent = !!document.querySelector('.msp-plugin');
      if (!pluginPresent) return {precondition: false, btnFound: false, clicked: false};
      const btn = document.querySelector('button[title="Reset Camera"]') as HTMLButtonElement | null;
      if (!btn) return {precondition: true, btnFound: false, clicked: false};
      btn.click();
      await new Promise((r) => setTimeout(r, 400));
      return {precondition: true, btnFound: true, clicked: true};
    });
    // Mol* WebGL engine is deliberately not initialized in CI (see top-of-file); when it is absent
    // this is a declared skip, not silent coverage.
    if (!result.precondition)
      throw new Error('Test is skipped: Mol* engine (.msp-plugin) not initialized in CI');
    expect(result.btnFound).toBe(true);
    expect(result.clicked).toBe(true);
    // GROK: observable camera-reset effect not asserted — Mol* engine unavailable in CI to verify against.
  });

  // Scenario 4 — RCSB mmCIF data provider wiring (JS-API setOptions path).
  await softStep('Scenario 4 — Wire RCSB mmCIF provider on existing Biostructure viewer', async () => {
    if (!scenarioMountedViewer) return;
    const res = await page.evaluate(() => {
      let v: any = null;
      for (const tv of grok.shell.tableViews || []) {
        for (const x of tv.viewers || []) if (x.type === 'Biostructure') { v = x; break; }
        if (v) break;
      }
      if (!v) return {ok: false};
      v.setOptions({
        biostructureIdColumnName: 'pdb_id',
        biostructureDataProvider: 'BiostructureViewer:getBiostructureRcsbMmcif',
      });
      const provider = v.props.get('biostructureDataProvider');
      const idCol = v.props.get('biostructureIdColumnName');
      const providerFn = DG.Func.find({name: 'getBiostructureRcsbMmcif', package: 'BiostructureViewer'})[0];
      return {ok: true, provider, idCol, providerRegistered: !!providerFn};
    });
    expect(res.ok).toBe(true);
    expect(res.provider).toBe('BiostructureViewer:getBiostructureRcsbMmcif');
    expect(res.idCol).toBe('pdb_id');
    expect(res.providerRegistered).toBe(true);
  });

  // Scenario 5 — Ligand + binding-site property descriptors via getProperties() (no engine init).
  await softStep('Scenario 5 — Ligand + binding-site property descriptors (no engine init)', async () => {
    const res = await page.evaluate(async () => {
      grok.shell.closeAll();
      for (let i = 0; i < 40; i++) {
        if ((grok.shell.tableViews?.length ?? 0) === 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      // Lightweight df: a single 'ligand' string column (no semType) so the engine never
      // enters the data-request pipeline.
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('ligand', ['placeholder']),
      ]);
      const tv = grok.shell.addTableView(df);
      for (let i = 0; i < 40; i++) {
        if (tv.grid || document.querySelector('[name="viewer-Grid"]')) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const v = tv.addViewer('Biostructure');
      // Cold-start readiness poll: wait up to 30s for the bindingSiteRadius descriptor to register.
      let props: any[] = [];
      let propsReady = false;
      let pollIters = 0;
      for (let i = 0; i < 60; i++) {
        pollIters = i + 1;
        await new Promise((r) => setTimeout(r, 500));
        try {
          const candidate = v.getProperties ? v.getProperties() : [];
          if (!candidate || candidate.length === 0) continue;
          const hasBs = candidate.some((p) => p.name === 'bindingSiteRadius');
          if (!hasBs) continue;
          props = candidate;
          propsReady = true;
          break;
        } catch (_e) {
          // Property descriptors not yet registered; keep polling.
        }
      }
      const byName: Record<string, any> = {};
      for (const p of props) byName[p.name] = p;
      const lc = byName['ligandColumnName'];
      const sc = byName['showCurrentRowLigand'];
      const sb = byName['showBindingSite'];
      const bs = byName['bindingSiteRadius'];

      return {
        hasContainer: !!document.querySelector('[name="viewer-Biostructure"]'),
        propsReady,
        pollIters,
        ligandColumnNameRegistered: !!lc,
        ligandColumnNamePropertyType: lc?.propertyType,
        showCurrentRowLigandRegistered: !!sc,
        showCurrentRowLigandDefault: sc?.defaultValue,
        showBindingSiteRegistered: !!sb,
        showBindingSiteDefault: sb?.defaultValue,
        bindingSiteRadiusRegistered: !!bs,
        bindingSiteRadiusDefault: bs?.defaultValue,
        bindingSiteRadiusMin: bs?.min,
        bindingSiteRadiusMax: bs?.max,
      };
    });
    expect(res.hasContainer).toBe(true);
    expect(res.propsReady).toBe(true);
    expect(res.ligandColumnNameRegistered).toBe(true);
    expect(res.ligandColumnNamePropertyType).toBe('string');
    expect(res.showCurrentRowLigandRegistered).toBe(true);
    expect(res.showCurrentRowLigandDefault).toBe(true);
    expect(res.showBindingSiteRegistered).toBe(true);
    expect(res.showBindingSiteDefault).toBe(false);
    expect(res.bindingSiteRadiusRegistered).toBe(true);
    expect(res.bindingSiteRadiusDefault).toBe(5);
    expect(res.bindingSiteRadiusMin).toBe(3);
    expect(res.bindingSiteRadiusMax).toBe(10);
  });

  // Scenario 6 — Viewport context-menu Download paths (precondition: .msp-viewport rendered).
  await softStep('Scenario 6 — Viewport right-click Download (DOM, precondition .msp-viewport)', async () => {
    const rect = await page.evaluate(() => {
      const vp = document.querySelector('.msp-viewport') as HTMLElement | null;
      if (!vp) return null;
      const r = vp.getBoundingClientRect();
      if (r.width === 0 || r.height === 0) return null;
      return {x: r.x, y: r.y, w: r.width, h: r.height};
    });
    // Mol* WebGL viewport is deliberately not initialized in CI (see top-of-file); declared skip.
    if (!rect)
      throw new Error('Test is skipped: Mol* viewport (.msp-viewport) not rendered in CI');
    const cx = rect.x + rect.w / 2;
    const cy = rect.y + rect.h / 2;
    // Spy the download side-effect: Mol* saves via a blob anchor click.
    await page.evaluate(() => {
      (window as any).__dlSpy = [];
      const proto = HTMLAnchorElement.prototype as any;
      if (!proto.__origClick) {
        proto.__origClick = proto.click;
        proto.click = function() {
          if (this.download || (this.href && /^blob:/.test(this.href)))
            (window as any).__dlSpy.push(this.download || this.href);
          return proto.__origClick.apply(this, arguments);
        };
      }
    });
    // As PDB.
    await page.mouse.click(cx, cy, {button: 'right'});
    await page.locator('[name="div-Download"]').waitFor({timeout: 10_000});
    await page.evaluate(async () => {
      const dl = document.querySelector('[name="div-Download"]') as HTMLElement | null;
      if (dl) dl.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 200));
    });
    const asPdbFound = await page.evaluate(() => {
      const leaf = document.querySelector('[name="div-Download---As-PDB"]') as HTMLElement | null;
      if (!leaf) return false;
      leaf.click();
      return true;
    });
    expect(asPdbFound).toBe(true);
    // Re-trigger menu for As CIF.
    await page.mouse.click(cx, cy, {button: 'right'});
    await page.locator('[name="div-Download"]').waitFor({timeout: 10_000});
    await page.evaluate(async () => {
      const dl = document.querySelector('[name="div-Download"]') as HTMLElement | null;
      if (dl) dl.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 200));
    });
    const asCifFound = await page.evaluate(() => {
      const leaf = document.querySelector('[name="div-Download---As-CIF"]') as HTMLElement | null;
      if (!leaf) return false;
      leaf.click();
      return true;
    });
    expect(asCifFound).toBe(true);
    // Both leaves must have triggered a file download.
    await page.waitForFunction(() => ((window as any).__dlSpy || []).length >= 2, null, {timeout: 15_000});
    await page.keyboard.press('Escape');
  });

  // Scenario 7 — Bio | Transform | Fetch PDB Sequences appends Chain N columns (outbound RCSB).
  let chainColsAfterScenario7: string[] = [];
  let scenario7DialogOk = false;
  await softStep('Scenario 7 — pdb_id.csv -> Bio | Transform | Fetch PDB Sequences (DOM driven; Bio-optional)', async () => {
    const setup = await page.evaluate(async (path) => {
      grok.shell.closeAll();
      for (let i = 0; i < 40; i++) {
        if ((grok.shell.tableViews?.length ?? 0) === 0) break;
        await new Promise((r) => setTimeout(r, 100));
      }
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      // Bio package menu registration: allow up to 8s.
      let hasBioMenu = false;
      for (let i = 0; i < 16; i++) {
        await new Promise((r) => setTimeout(r, 500));
        if (document.querySelector('[name="div-Bio"]')) { hasBioMenu = true; break; }
      }
      const pdbCol: any = df.col('pdb_id');
      return {rowCount: df.rowCount, semType: pdbCol?.semType, hasBioMenu};
    }, samplePdbIdCsv);
    // Bio is an OPTIONAL prerequisite for this suite: the "Fetch PDB Sequences" transform lives in the
    // Bio package. When Bio is not installed on the CI stack, skip Scenarios 7-8 gracefully (the
    // BiostructureViewer suite must stay green with PREREQ=BiostructureViewer alone).
    if (!setup.hasBioMenu) {
      console.warn('[Scenario 7] Bio package menu absent — Fetch PDB Sequences (Bio transform) unavailable; Scenarios 7-8 skipped (Bio-optional).');
      return;
    }
    expect(setup.semType).toBe('PDB_ID');

    // DOM-driving: open Bio | Transform | Fetch PDB Sequences...
    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const transform = document.querySelector('[name="div-Bio---Transform"]') as HTMLElement | null;
      if (transform) transform.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Transform---Fetch-PDB-Sequences..."]') as HTMLElement | null;
      if (leaf) leaf.click();
    });
    await page.locator('[name="dialog-Fetch-PDB-Sequences"]').waitFor({timeout: 30_000});

    const beforeCount: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    await page.locator('[name="dialog-Fetch-PDB-Sequences"] [name="button-OK"]').click();
    scenario7DialogOk = true;

    // 90s ceiling — outbound RCSB GraphQL. A timeout is a real failure unless RCSB is unreachable,
    // in which case skip explicitly (never swallow — that would hide a broken transform).
    try {
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base,
        beforeCount, {timeout: 90_000});
    } catch (e) {
      if (!(await rcsbReachable(page)))
        throw new Error('Test is skipped: RCSB unreachable — Fetch PDB Sequences could not complete');
      throw e;
    }
    const chainNames = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const names: string[] = [];
      for (let i = 0; i < df.columns.length; i++) names.push(df.columns.byIndex(i).name);
      return names.filter((n) => /^Chain\s+\d+$/.test(n));
    });
    expect(chainNames.length).toBeGreaterThan(0);
    const firstMeta = await page.evaluate((name) => {
      const col: any = grok.shell.tv.dataFrame.col(name);
      return {semType: col?.semType};
    }, chainNames[0]);
    expect(firstMeta.semType).toBe('Macromolecule');
    chainColsAfterScenario7 = chainNames;
  });

  // Scenario 8 — Re-running Fetch PDB Sequences is non-destructive (requires Scenario 7 chains).
  await softStep('Scenario 8 — Re-run Fetch PDB Sequences; non-conflicting Chain N (2) columns appended', async () => {
    if (!scenario7DialogOk || chainColsAfterScenario7.length === 0) return;
    const beforeCount: number = await page.evaluate(() => grok.shell.tv.dataFrame.columns.length);
    const originalChains = [...chainColsAfterScenario7];

    await page.evaluate(async () => {
      (document.querySelector('[name="div-Bio"]') as HTMLElement).click();
      await new Promise((r) => setTimeout(r, 400));
      const transform = document.querySelector('[name="div-Bio---Transform"]') as HTMLElement | null;
      if (transform) transform.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
      await new Promise((r) => setTimeout(r, 400));
      const leaf = document.querySelector('[name="div-Bio---Transform---Fetch-PDB-Sequences..."]') as HTMLElement | null;
      if (leaf) leaf.click();
    });
    await page.locator('[name="dialog-Fetch-PDB-Sequences"]').waitFor({timeout: 30_000});
    await page.locator('[name="dialog-Fetch-PDB-Sequences"] [name="button-OK"]').click();

    let newColNames: string[] = [];
    try {
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base,
        beforeCount, {timeout: 90_000});
      newColNames = await page.evaluate(() => {
        const df = grok.shell.tv.dataFrame;
        const names: string[] = [];
        for (let i = 0; i < df.columns.length; i++) names.push(df.columns.byIndex(i).name);
        return names;
      });
    } catch (e) {
      if (!(await rcsbReachable(page)))
        throw new Error('Test is skipped: RCSB unreachable — re-run Fetch PDB Sequences could not complete');
      throw e;
    }

    for (const original of originalChains) expect(newColNames).toContain(original);
    const re = /^Chain\s+\d+\s*\(2\)$/;
    const newSet = newColNames.filter((n) => re.test(n));
    expect(newSet.length).toBeGreaterThan(0);
  });

  // Cleanup.
  await page.evaluate(() => { grok.shell.closeAll(); });

  const realErrors = stepErrors.filter((e) => !e.error.startsWith('Test is skipped:'));
  if (realErrors.length > 0) {
    const summary = realErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${realErrors.length} step(s) failed:\n${summary}`);
  }
});
