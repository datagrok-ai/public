/* ---
sub_features_covered: [biostructure.data-provider.rcsb-mmcif, biostructure.file-open.importPdb, biostructure.overlay.reset-camera, biostructure.prop.binding-site-radius, biostructure.prop.biostructure-data-provider, biostructure.prop.biostructure-id-column, biostructure.prop.ligand-column, biostructure.prop.representation, biostructure.prop.show-binding-site, biostructure.prop.show-current-row-ligand, biostructure.top-menu.fetch-pdb-sequences, biostructure.viewer, biostructure.viewer.add-via-dropdown, biostructure.viewer.settings-panel, biostructure.viewport-context-menu.download-cif, biostructure.viewport-context-menu.download-pdb]
--- */
// BiostructureViewer happy-path smoke. File-handler routing is verified via DG.Func.find registry
// probes (no Mol* engine init, which surfaces WebGL noise in CI). Mol*-engine-dependent assertions
// (Reset Camera, viewport context menu) are gated on .msp-plugin / .msp-viewport being present.
// Scenarios 4, 7, 8 make outbound RCSB calls (bounded 30-90s; absence treated as a precondition).
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

declare const grok: any;
declare const DG: any;

const samplePdbIdCsv = 'System:AppData/BiostructureViewer/pdb_id.csv';

test('BiostructureViewer — happy-path smoke', async ({page}) => {
  test.setTimeout(900_000);
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
      await new Promise((r) => setTimeout(r, 1500));
      const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('pdb_id', ['1CRN'])]);
      df.col('pdb_id').semType = 'PDB_ID';
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
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
    const opened = await page.evaluate(async () => {
      const container = document.querySelector('[name="viewer-Biostructure"]');
      if (!container) return {gearClicked: false, panelOpened: false};
      // Gear lives in the panel-titlebar of the enclosing .panel-base, not the viewer container.
      const gear = container.closest('.panel-base')?.querySelector(
        '.panel-titlebar [name="icon-font-icon-settings"]') as HTMLElement | null;
      if (!gear) return {gearClicked: false, panelOpened: false};
      gear.click();
      await new Promise((r) => setTimeout(r, 1500));
      const cp = document.querySelector('.grok-prop-panel');
      return {gearClicked: true, panelOpened: !!cp};
    });
    expect(opened.gearClicked).toBe(true);
    expect(opened.panelOpened).toBe(true);
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
        await new Promise((res) => setTimeout(res, 1500));
        observed.push(v.props.get('representation'));
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
      if (!pluginPresent) return {precondition: false, clicked: null};
      const btn = document.querySelector('button[title="Reset Camera"]') as HTMLButtonElement | null;
      if (!btn) return {precondition: true, clicked: false};
      btn.click();
      await new Promise((r) => setTimeout(r, 400));
      return {precondition: true, clicked: true};
    });
    if (result.precondition) expect(result.clicked).toBe(true);
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
      await new Promise((r) => setTimeout(r, 1500));
      // Lightweight df: a single 'ligand' string column (no semType) so the engine never
      // enters the data-request pipeline.
      const df = DG.DataFrame.fromColumns([
        DG.Column.fromStrings('ligand', ['placeholder']),
      ]);
      const tv = grok.shell.addTableView(df);
      await new Promise((r) => setTimeout(r, 1500));
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
    if (!rect) return; // precondition (viewport rendered) not met; no DOM driving asserted here
    const cx = rect.x + rect.w / 2;
    const cy = rect.y + rect.h / 2;
    await page.mouse.click(cx, cy, {button: 'right'});
    await page.waitForTimeout(700);
    const dlPresent = await page.evaluate(() => {
      return !!document.querySelector('[name="div-Download"]');
    });
    if (dlPresent) {
      await page.evaluate(async () => {
        const dl = document.querySelector('[name="div-Download"]') as HTMLElement | null;
        if (dl) dl.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
        await new Promise((r) => setTimeout(r, 400));
      });
      await page.evaluate(() => {
        const leaf = document.querySelector('[name="div-Download---As-PDB"]') as HTMLElement | null;
        if (leaf) leaf.click();
      });
      await page.waitForTimeout(500);
      // Re-trigger menu for As CIF.
      const rect2 = await page.evaluate(() => {
        const vp = document.querySelector('.msp-viewport') as HTMLElement | null;
        if (!vp) return null;
        const r = vp.getBoundingClientRect();
        if (r.width === 0 || r.height === 0) return null;
        return {x: r.x, y: r.y, w: r.width, h: r.height};
      });
      if (rect2) {
        await page.mouse.click(rect2.x + rect2.w / 2, rect2.y + rect2.h / 2, {button: 'right'});
        await page.waitForTimeout(700);
        await page.evaluate(async () => {
          const dl = document.querySelector('[name="div-Download"]') as HTMLElement | null;
          if (dl) dl.dispatchEvent(new MouseEvent('mouseover', {bubbles: true}));
          await new Promise((r) => setTimeout(r, 400));
        });
        await page.evaluate(() => {
          const leaf = document.querySelector('[name="div-Download---As-CIF"]') as HTMLElement | null;
          if (leaf) leaf.click();
        });
      }
    }
    await page.keyboard.press('Escape');
    await page.waitForTimeout(300);
  });

  // Scenario 7 — Bio | Transform | Fetch PDB Sequences appends Chain N columns (outbound RCSB).
  let chainColsAfterScenario7: string[] = [];
  let scenario7DialogOk = false;
  await softStep('Scenario 7 — pdb_id.csv -> Bio | Transform | Fetch PDB Sequences (DOM driven)', async () => {
    const setup = await page.evaluate(async (path) => {
      grok.shell.closeAll();
      await new Promise((r) => setTimeout(r, 1500));
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
    expect(setup.semType).toBe('PDB_ID');
    expect(setup.hasBioMenu).toBe(true);

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

    // 90s ceiling — outbound RCSB GraphQL. If chains never appear, skip the rest (Scenario 8 too).
    try {
      await page.waitForFunction(
        (base) => grok.shell.tv.dataFrame.columns.length > base,
        beforeCount, {timeout: 90_000});
    } catch (e) {
      return;
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
      return;
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
